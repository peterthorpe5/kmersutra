#!/usr/bin/env python3
"""Run leakage-aware internal KmerSutra AI calibration validation.

The script validates a KmerSutra evidence-calibration model using held-out
sample, benchmark-family, panel and species splits. It trains only on
sample/species-level KmerSutra evidence features and deliberately excludes
benchmark-only leakage columns such as spike depth and true-target metadata.

The script assumes an AI-ready table produced by ``kmersutra-build-call-training``.
It does not replace KmerSutra's auditable rule layer; it evaluates whether the
optional calibration layer can reproduce or prioritise report-layer evidence
under harder held-out designs.
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import Sequence

from kmersutra_ai_common import (
    DEFAULT_SAFE_FEATURE_COLUMNS,
    add_common_parser_arguments,
    audit_feature_columns,
    configure_logging,
    ensure_output_dir,
    existing_file,
    find_present_features,
    label_counts,
    label_counts_text,
    safe_value,
    sanitise_name,
    stable_group_hash,
    write_json,
)

LABEL_COLUMN = "ml_report_label"
UNKNOWN_LABEL = "unknown_or_unresolved"


def import_kmersutra_helpers():
    """Import KmerSutra helpers lazily.

    Returns
    -------
    tuple
        KmerSutra helper functions.
    """
    from kmersutra.ai_calibration import evaluate_predictions
    from kmersutra.ml import predict_records, save_model, train_prototype_classifier
    from kmersutra.table_io import read_records_table, write_records_table

    return (
        evaluate_predictions,
        predict_records,
        save_model,
        train_prototype_classifier,
        read_records_table,
        write_records_table,
    )


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed arguments.
    """
    parser = argparse.ArgumentParser(
        description="Run leakage-aware KmerSutra AI holdout validation."
    )
    parser.add_argument(
        "--training_table",
        required=True,
        help="AI-ready training table from kmersutra-build-call-training.",
    )
    parser.add_argument(
        "--out_dir",
        required=True,
        help="Output directory for validation files.",
    )
    parser.add_argument(
        "--distance_quantile",
        type=float,
        default=0.95,
        help="Open-set novelty threshold quantile.",
    )
    parser.add_argument(
        "--sample_test_fraction",
        type=float,
        default=0.20,
        help="Deterministic held-out sample fraction for sample-hash split.",
    )
    parser.add_argument(
        "--label_column",
        default=LABEL_COLUMN,
        help="Training label column.",
    )
    parser.add_argument(
        "--unknown_label",
        default=UNKNOWN_LABEL,
        help="Open-set unresolved prediction label.",
    )
    parser.add_argument(
        "--feature_columns",
        nargs="*",
        default=DEFAULT_SAFE_FEATURE_COLUMNS,
        help="Explicit feature columns. Defaults exclude benchmark leakage columns.",
    )
    add_common_parser_arguments(parser)
    return parser.parse_args()


def make_sample_hash_split(
    *,
    records: list[dict[str, object]],
    test_fraction: float,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    """Create a deterministic sample-level hash split.

    Parameters
    ----------
    records : list of dict[str, object]
        Input records.
    test_fraction : float
        Test fraction between zero and one.

    Returns
    -------
    tuple[list[dict[str, object]], list[dict[str, object]]]
        Training and test records.
    """
    if not 0.0 < test_fraction < 1.0:
        raise ValueError("sample_test_fraction must be between 0 and 1")

    train_records: list[dict[str, object]] = []
    test_records: list[dict[str, object]] = []
    threshold = int(test_fraction * 10000)

    for record in records:
        sample_id = safe_value(record, "sample_id")
        bucket = stable_group_hash(sample_id) % 10000
        if bucket < threshold:
            test_records.append(record)
        else:
            train_records.append(record)

    return train_records, test_records


def make_leave_one_value_splits(
    *,
    records: list[dict[str, object]],
    column: str,
    prefix: str,
    label_column: str,
    only_expected_target_species: bool = False,
) -> list[tuple[str, str, str, list[dict[str, object]], list[dict[str, object]]]]:
    """Create leave-one-value-out validation splits.

    Parameters
    ----------
    records : list of dict[str, object]
        Input records.
    column : str
        Metadata column to hold out.
    prefix : str
        Validation split name prefix.
    label_column : str
        Label column.
    only_expected_target_species : bool, optional
        Restrict split values to species with expected-target rows.

    Returns
    -------
    list of tuple
        Split tuples.
    """
    values = sorted({safe_value(row, column) for row in records})
    expected_species = set()
    if only_expected_target_species:
        expected_species = {
            safe_value(row, column)
            for row in records
            if str(row.get(label_column, "")) == "expected_target"
        }

    splits = []
    for value in values:
        if value == "missing":
            continue
        if only_expected_target_species and value not in expected_species:
            continue

        train_records = [row for row in records if safe_value(row, column) != value]
        test_records = [row for row in records if safe_value(row, column) == value]
        validation_name = f"{prefix}_{sanitise_name(value)}"
        splits.append((validation_name, column, value, train_records, test_records))

    return splits


def should_skip_split(
    *,
    train_records: list[dict[str, object]],
    test_records: list[dict[str, object]],
    label_column: str,
) -> str:
    """Return a split skip reason, or an empty string if usable.

    Parameters
    ----------
    train_records : list of dict[str, object]
        Training rows.
    test_records : list of dict[str, object]
        Test rows.
    label_column : str
        Label column.

    Returns
    -------
    str
        Skip reason or empty string.
    """
    if not train_records:
        return "no_training_records"
    if not test_records:
        return "no_test_records"

    train_labels = label_counts(records=train_records, label_column=label_column)
    test_labels = label_counts(records=test_records, label_column=label_column)

    if len(train_labels) < 2:
        return "fewer_than_two_training_labels"
    if len(test_labels) < 1:
        return "no_test_labels"

    return ""


def build_validation_splits(
    *,
    records: list[dict[str, object]],
    label_column: str,
    sample_test_fraction: float,
) -> list[tuple[str, str, str, list[dict[str, object]], list[dict[str, object]]]]:
    """Build all configured validation splits.

    Parameters
    ----------
    records : list of dict[str, object]
        Input records.
    label_column : str
        Label column.
    sample_test_fraction : float
        Sample-hash split test fraction.

    Returns
    -------
    list of tuple
        Validation splits.
    """
    train_records, test_records = make_sample_hash_split(
        records=records,
        test_fraction=sample_test_fraction,
    )

    splits = [
        (
            "sample_group_hash_20pct",
            "sample_id_hash",
            "sample_group_hash_20pct",
            train_records,
            test_records,
        )
    ]

    for column, prefix in [
        ("benchmark_family", "leave_one_benchmark_family"),
        ("panel", "leave_one_panel"),
        ("species_name", "leave_one_species"),
    ]:
        splits.extend(
            make_leave_one_value_splits(
                records=records,
                column=column,
                prefix=prefix,
                label_column=label_column,
            )
        )

    splits.extend(
        make_leave_one_value_splits(
            records=records,
            column="species_name",
            prefix="leave_one_expected_target_species",
            label_column=label_column,
            only_expected_target_species=True,
        )
    )

    return splits


def summarise_metrics_by_split(
    metrics_rows: Sequence[dict[str, object]],
) -> list[dict[str, object]]:
    """Build a compact one-row-per-split summary table.

    Parameters
    ----------
    metrics_rows : sequence of dict[str, object]
        Per-label metric rows.

    Returns
    -------
    list[dict[str, object]]
        Summary rows.
    """
    by_split: dict[str, dict[str, dict[str, object]]] = {}
    for row in metrics_rows:
        split = str(row.get("validation_name", ""))
        label = str(row.get("label", ""))
        by_split.setdefault(split, {})[label] = dict(row)

    summary_rows = []
    for split, labels in sorted(by_split.items()):
        overall = labels.get("overall", {})
        expected = labels.get("expected_target", {})
        below = labels.get("observed_below_threshold", {})
        not_detected = labels.get("not_detected", {})
        off_target = labels.get("reportable_off_target_species", {})
        summary_rows.append(
            {
                "validation_name": split,
                "holdout_column": overall.get("holdout_column", ""),
                "holdout_value": overall.get("holdout_value", ""),
                "n_train": overall.get("n_train", ""),
                "n_test": overall.get("n_test", ""),
                "overall_accuracy": overall.get("f1", ""),
                "expected_target_n": expected.get("n", 0),
                "expected_target_precision": expected.get("precision", ""),
                "expected_target_recall": expected.get("recall", ""),
                "observed_below_threshold_precision": below.get("precision", ""),
                "observed_below_threshold_recall": below.get("recall", ""),
                "not_detected_precision": not_detected.get("precision", ""),
                "not_detected_recall": not_detected.get("recall", ""),
                "reportable_off_target_n": off_target.get("n", 0),
                "reportable_off_target_precision": off_target.get("precision", ""),
                "reportable_off_target_recall": off_target.get("recall", ""),
            }
        )
    return summary_rows


def run_validation(args: argparse.Namespace) -> None:
    """Run internal held-out validation.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed arguments.
    """
    training_table = existing_file(args.training_table, "training table")
    out_dir = ensure_output_dir(args.out_dir)
    logger = configure_logging(
        log_path=out_dir / "validate_kmersutra_ai_holdouts.log",
        verbose=args.verbose,
    )

    (
        evaluate_predictions,
        predict_records,
        save_model,
        train_prototype_classifier,
        read_records_table,
        write_records_table,
    ) = import_kmersutra_helpers()

    logger.info("Reading training table: %s", training_table)
    records = [dict(row) for row in read_records_table(input_path=training_table, logger=logger)]
    if not records:
        raise ValueError("Training table contained zero records")

    features = find_present_features(
        records=records,
        requested_features=args.feature_columns,
    )
    if not features:
        raise ValueError("No requested safe feature columns found in training table")

    logger.info("Records loaded: %d", len(records))
    logger.info("Features used: %s", "; ".join(features))
    logger.info(
        "Overall label counts: %s",
        label_counts_text(records=records, label_column=args.label_column),
    )

    write_records_table(
        records=[{"feature": feature} for feature in features],
        output_path=out_dir / "feature_columns_used.tsv",
        fieldnames=["feature"],
        logger=logger,
    )
    write_records_table(
        records=audit_feature_columns(features),
        output_path=out_dir / "feature_leakage_audit.tsv",
        fieldnames=["feature", "is_leakage_risk", "recommended_for_training"],
        logger=logger,
    )
    write_records_table(
        records=[
            {args.label_column: label, "n_records": count}
            for label, count in sorted(
                label_counts(records=records, label_column=args.label_column).items()
            )
        ],
        output_path=out_dir / "all_training_label_counts.tsv",
        fieldnames=[args.label_column, "n_records"],
        logger=logger,
    )

    validation_splits = build_validation_splits(
        records=records,
        label_column=args.label_column,
        sample_test_fraction=args.sample_test_fraction,
    )
    logger.info("Validation splits configured: %d", len(validation_splits))

    design_rows: list[dict[str, object]] = []
    skipped_rows: list[dict[str, object]] = []
    metrics_rows: list[dict[str, object]] = []
    prediction_rows: list[dict[str, object]] = []

    for validation_name, holdout_column, holdout_value, train_records, test_records in validation_splits:
        skip_reason = should_skip_split(
            train_records=train_records,
            test_records=test_records,
            label_column=args.label_column,
        )
        design_row = {
            "validation_name": validation_name,
            "holdout_column": holdout_column,
            "holdout_value": holdout_value,
            "n_train": len(train_records),
            "n_test": len(test_records),
            "train_label_counts": label_counts_text(
                records=train_records,
                label_column=args.label_column,
            ),
            "test_label_counts": label_counts_text(
                records=test_records,
                label_column=args.label_column,
            ),
            "skip_reason": skip_reason,
        }
        design_rows.append(design_row)

        if skip_reason:
            logger.info("Skipping %s: %s", validation_name, skip_reason)
            skipped_rows.append(design_row)
            continue

        logger.info(
            "Training split %s: n_train=%d n_test=%d",
            validation_name,
            len(train_records),
            len(test_records),
        )
        model = train_prototype_classifier(
            records=train_records,
            label_column=args.label_column,
            feature_columns=features,
            distance_quantile=args.distance_quantile,
            unknown_label=args.unknown_label,
            logger=logger,
        )
        predictions = predict_records(records=test_records, model=model, logger=logger)
        for row in predictions:
            row["validation_name"] = validation_name
            row["holdout_column"] = holdout_column
            row["holdout_value"] = holdout_value
            prediction_rows.append(row)

        split_metrics = evaluate_predictions(
            predictions=predictions,
            label_column=args.label_column,
        )
        for row in split_metrics:
            row["validation_name"] = validation_name
            row["holdout_column"] = holdout_column
            row["holdout_value"] = holdout_value
            row["n_train"] = len(train_records)
            row["n_test"] = len(test_records)
            metrics_rows.append(row)

        model_dir = out_dir / "models_by_split"
        model_dir.mkdir(parents=True, exist_ok=True)
        save_model(model=model, output_path=model_dir / f"{validation_name}.json")

    logger.info("Training final model on all records")
    final_model = train_prototype_classifier(
        records=records,
        label_column=args.label_column,
        feature_columns=features,
        distance_quantile=args.distance_quantile,
        unknown_label=args.unknown_label,
        logger=logger,
    )
    save_model(
        model=final_model,
        output_path=out_dir / "final_internal_plasmodium_calibrator_all_training.json",
    )

    final_model_summary = [
        {
            "label": label,
            "n_training_records": final_model.class_counts[label],
            "novelty_threshold": final_model.class_thresholds[label],
            "distance_quantile": final_model.distance_quantile,
        }
        for label in sorted(final_model.class_counts)
    ]
    write_records_table(
        records=final_model_summary,
        output_path=out_dir / "final_internal_model_training_summary.tsv",
        fieldnames=[
            "label",
            "n_training_records",
            "novelty_threshold",
            "distance_quantile",
        ],
        logger=logger,
    )

    metric_fields = [
        "validation_name",
        "holdout_column",
        "holdout_value",
        "n_train",
        "n_test",
        "label",
        "n",
        "tp",
        "fp",
        "fn",
        "precision",
        "recall",
        "f1",
    ]
    write_records_table(
        records=metrics_rows,
        output_path=out_dir / "holdout_metrics.tsv",
        fieldnames=metric_fields,
        logger=logger,
    )

    summary_fields = [
        "validation_name",
        "holdout_column",
        "holdout_value",
        "n_train",
        "n_test",
        "overall_accuracy",
        "expected_target_n",
        "expected_target_precision",
        "expected_target_recall",
        "observed_below_threshold_precision",
        "observed_below_threshold_recall",
        "not_detected_precision",
        "not_detected_recall",
        "reportable_off_target_n",
        "reportable_off_target_precision",
        "reportable_off_target_recall",
    ]
    write_records_table(
        records=summarise_metrics_by_split(metrics_rows),
        output_path=out_dir / "holdout_summary_by_split.tsv",
        fieldnames=summary_fields,
        logger=logger,
    )

    design_fields = [
        "validation_name",
        "holdout_column",
        "holdout_value",
        "n_train",
        "n_test",
        "train_label_counts",
        "test_label_counts",
        "skip_reason",
    ]
    write_records_table(
        records=design_rows,
        output_path=out_dir / "validation_design.tsv",
        fieldnames=design_fields,
        logger=logger,
    )
    write_records_table(
        records=skipped_rows,
        output_path=out_dir / "skipped_splits.tsv",
        fieldnames=design_fields,
        logger=logger,
    )

    if prediction_rows:
        write_records_table(
            records=prediction_rows,
            output_path=out_dir / "holdout_predictions.tsv.gz",
            fieldnames=list(prediction_rows[0].keys()),
            logger=logger,
        )

    manifest = {
        "training_table": str(training_table),
        "n_records": len(records),
        "feature_columns": features,
        "distance_quantile": args.distance_quantile,
        "sample_test_fraction": args.sample_test_fraction,
        "n_validation_splits": len(validation_splits),
        "n_completed_splits": len({row["validation_name"] for row in metrics_rows}),
        "n_skipped_splits": len(skipped_rows),
    }
    write_json(manifest, out_dir / "validation_manifest.json")
    logger.info("Completed validation: %s", out_dir)


def main() -> int:
    """Command-line entry point.

    Returns
    -------
    int
        Exit status.
    """
    args = parse_args()
    run_validation(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
