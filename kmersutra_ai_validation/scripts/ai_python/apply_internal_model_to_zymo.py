#!/usr/bin/env python3
"""Apply an internal KmerSutra AI calibrator to public Zymo outputs.

This script is an external cross-domain validation/stress test. It does not
retrain the model. It applies a model trained on internal Plasmodium
sample/species-level evidence to public Zymo species-detection rows and compares
predictions with an expected-organism truth mapping.

The output should be interpreted as evidence-calibrator generalisation, not as a
read-level classifier benchmark.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from kmersutra_ai_common import (
    add_common_parser_arguments,
    configure_logging,
    count_records_by_column,
    ensure_output_dir,
    existing_file,
    infer_public_truth_label,
    read_reference_label_map,
    write_json,
)


def import_kmersutra_helpers():
    """Import KmerSutra helper functions lazily.

    Returns
    -------
    tuple
        KmerSutra helper functions.
    """
    from kmersutra.ai_calibration import build_call_feature_record, evaluate_predictions
    from kmersutra.ml import load_model, predict_records
    from kmersutra.table_io import read_records_table, write_records_table

    return (
        build_call_feature_record,
        evaluate_predictions,
        load_model,
        predict_records,
        read_records_table,
        write_records_table,
    )


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Apply an internal KmerSutra AI model to public Zymo outputs."
    )
    parser.add_argument(
        "--model_json",
        required=True,
        help="Final internal model JSON from the Plasmodium validation run.",
    )
    parser.add_argument(
        "--calls_table",
        required=True,
        help="Public Zymo species_detection_calls.tsv or TSV.GZ table.",
    )
    parser.add_argument(
        "--reference_label_map",
        default="",
        help="Optional reference_label_map.tsv for expected target labels.",
    )
    parser.add_argument(
        "--out_dir",
        required=True,
        help="Output directory.",
    )
    add_common_parser_arguments(parser)
    return parser.parse_args()


def build_external_feature_records(
    *,
    calls: list[dict[str, object]],
    expected_reference_labels: set[str],
    expected_species_names: set[str],
) -> list[dict[str, object]]:
    """Build feature records and external truth labels.

    Parameters
    ----------
    calls : list of dict[str, object]
        Public Zymo species-detection call rows.
    expected_reference_labels : set[str]
        Expected official reference labels.
    expected_species_names : set[str]
        Expected species names.

    Returns
    -------
    list[dict[str, object]]
        Feature records compatible with KmerSutra model prediction.
    """
    build_call_feature_record, *_ = import_kmersutra_helpers()
    feature_records: list[dict[str, object]] = []

    for row in calls:
        feature_record = build_call_feature_record(record=row)
        feature_record["ml_report_label"] = infer_public_truth_label(
            row=row,
            expected_reference_labels=expected_reference_labels,
            expected_species_names=expected_species_names,
        )
        feature_record["external_dataset"] = "ERR5396170_ZymoBIOMICS_D6300_ONT_Q20"
        feature_record["public_call"] = row.get("call", row.get("call_status", ""))
        feature_record["public_species_name"] = row.get("species_name", "")
        feature_record["public_report_label"] = row.get("report_label", "")
        feature_record["public_reference_label"] = row.get("reference_label", "")
        feature_records.append(feature_record)

    return feature_records


def validate_feature_compatibility(
    *,
    feature_records: list[dict[str, object]],
    model,
) -> None:
    """Check that external records contain all model features.

    Parameters
    ----------
    feature_records : list of dict[str, object]
        External feature records.
    model : object
        KmerSutra prototype model.

    Raises
    ------
    ValueError
        If required model features are missing.
    """
    if not feature_records:
        raise ValueError("No external feature records were produced")

    missing = [
        column for column in model.feature_columns
        if column not in feature_records[0]
    ]
    if missing:
        raise ValueError(
            "External feature table is missing model feature columns: "
            + ", ".join(missing)
        )


def run_external_validation(args: argparse.Namespace) -> None:
    """Run external Zymo validation.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments.
    """
    model_json = existing_file(args.model_json, "model JSON")
    calls_table = existing_file(args.calls_table, "calls table")
    reference_map = Path(args.reference_label_map) if args.reference_label_map else None
    out_dir = ensure_output_dir(args.out_dir)

    logger = configure_logging(
        log_path=out_dir / "external_zymo_ai_validation.log",
        verbose=args.verbose,
    )

    (
        _,
        evaluate_predictions,
        load_model,
        predict_records,
        read_records_table,
        write_records_table,
    ) = import_kmersutra_helpers()

    logger.info("Loading model: %s", model_json)
    model = load_model(model_path=model_json)
    logger.info("Reading public calls table: %s", calls_table)
    calls = [dict(row) for row in read_records_table(input_path=calls_table, logger=logger)]

    expected_reference_labels, expected_species_names, reference_label_roles = (
        read_reference_label_map(
            reference_label_map=reference_map,
            logger=logger,
        )
    )

    feature_records = build_external_feature_records(
        calls=calls,
        expected_reference_labels=expected_reference_labels,
        expected_species_names=expected_species_names,
    )
    validate_feature_compatibility(feature_records=feature_records, model=model)

    logger.info("External feature records: %d", len(feature_records))
    logger.info("Applying model without retraining")
    predictions = predict_records(records=feature_records, model=model, logger=logger)
    metrics = evaluate_predictions(predictions=predictions, label_column=model.label_column)

    write_records_table(
        records=feature_records,
        output_path=out_dir / "external_zymo_feature_table.tsv.gz",
        fieldnames=list(feature_records[0].keys()),
        logger=logger,
    )
    write_records_table(
        records=predictions,
        output_path=out_dir / "external_zymo_predictions.tsv.gz",
        fieldnames=list(predictions[0].keys()) if predictions else [],
        logger=logger,
    )
    write_records_table(
        records=metrics,
        output_path=out_dir / "external_zymo_metrics.tsv",
        fieldnames=["label", "n", "tp", "fp", "fn", "precision", "recall", "f1"],
        logger=logger,
    )
    write_records_table(
        records=count_records_by_column(records=feature_records, column="ml_report_label"),
        output_path=out_dir / "external_zymo_truth_label_counts.tsv",
        fieldnames=["ml_report_label", "n_records"],
        logger=logger,
    )
    write_records_table(
        records=count_records_by_column(records=predictions, column="prediction"),
        output_path=out_dir / "external_zymo_prediction_counts.tsv",
        fieldnames=["prediction", "n_records"],
        logger=logger,
    )

    expected_rows = [
        row for row in predictions
        if str(row.get("ml_report_label", "")) == "expected_target"
    ]
    expected_fields = list(predictions[0].keys()) if predictions else []
    write_records_table(
        records=expected_rows,
        output_path=out_dir / "external_zymo_expected_target_predictions.tsv",
        fieldnames=expected_fields,
        logger=logger,
    )

    manifest = {
        "model_json": str(model_json),
        "calls_table": str(calls_table),
        "reference_label_map": str(reference_map or ""),
        "n_call_rows": len(calls),
        "n_feature_rows": len(feature_records),
        "model_feature_columns": model.feature_columns,
        "expected_reference_labels": sorted(expected_reference_labels),
        "expected_species_names": sorted(expected_species_names),
        "reference_label_roles_observed": reference_label_roles,
    }
    write_json(manifest, out_dir / "external_zymo_validation_manifest.json")
    logger.info("Completed external Zymo validation: %s", out_dir)


def main() -> int:
    """Command-line entry point.

    Returns
    -------
    int
        Exit status.
    """
    args = parse_args()
    run_external_validation(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
