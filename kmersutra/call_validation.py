"""Holdout validation helpers for KmerSutra call calibrators."""

from __future__ import annotations

import hashlib
import json
import logging
from collections import Counter
from pathlib import Path

from kmersutra.ai_calibration import (
    DEFAULT_LABEL_COLUMN,
    DEFAULT_UNKNOWN_LABEL,
    evaluate_predictions,
    infer_numeric_feature_columns,
)
from kmersutra.ml import predict_records, save_model, train_prototype_classifier
from kmersutra.table_io import read_records_table, write_records_table


def stable_text_hash(*, value: str) -> int:
    """Return a deterministic integer hash for text.

    Parameters
    ----------
    value : str
        Text value.

    Returns
    -------
    int
        Stable hash value.
    """
    digest = hashlib.sha256(value.encode("utf-8")).hexdigest()
    return int(digest[:12], 16)


def safe_value(*, record: dict[str, object], column: str) -> str:
    """Return a non-empty value for a split column.

    Parameters
    ----------
    record : dict[str, object]
        Input record.
    column : str
        Column name.

    Returns
    -------
    str
        Column value or ``missing``.
    """
    value = str(record.get(column, "")).strip()
    return value if value else "missing"


def label_count_summary(
    *,
    records: list[dict[str, object]],
    label_column: str = DEFAULT_LABEL_COLUMN,
) -> str:
    """Return a compact label-count string.

    Parameters
    ----------
    records : list of dict
        Input records.
    label_column : str, optional
        Label column.

    Returns
    -------
    str
        Semicolon-delimited label counts.
    """
    counts = Counter(str(row.get(label_column, "")) for row in records)
    return ";".join(f"{label}:{counts[label]}" for label in sorted(counts))


def label_count_records(
    *,
    records: list[dict[str, object]],
    label_column: str = DEFAULT_LABEL_COLUMN,
) -> list[dict[str, object]]:
    """Return label counts as table records.

    Parameters
    ----------
    records : list of dict
        Input records.
    label_column : str, optional
        Label column.

    Returns
    -------
    list[dict[str, object]]
        Count records.
    """
    counts = Counter(str(row.get(label_column, "")) for row in records)
    return [
        {label_column: label, "n_records": counts[label]}
        for label in sorted(counts)
    ]


def sample_hash_split(
    *,
    records: list[dict[str, object]],
    sample_column: str = "sample_id",
    test_fraction: float = 0.2,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    """Create a deterministic sample-hash validation split.

    Parameters
    ----------
    records : list of dict
        Input records.
    sample_column : str, optional
        Sample identifier column.
    test_fraction : float, optional
        Fraction of samples assigned to test.

    Returns
    -------
    tuple[list[dict], list[dict]]
        Train and test records.
    """
    if not 0.0 <= test_fraction < 1.0:
        raise ValueError("test_fraction must be >= 0 and < 1")
    train: list[dict[str, object]] = []
    test: list[dict[str, object]] = []
    threshold = int(test_fraction * 10000)
    for record in records:
        sample_id = safe_value(record=record, column=sample_column)
        bucket = stable_text_hash(value=sample_id) % 10000
        if bucket < threshold:
            test.append(dict(record))
        else:
            train.append(dict(record))
    return train, test


def leave_one_value_splits(
    *,
    records: list[dict[str, object]],
    column: str,
    prefix: str,
    label_column: str = DEFAULT_LABEL_COLUMN,
    only_expected_target_species: bool = False,
) -> list[tuple[str, str, str, list[dict[str, object]], list[dict[str, object]]]]:
    """Create leave-one-value-out validation splits.

    Parameters
    ----------
    records : list of dict
        Input records.
    column : str
        Holdout column.
    prefix : str
        Validation-name prefix.
    label_column : str, optional
        Label column.
    only_expected_target_species : bool, optional
        Restrict to species with expected-target records.

    Returns
    -------
    list of tuple
        Validation split tuples.
    """
    values = sorted({safe_value(record=row, column=column) for row in records})
    expected_values = set()
    if only_expected_target_species:
        expected_values = {
            safe_value(record=row, column=column)
            for row in records
            if str(row.get(label_column, "")) == "expected_target"
        }
    splits = []
    for value in values:
        if value == "missing":
            continue
        if only_expected_target_species and value not in expected_values:
            continue
        train = [row for row in records if safe_value(record=row, column=column) != value]
        test = [row for row in records if safe_value(record=row, column=column) == value]
        name = f"{prefix}_{value}".replace(" ", "_").replace("/", "_")
        splits.append((name, column, value, train, test))
    return splits


def split_skip_reason(
    *,
    train_records: list[dict[str, object]],
    test_records: list[dict[str, object]],
    label_column: str = DEFAULT_LABEL_COLUMN,
) -> str:
    """Return a skip reason for an invalid split.

    Parameters
    ----------
    train_records : list of dict
        Training records.
    test_records : list of dict
        Test records.
    label_column : str, optional
        Label column.

    Returns
    -------
    str
        Empty string if the split is valid, otherwise a reason.
    """
    if not train_records:
        return "no_training_records"
    if not test_records:
        return "no_test_records"
    train_labels = {
        str(row.get(label_column, ""))
        for row in train_records
        if str(row.get(label_column, ""))
    }
    test_labels = {
        str(row.get(label_column, ""))
        for row in test_records
        if str(row.get(label_column, ""))
    }
    if len(train_labels) < 2:
        return "fewer_than_two_training_labels"
    if not test_labels:
        return "no_test_labels"
    return ""


def build_validation_splits(
    *,
    records: list[dict[str, object]],
    sample_test_fraction: float = 0.2,
    label_column: str = DEFAULT_LABEL_COLUMN,
) -> list[tuple[str, str, str, list[dict[str, object]], list[dict[str, object]]]]:
    """Build the standard KmerSutra call-calibrator validation splits.

    Parameters
    ----------
    records : list of dict
        Training records.
    sample_test_fraction : float, optional
        Sample-hash test fraction.
    label_column : str, optional
        Label column.

    Returns
    -------
    list of tuple
        Validation split definitions.
    """
    train, test = sample_hash_split(
        records=records,
        test_fraction=sample_test_fraction,
    )
    splits = [
        (
            "sample_group_hash_20pct",
            "sample_id_hash",
            "sample_group_hash_20pct",
            train,
            test,
        )
    ]
    for column, prefix in [
        ("benchmark_family", "leave_one_benchmark_family"),
        ("panel", "leave_one_panel"),
        ("species_name", "leave_one_species"),
    ]:
        splits.extend(
            leave_one_value_splits(
                records=records,
                column=column,
                prefix=prefix,
                label_column=label_column,
            )
        )
    splits.extend(
        leave_one_value_splits(
            records=records,
            column="species_name",
            prefix="leave_one_expected_target_species",
            label_column=label_column,
            only_expected_target_species=True,
        )
    )
    return splits


def summarise_metrics_by_split(
    *,
    metrics: list[dict[str, object]],
) -> list[dict[str, object]]:
    """Extract one overall metric row per validation split.

    Parameters
    ----------
    metrics : list of dict
        Holdout metrics.

    Returns
    -------
    list[dict[str, object]]
        Overall split summaries.
    """
    summary = []
    for row in metrics:
        if str(row.get("label", "")) != "overall":
            continue
        summary.append(
            {
                "validation_name": row.get("validation_name", ""),
                "holdout_column": row.get("holdout_column", ""),
                "holdout_value": row.get("holdout_value", ""),
                "n_train": row.get("n_train", ""),
                "n_test": row.get("n_test", ""),
                "overall_accuracy": row.get("f1", ""),
            }
        )
    return summary


def write_json(*, data: object, output_path: str | Path) -> None:
    """Write JSON with stable formatting.

    Parameters
    ----------
    data : object
        JSON-serialisable data.
    output_path : str or pathlib.Path
        Output path.
    """
    Path(output_path).write_text(
        json.dumps(data, indent=2, sort_keys=True),
        encoding="utf-8",
    )


def validate_call_calibrator_holdouts(
    *,
    training_table: str | Path,
    out_dir: str | Path,
    feature_profile: str = "safe_transformed",
    feature_columns: list[str] | None = None,
    label_column: str = DEFAULT_LABEL_COLUMN,
    distance_quantile: float = 0.95,
    sample_test_fraction: float = 0.2,
    unknown_label: str = DEFAULT_UNKNOWN_LABEL,
    logger: logging.Logger | None = None,
) -> dict[str, object]:
    """Run full internal holdout validation and train a final model.

    Parameters
    ----------
    training_table : str or pathlib.Path
        AI call-training table.
    out_dir : str or pathlib.Path
        Output directory.
    feature_profile : str, optional
        Feature profile used when feature columns are not explicitly supplied.
    feature_columns : list[str] or None, optional
        Explicit feature columns.
    label_column : str, optional
        Label column.
    distance_quantile : float, optional
        Open-set threshold quantile.
    sample_test_fraction : float, optional
        Sample-hash validation test fraction.
    unknown_label : str, optional
        Unknown/open-set label.
    logger : logging.Logger or None, optional
        Logger.

    Returns
    -------
    dict[str, object]
        Validation manifest.
    """
    output_dir = Path(out_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    records = read_records_table(input_path=training_table, logger=logger)
    if not records:
        raise ValueError("Cannot validate from zero training records")
    features = feature_columns or infer_numeric_feature_columns(
        records=records,
        feature_profile=feature_profile,
    )
    if not features:
        raise ValueError("No feature columns available for validation")

    if logger:
        logger.info("Loaded %d training records", len(records))
        logger.info("Feature profile: %s", feature_profile)
        logger.info("Feature columns: %s", "; ".join(features))

    splits = build_validation_splits(
        records=records,
        sample_test_fraction=sample_test_fraction,
        label_column=label_column,
    )

    design_rows: list[dict[str, object]] = []
    skipped_rows: list[dict[str, object]] = []
    metrics_rows: list[dict[str, object]] = []
    prediction_rows: list[dict[str, object]] = []
    model_dir = output_dir / "models_by_split"
    model_dir.mkdir(parents=True, exist_ok=True)

    for name, column, value, train_records, test_records in splits:
        reason = split_skip_reason(
            train_records=train_records,
            test_records=test_records,
            label_column=label_column,
        )
        design_row = {
            "validation_name": name,
            "holdout_column": column,
            "holdout_value": value,
            "n_train": len(train_records),
            "n_test": len(test_records),
            "train_label_counts": label_count_summary(
                records=train_records,
                label_column=label_column,
            ),
            "test_label_counts": label_count_summary(
                records=test_records,
                label_column=label_column,
            ),
            "skip_reason": reason,
        }
        design_rows.append(design_row)
        if reason:
            skipped_rows.append(design_row)
            if logger:
                logger.info("Skipping %s: %s", name, reason)
            continue

        if logger:
            logger.info(
                "Training split %s with %d train and %d test rows",
                name,
                len(train_records),
                len(test_records),
            )
        model = train_prototype_classifier(
            records=train_records,
            label_column=label_column,
            feature_columns=features,
            distance_quantile=distance_quantile,
            unknown_label=unknown_label,
            logger=logger,
        )
        predictions = predict_records(
            records=test_records,
            model=model,
            logger=logger,
        )
        for row in predictions:
            row["validation_name"] = name
            row["holdout_column"] = column
            row["holdout_value"] = value
            prediction_rows.append(row)
        split_metrics = evaluate_predictions(
            predictions=predictions,
            label_column=label_column,
        )
        for row in split_metrics:
            row["validation_name"] = name
            row["holdout_column"] = column
            row["holdout_value"] = value
            row["n_train"] = len(train_records)
            row["n_test"] = len(test_records)
            metrics_rows.append(row)
        save_model(model=model, output_path=model_dir / f"{name}.json")

    final_model = train_prototype_classifier(
        records=records,
        label_column=label_column,
        feature_columns=features,
        distance_quantile=distance_quantile,
        unknown_label=unknown_label,
        logger=logger,
    )
    save_model(
        model=final_model,
        output_path=output_dir / "final_internal_calibrator_all_training.json",
    )

    final_summary = [
        {
            "label": label,
            "n_training_records": final_model.class_counts[label],
            "novelty_threshold": final_model.class_thresholds[label],
            "distance_quantile": final_model.distance_quantile,
            "feature_profile": feature_profile,
        }
        for label in sorted(final_model.class_counts)
    ]

    write_records_table(
        records=[{"feature": feature} for feature in features],
        output_path=output_dir / "feature_columns_used.tsv",
        fieldnames=["feature"],
        logger=logger,
    )
    write_records_table(
        records=label_count_records(records=records, label_column=label_column),
        output_path=output_dir / "all_training_label_counts.tsv",
        fieldnames=[label_column, "n_records"],
        logger=logger,
    )
    write_records_table(
        records=design_rows,
        output_path=output_dir / "validation_design.tsv",
        fieldnames=[
            "validation_name",
            "holdout_column",
            "holdout_value",
            "n_train",
            "n_test",
            "train_label_counts",
            "test_label_counts",
            "skip_reason",
        ],
        logger=logger,
    )
    write_records_table(
        records=skipped_rows,
        output_path=output_dir / "skipped_splits.tsv",
        fieldnames=[
            "validation_name",
            "holdout_column",
            "holdout_value",
            "n_train",
            "n_test",
            "train_label_counts",
            "test_label_counts",
            "skip_reason",
        ],
        logger=logger,
    )
    write_records_table(
        records=metrics_rows,
        output_path=output_dir / "holdout_metrics.tsv",
        fieldnames=[
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
        ],
        logger=logger,
    )
    if prediction_rows:
        write_records_table(
            records=prediction_rows,
            output_path=output_dir / "holdout_predictions.tsv.gz",
            fieldnames=list(prediction_rows[0].keys()),
            logger=logger,
        )
    write_records_table(
        records=summarise_metrics_by_split(metrics=metrics_rows),
        output_path=output_dir / "holdout_summary_by_split.tsv",
        fieldnames=[
            "validation_name",
            "holdout_column",
            "holdout_value",
            "n_train",
            "n_test",
            "overall_accuracy",
        ],
        logger=logger,
    )
    write_records_table(
        records=final_summary,
        output_path=output_dir / "final_internal_model_training_summary.tsv",
        fieldnames=[
            "label",
            "n_training_records",
            "novelty_threshold",
            "distance_quantile",
            "feature_profile",
        ],
        logger=logger,
    )

    manifest = {
        "training_table": str(training_table),
        "n_records": len(records),
        "feature_profile": feature_profile,
        "feature_columns": features,
        "distance_quantile": distance_quantile,
        "sample_test_fraction": sample_test_fraction,
        "n_validation_splits": len(splits),
        "n_completed_splits": len({row["validation_name"] for row in metrics_rows}),
        "n_skipped_splits": len(skipped_rows),
    }
    write_json(data=manifest, output_path=output_dir / "validation_manifest.json")
    if logger:
        logger.info("Completed call-calibrator validation")
    return manifest
