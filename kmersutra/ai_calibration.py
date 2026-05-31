"""AI-ready call calibration helpers for KmerSutra benchmark outputs.

The initial KmerSutra AI layer is deliberately interpretable. It converts
summary/evidence rows into a training table and reuses the lightweight
open-set prototype classifier from :mod:`kmersutra.ml`. The aim is to calibrate
report-layer evidence rather than replace the auditable rule-based calls.
"""

from __future__ import annotations

import hashlib
import logging
from collections import Counter, defaultdict
from collections.abc import Iterable
from pathlib import Path

from kmersutra.table_io import read_records_table, write_records_table
from kmersutra.ml import (
    predict_records,
    save_model,
    train_prototype_classifier,
)

CALL_FEATURE_COLUMNS = [
    "n_hits",
    "n_unique_kmers",
    "n_positive_sequences",
    "n_k_values_positive",
    "best_k",
    "n_exact_hits",
    "n_fuzzy_hits",
    "conflicting_unique_kmers",
    "conflict_ratio",
    "reportable_conflicting_unique_kmers",
    "reportable_conflict_ratio",
    "mixed_species_support_fraction",
    "confidence_score",
    "signal_confidence_score",
    "spike_n",
    "spike_n_per_genome",
    "total_spike_n",
]

DEFAULT_LABEL_COLUMN = "ml_report_label"
DEFAULT_UNKNOWN_LABEL = "unknown_or_unresolved"


def normalise_bool(value: object) -> bool:
    """Return a robust boolean interpretation of common TSV values.

    Parameters
    ----------
    value : object
        Input value from a parsed TSV row.

    Returns
    -------
    bool
        Boolean interpretation.
    """
    return str(value).strip().lower() in {"1", "true", "t", "yes", "y"}


def normalise_float(value: object) -> float:
    """Return a float value with defensive missing-value handling.

    Parameters
    ----------
    value : object
        Input value.

    Returns
    -------
    float
        Parsed float, or zero for missing/non-numeric values.
    """
    if value in {None, ""}:
        return 0.0
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def infer_report_label(*, record: dict[str, object]) -> str:
    """Infer an AI training label from a KmerSutra call/evidence row.

    Parameters
    ----------
    record : dict[str, object]
        Combined detection-call row from a comparable summary.

    Returns
    -------
    str
        Training label for call-calibration models.
    """
    report_layer = str(
        record.get("report_layer")
        or record.get("interpretation_layer")
        or record.get("call_layer")
        or ""
    ).strip().lower()
    if report_layer in {
        "background_candidate_signal",
        "background_candidate",
        "empirical_background_signal",
    }:
        return "background_candidate_signal"
    if report_layer in {"neighbour_lineage_evidence", "neighbor_lineage_evidence"}:
        return "neighbour_lineage_evidence"

    if normalise_bool(record.get("is_positive_expected", False)):
        return "expected_target"
    if normalise_bool(record.get("is_background_candidate", False)):
        return "background_candidate_signal"
    if normalise_bool(record.get("is_neighbour_lineage", False)):
        return "neighbour_lineage_evidence"
    if normalise_bool(record.get("is_positive_off_target", False)):
        return "reportable_off_target_species"
    if normalise_bool(record.get("is_positive_call", False)):
        return "positive_unclassified"

    call = str(record.get("call", "")).strip().lower()
    if call in {"present", "species_detected", "detected", "positive"}:
        return "positive_unclassified"
    if call in {"observed_below_threshold", "low_evidence"}:
        return "observed_below_threshold"
    return "not_detected"


def build_call_feature_record(
    *,
    record: dict[str, object],
    label_column: str = DEFAULT_LABEL_COLUMN,
) -> dict[str, object]:
    """Convert one summary row into an AI-ready feature record.

    Parameters
    ----------
    record : dict[str, object]
        Input detection/evidence row.
    label_column : str, optional
        Output label column name.

    Returns
    -------
    dict[str, object]
        Feature record with numeric features and metadata.
    """
    feature_record: dict[str, object] = {
        "sample_id": record.get("sample_id", ""),
        "species_name": record.get("species_name", ""),
        "clade": record.get("clade", ""),
        "benchmark_family": record.get("benchmark_family", ""),
        "panel": record.get("panel", ""),
        "replicate": record.get("replicate", ""),
        "is_negative": str(record.get("is_negative", "")),
        "is_shuffled_control": str(record.get("is_shuffled_control", "")),
        label_column: infer_report_label(record=record),
    }
    for column in CALL_FEATURE_COLUMNS:
        feature_record[column] = normalise_float(record.get(column, 0.0))
    feature_record["has_long_k_support"] = 1.0 if feature_record["best_k"] >= 101 else 0.0
    feature_record["has_multi_k_support"] = (
        1.0 if feature_record["n_k_values_positive"] >= 2 else 0.0
    )
    feature_record["exact_hit_fraction"] = (
        feature_record["n_exact_hits"] / feature_record["n_hits"]
        if feature_record["n_hits"] > 0
        else 0.0
    )
    return feature_record


def build_call_training_table(
    *,
    records: Iterable[dict[str, object]],
    label_column: str = DEFAULT_LABEL_COLUMN,
    include_not_detected: bool = True,
    max_not_detected_per_label: int | None = None,
    logger: logging.Logger | None = None,
) -> list[dict[str, object]]:
    """Build an AI-ready training table from detection-call rows.

    Parameters
    ----------
    records : iterable of dict[str, object]
        Detection-call records.
    label_column : str, optional
        Output label column name.
    include_not_detected : bool, optional
        Whether to retain not-detected rows.
    max_not_detected_per_label : int or None, optional
        Optional cap on not-detected rows to reduce class imbalance.
    logger : logging.Logger or None, optional
        Logger.

    Returns
    -------
    list[dict[str, object]]
        Feature records.
    """
    output: list[dict[str, object]] = []
    not_detected_count = 0
    for record in records:
        feature_record = build_call_feature_record(
            record=record,
            label_column=label_column,
        )
        label = str(feature_record[label_column])
        if label == "not_detected":
            if not include_not_detected:
                continue
            if (
                max_not_detected_per_label is not None
                and not_detected_count >= max_not_detected_per_label
            ):
                continue
            not_detected_count += 1
        output.append(feature_record)
    if logger:
        counts = Counter(str(row[label_column]) for row in output)
        for label, count in sorted(counts.items()):
            logger.info("AI training label %s: %d", label, count)
    return output


def stable_group_hash(*, record: dict[str, object], group_columns: list[str]) -> int:
    """Return a deterministic hash for grouped train/test splitting.

    Parameters
    ----------
    record : dict[str, object]
        Input record.
    group_columns : list[str]
        Columns defining the split group.

    Returns
    -------
    int
        Stable integer hash.
    """
    payload = "\t".join(str(record.get(column, "")) for column in group_columns)
    digest = hashlib.sha256(payload.encode("utf-8")).hexdigest()
    return int(digest[:12], 16)


def split_records_by_group(
    *,
    records: Iterable[dict[str, object]],
    group_columns: list[str] | None = None,
    test_fraction: float = 0.2,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    """Split records deterministically by sample/run group.

    Parameters
    ----------
    records : iterable of dict[str, object]
        Feature records.
    group_columns : list[str] or None, optional
        Grouping columns. Defaults to ``sample_id``.
    test_fraction : float, optional
        Approximate fraction assigned to the test split.

    Returns
    -------
    tuple[list[dict[str, object]], list[dict[str, object]]]
        Train and test records.
    """
    if not 0.0 <= test_fraction < 1.0:
        raise ValueError("test_fraction must be >= 0 and < 1")
    groups = group_columns or ["sample_id"]
    train: list[dict[str, object]] = []
    test: list[dict[str, object]] = []
    for record in records:
        bucket = stable_group_hash(record=record, group_columns=groups) % 10000
        if bucket < int(test_fraction * 10000):
            test.append(dict(record))
        else:
            train.append(dict(record))
    if not train and test:
        train.append(test.pop())
    return train, test


def evaluate_predictions(
    *,
    predictions: Iterable[dict[str, object]],
    label_column: str = DEFAULT_LABEL_COLUMN,
    prediction_column: str = "prediction",
) -> list[dict[str, object]]:
    """Evaluate prediction performance by class.

    Parameters
    ----------
    predictions : iterable of dict[str, object]
        Prediction rows containing true and predicted labels.
    label_column : str, optional
        True label column.
    prediction_column : str, optional
        Prediction column.

    Returns
    -------
    list[dict[str, object]]
        Per-label and overall metrics.
    """
    rows = list(predictions)
    labels = sorted(
        {
            str(row.get(label_column, ""))
            for row in rows
            if str(row.get(label_column, ""))
        }
        | {
            str(row.get(prediction_column, ""))
            for row in rows
            if str(row.get(prediction_column, ""))
        }
    )
    metrics: list[dict[str, object]] = []
    total_correct = sum(
        1 for row in rows
        if str(row.get(label_column, "")) == str(row.get(prediction_column, ""))
    )
    for label in labels:
        tp = sum(
            1 for row in rows
            if str(row.get(label_column, "")) == label
            and str(row.get(prediction_column, "")) == label
        )
        fp = sum(
            1 for row in rows
            if str(row.get(label_column, "")) != label
            and str(row.get(prediction_column, "")) == label
        )
        fn = sum(
            1 for row in rows
            if str(row.get(label_column, "")) == label
            and str(row.get(prediction_column, "")) != label
        )
        precision = tp / (tp + fp) if tp + fp else 0.0
        recall = tp / (tp + fn) if tp + fn else 0.0
        f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
        metrics.append(
            {
                "label": label,
                "n": tp + fn,
                "tp": tp,
                "fp": fp,
                "fn": fn,
                "precision": round(precision, 6),
                "recall": round(recall, 6),
                "f1": round(f1, 6),
            }
        )
    metrics.append(
        {
            "label": "overall",
            "n": len(rows),
            "tp": total_correct,
            "fp": "",
            "fn": "",
            "precision": "",
            "recall": "",
            "f1": round(total_correct / len(rows), 6) if rows else 0.0,
        }
    )
    return metrics


def write_call_training_table_from_table(
    *,
    calls_table: str | Path,
    output_table: str | Path,
    label_column: str = DEFAULT_LABEL_COLUMN,
    include_not_detected: bool = True,
    max_not_detected: int | None = None,
    logger: logging.Logger | None = None,
) -> list[dict[str, object]]:
    """Read calls and write an AI-ready training table.

    Parameters
    ----------
    calls_table : str or pathlib.Path
        Combined KmerSutra detection-call table. Supported suffixes are
        ``.tsv``, ``.tsv.gz`` and ``.parquet``.
    output_table : str or pathlib.Path
        Output training table. Supported suffixes are ``.tsv``, ``.tsv.gz``
        and ``.parquet``.
    label_column : str, optional
        Label column name.
    include_not_detected : bool, optional
        Whether to keep not-detected rows.
    max_not_detected : int or None, optional
        Optional cap on not-detected rows.
    logger : logging.Logger or None, optional
        Logger.

    Returns
    -------
    list[dict[str, object]]
        Written feature records.
    """
    records = read_records_table(input_path=calls_table, logger=logger)
    feature_records = build_call_training_table(
        records=records,
        label_column=label_column,
        include_not_detected=include_not_detected,
        max_not_detected_per_label=max_not_detected,
        logger=logger,
    )
    if not feature_records:
        raise ValueError("No AI training records were generated")
    fieldnames = list(feature_records[0].keys())
    write_records_table(
        records=feature_records,
        output_path=output_table,
        fieldnames=fieldnames,
        logger=logger,
    )
    return feature_records


def write_call_training_table_from_tsv(
    *,
    calls_tsv: str | Path,
    output_tsv: str | Path,
    label_column: str = DEFAULT_LABEL_COLUMN,
    include_not_detected: bool = True,
    max_not_detected: int | None = None,
    logger: logging.Logger | None = None,
) -> list[dict[str, object]]:
    """Backward-compatible wrapper for TSV-style call-training output.

    Parameters
    ----------
    calls_tsv : str or pathlib.Path
        Input table path. Despite the historical argument name, ``.tsv``,
        ``.tsv.gz`` and ``.parquet`` are supported.
    output_tsv : str or pathlib.Path
        Output table path. Despite the historical argument name, ``.tsv``,
        ``.tsv.gz`` and ``.parquet`` are supported.
    label_column : str, optional
        Label column name.
    include_not_detected : bool, optional
        Whether to keep not-detected rows.
    max_not_detected : int or None, optional
        Optional cap on not-detected rows.
    logger : logging.Logger or None, optional
        Logger.

    Returns
    -------
    list[dict[str, object]]
        Written feature records.
    """
    return write_call_training_table_from_table(
        calls_table=calls_tsv,
        output_table=output_tsv,
        label_column=label_column,
        include_not_detected=include_not_detected,
        max_not_detected=max_not_detected,
        logger=logger,
    )


def train_evaluate_call_calibrator(
    *,
    training_tsv: str | Path | None = None,
    model_json: str | Path,
    summary_tsv: str | Path | None = None,
    evaluation_tsv: str | Path | None = None,
    training_table: str | Path | None = None,
    summary_table: str | Path | None = None,
    evaluation_table: str | Path | None = None,
    label_column: str = DEFAULT_LABEL_COLUMN,
    feature_columns: list[str] | None = None,
    test_fraction: float = 0.2,
    group_columns: list[str] | None = None,
    distance_quantile: float = 0.95,
    unknown_label: str = DEFAULT_UNKNOWN_LABEL,
    logger: logging.Logger | None = None,
):
    """Train and evaluate an interpretable call-calibration model.

    Parameters
    ----------
    training_tsv : str or pathlib.Path or None, optional
        Backward-compatible training-table path. Despite the historical name,
        ``.tsv``, ``.tsv.gz`` and ``.parquet`` are supported.
    training_table : str or pathlib.Path or None, optional
        Preferred generic training-table path.
    model_json : str or pathlib.Path
        Output model JSON.
    summary_tsv : str or pathlib.Path or None, optional
        Backward-compatible training summary path.
    evaluation_tsv : str or pathlib.Path or None, optional
        Backward-compatible evaluation path.
    summary_table : str or pathlib.Path or None, optional
        Preferred generic training summary path.
    evaluation_table : str or pathlib.Path or None, optional
        Preferred generic evaluation path.
    label_column : str, optional
        True label column.
    feature_columns : list[str] or None, optional
        Explicit feature columns. Defaults to AI numeric features.
    test_fraction : float, optional
        Grouped test fraction.
    group_columns : list[str] or None, optional
        Grouping columns for split.
    distance_quantile : float, optional
        Open-set threshold quantile.
    unknown_label : str, optional
        Label for unresolved predictions.
    logger : logging.Logger or None, optional
        Logger.

    Returns
    -------
    tuple[PrototypeClassifier, list[dict[str, object]], list[dict[str, object]]]
        Model, test predictions and evaluation metrics.
    """
    resolved_training_table = training_table or training_tsv
    resolved_summary_table = summary_table or summary_tsv
    resolved_evaluation_table = evaluation_table or evaluation_tsv
    if resolved_training_table is None:
        raise ValueError("training_table or training_tsv is required")
    if resolved_summary_table is None:
        raise ValueError("summary_table or summary_tsv is required")
    if resolved_evaluation_table is None:
        raise ValueError("evaluation_table or evaluation_tsv is required")
    records = read_records_table(input_path=resolved_training_table, logger=logger)
    if not records:
        raise ValueError("Cannot train call calibrator from zero records")
    features = feature_columns or [
        column for column in CALL_FEATURE_COLUMNS
        if column in records[0]
    ] + ["has_long_k_support", "has_multi_k_support", "exact_hit_fraction"]
    train_records, test_records = split_records_by_group(
        records=records,
        group_columns=group_columns,
        test_fraction=test_fraction,
    )
    model = train_prototype_classifier(
        records=train_records,
        label_column=label_column,
        feature_columns=features,
        distance_quantile=distance_quantile,
        unknown_label=unknown_label,
        logger=logger,
    )
    save_model(model=model, output_path=model_json)
    summary = [
        {
            "label": label,
            "n_training_records": model.class_counts[label],
            "novelty_threshold": model.class_thresholds[label],
            "distance_quantile": model.distance_quantile,
        }
        for label in sorted(model.class_counts)
    ]
    write_records_table(
        records=summary,
        output_path=resolved_summary_table,
        fieldnames=["label", "n_training_records", "novelty_threshold", "distance_quantile"],
        logger=logger,
    )
    predictions = predict_records(records=test_records, model=model, logger=logger)
    metrics = evaluate_predictions(predictions=predictions, label_column=label_column)
    write_records_table(
        records=metrics,
        output_path=resolved_evaluation_table,
        fieldnames=["label", "n", "tp", "fp", "fn", "precision", "recall", "f1"],
        logger=logger,
    )
    return model, predictions, metrics
