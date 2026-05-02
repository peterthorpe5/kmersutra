"""Lightweight open-set classifier for KmerSutra feature tables.

The first KmerSutra-ML model is intentionally simple and interpretable. It
uses per-class feature centroids and a distance-to-centroid novelty threshold.
This avoids adding heavy machine-learning dependencies while providing a
reproducible baseline for clade/species classification and unresolved calls.
"""

from __future__ import annotations

import logging
import math
from collections import Counter, defaultdict
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

from kmersutra.features import infer_numeric_feature_columns
from kmersutra.io import read_tsv, write_json, write_tsv
import json


@dataclass(frozen=True)
class PrototypeClassifier:
    """Centroid-based open-set classifier.

    Attributes
    ----------
    label_column : str
        Training label column.
    feature_columns : list[str]
        Numeric feature columns used by the model.
    class_centroids : dict[str, dict[str, float]]
        Mean feature values by class.
    class_thresholds : dict[str, float]
        Class-specific novelty thresholds.
    global_means : dict[str, float]
        Global feature means used for standardisation.
    global_stds : dict[str, float]
        Global feature standard deviations used for standardisation.
    class_counts : dict[str, int]
        Number of training records by class.
    unknown_label : str
        Label returned when evidence is outside class thresholds.
    distance_quantile : float
        Quantile used when estimating novelty thresholds.
    """

    label_column: str
    feature_columns: list[str]
    class_centroids: dict[str, dict[str, float]]
    class_thresholds: dict[str, float]
    global_means: dict[str, float]
    global_stds: dict[str, float]
    class_counts: dict[str, int]
    unknown_label: str = "unknown_or_unresolved"
    distance_quantile: float = 0.95

    def to_dict(self) -> dict[str, object]:
        """Return a JSON-serialisable model dictionary.

        Returns
        -------
        dict[str, object]
            Model dictionary.
        """
        return {
            "model_type": "prototype_open_set",
            "label_column": self.label_column,
            "feature_columns": self.feature_columns,
            "class_centroids": self.class_centroids,
            "class_thresholds": self.class_thresholds,
            "global_means": self.global_means,
            "global_stds": self.global_stds,
            "class_counts": self.class_counts,
            "unknown_label": self.unknown_label,
            "distance_quantile": self.distance_quantile,
        }

    @classmethod
    def from_dict(cls, data: dict[str, object]) -> "PrototypeClassifier":
        """Create a classifier from a model dictionary.

        Parameters
        ----------
        data : dict[str, object]
            Model dictionary.

        Returns
        -------
        PrototypeClassifier
            Loaded classifier.
        """
        if data.get("model_type") != "prototype_open_set":
            raise ValueError("Unsupported model type")
        return cls(
            label_column=str(data["label_column"]),
            feature_columns=list(data["feature_columns"]),
            class_centroids={
                str(label): {str(key): float(value) for key, value in values.items()}
                for label, values in dict(data["class_centroids"]).items()
            },
            class_thresholds={
                str(label): float(value)
                for label, value in dict(data["class_thresholds"]).items()
            },
            global_means={
                str(key): float(value)
                for key, value in dict(data["global_means"]).items()
            },
            global_stds={
                str(key): float(value)
                for key, value in dict(data["global_stds"]).items()
            },
            class_counts={
                str(label): int(value)
                for label, value in dict(data["class_counts"]).items()
            },
            unknown_label=str(data.get("unknown_label", "unknown_or_unresolved")),
            distance_quantile=float(data.get("distance_quantile", 0.95)),
        )


def _to_float(value: object) -> float:
    """Convert values to float, treating missing values as zero."""
    if value in {None, ""}:
        return 0.0
    return float(value)


def _quantile(values: list[float], quantile: float) -> float:
    """Return a deterministic nearest-rank quantile.

    Parameters
    ----------
    values : list[float]
        Values to summarise.
    quantile : float
        Quantile between zero and one.

    Returns
    -------
    float
        Quantile value.
    """
    if not values:
        return 0.0
    if not 0.0 <= quantile <= 1.0:
        raise ValueError("quantile must be between 0 and 1")
    ordered = sorted(values)
    index = math.ceil(quantile * len(ordered)) - 1
    index = max(0, min(index, len(ordered) - 1))
    return ordered[index]


def calculate_global_scaling(
    *,
    records: list[dict[str, object]],
    feature_columns: list[str],
) -> tuple[dict[str, float], dict[str, float]]:
    """Calculate global means and standard deviations.

    Parameters
    ----------
    records : list[dict[str, object]]
        Training feature records.
    feature_columns : list[str]
        Numeric feature columns.

    Returns
    -------
    tuple[dict[str, float], dict[str, float]]
        Means and standard deviations.
    """
    means: dict[str, float] = {}
    stds: dict[str, float] = {}
    for column in feature_columns:
        values = [_to_float(row.get(column, 0.0)) for row in records]
        mean = sum(values) / len(values) if values else 0.0
        variance = (
            sum((value - mean) ** 2 for value in values) / len(values)
            if values
            else 0.0
        )
        std = math.sqrt(variance) or 1.0
        means[column] = mean
        stds[column] = std
    return means, stds


def standardise_record(
    *,
    record: dict[str, object],
    feature_columns: list[str],
    means: dict[str, float],
    stds: dict[str, float],
) -> dict[str, float]:
    """Standardise one feature record.

    Parameters
    ----------
    record : dict[str, object]
        Feature row.
    feature_columns : list[str]
        Numeric feature columns.
    means : dict[str, float]
        Global feature means.
    stds : dict[str, float]
        Global feature standard deviations.

    Returns
    -------
    dict[str, float]
        Standardised feature values.
    """
    return {
        column: (_to_float(record.get(column, 0.0)) - means[column]) / stds[column]
        for column in feature_columns
    }


def euclidean_distance(
    *,
    left: dict[str, float],
    right: dict[str, float],
    feature_columns: list[str],
) -> float:
    """Calculate Euclidean distance between standardised records.

    Parameters
    ----------
    left : dict[str, float]
        First feature vector.
    right : dict[str, float]
        Second feature vector.
    feature_columns : list[str]
        Feature columns to compare.

    Returns
    -------
    float
        Euclidean distance.
    """
    return math.sqrt(sum((left[column] - right[column]) ** 2 for column in feature_columns))


def train_prototype_classifier(
    *,
    records: Iterable[dict[str, object]],
    label_column: str,
    feature_columns: list[str] | None = None,
    distance_quantile: float = 0.95,
    unknown_label: str = "unknown_or_unresolved",
    logger: logging.Logger | None = None,
) -> PrototypeClassifier:
    """Train a centroid-based open-set classifier.

    Parameters
    ----------
    records : iterable of dict[str, object]
        Feature table rows containing labels.
    label_column : str
        Column containing the training label.
    feature_columns : list[str] | None, optional
        Explicit numeric feature columns. If omitted, numeric columns are
        inferred.
    distance_quantile : float, optional
        Quantile of within-class distances used as novelty threshold.
    unknown_label : str, optional
        Output label for unresolved/open-set predictions.
    logger : logging.Logger | None, optional
        Logger for progress messages.

    Returns
    -------
    PrototypeClassifier
        Trained classifier.
    """
    rows = [dict(row) for row in records]
    if not rows:
        raise ValueError("Cannot train classifier from zero records")
    if label_column not in rows[0]:
        raise ValueError(f"Training label column not found: {label_column}")
    if not 0.0 <= distance_quantile <= 1.0:
        raise ValueError("distance_quantile must be between 0 and 1")

    features = feature_columns or infer_numeric_feature_columns(
        records=rows,
        label_column=label_column,
    )
    if not features:
        raise ValueError("No numeric feature columns are available for training")

    labelled_rows = [row for row in rows if str(row.get(label_column, ""))]
    if not labelled_rows:
        raise ValueError("No labelled training records found")

    means, stds = calculate_global_scaling(records=labelled_rows, feature_columns=features)
    by_label: dict[str, list[dict[str, float]]] = defaultdict(list)
    for row in labelled_rows:
        label = str(row[label_column])
        by_label[label].append(
            standardise_record(
                record=row,
                feature_columns=features,
                means=means,
                stds=stds,
            )
        )

    if logger:
        logger.info(
            "Training open-set prototype classifier with %d labels and %d features",
            len(by_label),
            len(features),
        )

    centroids: dict[str, dict[str, float]] = {}
    thresholds: dict[str, float] = {}
    class_counts: dict[str, int] = {}
    for label, vectors in sorted(by_label.items()):
        class_counts[label] = len(vectors)
        centroid = {
            column: sum(vector[column] for vector in vectors) / len(vectors)
            for column in features
        }
        centroids[label] = centroid
        distances = [
            euclidean_distance(left=vector, right=centroid, feature_columns=features)
            for vector in vectors
        ]
        threshold = _quantile(distances, distance_quantile)
        thresholds[label] = max(threshold, 1e-9)
        if logger:
            logger.info(
                "Label %s: n=%d threshold=%.4f",
                label,
                len(vectors),
                thresholds[label],
            )

    return PrototypeClassifier(
        label_column=label_column,
        feature_columns=features,
        class_centroids=centroids,
        class_thresholds=thresholds,
        global_means=means,
        global_stds=stds,
        class_counts=class_counts,
        unknown_label=unknown_label,
        distance_quantile=distance_quantile,
    )


def predict_record(
    *,
    record: dict[str, object],
    model: PrototypeClassifier,
    novelty_scale: float = 1.0,
) -> dict[str, object]:
    """Predict a label for one feature record.

    Parameters
    ----------
    record : dict[str, object]
        Feature row.
    model : PrototypeClassifier
        Trained classifier.
    novelty_scale : float, optional
        Multiplier applied to class novelty thresholds.

    Returns
    -------
    dict[str, object]
        Prediction row.
    """
    if novelty_scale <= 0:
        raise ValueError("novelty_scale must be positive")
    vector = standardise_record(
        record=record,
        feature_columns=model.feature_columns,
        means=model.global_means,
        stds=model.global_stds,
    )
    distances = {
        label: euclidean_distance(
            left=vector,
            right=centroid,
            feature_columns=model.feature_columns,
        )
        for label, centroid in model.class_centroids.items()
    }
    ranked = sorted(distances.items(), key=lambda item: (item[1], item[0]))
    best_label, best_distance = ranked[0]
    second_label, second_distance = ranked[1] if len(ranked) > 1 else ("", math.inf)
    threshold = model.class_thresholds[best_label] * novelty_scale
    in_distribution = best_distance <= threshold
    prediction = best_label if in_distribution else model.unknown_label
    confidence = 1.0 / (1.0 + (best_distance / (threshold or 1e-9)))
    margin = second_distance - best_distance if math.isfinite(second_distance) else best_distance
    return {
        **record,
        "prediction": prediction,
        "best_label": best_label,
        "best_distance": round(best_distance, 6),
        "best_threshold": round(threshold, 6),
        "second_label": second_label,
        "second_distance": round(second_distance, 6) if math.isfinite(second_distance) else "",
        "distance_margin": round(margin, 6),
        "open_set_status": "known_like" if in_distribution else "unknown_or_unresolved",
        "ml_confidence_score": round(max(0.0, min(1.0, confidence)), 6),
    }


def predict_records(
    *,
    records: Iterable[dict[str, object]],
    model: PrototypeClassifier,
    novelty_scale: float = 1.0,
    logger: logging.Logger | None = None,
) -> list[dict[str, object]]:
    """Predict labels for feature records.

    Parameters
    ----------
    records : iterable of dict[str, object]
        Feature rows.
    model : PrototypeClassifier
        Trained model.
    novelty_scale : float, optional
        Multiplier for novelty thresholds.
    logger : logging.Logger | None, optional
        Logger for progress messages.

    Returns
    -------
    list[dict[str, object]]
        Prediction rows.
    """
    rows = list(records)
    if logger:
        logger.info("Predicting labels for %d feature records", len(rows))
    predictions = [
        predict_record(record=row, model=model, novelty_scale=novelty_scale)
        for row in rows
    ]
    if logger:
        counts = Counter(str(row["prediction"]) for row in predictions)
        for label, count in sorted(counts.items()):
            logger.info("Prediction %s: %d", label, count)
    return predictions


def save_model(*, model: PrototypeClassifier, output_path: str | Path) -> None:
    """Write a trained model to JSON.

    Parameters
    ----------
    model : PrototypeClassifier
        Model to save.
    output_path : str or pathlib.Path
        Output JSON path.
    """
    write_json(data=model.to_dict(), output_path=output_path)


def load_model(*, model_path: str | Path) -> PrototypeClassifier:
    """Load a trained model from JSON.

    Parameters
    ----------
    model_path : str or pathlib.Path
        Model JSON path.

    Returns
    -------
    PrototypeClassifier
        Loaded model.
    """
    data = json.loads(Path(model_path).read_text(encoding="utf-8"))
    return PrototypeClassifier.from_dict(data)


def train_model_from_tsv(
    *,
    features_tsv: str | Path,
    label_column: str,
    model_json: str | Path,
    summary_tsv: str | Path,
    distance_quantile: float = 0.95,
    unknown_label: str = "unknown_or_unresolved",
    logger: logging.Logger | None = None,
) -> PrototypeClassifier:
    """Train and save a prototype classifier from a TSV file.

    Parameters
    ----------
    features_tsv : str or pathlib.Path
        Input feature TSV containing labels.
    label_column : str
        Label column.
    model_json : str or pathlib.Path
        Output model JSON.
    summary_tsv : str or pathlib.Path
        Output training summary TSV.
    distance_quantile : float, optional
        Novelty threshold quantile.
    unknown_label : str, optional
        Unknown label.
    logger : logging.Logger | None, optional
        Logger.

    Returns
    -------
    PrototypeClassifier
        Trained model.
    """
    records = read_tsv(input_path=features_tsv)
    if logger:
        logger.info("Loaded %d training feature records from %s", len(records), features_tsv)
    model = train_prototype_classifier(
        records=records,
        label_column=label_column,
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
    write_tsv(
        records=summary,
        output_path=summary_tsv,
        fieldnames=["label", "n_training_records", "novelty_threshold", "distance_quantile"],
    )
    return model


def predict_tsv(
    *,
    features_tsv: str | Path,
    model_json: str | Path,
    output_tsv: str | Path,
    novelty_scale: float = 1.0,
    logger: logging.Logger | None = None,
) -> list[dict[str, object]]:
    """Predict labels for a feature TSV.

    Parameters
    ----------
    features_tsv : str or pathlib.Path
        Input feature table.
    model_json : str or pathlib.Path
        Trained model JSON.
    output_tsv : str or pathlib.Path
        Output prediction TSV.
    novelty_scale : float, optional
        Multiplier for open-set thresholds.
    logger : logging.Logger | None, optional
        Logger.

    Returns
    -------
    list[dict[str, object]]
        Prediction records.
    """
    model = load_model(model_path=model_json)
    records = read_tsv(input_path=features_tsv)
    predictions = predict_records(
        records=records,
        model=model,
        novelty_scale=novelty_scale,
        logger=logger,
    )
    prediction_columns = [
        "prediction",
        "best_label",
        "best_distance",
        "best_threshold",
        "second_label",
        "second_distance",
        "distance_margin",
        "open_set_status",
        "ml_confidence_score",
    ]
    original_columns = list(records[0].keys()) if records else []
    fieldnames = original_columns + [column for column in prediction_columns if column not in original_columns]
    write_tsv(records=predictions, output_path=output_tsv, fieldnames=fieldnames)
    return predictions
