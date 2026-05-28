"""Tests for the lightweight KmerSutra-ML classifier."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from kmersutra.io import read_tsv, write_tsv
from kmersutra.ml import (
    calculate_global_scaling,
    euclidean_distance,
    load_model,
    predict_record,
    predict_records,
    predict_tsv,
    save_model,
    standardise_record,
    train_model_from_tsv,
    train_prototype_classifier,
)


TRAINING_RECORDS = [
    {
        "sequence_id": "a1",
        "true_species": "Species A",
        "n_total_hits": 10,
        "n_unique_species_kmers": 8,
        "best_species_hits": 9,
        "conflict_ratio": 0.0,
        "longest_k": 101,
    },
    {
        "sequence_id": "a2",
        "true_species": "Species A",
        "n_total_hits": 11,
        "n_unique_species_kmers": 9,
        "best_species_hits": 10,
        "conflict_ratio": 0.0,
        "longest_k": 101,
    },
    {
        "sequence_id": "b1",
        "true_species": "Species B",
        "n_total_hits": 3,
        "n_unique_species_kmers": 2,
        "best_species_hits": 3,
        "conflict_ratio": 0.2,
        "longest_k": 71,
    },
    {
        "sequence_id": "b2",
        "true_species": "Species B",
        "n_total_hits": 4,
        "n_unique_species_kmers": 2,
        "best_species_hits": 4,
        "conflict_ratio": 0.1,
        "longest_k": 71,
    },
]

FEATURE_COLUMNS = [
    "n_total_hits",
    "n_unique_species_kmers",
    "best_species_hits",
    "conflict_ratio",
    "longest_k",
]


class TestMachineLearning(unittest.TestCase):
    """Test training, saving, loading and prediction behaviour."""

    def test_calculate_global_scaling_outputs_nonzero_stds(self):
        """Calculate scaling parameters with protected non-zero standard deviations."""
        means, stds = calculate_global_scaling(
            records=TRAINING_RECORDS,
            feature_columns=FEATURE_COLUMNS,
        )
        self.assertIn("n_total_hits", means)
        self.assertGreater(stds["n_total_hits"], 0.0)

    def test_standardise_record_centres_values(self):
        """Standardise a feature record using supplied means and standard deviations."""
        means, stds = calculate_global_scaling(
            records=TRAINING_RECORDS,
            feature_columns=FEATURE_COLUMNS,
        )
        vector = standardise_record(
            record=TRAINING_RECORDS[0],
            feature_columns=FEATURE_COLUMNS,
            means=means,
            stds=stds,
        )
        self.assertEqual(sorted(vector), sorted(FEATURE_COLUMNS))

    def test_euclidean_distance_is_zero_for_identical_vectors(self):
        """Return zero distance when feature vectors are identical."""
        vector = {column: 1.0 for column in FEATURE_COLUMNS}
        distance = euclidean_distance(
            left=vector,
            right=vector,
            feature_columns=FEATURE_COLUMNS,
        )
        self.assertEqual(distance, 0.0)

    def test_train_prototype_classifier_learns_labels(self):
        """Train a prototype classifier and retain class counts."""
        model = train_prototype_classifier(
            records=TRAINING_RECORDS,
            label_column="true_species",
            feature_columns=FEATURE_COLUMNS,
            distance_quantile=1.0,
        )
        self.assertEqual(set(model.class_counts), {"Species A", "Species B"})
        self.assertEqual(model.class_counts["Species A"], 2)

    def test_predict_known_like_record(self):
        """Predict a training-like record as the closest known species."""
        model = train_prototype_classifier(
            records=TRAINING_RECORDS,
            label_column="true_species",
            feature_columns=FEATURE_COLUMNS,
            distance_quantile=1.0,
        )
        prediction = predict_record(record=TRAINING_RECORDS[0], model=model)
        self.assertEqual(prediction["prediction"], "Species A")
        self.assertEqual(prediction["open_set_status"], "known_like")

    def test_predict_unknown_record(self):
        """Flag a distant record as unknown or unresolved."""
        model = train_prototype_classifier(
            records=TRAINING_RECORDS,
            label_column="true_species",
            feature_columns=FEATURE_COLUMNS,
            distance_quantile=0.5,
        )
        distant_record = {
            "sequence_id": "x1",
            "n_total_hits": 100,
            "n_unique_species_kmers": 0,
            "best_species_hits": 0,
            "conflict_ratio": 0.9,
            "longest_k": 31,
        }
        prediction = predict_record(record=distant_record, model=model)
        self.assertEqual(prediction["prediction"], "unknown_or_unresolved")
        self.assertEqual(prediction["open_set_status"], "unknown_or_unresolved")

    def test_predict_records_counts_multiple_rows(self):
        """Predict multiple feature rows in one call."""
        model = train_prototype_classifier(
            records=TRAINING_RECORDS,
            label_column="true_species",
            feature_columns=FEATURE_COLUMNS,
            distance_quantile=1.0,
        )
        predictions = predict_records(records=TRAINING_RECORDS[:2], model=model)
        self.assertEqual(len(predictions), 2)
        self.assertTrue(all("ml_confidence_score" in row for row in predictions))

    def test_save_and_load_model_round_trip(self):
        """Save a trained model to JSON and load it back."""
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "model.json"
            model = train_prototype_classifier(
                records=TRAINING_RECORDS,
                label_column="true_species",
                feature_columns=FEATURE_COLUMNS,
                distance_quantile=1.0,
            )
            save_model(model=model, output_path=path)
            loaded = load_model(model_path=path)
            self.assertEqual(loaded.feature_columns, FEATURE_COLUMNS)
            self.assertEqual(loaded.class_counts, model.class_counts)

    def test_train_model_from_tsv_writes_model_and_summary(self):
        """Train from TSV and write both model JSON and summary TSV."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            features_tsv = tmp_path / "features.tsv"
            model_json = tmp_path / "model.json"
            summary_tsv = tmp_path / "summary.tsv"
            write_tsv(
                records=TRAINING_RECORDS,
                output_path=features_tsv,
                fieldnames=list(TRAINING_RECORDS[0].keys()),
            )
            model = train_model_from_tsv(
                features_tsv=features_tsv,
                label_column="true_species",
                model_json=model_json,
                summary_tsv=summary_tsv,
                distance_quantile=1.0,
            )
            self.assertTrue(model_json.exists())
            self.assertTrue(summary_tsv.exists())
            self.assertEqual(set(model.class_counts), {"Species A", "Species B"})

    def test_predict_tsv_writes_predictions(self):
        """Predict from a TSV feature table and write prediction output."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            features_tsv = tmp_path / "features.tsv"
            model_json = tmp_path / "model.json"
            predictions_tsv = tmp_path / "predictions.tsv"
            write_tsv(
                records=TRAINING_RECORDS,
                output_path=features_tsv,
                fieldnames=list(TRAINING_RECORDS[0].keys()),
            )
            model = train_prototype_classifier(
                records=TRAINING_RECORDS,
                label_column="true_species",
                feature_columns=FEATURE_COLUMNS,
                distance_quantile=1.0,
            )
            save_model(model=model, output_path=model_json)
            predictions = predict_tsv(
                features_tsv=features_tsv,
                model_json=model_json,
                output_tsv=predictions_tsv,
            )
            loaded = read_tsv(input_path=predictions_tsv)
            self.assertEqual(len(predictions), 4)
            self.assertEqual(len(loaded), 4)
            self.assertIn("prediction", loaded[0])

    def test_train_classifier_rejects_missing_label(self):
        """Raise a useful error when the requested label column is absent."""
        with self.assertRaises(ValueError):
            train_prototype_classifier(
                records=TRAINING_RECORDS,
                label_column="missing_label",
                feature_columns=FEATURE_COLUMNS,
            )

    def test_predict_rejects_bad_novelty_scale(self):
        """Reject non-positive novelty scale values."""
        model = train_prototype_classifier(
            records=TRAINING_RECORDS,
            label_column="true_species",
            feature_columns=FEATURE_COLUMNS,
        )
        with self.assertRaises(ValueError):
            predict_record(record=TRAINING_RECORDS[0], model=model, novelty_scale=0.0)


if __name__ == "__main__":
    unittest.main()
