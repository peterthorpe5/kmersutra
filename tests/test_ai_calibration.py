"""Tests for AI-ready call calibration helpers."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from kmersutra.ai_calibration import (
    build_call_feature_record,
    build_call_training_table,
    evaluate_predictions,
    infer_report_label,
    split_records_by_group,
    train_evaluate_call_calibrator,
    write_call_training_table_from_table,
    write_call_training_table_from_tsv,
)
from kmersutra.io import read_tsv, write_tsv
from kmersutra.table_io import read_records_table


CALL_ROWS = [
    {
        "sample_id": "s1",
        "species_name": "Plasmodium vivax",
        "benchmark_family": "single_genome",
        "panel": "panel1",
        "is_positive_expected": "True",
        "is_positive_call": "True",
        "n_hits": "100",
        "n_unique_kmers": "50",
        "n_positive_sequences": "20",
        "n_k_values_positive": "3",
        "best_k": "151",
        "n_exact_hits": "100",
        "conflict_ratio": "0.0",
        "confidence_score": "0.95",
    },
    {
        "sample_id": "s2",
        "species_name": "Hammondia hammondi",
        "report_layer": "background_candidate_signal",
        "is_positive_call": "True",
        "n_hits": "80",
        "n_unique_kmers": "28",
        "n_positive_sequences": "73",
        "n_k_values_positive": "4",
        "best_k": "151",
        "n_exact_hits": "80",
        "conflict_ratio": "0.1",
        "confidence_score": "0.7",
    },
    {
        "sample_id": "s3",
        "species_name": "Plasmodium simium",
        "report_layer": "neighbour_lineage_evidence",
        "is_positive_call": "True",
        "n_hits": "50",
        "n_unique_kmers": "20",
        "n_positive_sequences": "10",
        "n_k_values_positive": "2",
        "best_k": "101",
        "n_exact_hits": "50",
        "conflict_ratio": "0.2",
        "confidence_score": "0.6",
    },
    {
        "sample_id": "s4",
        "species_name": "Toxoplasma gondii",
        "is_positive_off_target": "True",
        "is_positive_call": "True",
        "n_hits": "40",
        "n_unique_kmers": "10",
        "n_positive_sequences": "5",
        "n_k_values_positive": "1",
        "best_k": "77",
        "n_exact_hits": "40",
        "conflict_ratio": "0.3",
        "confidence_score": "0.5",
    },
    {
        "sample_id": "s5",
        "species_name": "Babesia bigemina",
        "call": "not_detected",
        "is_positive_call": "False",
        "n_hits": "0",
        "n_unique_kmers": "0",
        "n_positive_sequences": "0",
        "n_k_values_positive": "0",
        "best_k": "0",
        "n_exact_hits": "0",
        "conflict_ratio": "0",
        "confidence_score": "0",
    },
]


class TestAICalibration(unittest.TestCase):
    """Test AI-ready calibration data preparation and modelling."""

    def test_infer_report_label_uses_report_layers(self) -> None:
        """Training labels should follow report-layer categories."""
        self.assertEqual(infer_report_label(record=CALL_ROWS[0]), "expected_target")
        self.assertEqual(
            infer_report_label(record=CALL_ROWS[1]),
            "background_candidate_signal",
        )
        self.assertEqual(
            infer_report_label(record=CALL_ROWS[2]),
            "neighbour_lineage_evidence",
        )
        self.assertEqual(
            infer_report_label(record=CALL_ROWS[3]),
            "reportable_off_target_species",
        )
        self.assertEqual(infer_report_label(record=CALL_ROWS[4]), "not_detected")

    def test_build_call_feature_record_adds_derived_features(self) -> None:
        """Feature records should include numeric and derived AI features."""
        record = build_call_feature_record(record=CALL_ROWS[0])
        self.assertEqual(record["ml_report_label"], "expected_target")
        self.assertEqual(record["has_long_k_support"], 1.0)
        self.assertEqual(record["has_multi_k_support"], 1.0)
        self.assertEqual(record["exact_hit_fraction"], 1.0)

    def test_build_call_training_table_can_drop_not_detected(self) -> None:
        """Training construction should optionally remove not-detected rows."""
        records = build_call_training_table(
            records=CALL_ROWS,
            include_not_detected=False,
        )
        labels = {row["ml_report_label"] for row in records}
        self.assertNotIn("not_detected", labels)

    def test_split_records_by_group_is_deterministic(self) -> None:
        """Grouped splitting should be repeatable."""
        first = split_records_by_group(records=CALL_ROWS, test_fraction=0.4)
        second = split_records_by_group(records=CALL_ROWS, test_fraction=0.4)
        self.assertEqual(first, second)

    def test_evaluate_predictions_outputs_overall_metric(self) -> None:
        """Evaluation should include per-label and overall rows."""
        predictions = [
            {"ml_report_label": "a", "prediction": "a"},
            {"ml_report_label": "b", "prediction": "a"},
        ]
        metrics = evaluate_predictions(predictions=predictions)
        labels = {row["label"] for row in metrics}
        self.assertIn("overall", labels)


    def test_write_training_table_accepts_compressed_tsv(self) -> None:
        """AI table construction should support TSV.GZ inputs and outputs."""
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            calls_tsv = root / "calls.tsv.gz"
            training_tsv = root / "training.tsv.gz"
            write_tsv(
                records=CALL_ROWS,
                output_path=calls_tsv,
                fieldnames=list(CALL_ROWS[0].keys()),
            )
            features = write_call_training_table_from_table(
                calls_table=calls_tsv,
                output_table=training_tsv,
            )
            observed = read_records_table(input_path=training_tsv)
        self.assertTrue(features)
        self.assertEqual(len(observed), len(features))

    def test_train_calibrator_accepts_compressed_tables(self) -> None:
        """Call-calibrator training should support compressed table paths."""
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            calls_tsv = root / "calls.tsv.gz"
            training_tsv = root / "training.tsv.gz"
            model_json = root / "model.json"
            summary_tsv = root / "summary.tsv.gz"
            evaluation_tsv = root / "evaluation.tsv.gz"
            repeated_rows = []
            for index in range(6):
                for row in CALL_ROWS:
                    copied = dict(row)
                    copied["sample_id"] = f"{row['sample_id']}_{index}"
                    repeated_rows.append(copied)
            write_tsv(
                records=repeated_rows,
                output_path=calls_tsv,
                fieldnames=list(repeated_rows[0].keys()),
            )
            write_call_training_table_from_table(
                calls_table=calls_tsv,
                output_table=training_tsv,
            )
            _model, predictions, metrics = train_evaluate_call_calibrator(
                training_table=training_tsv,
                model_json=model_json,
                summary_table=summary_tsv,
                evaluation_table=evaluation_tsv,
                test_fraction=0.25,
                distance_quantile=1.0,
            )
            summary_records = read_records_table(input_path=summary_tsv)
            evaluation_records = read_records_table(input_path=evaluation_tsv)
        self.assertTrue(predictions)
        self.assertTrue(metrics)
        self.assertTrue(summary_records)
        self.assertTrue(evaluation_records)

    def test_write_training_table_and_train_calibrator(self) -> None:
        """CLI-facing helper functions should write model and evaluation files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            calls_tsv = root / "calls.tsv"
            training_tsv = root / "training.tsv"
            model_json = root / "model.json"
            summary_tsv = root / "summary.tsv"
            evaluation_tsv = root / "evaluation.tsv"
            repeated_rows = []
            for index in range(6):
                for row in CALL_ROWS:
                    copied = dict(row)
                    copied["sample_id"] = f"{row['sample_id']}_{index}"
                    repeated_rows.append(copied)
            write_tsv(
                records=repeated_rows,
                output_path=calls_tsv,
                fieldnames=list(repeated_rows[0].keys()),
            )
            features = write_call_training_table_from_tsv(
                calls_tsv=calls_tsv,
                output_tsv=training_tsv,
            )
            model, predictions, metrics = train_evaluate_call_calibrator(
                training_tsv=training_tsv,
                model_json=model_json,
                summary_tsv=summary_tsv,
                evaluation_tsv=evaluation_tsv,
                test_fraction=0.25,
                distance_quantile=1.0,
            )
            self.assertTrue(features)
            self.assertTrue(model_json.exists())
            self.assertTrue(summary_tsv.exists())
            self.assertTrue(evaluation_tsv.exists())
            self.assertTrue(predictions)
            self.assertTrue(metrics)
            self.assertIn("expected_target", model.class_counts)
            self.assertTrue(read_tsv(input_path=evaluation_tsv))


if __name__ == "__main__":
    unittest.main()

class TestAICalibrationLcaFeatures(unittest.TestCase):
    """Test optional LCA-derived AI features."""

    def test_write_training_table_can_merge_lca_features(self) -> None:
        """AI training output should include numeric LCA features when supplied."""
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            calls_tsv = root / "calls.tsv.gz"
            lca_tsv = root / "lca.tsv.gz"
            training_tsv = root / "training.tsv.gz"
            lca_rows = [
                {
                    "sample_id": "s1",
                    "lca_scope": "dominant_lineage",
                    "lca_taxid": "11",
                    "lca_rank": "species",
                    "n_taxa": "1",
                    "total_unique_kmers": "50",
                    "total_positive_sequences": "20",
                    "max_best_k": "151",
                    "max_k_values_positive": "3",
                    "top_unique_kmers": "50",
                    "top_positive_sequences": "20",
                    "top_best_k": "151",
                    "top_score": "1.0",
                },
                {
                    "sample_id": "s2",
                    "lca_scope": "all_supported_evidence",
                    "lca_taxid": "10",
                    "lca_rank": "genus",
                    "n_taxa": "2",
                    "total_unique_kmers": "28",
                    "total_positive_sequences": "73",
                    "max_best_k": "151",
                    "max_k_values_positive": "4",
                    "top_unique_kmers": "28",
                    "top_positive_sequences": "73",
                    "top_best_k": "151",
                    "top_score": "1.0",
                },
            ]
            write_tsv(
                records=CALL_ROWS,
                output_path=calls_tsv,
                fieldnames=list(CALL_ROWS[0].keys()),
            )
            write_tsv(
                records=lca_rows,
                output_path=lca_tsv,
                fieldnames=list(lca_rows[0].keys()),
            )
            write_call_training_table_from_table(
                calls_table=calls_tsv,
                output_table=training_tsv,
                lca_table=lca_tsv,
            )
            observed = read_records_table(input_path=training_tsv)
        s1 = next(row for row in observed if row["sample_id"] == "s1")
        s2 = next(row for row in observed if row["sample_id"] == "s2")
        self.assertEqual(float(s1["lca_dominant_lineage_rank_depth"]), 7.0)
        self.assertEqual(float(s1["lca_dominant_lineage_is_species"]), 1.0)
        self.assertEqual(float(s2["lca_all_supported_evidence_rank_depth"]), 6.0)

    def test_train_calibrator_infers_lca_feature_columns(self) -> None:
        """Default training should include lca_* features when present."""
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            training_tsv = root / "training.tsv"
            model_json = root / "model.json"
            summary_tsv = root / "summary.tsv"
            evaluation_tsv = root / "evaluation.tsv"
            records = build_call_training_table(records=CALL_ROWS * 4)
            for index, row in enumerate(records):
                row["sample_id"] = f"{row['sample_id']}_{index}"
                row["lca_dominant_lineage_rank_depth"] = 7.0 if index % 2 else 6.0
            write_tsv(
                records=records,
                output_path=training_tsv,
                fieldnames=list(records[0].keys()),
            )
            model, _predictions, _metrics = train_evaluate_call_calibrator(
                training_table=training_tsv,
                model_json=model_json,
                summary_table=summary_tsv,
                evaluation_table=evaluation_tsv,
                test_fraction=0.25,
                distance_quantile=1.0,
            )
        self.assertIn("lca_dominant_lineage_rank_depth", model.feature_columns)
