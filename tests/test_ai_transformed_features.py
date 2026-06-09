"""Tests for transformed AI call-calibration features."""

from __future__ import annotations

import math
import tempfile
import unittest
from pathlib import Path

from kmersutra.ai_calibration import (
    get_call_feature_columns,
    build_call_feature_record,
    train_evaluate_call_calibrator,
    write_call_training_table_from_table,
)
from kmersutra.table_io import read_records_table, write_records_table


CALL_ROWS = [
    {
        "sample_id": "s1",
        "species_name": "Plasmodium vivax",
        "is_positive_expected": "True",
        "n_hits": "1000",
        "n_unique_kmers": "100",
        "n_positive_sequences": "50",
        "n_k_values_positive": "4",
        "best_k": "151",
        "n_exact_hits": "1000",
        "n_fuzzy_hits": "0",
        "conflicting_unique_kmers": "20",
        "reportable_conflicting_unique_kmers": "10",
        "mixed_species_support_fraction": "0.5",
        "confidence_score": "0.9",
        "signal_confidence_score": "1.0",
    },
    {
        "sample_id": "s2",
        "species_name": "Plasmodium falciparum",
        "call": "observed_below_threshold",
        "n_hits": "8",
        "n_unique_kmers": "4",
        "n_positive_sequences": "2",
        "n_k_values_positive": "1",
        "best_k": "51",
        "n_exact_hits": "8",
        "n_fuzzy_hits": "0",
        "conflicting_unique_kmers": "40",
        "reportable_conflicting_unique_kmers": "30",
        "mixed_species_support_fraction": "0.02",
        "confidence_score": "0.1",
        "signal_confidence_score": "0.4",
    },
    {
        "sample_id": "s3",
        "species_name": "Toxoplasma gondii",
        "call": "not_detected",
        "n_hits": "0",
        "n_unique_kmers": "0",
        "n_positive_sequences": "0",
        "n_k_values_positive": "0",
        "best_k": "0",
        "n_exact_hits": "0",
        "n_fuzzy_hits": "0",
        "conflicting_unique_kmers": "0",
        "reportable_conflicting_unique_kmers": "0",
        "mixed_species_support_fraction": "0",
        "confidence_score": "0",
        "signal_confidence_score": "0",
    },
]


class TestAITransformedFeatures(unittest.TestCase):
    """Test log1p and leakage-safe feature profiles."""

    def test_safe_transformed_profile_excludes_spike_columns(self) -> None:
        """The transformed profile should not contain benchmark truth fields."""
        features = get_call_feature_columns(profile="safe_transformed")
        self.assertIn("log1p_n_hits", features)
        self.assertIn("positive_sequences_per_unique_kmer", features)
        self.assertNotIn("spike_n", features)
        self.assertNotIn("spike_n_per_genome", features)
        self.assertNotIn("total_spike_n", features)

    def test_bounded_transformed_profile_excludes_unbounded_ratios(self) -> None:
        """The bounded profile should drop high-magnitude ratio features."""
        features = get_call_feature_columns(profile="safe_transformed_bounded")
        self.assertIn("log1p_n_hits", features)
        self.assertIn("exact_hit_fraction", features)
        self.assertNotIn("positive_sequences_per_unique_kmer", features)
        self.assertNotIn("unique_kmers_per_positive_sequence", features)
        self.assertNotIn("exact_hits_per_unique_kmer", features)
        self.assertNotIn("spike_n", features)

    def test_build_feature_record_adds_log_and_ratio_features(self) -> None:
        """Feature construction should add log1p and ratio features."""
        record = build_call_feature_record(record=CALL_ROWS[0])
        self.assertAlmostEqual(record["log1p_n_hits"], math.log1p(1000.0))
        self.assertAlmostEqual(record["positive_sequences_per_unique_kmer"], 0.5)
        self.assertAlmostEqual(record["unique_kmers_per_positive_sequence"], 2.0)
        self.assertAlmostEqual(record["conflicting_unique_kmer_fraction"], 20.0 / 120.0)

    def test_training_cli_helper_accepts_safe_transformed_profile(self) -> None:
        """Call calibrator helper should train using transformed features."""
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            calls = root / "calls.tsv.gz"
            training = root / "training.tsv.gz"
            model = root / "model.json"
            summary = root / "summary.tsv"
            evaluation = root / "evaluation.tsv"
            rows = []
            for index in range(8):
                for row in CALL_ROWS:
                    copied = dict(row)
                    copied["sample_id"] = f"{row['sample_id']}_{index}"
                    rows.append(copied)
            write_records_table(
                records=rows,
                output_path=calls,
                fieldnames=list(rows[0].keys()),
            )
            write_call_training_table_from_table(
                calls_table=calls,
                output_table=training,
            )
            trained, predictions, metrics = train_evaluate_call_calibrator(
                training_table=training,
                model_json=model,
                summary_table=summary,
                evaluation_table=evaluation,
                feature_profile="safe_transformed_bounded",
                test_fraction=0.25,
                distance_quantile=1.0,
            )
            observed_summary = read_records_table(input_path=summary)
        self.assertIn("log1p_n_hits", trained.feature_columns)
        self.assertNotIn("spike_n", trained.feature_columns)
        self.assertTrue(predictions)
        self.assertTrue(metrics)
        self.assertTrue(observed_summary)


if __name__ == "__main__":
    unittest.main()
