"""Tests for call-calibrator holdout validation."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from kmersutra.call_validation import validate_call_calibrator_holdouts
from kmersutra.table_io import read_records_table, write_records_table


def feature_row(index: int, label: str, species: str, panel: str) -> dict[str, object]:
    """Build a compact transformed feature row."""
    value = float(index + 1)
    return {
        "sample_id": f"sample_{index}",
        "species_name": species,
        "benchmark_family": panel,
        "panel": panel,
        "ml_report_label": label,
        "log1p_n_hits": value if label != "not_detected" else 0.0,
        "log1p_n_unique_kmers": value / 2.0 if label != "not_detected" else 0.0,
        "log1p_n_positive_sequences": value / 3.0 if label != "not_detected" else 0.0,
        "n_k_values_positive": 2.0 if label != "not_detected" else 0.0,
        "log1p_best_k": 4.0 if label != "not_detected" else 0.0,
        "log1p_n_exact_hits": value if label != "not_detected" else 0.0,
        "log1p_n_fuzzy_hits": 0.0,
        "log1p_conflicting_unique_kmers": 0.0,
        "conflict_ratio": 0.0,
        "log1p_reportable_conflicting_unique_kmers": 0.0,
        "reportable_conflict_ratio": 0.0,
        "mixed_species_support_fraction": 0.5 if label == "expected_target" else 0.0,
        "confidence_score": 0.8 if label == "expected_target" else 0.1,
        "signal_confidence_score": 1.0 if label != "not_detected" else 0.0,
        "has_long_k_support": 1.0 if label != "not_detected" else 0.0,
        "has_multi_k_support": 1.0 if label != "not_detected" else 0.0,
        "exact_hit_fraction": 1.0 if label != "not_detected" else 0.0,
        "fuzzy_hit_fraction": 0.0,
        "positive_sequences_per_unique_kmer": 1.0,
        "unique_kmers_per_positive_sequence": 1.0,
        "exact_hits_per_unique_kmer": 1.0,
        "conflicting_unique_kmer_fraction": 0.0,
        "reportable_conflicting_unique_kmer_fraction": 0.0,
    }


class TestCallValidation(unittest.TestCase):
    """Test packaged holdout validation helper."""

    def test_validate_call_calibrator_holdouts_writes_outputs(self) -> None:
        """Full validation helper should create key output files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            table = root / "training.tsv.gz"
            out_dir = root / "out"
            rows = []
            labels = ["expected_target", "observed_below_threshold", "not_detected"]
            for index in range(30):
                label = labels[index % len(labels)]
                species = "Plasmodium vivax" if index % 2 else "Plasmodium falciparum"
                panel = "panel1" if index % 3 else "panel2"
                rows.append(feature_row(index, label, species, panel))
            write_records_table(
                records=rows,
                output_path=table,
                fieldnames=list(rows[0].keys()),
            )
            manifest = validate_call_calibrator_holdouts(
                training_table=table,
                out_dir=out_dir,
                feature_profile="safe_transformed",
                distance_quantile=1.0,
            )
            metrics = read_records_table(input_path=out_dir / "holdout_metrics.tsv")
            model_exists = (
                out_dir / "final_internal_calibrator_all_training.json"
            ).exists()
        self.assertEqual(manifest["n_records"], 30)
        self.assertTrue(model_exists)
        self.assertTrue(metrics)


if __name__ == "__main__":
    unittest.main()
