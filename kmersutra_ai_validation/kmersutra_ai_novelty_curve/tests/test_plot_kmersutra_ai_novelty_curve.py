"""Unit tests for plot_kmersutra_ai_novelty_curve.py."""

from __future__ import annotations

import csv
import gzip
import tempfile
import unittest
from pathlib import Path

import sys

SCRIPT_DIR = Path(__file__).resolve().parents[1] / "scripts" / "ai_python"
sys.path.insert(0, str(SCRIPT_DIR))

import plot_kmersutra_ai_novelty_curve as plotter  # noqa: E402


class TestKmerSutraAiNoveltyCurve(unittest.TestCase):
    """Test novelty-sweep summarisation helpers."""

    def test_summarise_prediction_table(self) -> None:
        """A small prediction table is summarised correctly."""
        with tempfile.TemporaryDirectory() as temp_dir:
            out_dir = Path(temp_dir) / "run" / "outputs"
            out_dir.mkdir(parents=True)
            pred_path = out_dir / "external_zymo_predictions.tsv.gz"
            log_path = out_dir / "external_zymo_predictions.tsv.log"

            rows = [
                {
                    "zymo_truth_category": "expected_species",
                    "prediction": "expected_target",
                },
                {
                    "zymo_truth_category": "expected_reference_label",
                    "prediction": "unknown_or_unresolved",
                },
                {
                    "zymo_truth_category": "not_detected",
                    "prediction": "not_detected",
                },
                {
                    "zymo_truth_category": "near_neighbour_evidence",
                    "prediction": "expected_target",
                },
                {
                    "zymo_truth_category": "same_species_compatible_reference",
                    "prediction": "expected_target",
                },
            ]

            with gzip.open(pred_path, "wt", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    delimiter="\t",
                    fieldnames=["zymo_truth_category", "prediction"],
                    lineterminator="\n",
                )
                writer.writeheader()
                writer.writerows(rows)

            log_path.write_text("INFO Novelty scale: 2.500\n")

            summary = plotter.summarise_prediction_table(pred_path)

            self.assertEqual(summary.n_all_rows, 5)
            self.assertEqual(summary.n_strict_rows, 4)
            self.assertEqual(summary.n_strict_expected, 2)
            self.assertEqual(summary.n_strict_expected_tp, 1)
            self.assertEqual(summary.n_strict_expected_fp, 1)
            self.assertEqual(summary.n_strict_not_detected_tp, 1)
            self.assertEqual(summary.n_neighbour_expected_overpromoted, 1)
            self.assertAlmostEqual(summary.novelty_scale, 2.5)
            self.assertAlmostEqual(summary.expected_precision, 0.5)
            self.assertAlmostEqual(summary.expected_recall, 0.5)

    def test_extract_novelty_scale_from_settings(self) -> None:
        """The novelty scale can be read from settings when no log exists."""
        with tempfile.TemporaryDirectory() as temp_dir:
            out_dir = Path(temp_dir) / "run" / "outputs"
            out_dir.mkdir(parents=True)
            pred_path = out_dir / "external_zymo_predictions.tsv.gz"
            settings_path = out_dir / "run_submission_settings.tsv"
            pred_path.write_bytes(b"")
            settings_path.write_text("setting\tvalue\nnovelty_scale\t3.0\n")

            scale = plotter.extract_novelty_scale(pred_path)
            self.assertAlmostEqual(scale, 3.0)


if __name__ == "__main__":
    unittest.main()
