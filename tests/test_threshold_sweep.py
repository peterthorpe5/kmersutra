"""Tests for KmerSutra threshold-sweep diagnostics."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import pandas as pd

from kmersutra.threshold_sweep import (
    annotate_threshold,
    parse_threshold_spec,
    read_calls_table,
    run_threshold_sweep,
    summarise_targets,
)


class TestThresholdSweep(unittest.TestCase):
    """Test threshold-sweep parsing and benchmark reinterpretation."""

    def test_parse_threshold_spec(self) -> None:
        """Compact threshold specifications should parse predictably."""
        observed = parse_threshold_spec(
            "final_stress:k1,u8,p1,e8,best101,conf0.50,conflict0.10"
        )
        self.assertEqual(observed.name, "final_stress")
        self.assertEqual(observed.min_k_values_positive, 1)
        self.assertEqual(observed.min_unique_kmers, 8)
        self.assertEqual(observed.min_positive_sequences, 1)
        self.assertEqual(observed.min_exact_hits, 8)
        self.assertEqual(observed.min_best_k, 101)
        self.assertAlmostEqual(observed.min_confidence_score, 0.50)
        self.assertAlmostEqual(observed.max_reportable_conflict_ratio, 0.10)

    def make_calls(self, path: Path) -> None:
        """Write a minimal long-call table for testing.

        Parameters
        ----------
        path : pathlib.Path
            Output table path.
        """
        rows = [
            {
                "sample_id": "pos1",
                "benchmark_family": "single_genome",
                "panel": "panel1",
                "spike_n": "100",
                "is_negative": "False",
                "expected_targets": "Plasmodium vivax",
                "report_label": "Plasmodium vivax",
                "is_species_level": "True",
                "is_expected_target": "True",
                "n_unique_kmers": "9",
                "n_positive_sequences": "1",
                "n_k_values_positive": "1",
                "best_k": "151",
                "n_exact_hits": "9",
                "reportable_conflict_ratio": "0.0",
                "confidence_score": "0.80",
                "is_expected_genus_neighbour": "False",
                "is_background_candidate_taxon": "False",
                "is_background_candidate_signal": "False",
            },
            {
                "sample_id": "pos1",
                "benchmark_family": "single_genome",
                "panel": "panel1",
                "spike_n": "100",
                "is_negative": "False",
                "expected_targets": "Plasmodium vivax",
                "report_label": "Plasmodium simium",
                "is_species_level": "True",
                "is_expected_target": "False",
                "n_unique_kmers": "9",
                "n_positive_sequences": "1",
                "n_k_values_positive": "1",
                "best_k": "151",
                "n_exact_hits": "9",
                "reportable_conflict_ratio": "0.0",
                "confidence_score": "0.80",
                "is_expected_genus_neighbour": "True",
                "is_background_candidate_taxon": "False",
                "is_background_candidate_signal": "False",
            },
            {
                "sample_id": "neg1",
                "benchmark_family": "single_genome",
                "panel": "panel1",
                "spike_n": "0",
                "is_negative": "True",
                "expected_targets": "Plasmodium vivax",
                "report_label": "Plasmodium vivax",
                "is_species_level": "True",
                "is_expected_target": "True",
                "n_unique_kmers": "0",
                "n_positive_sequences": "0",
                "n_k_values_positive": "0",
                "best_k": "0",
                "n_exact_hits": "0",
                "reportable_conflict_ratio": "1.0",
                "confidence_score": "0.0",
                "is_expected_genus_neighbour": "False",
                "is_background_candidate_taxon": "False",
                "is_background_candidate_signal": "False",
            },
        ]
        pd.DataFrame(rows).to_csv(path, sep="\t", index=False)

    def test_threshold_sweep_promotes_only_relaxed_evidence(self) -> None:
        """Relaxed thresholds should promote evidence that strict thresholds reject."""
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "calls.tsv"
            self.make_calls(path)
            calls = read_calls_table(path=path)
            strict = parse_threshold_spec("strict:k1,u10,p2,e10,best101,conf0.50,conflict0.10")
            relaxed = parse_threshold_spec("relaxed:k1,u8,p1,e8,best101,conf0.50,conflict0.10")
            strict_annotated = annotate_threshold(calls=calls, threshold=strict)
            relaxed_annotated = annotate_threshold(calls=calls, threshold=relaxed)
            self.assertEqual(int(strict_annotated["sweep_positive"].sum()), 0)
            self.assertEqual(int(relaxed_annotated["sweep_expected_positive"].sum()), 1)
            self.assertEqual(int(relaxed_annotated["sweep_neighbour_lineage"].sum()), 1)
            target_summary = summarise_targets(annotated=relaxed_annotated)
            self.assertEqual(float(target_summary.iloc[0]["sensitivity"]), 1.0)
            self.assertEqual(float(target_summary.iloc[0]["specificity"]), 1.0)

    def test_run_threshold_sweep_writes_expected_outputs(self) -> None:
        """The command should write TSV, Excel and HTML outputs."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            calls = root / "calls.tsv"
            out_dir = root / "sweep"
            self.make_calls(calls)
            paths = run_threshold_sweep(
                calls_table=calls,
                out_dir=out_dir,
                threshold_specs=["relaxed:k1,u8,p1,e8,best101,conf0.50,conflict0.10"],
                report_name="unit_sweep",
            )
            self.assertTrue(paths.target_performance.exists())
            self.assertTrue(paths.real_world_summary.exists())
            self.assertTrue(paths.neighbour_lineage_summary.exists())
            self.assertTrue(paths.workbook.exists())
            self.assertTrue(paths.html_report.exists())


if __name__ == "__main__":
    unittest.main()
