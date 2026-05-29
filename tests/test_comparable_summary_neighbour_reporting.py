"""Tests for expected-genus neighbour reporting in comparable summaries."""

from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path

import pandas as pd


MODULE_PATH = Path(__file__).resolve().parents[1] / "summarise_kmersutra_comparable_benchmark.py"
SPEC = importlib.util.spec_from_file_location("kmersutra_comparable_summary", MODULE_PATH)
SUMMARY = importlib.util.module_from_spec(SPEC)
assert SPEC is not None and SPEC.loader is not None
sys.modules["kmersutra_comparable_summary"] = SUMMARY
SPEC.loader.exec_module(SUMMARY)


class TestComparableSummaryNeighbourReporting(unittest.TestCase):
    """Test summary-level demotion of expected-genus neighbours."""

    def test_add_call_metadata_demotes_positive_expected_genus_neighbour(self) -> None:
        """A same-genus non-target in a positive sample should not be off-target."""
        calls = pd.DataFrame(
            [
                {
                    "sample_id": "s1",
                    "species_name": "Plasmodium vivax",
                    "call": "present_high_confidence",
                    "n_unique_kmers": "100",
                },
                {
                    "sample_id": "s1",
                    "species_name": "Plasmodium simium",
                    "call": "present_high_confidence",
                    "n_unique_kmers": "80",
                },
                {
                    "sample_id": "s1",
                    "species_name": "Toxoplasma gondii",
                    "call": "present_high_confidence",
                    "n_unique_kmers": "20",
                },
            ]
        )
        row = {
            "sample_id": "s1",
            "benchmark_family": "single_genome",
            "panel": "panel1",
            "replicate": "1",
            "spike_n": "1000",
            "input_fastq": "/tmp/s1.fastq.gz",
            "source_run_dir": "/tmp",
            "source_relative_dir": "s1",
        }
        observed = SUMMARY.add_call_metadata(
            calls_df=calls,
            row=row,
            out_root=Path("/tmp/out"),
            expected_targets=["Plasmodium vivax"],
            positive_calls={"present_high_confidence"},
            background_candidate_calls=set(),
            background_candidate_taxa=set(),
            demote_expected_genus_neighbours=True,
        )
        by_label = {record["report_label"]: record for record in observed.to_dict("records")}
        self.assertTrue(by_label["Plasmodium simium"]["is_positive_neighbour_lineage"])
        self.assertFalse(by_label["Plasmodium simium"]["is_positive_off_target"])
        self.assertTrue(by_label["Toxoplasma gondii"]["is_positive_off_target"])

    def test_negative_sample_does_not_demote_same_genus_signal(self) -> None:
        """Unexpected same-genus signal in a no-spike sample remains off-target."""
        calls = pd.DataFrame(
            [
                {
                    "sample_id": "s1",
                    "species_name": "Plasmodium simium",
                    "call": "present_high_confidence",
                    "n_unique_kmers": "80",
                },
            ]
        )
        row = {
            "sample_id": "s1",
            "benchmark_family": "single_genome",
            "panel": "panel1",
            "replicate": "1",
            "spike_n": "0",
            "input_fastq": "/tmp/s1.fastq.gz",
            "source_run_dir": "/tmp",
            "source_relative_dir": "s1",
        }
        observed = SUMMARY.add_call_metadata(
            calls_df=calls,
            row=row,
            out_root=Path("/tmp/out"),
            expected_targets=["Plasmodium vivax"],
            positive_calls={"present_high_confidence"},
            demote_expected_genus_neighbours=True,
        )
        self.assertFalse(bool(observed.iloc[0]["is_positive_neighbour_lineage"]))
        self.assertTrue(bool(observed.iloc[0]["is_positive_off_target"]))

    def test_sample_summary_counts_neighbour_lineage_separately(self) -> None:
        """Sample summaries should include neighbour-lineage labels and counts."""
        status = pd.DataFrame(
            [
                {
                    "sample_id": "s1",
                    "benchmark_family": "single_genome",
                    "panel": "panel1",
                    "replicate": 1,
                    "spike_n": 1000,
                    "spike_n_per_genome": 1000,
                    "n_expected_genomes": 1,
                    "total_spike_n": 1000,
                    "is_shuffled_control": False,
                    "is_negative": False,
                    "expected_targets": "Plasmodium vivax",
                    "input_fastq": "/tmp/s1.fastq.gz",
                    "source_run_dir": "/tmp",
                    "source_relative_dir": "s1",
                    "kmersutra_sample_out": "/tmp/out/s1",
                    "screen_status": "ok",
                }
            ]
        )
        calls = pd.DataFrame(
            [
                {
                    "sample_id": "s1",
                    "report_label": "Plasmodium vivax",
                    "is_positive_call": True,
                    "is_species_level": True,
                    "is_positive_expected": True,
                    "is_positive_off_target": False,
                    "is_positive_plasmodium_off_target": False,
                    "is_positive_background_candidate": False,
                    "is_positive_neighbour_lineage": False,
                    "is_positive_off_target_raw": False,
                },
                {
                    "sample_id": "s1",
                    "report_label": "Plasmodium simium",
                    "is_positive_call": True,
                    "is_species_level": True,
                    "is_positive_expected": False,
                    "is_positive_off_target": False,
                    "is_positive_plasmodium_off_target": False,
                    "is_positive_background_candidate": False,
                    "is_positive_neighbour_lineage": True,
                    "is_positive_off_target_raw": True,
                },
            ]
        )
        summary = SUMMARY.build_sample_summary(status_df=status, calls_long=calls)
        self.assertEqual(int(summary.iloc[0]["n_off_target_species"]), 0)
        self.assertEqual(int(summary.iloc[0]["n_neighbour_lineage_species"]), 1)
        self.assertEqual(summary.iloc[0]["neighbour_lineage_labels"], "Plasmodium simium")


if __name__ == "__main__":
    unittest.main()
