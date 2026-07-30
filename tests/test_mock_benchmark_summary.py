"""Tests for mock-community benchmark summaries."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from kmersutra.mock_benchmark_summary import (
    summarise_mock_benchmark,
    summarise_one_sample,
)
from kmersutra.mock_community import load_truth_manifest
from kmersutra.table_io import read_records_table


class TestMockBenchmarkSummary(unittest.TestCase):
    """Test raw/reportable and tier-level benchmark outcomes."""

    def write_truth(self, root: Path) -> Path:
        """Write a two-species truth manifest."""
        path = root / "truth.tsv"
        path.write_text(
            "benchmark_id\torganism_id\tspecies_name\texpected\t"
            "expected_abundance_fraction\tabundance_tier\n"
            "example\tone\tSpecies alpha\ttrue\t0.9\thigh\n"
            "example\ttwo\tSpecies beta\ttrue\t0.1\tlow\n",
            encoding="utf-8",
        )
        return path

    def call_rows(self) -> list[dict[str, str]]:
        """Return expected, below-threshold and off-target calls."""
        return [
            {
                "species_name": "Species alpha",
                "call": "present",
                "n_unique_kmers": "20",
                "n_positive_sequences": "5",
                "best_k": "151",
            },
            {
                "species_name": "Species beta",
                "call": "observed_below_threshold",
                "n_unique_kmers": "2",
                "n_positive_sequences": "1",
                "best_k": "51",
            },
            {
                "species_name": "Species gamma",
                "call": "present",
                "n_unique_kmers": "5",
                "n_positive_sequences": "2",
                "best_k": "77",
            },
        ]

    def test_one_sample_separates_raw_and_reportable_detection(self) -> None:
        """Raw evidence should not be conflated with reportable species calls."""
        with tempfile.TemporaryDirectory() as temporary:
            truth = load_truth_manifest(
                manifest_path=self.write_truth(Path(temporary))
            )
        expected, off_targets, summary = summarise_one_sample(
            task={
                "sample_id": "sample",
                "analysis_type": "full",
                "fraction": "1.0",
                "seed": "",
                "k_value": "multi_k",
            },
            truth_records=truth,
            call_rows=self.call_rows(),
        )
        self.assertEqual(len(expected), 2)
        self.assertEqual(summary["n_expected_raw_detected"], 2)
        self.assertEqual(summary["n_expected_reportable"], 1)
        self.assertEqual(summary["n_reportable_off_target_species"], 1)
        self.assertEqual(len(off_targets), 1)

    def test_summary_workflow_writes_all_tables(self) -> None:
        """The summary workflow should write sample, tier and row-level tables."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            truth_path = self.write_truth(root)
            calls_path = root / "calls.tsv"
            header = [
                "species_name",
                "call",
                "n_unique_kmers",
                "n_positive_sequences",
                "best_k",
            ]
            with calls_path.open("w", encoding="utf-8") as handle:
                handle.write("\t".join(header) + "\n")
                for row in self.call_rows():
                    handle.write("\t".join(row[column] for column in header) + "\n")
            tasks = root / "tasks.tsv"
            tasks.write_text(
                "sample_id\tanalysis_type\tfraction\tseed\tk_value\tcalls_table\n"
                f"sample\tfull\t1.0\t\tmulti_k\t{calls_path}\n",
                encoding="utf-8",
            )
            paths = summarise_mock_benchmark(
                task_manifest=tasks,
                truth_manifest=truth_path,
                output_dir=root / "summary",
            )
            sample_rows = read_records_table(input_path=paths["sample_summary"])
            tier_rows = read_records_table(input_path=paths["tier_summary"])
        self.assertEqual(len(paths), 4)
        self.assertEqual(sample_rows[0]["n_expected_raw_detected"], "2")
        self.assertEqual(len(tier_rows), 2)


if __name__ == "__main__":
    unittest.main()
