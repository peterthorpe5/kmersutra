"""Tests for the registered ATCC MSA-1003 benchmark specification."""

from __future__ import annotations

import json
import unittest
from collections import Counter
from pathlib import Path

from kmersutra.benchmark_workflow import validate_locked_atcc_config
from kmersutra.mock_community import load_truth_manifest


BENCHMARK_ROOT = (
    Path(__file__).resolve().parents[1] / "benchmarks" / "atcc_msa1003_hifi"
)


class TestAtccBenchmarkSpecification(unittest.TestCase):
    """Protect the frozen dataset, composition and primary settings."""

    def test_truth_manifest_has_twenty_expected_organisms_and_four_tiers(self) -> None:
        """The ATCC truth manifest should retain the registered composition."""
        records = load_truth_manifest(
            manifest_path=BENCHMARK_ROOT / "truth_manifest.tsv"
        )
        self.assertEqual(len(records), 20)
        self.assertTrue(all(record.expected for record in records))
        self.assertAlmostEqual(
            sum(record.expected_abundance_fraction for record in records),
            1.0,
        )
        self.assertEqual(
            Counter(record.expected_abundance_fraction for record in records),
            Counter({0.18: 5, 0.018: 5, 0.0018: 5, 0.0002: 5}),
        )
        self.assertEqual(
            len(
                {
                    accession
                    for record in records
                    for accession in record.truth_sequence_accessions
                }
            ),
            20,
        )

    def test_example_config_obeys_locked_primary_design(self) -> None:
        """The distributed example must not drift from the registered design."""
        with (BENCHMARK_ROOT / "config.example.json").open(
            "r",
            encoding="utf-8",
        ) as handle:
            config = json.load(handle)
        validate_locked_atcc_config(config=config)

    def test_locked_specification_matches_example_config(self) -> None:
        """The human-readable registration JSON should match executable config."""
        with (BENCHMARK_ROOT / "config.example.json").open(
            "r",
            encoding="utf-8",
        ) as handle:
            config = json.load(handle)
        with (BENCHMARK_ROOT / "locked_analysis_specification.json").open(
            "r",
            encoding="utf-8",
        ) as handle:
            specification = json.load(handle)
        self.assertEqual(
            config["benchmark_id"],
            specification["benchmark_id"],
        )
        self.assertEqual(config["dataset"], {
            **specification["dataset"],
            "input_reads": "",
        })
        self.assertEqual(config["screen"]["k_values"], [51, 77, 101, 151])
        self.assertEqual(config["ai"]["novelty_scale"], 2.9)
        self.assertFalse(config["allow_primary_threshold_tuning"])


if __name__ == "__main__":
    unittest.main()
