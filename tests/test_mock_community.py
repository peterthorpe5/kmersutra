"""Tests for generic mock-community truth handling."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from kmersutra.mock_community import (
    build_mock_ai_feature_records,
    classify_mock_truth_category,
    load_truth_manifest,
    normalise_species_key,
    truth_by_species,
    write_mock_ai_feature_table,
)
from kmersutra.table_io import read_records_table


class TestMockCommunityTruth(unittest.TestCase):
    """Test truth-manifest validation and call labelling."""

    def write_truth(self, directory: Path) -> Path:
        """Write a minimal truth manifest."""
        path = directory / "truth.tsv"
        path.write_text(
            "benchmark_id\torganism_id\tspecies_name\t"
            "accepted_species_names\texpected\texpected_abundance_fraction\t"
            "abundance_tier\ttruth_sequence_accessions\n"
            "example\tone\tSpecies alpha\tOld_alpha\ttrue\t0.75\t"
            "high\tNC_000001.1\n"
            "example\ttwo\tSpecies beta\t\ttrue\t0.25\tlow\tNC_000002.1\n",
            encoding="utf-8",
        )
        return path

    def test_load_truth_manifest_validates_and_normalises(self) -> None:
        """Truth records should retain aliases and normalised accessions."""
        with tempfile.TemporaryDirectory() as temporary:
            truth = load_truth_manifest(
                manifest_path=self.write_truth(Path(temporary))
            )
        self.assertEqual(len(truth), 2)
        self.assertIn("old alpha", truth[0].accepted_species_names)
        self.assertEqual(truth[0].truth_sequence_accessions, {"NC_000001"})

    def test_truth_manifest_rejects_fraction_sum_above_one(self) -> None:
        """Expected fractions above one should fail clearly."""
        with tempfile.TemporaryDirectory() as temporary:
            path = self.write_truth(Path(temporary))
            text = path.read_text(encoding="utf-8").replace("\t0.75\t", "\t0.85\t")
            path.write_text(text, encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "exceed one"):
                load_truth_manifest(manifest_path=path)

    def test_classify_expected_reportable_species(self) -> None:
        """A reportable expected species should become expected_target."""
        with tempfile.TemporaryDirectory() as temporary:
            truth = load_truth_manifest(
                manifest_path=self.write_truth(Path(temporary))
            )
        mapping = truth_by_species(truth_records=truth)
        labelled = classify_mock_truth_category(
            record={
                "species_name": "Old_alpha",
                "call": "present",
                "n_unique_kmers": "10",
            },
            truth_mapping=mapping,
        )
        self.assertEqual(labelled["mock_truth_category"], "expected_species_reportable")
        self.assertEqual(labelled["ml_report_label"], "expected_target")

    def test_classify_unexpected_reportable_species(self) -> None:
        """A reportable unexpected species should be an off-target."""
        with tempfile.TemporaryDirectory() as temporary:
            truth = load_truth_manifest(
                manifest_path=self.write_truth(Path(temporary))
            )
        labelled = classify_mock_truth_category(
            record={
                "species_name": "Species gamma",
                "call": "present",
                "n_unique_kmers": "10",
            },
            truth_mapping=truth_by_species(truth_records=truth),
        )
        self.assertEqual(
            labelled["ml_report_label"],
            "reportable_off_target_species",
        )

    def test_build_feature_records_preserves_ai_features(self) -> None:
        """Generic labelling should retain numeric AI feature construction."""
        with tempfile.TemporaryDirectory() as temporary:
            truth = load_truth_manifest(
                manifest_path=self.write_truth(Path(temporary))
            )
        records = build_mock_ai_feature_records(
            records=[
                {
                    "species_name": "Species alpha",
                    "call": "present",
                    "n_hits": "20",
                    "n_unique_kmers": "10",
                    "n_positive_sequences": "5",
                    "n_k_values_positive": "2",
                    "best_k": "101",
                }
            ],
            truth_records=truth,
        )
        self.assertEqual(len(records), 1)
        self.assertIn("log1p_n_unique_kmers", records[0])
        self.assertEqual(records[0]["ml_report_label"], "expected_target")

    def test_write_mock_feature_table_outputs_counts(self) -> None:
        """The file workflow should write labelled rows and both count tables."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            truth_path = self.write_truth(root)
            calls_path = root / "calls.tsv"
            calls_path.write_text(
                "species_name\tcall\tn_hits\tn_unique_kmers\t"
                "n_positive_sequences\tn_k_values_positive\tbest_k\n"
                "Species alpha\tpresent\t20\t10\t5\t2\t101\n",
                encoding="utf-8",
            )
            output = root / "labelled.tsv.gz"
            category_counts = root / "category.tsv"
            coarse_counts = root / "coarse.tsv"
            write_mock_ai_feature_table(
                calls_table=calls_path,
                truth_manifest=truth_path,
                output_table=output,
                category_counts_table=category_counts,
                coarse_label_counts_table=coarse_counts,
            )
            labelled = read_records_table(input_path=output)
            categories = read_records_table(input_path=category_counts)
        self.assertEqual(len(labelled), 1)
        self.assertEqual(categories[0]["n_records"], "1")

    def test_registered_atcc_manifest_has_twenty_expected_species(self) -> None:
        """The registered ATCC manifest should sum to a complete community."""
        path = (
            Path(__file__).resolve().parents[1]
            / "benchmarks"
            / "atcc_msa1003_hifi"
            / "truth_manifest.tsv"
        )
        truth = load_truth_manifest(manifest_path=path)
        self.assertEqual(len(truth), 20)
        self.assertAlmostEqual(
            sum(record.expected_abundance_fraction for record in truth),
            1.0,
        )
        tiers = {record.abundance_tier for record in truth}
        self.assertEqual(len(tiers), 4)
        self.assertEqual(
            normalise_species_key("Propionibacterium_acnes"),
            "propionibacterium acnes",
        )


if __name__ == "__main__":
    unittest.main()
