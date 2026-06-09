#!/usr/bin/env python3
"""Unit tests for common KmerSutra AI validation helpers."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from kmersutra_ai_common import (
    audit_feature_columns,
    infer_public_truth_label,
    label_counts_text,
    read_reference_label_map,
    sanitise_name,
    stable_group_hash,
    write_tsv_simple,
)


class TestCommonHelpers(unittest.TestCase):
    """Tests for common helper functions."""

    def test_stable_group_hash_is_deterministic(self) -> None:
        """Hashing should be deterministic for grouped validation."""
        self.assertEqual(
            stable_group_hash("sample_001"),
            stable_group_hash("sample_001"),
        )
        self.assertNotEqual(
            stable_group_hash("sample_001"),
            stable_group_hash("sample_002"),
        )

    def test_sanitise_name_removes_unsafe_characters(self) -> None:
        """Validation split names should be filesystem safe."""
        observed = sanitise_name("Plasmodium sp. / panel:1")
        self.assertNotIn("/", observed)
        self.assertNotIn(":", observed)
        self.assertIn("Plasmodium", observed)

    def test_label_counts_text(self) -> None:
        """Label counts should be sorted and compact."""
        records = [
            {"label": "b"},
            {"label": "a"},
            {"label": "b"},
        ]
        self.assertEqual(
            label_counts_text(records=records, label_column="label"),
            "a:1;b:2",
        )

    def test_audit_feature_columns_flags_leakage(self) -> None:
        """Leakage-prone columns should not be recommended."""
        rows = audit_feature_columns(["n_hits", "spike_n", "sample_id"])
        by_feature = {row["feature"]: row for row in rows}
        self.assertEqual(by_feature["n_hits"]["recommended_for_training"], "yes")
        self.assertEqual(by_feature["spike_n"]["recommended_for_training"], "no")
        self.assertEqual(by_feature["sample_id"]["is_leakage_risk"], "yes")

    def test_read_reference_label_map(self) -> None:
        """Expected target labels should be read from role target_species."""
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "reference_label_map.tsv"
            write_tsv_simple(
                records=[
                    {
                        "reference_label": "Ecoli_ref",
                        "species_name": "Escherichia coli",
                        "original_species_name": "Escherichia coli",
                        "role": "target_species",
                    },
                    {
                        "reference_label": "Other_ref",
                        "species_name": "Other species",
                        "original_species_name": "Other species",
                        "role": "near_neighbour",
                    },
                ],
                output_path=path,
                fieldnames=[
                    "reference_label",
                    "species_name",
                    "original_species_name",
                    "role",
                ],
            )
            expected_refs, expected_species, roles = read_reference_label_map(
                reference_label_map=path,
            )
        self.assertIn("Ecoli_ref", expected_refs)
        self.assertIn("Escherichia coli", expected_species)
        self.assertNotIn("Other_ref", expected_refs)
        self.assertEqual(roles["Other_ref"], "near_neighbour")

    def test_infer_public_truth_label(self) -> None:
        """Public truth labels should distinguish expected and off-target rows."""
        expected_refs = {"Ecoli_ref"}
        expected_species = {"Escherichia coli"}

        self.assertEqual(
            infer_public_truth_label(
                row={"reference_label": "Ecoli_ref", "call": "present_in_mixed_sample"},
                expected_reference_labels=expected_refs,
                expected_species_names=expected_species,
            ),
            "expected_target",
        )
        self.assertEqual(
            infer_public_truth_label(
                row={"species_name": "Other species", "call": "present_in_mixed_sample"},
                expected_reference_labels=expected_refs,
                expected_species_names=expected_species,
            ),
            "reportable_off_target_species",
        )
        self.assertEqual(
            infer_public_truth_label(
                row={"species_name": "Other species", "n_unique_kmers": "3"},
                expected_reference_labels=expected_refs,
                expected_species_names=expected_species,
            ),
            "observed_below_threshold",
        )
        self.assertEqual(
            infer_public_truth_label(
                row={"species_name": "Other species", "n_unique_kmers": "0"},
                expected_reference_labels=expected_refs,
                expected_species_names=expected_species,
            ),
            "not_detected",
        )


if __name__ == "__main__":
    unittest.main()
