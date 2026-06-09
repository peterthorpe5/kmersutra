"""Tests for clean Zymo truth labelling."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from kmersutra.table_io import read_records_table, write_records_table
from kmersutra.zymo_truth import (
    classify_zymo_truth_category,
    read_reference_label_map,
    write_zymo_ai_feature_table,
)


REFERENCE_MAP_ROWS = [
    {
        "reference_label": "Escherichia_coli__NRRL_B-1109",
        "species_name": "Escherichia coli",
        "role": "target_species",
    },
    {
        "reference_label": "Escherichia_coli__ALT_1",
        "species_name": "Escherichia coli",
        "role": "same_species_competitor",
    },
    {
        "reference_label": "Shigella_flexneri__ALT_1",
        "species_name": "Shigella flexneri",
        "role": "near_neighbour",
    },
]


def row_for(
    *,
    species_name: str,
    call: str,
    n_unique_kmers: str = "10",
) -> dict[str, object]:
    """Build a compact Zymo-like call row."""
    return {
        "sample_id": "ERR5396170",
        "species_name": species_name,
        "clade": species_name if "__" not in species_name else species_name.split("__", 1)[0],
        "call": call,
        "n_hits": "20" if n_unique_kmers != "0" else "0",
        "n_unique_kmers": n_unique_kmers,
        "n_positive_sequences": "5" if n_unique_kmers != "0" else "0",
        "n_k_values_positive": "2" if n_unique_kmers != "0" else "0",
        "best_k": "101" if n_unique_kmers != "0" else "0",
        "n_exact_hits": "20" if n_unique_kmers != "0" else "0",
        "n_fuzzy_hits": "0",
        "conflicting_unique_kmers": "0",
        "reportable_conflicting_unique_kmers": "0",
        "mixed_species_support_fraction": "0.1",
        "confidence_score": "0.5",
        "signal_confidence_score": "0.8",
    }


class TestZymoTruth(unittest.TestCase):
    """Test Zymo fine-category and coarse ML labels."""

    def test_classify_expected_reference_species_and_neighbour(self) -> None:
        """Rows should be separated into fine truth categories."""
        roles = {
            "Escherichia_coli__NRRL_B-1109": "target_species",
            "Escherichia_coli__ALT_1": "same_species_competitor",
            "Shigella_flexneri__ALT_1": "near_neighbour",
        }
        expected_labels = {"Escherichia_coli__NRRL_B-1109"}
        expected_species = {"Escherichia coli"}

        expected_reference = classify_zymo_truth_category(
            record=row_for(species_name="Escherichia_coli__NRRL_B-1109", call="present"),
            reference_label_roles=roles,
            expected_reference_labels=expected_labels,
            expected_species_names=expected_species,
        )
        same_species = classify_zymo_truth_category(
            record=row_for(species_name="Escherichia_coli__ALT_1", call="observed_below_threshold"),
            reference_label_roles=roles,
            expected_reference_labels=expected_labels,
            expected_species_names=expected_species,
        )
        near_neighbour = classify_zymo_truth_category(
            record=row_for(species_name="Shigella_flexneri__ALT_1", call="observed_below_threshold"),
            reference_label_roles=roles,
            expected_reference_labels=expected_labels,
            expected_species_names=expected_species,
        )
        off_target = classify_zymo_truth_category(
            record=row_for(species_name="Shigella_flexneri__ALT_1", call="present"),
            reference_label_roles=roles,
            expected_reference_labels=expected_labels,
            expected_species_names=expected_species,
        )
        absent = classify_zymo_truth_category(
            record=row_for(species_name="Shigella_flexneri__ALT_1", call="not_detected", n_unique_kmers="0"),
            reference_label_roles=roles,
            expected_reference_labels=expected_labels,
            expected_species_names=expected_species,
        )

        self.assertEqual(expected_reference["zymo_truth_category"], "expected_reference_label")
        self.assertEqual(same_species["zymo_truth_category"], "same_species_compatible_reference")
        self.assertEqual(near_neighbour["zymo_truth_category"], "near_neighbour_evidence")
        self.assertEqual(off_target["zymo_truth_category"], "true_off_target_reportable")
        self.assertEqual(absent["zymo_truth_category"], "not_detected")

    def test_write_zymo_ai_feature_table_writes_counts(self) -> None:
        """Zymo table writer should emit feature and count tables."""
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            ref_map = root / "reference_label_map.tsv"
            calls = root / "species_detection_calls.tsv"
            features = root / "features.tsv.gz"
            categories = root / "categories.tsv"
            coarse = root / "coarse.tsv"
            write_records_table(
                records=REFERENCE_MAP_ROWS,
                output_path=ref_map,
                fieldnames=list(REFERENCE_MAP_ROWS[0].keys()),
            )
            calls_rows = [
                row_for(species_name="Escherichia_coli__NRRL_B-1109", call="present"),
                row_for(species_name="Escherichia_coli__ALT_1", call="observed_below_threshold"),
                row_for(species_name="Shigella_flexneri__ALT_1", call="observed_below_threshold"),
            ]
            write_records_table(
                records=calls_rows,
                output_path=calls,
                fieldnames=list(calls_rows[0].keys()),
            )
            written = write_zymo_ai_feature_table(
                calls_table=calls,
                output_table=features,
                reference_label_map=ref_map,
                category_counts_table=categories,
                coarse_label_counts_table=coarse,
            )
            observed = read_records_table(input_path=features)
            category_counts = read_records_table(input_path=categories)
        self.assertEqual(len(written), 3)
        self.assertIn("zymo_truth_category", observed[0])
        self.assertTrue(category_counts)


if __name__ == "__main__":
    unittest.main()
