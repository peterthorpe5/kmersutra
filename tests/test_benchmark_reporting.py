"""Tests for benchmark-aware KmerSutra reporting helpers."""

from __future__ import annotations

import unittest

from kmersutra.benchmark_reporting import (
    expected_genera_from_targets,
    extract_genus,
    is_expected_genus_neighbour,
    normalise_taxon_name,
    normalise_taxa,
    reporting_layer_for_call,
)


class TestBenchmarkReporting(unittest.TestCase):
    """Test benchmark-aware report-layer decisions."""

    def test_normalise_taxon_name_handles_case_underscores_and_spaces(self) -> None:
        """Taxon names should be normalised before comparison."""
        self.assertEqual(
            normalise_taxon_name(value="  Plasmodium_simium  "),
            "plasmodium simium",
        )

    def test_extract_genus_returns_first_token(self) -> None:
        """Genus extraction should return the first normalised token."""
        self.assertEqual(extract_genus(taxon_name="Plasmodium vivax"), "plasmodium")
        self.assertEqual(extract_genus(taxon_name=""), "")

    def test_normalise_taxa_drops_empty_values(self) -> None:
        """A taxon set should not retain blank labels."""
        self.assertEqual(
            normalise_taxa(taxa=["Hammondia hammondi", "", None]),
            {"hammondia hammondi"},
        )

    def test_expected_genera_from_targets_handles_multiple_targets(self) -> None:
        """Expected target genera should be inferred from all target labels."""
        self.assertEqual(
            expected_genera_from_targets(
                expected_targets=["Plasmodium vivax", "Babesia microti"],
            ),
            {"plasmodium", "babesia"},
        )

    def test_positive_same_genus_non_target_is_neighbour(self) -> None:
        """A non-target same-genus species in a positive sample is neighbour evidence."""
        self.assertTrue(
            is_expected_genus_neighbour(
                report_label="Plasmodium simium",
                expected_targets=["Plasmodium vivax"],
                is_expected_target=False,
                is_negative_sample=False,
                is_background_candidate=False,
            )
        )

    def test_negative_same_genus_non_target_is_not_demoted(self) -> None:
        """No-spike negatives should not demote unexpected same-genus signal."""
        self.assertFalse(
            is_expected_genus_neighbour(
                report_label="Plasmodium simium",
                expected_targets=["Plasmodium vivax"],
                is_expected_target=False,
                is_negative_sample=True,
                is_background_candidate=False,
            )
        )

    def test_background_candidate_is_not_neighbour(self) -> None:
        """Background candidates should keep the background reporting layer."""
        self.assertFalse(
            is_expected_genus_neighbour(
                report_label="Plasmodium simium",
                expected_targets=["Plasmodium vivax"],
                is_expected_target=False,
                is_negative_sample=False,
                is_background_candidate=True,
            )
        )

    def test_reporting_layer_priority(self) -> None:
        """Report layer assignment should follow expected priority order."""
        self.assertEqual(
            reporting_layer_for_call(
                is_positive_call=True,
                is_species_level=True,
                is_expected_target=False,
                is_background_candidate=False,
                is_expected_genus_neighbour_call=True,
            ),
            "expected_genus_neighbour",
        )
        self.assertEqual(
            reporting_layer_for_call(
                is_positive_call=True,
                is_species_level=True,
                is_expected_target=False,
                is_background_candidate=False,
                is_expected_genus_neighbour_call=False,
            ),
            "reportable_off_target",
        )


if __name__ == "__main__":
    unittest.main()
