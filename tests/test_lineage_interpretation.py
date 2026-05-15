"""Tests for unresolved and possible-novel lineage interpretation."""

from __future__ import annotations

import unittest

from kmersutra.lineage_interpretation import (
    calculate_taxonomic_confidence_score,
    interpret_lineage_evidence,
)


class TestLineageInterpretation(unittest.TestCase):
    """Tests for sample-level lineage interpretation."""

    def test_taxonomic_confidence_score_is_bounded(self) -> None:
        """Taxonomic confidence scores should be bounded between zero and one."""
        score = calculate_taxonomic_confidence_score(
            n_unique_kmers=500,
            n_positive_sequences=100,
            n_k_values_positive=4,
            best_k=201,
        )
        self.assertEqual(score, 1.0)

    def test_no_evidence_reports_no_supported_signal(self) -> None:
        """Zero-evidence samples should report no supported signal."""
        records = interpret_lineage_evidence(
            species_calls=[
                {
                    "sample_id": "s1",
                    "species_name": "Alpha",
                    "call": "not_detected",
                    "n_unique_kmers": 0,
                    "n_positive_sequences": 0,
                    "best_k": 0,
                }
            ],
            taxonomic_evidence=[],
        )
        self.assertEqual(records[0]["lineage_call"], "no_supported_signal")
        self.assertFalse(records[0]["possible_novel_lineage"])

    def test_single_reportable_species_is_not_novel(self) -> None:
        """A confident species call should remain a species-level report."""
        records = interpret_lineage_evidence(
            species_calls=[
                {
                    "sample_id": "s1",
                    "species_name": "Alpha",
                    "call": "present_high_confidence",
                    "n_unique_kmers": 100,
                    "n_positive_sequences": 12,
                    "best_k": 101,
                }
            ],
            taxonomic_evidence=[
                {
                    "sample_id": "s1",
                    "evidence_rank": "genus",
                    "evidence_name": "Demo",
                    "evidence_taxid": "123",
                    "n_unique_kmers": 200,
                    "n_positive_sequences": 20,
                    "n_k_values_positive": 2,
                    "best_k": 101,
                }
            ],
        )
        self.assertEqual(records[0]["lineage_call"], "species_detected")
        self.assertEqual(records[0]["report_rank"], "species")
        self.assertEqual(records[0]["report_name"], "Alpha")
        self.assertFalse(records[0]["possible_novel_lineage"])

    def test_multiple_reportable_species_are_mixed_not_novel(self) -> None:
        """True mixed-species evidence should not be flagged as novel."""
        records = interpret_lineage_evidence(
            species_calls=[
                {
                    "sample_id": "s1",
                    "species_name": "Alpha",
                    "call": "present_in_mixed_sample",
                    "n_unique_kmers": 100,
                    "n_positive_sequences": 12,
                    "best_k": 101,
                },
                {
                    "sample_id": "s1",
                    "species_name": "Beta",
                    "call": "present_in_mixed_sample",
                    "n_unique_kmers": 80,
                    "n_positive_sequences": 10,
                    "best_k": 101,
                },
            ],
            taxonomic_evidence=[],
        )
        self.assertEqual(records[0]["lineage_call"], "mixed_species_detected")
        self.assertEqual(records[0]["n_reportable_species"], 2)
        self.assertFalse(records[0]["possible_novel_lineage"])

    def test_strong_genus_signal_with_weak_neighbours_flags_possible_novelty(self) -> None:
        """Strong genus signal plus weak species-neighbour evidence should flag novelty."""
        records = interpret_lineage_evidence(
            species_calls=[
                {
                    "sample_id": "s1",
                    "species_name": "Alpha",
                    "call": "observed_below_threshold",
                    "n_unique_kmers": 30,
                    "n_positive_sequences": 6,
                    "best_k": 101,
                },
                {
                    "sample_id": "s1",
                    "species_name": "Beta",
                    "call": "observed_below_threshold",
                    "n_unique_kmers": 25,
                    "n_positive_sequences": 6,
                    "best_k": 101,
                },
            ],
            taxonomic_evidence=[
                {
                    "sample_id": "s1",
                    "evidence_rank": "genus",
                    "evidence_name": "Demo",
                    "evidence_taxid": "123",
                    "n_unique_kmers": 80,
                    "n_positive_sequences": 14,
                    "n_k_values_positive": 2,
                    "best_k": 101,
                }
            ],
        )
        self.assertEqual(
            records[0]["lineage_call"],
            "possible_novel_or_unsampled_lineage",
        )
        self.assertTrue(records[0]["possible_novel_lineage"])
        self.assertEqual(records[0]["report_rank"], "genus")
        self.assertEqual(records[0]["report_name"], "Demo")

    def test_strong_genus_signal_without_multiple_neighbours_is_unresolved(self) -> None:
        """Strong taxonomic evidence alone should be unresolved rather than novel."""
        records = interpret_lineage_evidence(
            species_calls=[
                {
                    "sample_id": "s1",
                    "species_name": "Alpha",
                    "call": "observed_below_threshold",
                    "n_unique_kmers": 30,
                    "n_positive_sequences": 6,
                    "best_k": 101,
                }
            ],
            taxonomic_evidence=[
                {
                    "sample_id": "s1",
                    "evidence_rank": "genus",
                    "evidence_name": "Demo",
                    "evidence_taxid": "123",
                    "n_unique_kmers": 80,
                    "n_positive_sequences": 14,
                    "n_k_values_positive": 2,
                    "best_k": 101,
                }
            ],
        )
        self.assertEqual(records[0]["lineage_call"], "unresolved_taxonomic_signal")
        self.assertFalse(records[0]["possible_novel_lineage"])

    def test_weak_neighbour_evidence_below_taxonomic_threshold_stays_weak(self) -> None:
        """Weak species evidence should not become a novel-lineage call by itself."""
        records = interpret_lineage_evidence(
            species_calls=[
                {
                    "sample_id": "s1",
                    "species_name": "Alpha",
                    "call": "observed_below_threshold",
                    "n_unique_kmers": 3,
                    "n_positive_sequences": 1,
                    "best_k": 77,
                }
            ],
            taxonomic_evidence=[
                {
                    "sample_id": "s1",
                    "evidence_rank": "genus",
                    "evidence_name": "Demo",
                    "evidence_taxid": "123",
                    "n_unique_kmers": 4,
                    "n_positive_sequences": 1,
                    "n_k_values_positive": 1,
                    "best_k": 77,
                }
            ],
        )
        self.assertEqual(
            records[0]["lineage_call"],
            "weak_unresolved_neighbour_signal",
        )
        self.assertFalse(records[0]["possible_novel_lineage"])

    def test_species_rank_taxonomic_evidence_is_not_used_for_unresolved_report(self) -> None:
        """Weak species-level rows should not be promoted to unresolved reports."""
        records = interpret_lineage_evidence(
            species_calls=[
                {
                    "sample_id": "s1",
                    "species_name": "Alpha",
                    "call": "observed_below_threshold",
                    "n_unique_kmers": 30,
                    "n_positive_sequences": 6,
                    "best_k": 101,
                },
                {
                    "sample_id": "s1",
                    "species_name": "Beta",
                    "call": "observed_below_threshold",
                    "n_unique_kmers": 24,
                    "n_positive_sequences": 5,
                    "best_k": 101,
                },
            ],
            taxonomic_evidence=[
                {
                    "sample_id": "s1",
                    "evidence_rank": "species",
                    "evidence_name": "Alpha",
                    "evidence_taxid": "1",
                    "n_unique_kmers": 1000,
                    "n_positive_sequences": 50,
                    "n_k_values_positive": 2,
                    "best_k": 101,
                }
            ],
        )
        self.assertEqual(
            records[0]["lineage_call"],
            "weak_unresolved_neighbour_signal",
        )
        self.assertEqual(records[0]["report_rank"], "unresolved")

    def test_family_level_evidence_can_be_reported_when_genus_absent(self) -> None:
        """Family-level evidence can support unresolved broader lineage calls."""
        records = interpret_lineage_evidence(
            species_calls=[],
            taxonomic_evidence=[
                {
                    "sample_id": "s1",
                    "evidence_rank": "family",
                    "evidence_name": "Demoaceae",
                    "evidence_taxid": "555",
                    "n_unique_kmers": 50,
                    "n_positive_sequences": 8,
                    "n_k_values_positive": 1,
                    "best_k": 77,
                }
            ],
        )
        self.assertEqual(records[0]["lineage_call"], "unresolved_taxonomic_signal")
        self.assertEqual(records[0]["report_rank"], "family")
        self.assertEqual(records[0]["report_name"], "Demoaceae")

    def test_best_species_context_reports_margin_and_ratio(self) -> None:
        """Lineage output should include best and second-best species context."""
        records = interpret_lineage_evidence(
            species_calls=[
                {
                    "sample_id": "s1",
                    "species_name": "Alpha",
                    "call": "observed_below_threshold",
                    "n_unique_kmers": 60,
                    "n_positive_sequences": 6,
                    "best_k": 101,
                },
                {
                    "sample_id": "s1",
                    "species_name": "Beta",
                    "call": "observed_below_threshold",
                    "n_unique_kmers": 30,
                    "n_positive_sequences": 5,
                    "best_k": 101,
                },
            ],
            taxonomic_evidence=[],
        )
        self.assertEqual(records[0]["best_species"], "Alpha")
        self.assertEqual(records[0]["second_species"], "Beta")
        self.assertEqual(records[0]["best_species_margin"], 30)
        self.assertEqual(records[0]["best_species_ratio"], 2.0)


if __name__ == "__main__":
    unittest.main()


class TestLineageAwareNeighbourEvidence(unittest.TestCase):
    """Tests for neighbour-lineage evidence interpretation."""

    def test_neighbour_lineage_evidence_can_support_possible_novelty(self) -> None:
        """Strong genus evidence plus demoted neighbours should flag novelty."""
        records = interpret_lineage_evidence(
            species_calls=[
                {
                    "sample_id": "s1",
                    "species_name": "Near A",
                    "call": "neighbour_lineage_evidence",
                    "n_unique_kmers": 30,
                    "n_positive_sequences": 6,
                    "best_k": 101,
                },
                {
                    "sample_id": "s1",
                    "species_name": "Near B",
                    "call": "neighbour_lineage_evidence",
                    "n_unique_kmers": 28,
                    "n_positive_sequences": 6,
                    "best_k": 101,
                },
            ],
            taxonomic_evidence=[
                {
                    "sample_id": "s1",
                    "evidence_rank": "genus",
                    "evidence_name": "Plasmodium",
                    "evidence_taxid": "5820",
                    "n_unique_kmers": 90,
                    "n_positive_sequences": 20,
                    "n_k_values_positive": 2,
                    "best_k": 101,
                }
            ],
        )
        self.assertEqual(
            records[0]["lineage_call"],
            "possible_novel_or_unsampled_lineage",
        )
        self.assertEqual(records[0]["n_neighbour_lineage_species"], 2)
        self.assertTrue(records[0]["possible_novel_lineage"])

    def test_neighbour_lineage_evidence_is_not_reportable_species(self) -> None:
        """Demoted neighbours should not count as reportable species."""
        records = interpret_lineage_evidence(
            species_calls=[
                {
                    "sample_id": "s1",
                    "species_name": "Near A",
                    "call": "neighbour_lineage_evidence",
                    "n_unique_kmers": 30,
                    "n_positive_sequences": 6,
                    "best_k": 101,
                }
            ],
            taxonomic_evidence=[],
        )
        self.assertEqual(records[0]["n_reportable_species"], 0)
        self.assertEqual(records[0]["lineage_call"], "weak_unresolved_neighbour_signal")
