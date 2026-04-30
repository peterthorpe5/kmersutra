"""Tests for detection thresholds and confidence scoring."""

import unittest

from kmersutra.thresholds import (
    call_species_presence,
    calculate_confidence_score,
    calculate_conflict_ratio,
)


class TestThresholds(unittest.TestCase):
    """Tests for detection-call logic."""

    def test_calculate_conflict_ratio(self) -> None:
        """Conflict ratio should use target plus conflict evidence."""
        self.assertEqual(calculate_conflict_ratio(target_unique_kmers=8, conflicting_unique_kmers=2), 0.2)

    def test_calculate_conflict_ratio_handles_zero_total(self) -> None:
        """Conflict ratio should be zero when no evidence is present."""
        self.assertEqual(calculate_conflict_ratio(target_unique_kmers=0, conflicting_unique_kmers=0), 0.0)

    def test_confidence_score_is_bounded(self) -> None:
        """Confidence score should be bounded between zero and one."""
        score = calculate_confidence_score(
            n_unique_kmers=100,
            n_positive_sequences=100,
            n_k_values_positive=5,
            best_k=151,
            conflict_ratio=0,
        )
        self.assertEqual(score, 1.0)

    def test_confidence_score_decreases_with_conflict(self) -> None:
        """Confidence score should decrease when conflict ratio increases."""
        clean = calculate_confidence_score(
            n_unique_kmers=10,
            n_positive_sequences=5,
            n_k_values_positive=2,
            best_k=101,
            conflict_ratio=0,
        )
        conflicted = calculate_confidence_score(
            n_unique_kmers=10,
            n_positive_sequences=5,
            n_k_values_positive=2,
            best_k=101,
            conflict_ratio=0.5,
        )
        self.assertLess(conflicted, clean)

    def test_call_species_presence_high_confidence(self) -> None:
        """Sufficient clean evidence should call high confidence presence."""
        evidence = [
            {
                "sample_id": "s1",
                "species_name": "Alpha",
                "clade": "Demo",
                "n_hits": 10,
                "n_unique_kmers": 5,
                "n_positive_sequences": 3,
                "n_k_values_positive": 2,
                "best_k": 71,
                "n_exact_hits": 10,
                "n_fuzzy_hits": 0,
            }
        ]
        calls = call_species_presence(evidence_records=evidence)
        self.assertEqual(calls[0]["call"], "present_high_confidence")

    def test_call_species_presence_low_confidence(self) -> None:
        """Sub-threshold evidence should call low-confidence presence."""
        evidence = [
            {
                "sample_id": "s1",
                "species_name": "Alpha",
                "clade": "Demo",
                "n_hits": 1,
                "n_unique_kmers": 1,
                "n_positive_sequences": 1,
                "n_k_values_positive": 1,
                "best_k": 71,
                "n_exact_hits": 1,
                "n_fuzzy_hits": 0,
            }
        ]
        calls = call_species_presence(evidence_records=evidence)
        self.assertEqual(calls[0]["call"], "present_low_confidence")

    def test_call_species_presence_ambiguous(self) -> None:
        """Conflicting evidence should produce an ambiguous call."""
        evidence = [
            {
                "sample_id": "s1",
                "species_name": "Alpha",
                "clade": "Demo",
                "n_hits": 10,
                "n_unique_kmers": 5,
                "n_positive_sequences": 3,
                "n_k_values_positive": 2,
                "best_k": 71,
                "n_exact_hits": 10,
                "n_fuzzy_hits": 0,
            },
            {
                "sample_id": "s1",
                "species_name": "Beta",
                "clade": "Demo",
                "n_hits": 8,
                "n_unique_kmers": 4,
                "n_positive_sequences": 3,
                "n_k_values_positive": 2,
                "best_k": 71,
                "n_exact_hits": 8,
                "n_fuzzy_hits": 0,
            },
        ]
        calls = call_species_presence(evidence_records=evidence)
        self.assertEqual({call["call"] for call in calls}, {"ambiguous_mixed_signal"})

    def test_call_species_presence_handles_empty_evidence(self) -> None:
        """Empty evidence should produce no calls."""
        self.assertEqual(call_species_presence(evidence_records=[]), [])
