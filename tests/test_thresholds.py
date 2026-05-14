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
        """Sufficient clean evidence should call high-confidence presence."""
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

    def test_call_species_presence_identifies_true_mixed_samples(self) -> None:
        """Multiple independently supported species should be mixed-present calls."""
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
        self.assertEqual({call["call"] for call in calls}, {"present_in_mixed_sample"})
        self.assertTrue(all(float(call["confidence_score"]) > 0.5 for call in calls))

    def test_call_species_presence_can_disallow_mixed_samples(self) -> None:
        """Mixed evidence can still be treated as conflicting when requested."""
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
        calls = call_species_presence(evidence_records=evidence, allow_mixed_species=False)
        self.assertEqual({call["call"] for call in calls}, {"ambiguous_conflicting_signal"})

    def test_call_species_presence_reports_not_detected(self) -> None:
        """Explicit zero-evidence rows should be called not detected."""
        evidence = [
            {
                "sample_id": "s1",
                "species_name": "Alpha",
                "clade": "Demo",
                "n_hits": 0,
                "n_unique_kmers": 0,
                "n_positive_sequences": 0,
                "n_k_values_positive": 0,
                "best_k": 0,
                "n_exact_hits": 0,
                "n_fuzzy_hits": 0,
            }
        ]
        calls = call_species_presence(evidence_records=evidence)
        self.assertEqual(calls[0]["call"], "not_detected")

    def test_call_species_presence_handles_empty_evidence(self) -> None:
        """Empty evidence should produce no calls."""
        self.assertEqual(call_species_presence(evidence_records=[]), [])


class TestConservativeThresholds(unittest.TestCase):
    """Tests for conservative v0.15 species-call behaviour."""

    def test_call_preset_returns_conservative_settings(self) -> None:
        """Conservative preset should require multi-k and long-k evidence."""
        from kmersutra.thresholds import apply_species_call_preset

        settings = apply_species_call_preset(preset_name="conservative")
        self.assertEqual(settings["min_k_values_positive"], 2)
        self.assertEqual(settings["min_best_k"], 101)
        self.assertEqual(settings["low_evidence_call"], "observed_below_threshold")

    def test_conservative_threshold_labels_weak_signal_below_threshold(self) -> None:
        """Weak off-target-like evidence should not become low-confidence presence."""
        evidence = [
            {
                "sample_id": "s1",
                "species_name": "Weak off-target",
                "clade": "Demo",
                "n_hits": 4,
                "n_unique_kmers": 4,
                "n_positive_sequences": 2,
                "n_k_values_positive": 1,
                "best_k": 77,
                "n_exact_hits": 4,
                "n_fuzzy_hits": 0,
            }
        ]
        calls = call_species_presence(
            evidence_records=evidence,
            min_unique_kmers=20,
            min_positive_sequences=5,
            min_k_values_positive=2,
            min_best_k=101,
            min_exact_hits=20,
            min_confidence_score=0.5,
            low_evidence_call="observed_below_threshold",
        )
        self.assertEqual(calls[0]["call"], "observed_below_threshold")

    def test_conservative_threshold_accepts_strong_multi_k_signal(self) -> None:
        """Strong exact evidence across long k values should pass."""
        evidence = [
            {
                "sample_id": "s1",
                "species_name": "Alpha",
                "clade": "Demo",
                "n_hits": 120,
                "n_unique_kmers": 60,
                "n_positive_sequences": 15,
                "n_k_values_positive": 2,
                "best_k": 101,
                "n_exact_hits": 120,
                "n_fuzzy_hits": 0,
            }
        ]
        calls = call_species_presence(
            evidence_records=evidence,
            min_unique_kmers=20,
            min_positive_sequences=5,
            min_k_values_positive=2,
            min_best_k=101,
            min_exact_hits=20,
            min_confidence_score=0.5,
            low_evidence_call="observed_below_threshold",
        )
        self.assertEqual(calls[0]["call"], "present_high_confidence")

    def test_relative_margin_can_filter_second_best_neighbour(self) -> None:
        """A required margin should prevent near-tied neighbour calls."""
        evidence = [
            {
                "sample_id": "s1",
                "species_name": "Alpha",
                "clade": "Demo",
                "n_hits": 100,
                "n_unique_kmers": 60,
                "n_positive_sequences": 20,
                "n_k_values_positive": 2,
                "best_k": 101,
                "n_exact_hits": 100,
                "n_fuzzy_hits": 0,
            },
            {
                "sample_id": "s1",
                "species_name": "Close neighbour",
                "clade": "Demo",
                "n_hits": 95,
                "n_unique_kmers": 55,
                "n_positive_sequences": 18,
                "n_k_values_positive": 2,
                "best_k": 101,
                "n_exact_hits": 95,
                "n_fuzzy_hits": 0,
            },
        ]
        calls = call_species_presence(
            evidence_records=evidence,
            min_unique_kmers=20,
            min_positive_sequences=5,
            min_k_values_positive=2,
            min_best_k=101,
            min_exact_hits=20,
            low_evidence_call="observed_below_threshold",
            min_unique_kmer_margin=10,
        )
        call_map = {call["species_name"]: call["call"] for call in calls}
        self.assertEqual(call_map["Alpha"], "observed_below_threshold")
        self.assertEqual(call_map["Close neighbour"], "observed_below_threshold")
