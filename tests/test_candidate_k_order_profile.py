"""Tests for raw-ONT candidate k-order profile handling."""

from __future__ import annotations

import unittest

from kmersutra.cli.build_clade_kmer_panel import resolve_candidate_k_order
from kmersutra.global_candidate_evidence import order_candidate_k_values


class TestCandidateKOrder(unittest.TestCase):
    """Test candidate k order resolution and validation."""

    def test_default_profile_preserves_long_to_short(self) -> None:
        """The default profile should preserve historical long-k-first sampling."""
        self.assertEqual(
            resolve_candidate_k_order(
                marker_profile="default",
                candidate_k_order="auto",
            ),
            "long_to_short",
        )

    def test_raw_ont_profile_uses_short_to_long(self) -> None:
        """The raw ONT profile should prefer shorter markers during sampling."""
        self.assertEqual(
            resolve_candidate_k_order(
                marker_profile="raw_ont_balanced",
                candidate_k_order="auto",
            ),
            "short_to_long",
        )

    def test_explicit_order_overrides_profile(self) -> None:
        """A user-specified order should override profile defaults."""
        self.assertEqual(
            resolve_candidate_k_order(
                marker_profile="raw_ont_balanced",
                candidate_k_order="long_to_short",
            ),
            "long_to_short",
        )

    def test_order_candidate_k_values(self) -> None:
        """Candidate k values should be ordered deterministically."""
        values = [101, 51, 151, 77]
        self.assertEqual(
            order_candidate_k_values(
                k_values=values,
                candidate_k_order="long_to_short",
            ),
            [151, 101, 77, 51],
        )
        self.assertEqual(
            order_candidate_k_values(
                k_values=values,
                candidate_k_order="short_to_long",
            ),
            [51, 77, 101, 151],
        )
        self.assertEqual(
            order_candidate_k_values(k_values=values, candidate_k_order="input"),
            values,
        )

    def test_invalid_order_raises(self) -> None:
        """Unsupported candidate k orders should fail fast."""
        with self.assertRaises(ValueError):
            order_candidate_k_values(k_values=[51, 77], candidate_k_order="bad")


if __name__ == "__main__":
    unittest.main()
