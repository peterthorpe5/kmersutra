"""Tests for k-mer helpers."""

import unittest

from kmersutra.kmers import (
    canonical_kmer,
    hamming_distance,
    is_valid_kmer,
    iter_kmers,
    reverse_complement,
)


class TestKmers(unittest.TestCase):
    """Tests for k-mer functions."""

    def test_reverse_complement(self) -> None:
        """Reverse complement should be correct."""
        self.assertEqual(reverse_complement("ACGTTA"), "TAACGT")

    def test_canonical_kmer(self) -> None:
        """Canonical k-mer should be lexicographically minimal strand."""
        self.assertEqual(canonical_kmer("TTTT"), "AAAA")

    def test_iter_kmers_skips_ambiguous(self) -> None:
        """Ambiguous k-mers should be skipped by default."""
        self.assertEqual(list(iter_kmers(sequence="ACNGT", k=3)), [])

    def test_iter_kmers_positions(self) -> None:
        """K-mer iteration should preserve positions."""
        observed = list(iter_kmers(sequence="AACCG", k=3, canonical=False))
        self.assertEqual(observed, [(0, "AAC"), (1, "ACC"), (2, "CCG")])

    def test_iter_kmers_rejects_non_positive_k(self) -> None:
        """K-mer iteration should reject non-positive k values."""
        with self.assertRaises(ValueError):
            list(iter_kmers(sequence="ACGT", k=0))

    def test_is_valid_kmer(self) -> None:
        """K-mer validator should reject ambiguous and empty k-mers."""
        self.assertTrue(is_valid_kmer("ACGT"))
        self.assertFalse(is_valid_kmer("ACGN"))
        self.assertFalse(is_valid_kmer(""))

    def test_hamming_distance(self) -> None:
        """Hamming distance should count mismatches."""
        self.assertEqual(hamming_distance(left="AAAA", right="AATA"), 1)

    def test_hamming_distance_rejects_unequal_lengths(self) -> None:
        """Hamming distance should require equal-length strings."""
        with self.assertRaises(ValueError):
            hamming_distance(left="AAAA", right="AAA")
