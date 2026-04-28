"""Tests for k-mer helpers."""

import unittest

from kmersutra.kmers import canonical_kmer, hamming_distance, iter_kmers, reverse_complement


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

    def test_hamming_distance(self) -> None:
        """Hamming distance should count mismatches."""
        self.assertEqual(hamming_distance(left="AAAA", right="AATA"), 1)
