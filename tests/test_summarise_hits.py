"""Tests for hit summarisation."""

import unittest

from kmersutra.screen_reads import KmerHit
from kmersutra.summarise_hits import summarise_sample_species_evidence, summarise_species_hits


class TestSummariseHits(unittest.TestCase):
    """Tests for hit summary functions."""

    def test_summarise_species_hits(self) -> None:
        """Hit summarisation should count unique k-mers and sequences."""
        hits = [
            KmerHit(sample_id="s1", sequence_id="r1", sequence_type="read", k=5, query_position=0, matched_kmer="AAAAA", query_kmer="AAAAA", mismatches=0, panel_type="species_unique", species_name="Alpha", clade="Demo"),
            KmerHit(sample_id="s1", sequence_id="r2", sequence_type="read", k=5, query_position=1, matched_kmer="AAAAC", query_kmer="AAAAC", mismatches=0, panel_type="species_unique", species_name="Alpha", clade="Demo"),
        ]
        summary = summarise_species_hits(hits=hits)
        self.assertEqual(summary[0]["n_unique_kmers"], 2)
        self.assertEqual(summary[0]["n_positive_sequences"], 2)
        evidence = summarise_sample_species_evidence(species_summary=summary)
        self.assertEqual(evidence[0]["species_name"], "Alpha")
