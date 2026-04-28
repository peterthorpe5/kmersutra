"""Tests for sequence screening."""

import unittest

from kmersutra.build_panel import DiagnosticKmer
from kmersutra.fasta import SequenceRecord
from kmersutra.screen_reads import screen_sequence_for_kmers


class TestScreenReads(unittest.TestCase):
    """Tests for screening sequences against a panel."""

    def test_screen_sequence_exact_hit(self) -> None:
        """Screening should identify exact species-unique k-mer hits."""
        diagnostic = DiagnosticKmer(
            kmer="AAAAA",
            k=5,
            panel_type="species_unique",
            species_name="Alpha",
            clade="Demo",
            source_genomes="g1",
            source_contigs="c1",
            example_position=0,
        )
        panel = {5: {"AAAAA": [diagnostic]}}
        record = SequenceRecord(identifier="read1", description="read1", sequence="GGGAAAAATTT")
        hits = screen_sequence_for_kmers(
            sequence_record=record,
            panel_index=panel,
            sample_id="sample1",
            sequence_type="read",
        )
        self.assertEqual(len(hits), 1)
        self.assertEqual(hits[0].species_name, "Alpha")

    def test_screen_sequence_fuzzy_hit(self) -> None:
        """Screening should optionally identify fuzzy long-k hits."""
        diagnostic = DiagnosticKmer(
            kmer="A" * 71,
            k=71,
            panel_type="species_unique",
            species_name="Alpha",
            clade="Demo",
            source_genomes="g1",
            source_contigs="c1",
            example_position=0,
        )
        panel = {71: {"A" * 71: [diagnostic]}}
        record = SequenceRecord(identifier="read1", description="read1", sequence=("A" * 70) + "C")
        hits = screen_sequence_for_kmers(
            sequence_record=record,
            panel_index=panel,
            sample_id="sample1",
            sequence_type="read",
            max_mismatches=1,
            fuzzy_min_k=71,
        )
        self.assertEqual(len(hits), 1)
        self.assertEqual(hits[0].mismatches, 1)
