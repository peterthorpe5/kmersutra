"""Tests for sequence screening."""

import unittest

from kmersutra.build_panel import DiagnosticKmer
from kmersutra.fasta import SequenceRecord
from kmersutra.screen_reads import (
    iter_mismatch_neighbourhood,
    screen_records_for_species_kmers,
    screen_sequence_for_kmers,
)


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

    def test_screen_sequence_fuzzy_respects_min_k(self) -> None:
        """Fuzzy matching should not run below the configured minimum k."""
        diagnostic = DiagnosticKmer(
            kmer="A" * 10,
            k=10,
            panel_type="species_unique",
            species_name="Alpha",
            clade="Demo",
            source_genomes="g1",
            source_contigs="c1",
            example_position=0,
        )
        panel = {10: {"A" * 10: [diagnostic]}}
        record = SequenceRecord(identifier="read1", description="read1", sequence=("A" * 9) + "C")
        hits = screen_sequence_for_kmers(
            sequence_record=record,
            panel_index=panel,
            sample_id="sample1",
            sequence_type="read",
            max_mismatches=1,
            fuzzy_min_k=71,
        )
        self.assertEqual(hits, [])

    def test_screen_sequence_rejects_unsupported_mismatches(self) -> None:
        """Screening should reject fuzzy settings above two mismatches."""
        panel = {}
        record = SequenceRecord(identifier="read1", description="read1", sequence="AAAAA")
        with self.assertRaises(ValueError):
            screen_sequence_for_kmers(
                sequence_record=record,
                panel_index=panel,
                sample_id="sample1",
                sequence_type="read",
                max_mismatches=3,
            )

    def test_mismatch_neighbourhood_one_mismatch(self) -> None:
        """One-mismatch neighbourhood should generate expected substitutions."""
        neighbours = set(iter_mismatch_neighbourhood(kmer="AAA", max_mismatches=1))
        self.assertIn("CAA", neighbours)
        self.assertIn("ACA", neighbours)
        self.assertIn("AAC", neighbours)
        self.assertEqual(len(neighbours), 9)

    def test_mismatch_neighbourhood_rejects_above_two(self) -> None:
        """Mismatch neighbourhood should reject settings above two mismatches."""
        with self.assertRaises(ValueError):
            list(iter_mismatch_neighbourhood(kmer="AAA", max_mismatches=3))

    def test_parallel_screening_matches_single_worker(self) -> None:
        """Parallel screening should match single-worker screening."""
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
        records = [
            SequenceRecord(identifier="read1", description="read1", sequence="GGGAAAAATTT"),
            SequenceRecord(identifier="read2", description="read2", sequence="CCCCCAAAAAG"),
        ]
        serial = screen_records_for_species_kmers(
            records=records,
            panel_index=panel,
            sample_id="sample1",
            sequence_type="read",
            threads=1,
            chunk_size=1,
        )
        parallel = screen_records_for_species_kmers(
            records=records,
            panel_index=panel,
            sample_id="sample1",
            sequence_type="read",
            threads=2,
            chunk_size=1,
        )
        serial_keys = {(hit.sequence_id, hit.query_position, hit.matched_kmer) for hit in serial}
        parallel_keys = {(hit.sequence_id, hit.query_position, hit.matched_kmer) for hit in parallel}
        self.assertEqual(serial_keys, parallel_keys)

    def test_parallel_screening_rejects_bad_thread_count(self) -> None:
        """Parallel screening should reject non-positive thread counts."""
        with self.assertRaises(ValueError):
            screen_records_for_species_kmers(
                records=[],
                panel_index={},
                sample_id="sample1",
                sequence_type="read",
                threads=0,
            )

    def test_parallel_screening_rejects_bad_chunk_size(self) -> None:
        """Parallel screening should reject non-positive chunk sizes."""
        with self.assertRaises(ValueError):
            screen_records_for_species_kmers(
                records=[],
                panel_index={},
                sample_id="sample1",
                sequence_type="read",
                chunk_size=0,
            )

    def test_streaming_parallel_screening_rejects_bad_pending_count(self) -> None:
        """Streaming scheduler should reject non-positive pending chunk limits."""
        with self.assertRaises(ValueError):
            screen_records_for_species_kmers(
                records=[],
                panel_index={},
                sample_id="sample1",
                sequence_type="read",
                threads=2,
                max_pending_chunks=0,
            )

    def test_streaming_parallel_screening_matches_single_worker(self) -> None:
        """Streaming scheduler should preserve hits when pending chunks are limited."""
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
        records = [
            SequenceRecord(identifier=f"read{i}", description=f"read{i}", sequence="GGGAAAAATTT")
            for i in range(6)
        ]
        serial = screen_records_for_species_kmers(
            records=records,
            panel_index=panel,
            sample_id="sample1",
            sequence_type="read",
            threads=1,
            chunk_size=2,
        )
        parallel = screen_records_for_species_kmers(
            records=records,
            panel_index=panel,
            sample_id="sample1",
            sequence_type="read",
            threads=2,
            chunk_size=2,
            max_pending_chunks=1,
        )
        serial_ids = sorted(hit.sequence_id for hit in serial)
        parallel_ids = sorted(hit.sequence_id for hit in parallel)
        self.assertEqual(serial_ids, parallel_ids)
