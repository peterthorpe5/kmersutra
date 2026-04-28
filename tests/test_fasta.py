"""Tests for FASTA and FASTQ parsing."""

import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from kmersutra.fasta import normalise_sequence, read_fasta_records, read_fastq_records


class TestFasta(unittest.TestCase):
    """Tests for sequence parsers."""

    def test_normalise_sequence_removes_whitespace_and_converts_u(self) -> None:
        """Normalisation should uppercase, remove whitespace and convert U."""
        self.assertEqual(normalise_sequence(" acg u\n t "), "ACGTT")

    def test_read_fasta_records_multiline(self) -> None:
        """FASTA parser should join multi-line records."""
        with TemporaryDirectory() as tmpdir:
            fasta_path = Path(tmpdir) / "input.fna"
            fasta_path.write_text(">seq1 description\nACG\nTTT\n>seq2\nNNN\n", encoding="utf-8")
            records = list(read_fasta_records(fasta_path=fasta_path))
        self.assertEqual(records[0].identifier, "seq1")
        self.assertEqual(records[0].sequence, "ACGTTT")
        self.assertEqual(records[1].identifier, "seq2")
        self.assertEqual(records[1].sequence, "NNN")

    def test_read_fastq_records(self) -> None:
        """FASTQ parser should return identifiers and normalised sequences."""
        with TemporaryDirectory() as tmpdir:
            fastq_path = Path(tmpdir) / "input.fastq"
            fastq_path.write_text("@read1 comment\nacgt\n+\n!!!!\n", encoding="utf-8")
            records = list(read_fastq_records(fastq_path=fastq_path))
        self.assertEqual(records[0].identifier, "read1")
        self.assertEqual(records[0].sequence, "ACGT")
