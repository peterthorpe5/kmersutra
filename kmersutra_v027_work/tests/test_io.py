"""Tests for input/output helpers."""

import gzip
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from kmersutra.io import read_tsv, write_tsv


class TestIo(unittest.TestCase):
    """Tests for TSV reading and writing."""

    def test_read_tsv_pads_missing_values(self) -> None:
        """TSV reader should pad missing trailing values with blanks."""
        with TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "input.tsv"
            path.write_text("a\tb\tc\n1\t2\n", encoding="utf-8")
            records = read_tsv(input_path=path)
        self.assertEqual(records, [{"a": "1", "b": "2", "c": ""}])

    def test_read_tsv_truncates_extra_values(self) -> None:
        """TSV reader should ignore extra fields beyond the header."""
        with TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "input.tsv"
            path.write_text("a\tb\n1\t2\t3\n", encoding="utf-8")
            records = read_tsv(input_path=path)
        self.assertEqual(records, [{"a": "1", "b": "2"}])

    def test_write_tsv_empty_records_with_fieldnames(self) -> None:
        """TSV writer should allow header-only files when columns are known."""
        with TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "empty.tsv"
            write_tsv(records=[], output_path=path, fieldnames=["a", "b"])
            self.assertEqual(path.read_text(encoding="utf-8"), "a\tb\n")

    def test_write_tsv_gzip_output(self) -> None:
        """TSV writer should support gzip-compressed output."""
        with TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "records.tsv.gz"
            write_tsv(
                records=[{"a": "1", "b": "2"}],
                output_path=path,
                fieldnames=["a", "b"],
            )
            with gzip.open(path, "rt", encoding="utf-8") as handle:
                text = handle.read()
        self.assertEqual(text, "a\tb\n1\t2\n")
