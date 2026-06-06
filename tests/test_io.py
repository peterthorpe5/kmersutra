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


class TestInputDecompression(unittest.TestCase):
    """Tests for optional input decompression helpers."""

    def test_open_text_reader_reads_python_gzip(self) -> None:
        """Python gzip decompression should preserve input text."""
        from kmersutra.io import open_text_reader

        with TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "input.txt.gz"
            with gzip.open(path, "wt", encoding="utf-8") as handle:
                handle.write("a\nb\n")
            with open_text_reader(path, decompressor="python") as handle:
                text = handle.read()
        self.assertEqual(text, "a\nb\n")

    def test_auto_decompressor_falls_back_without_pigz(self) -> None:
        """Auto mode should use Python gzip when pigz is unavailable."""
        from unittest.mock import patch

        from kmersutra.io import resolve_decompressor

        with patch("kmersutra.io.shutil.which", return_value=None):
            resolved = resolve_decompressor(
                path="reads.fastq.gz",
                decompressor="auto",
            )
        self.assertEqual(resolved, "python")

    def test_explicit_pigz_requires_executable(self) -> None:
        """Explicit pigz mode should fail clearly when pigz is unavailable."""
        from unittest.mock import patch

        from kmersutra.io import resolve_decompressor

        with patch("kmersutra.io.shutil.which", return_value=None):
            with self.assertRaises(RuntimeError):
                resolve_decompressor(path="reads.fastq.gz", decompressor="pigz")

    def test_plain_input_ignores_decompressor(self) -> None:
        """Plain input should resolve to a normal text reader."""
        from kmersutra.io import resolve_decompressor

        self.assertEqual(
            resolve_decompressor(path="reads.fastq", decompressor="pigz"),
            "plain",
        )

    def test_pigz_decompression_matches_python_when_available(self) -> None:
        """pigz decompression should match Python gzip when pigz is installed."""
        import shutil

        from kmersutra.io import open_text_reader

        if shutil.which("pigz") is None:
            self.skipTest("pigz is not installed on this test host")
        with TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "input.txt.gz"
            with gzip.open(path, "wt", encoding="utf-8") as handle:
                handle.write("a\nb\n")
            with open_text_reader(path, decompressor="python") as handle:
                python_text = handle.read()
            with open_text_reader(path, decompressor="pigz") as handle:
                pigz_text = handle.read()
        self.assertEqual(pigz_text, python_text)
