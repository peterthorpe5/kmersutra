"""Tests for input/output helpers."""

import gzip
import io
import logging
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import Mock, patch

from kmersutra.io import (
    open_text,
    open_text_reader,
    read_tsv,
    resolve_decompressor,
    write_json,
    write_tsv,
    write_tsv_stream,
)


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

    def test_write_tsv_requires_columns_for_empty_records(self) -> None:
        """An empty iterable without columns should fail clearly."""
        with TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "empty.tsv"
            with self.assertRaisesRegex(ValueError, "Cannot infer TSV"):
                write_tsv(records=[], output_path=path)

    def test_read_tsv_empty_file_returns_no_records(self) -> None:
        """An empty TSV should produce an empty record list."""
        with TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "empty.tsv"
            path.write_text("", encoding="utf-8")
            self.assertEqual(read_tsv(input_path=path), [])

    def test_open_text_rejects_binary_mode(self) -> None:
        """The text-only helper should reject binary modes."""
        with self.assertRaisesRegex(ValueError, "text modes"):
            open_text("unused.tsv", "rb")

    def test_write_json_creates_parent_and_stable_text(self) -> None:
        """JSON output should create parents and use sorted indentation."""
        with TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "nested" / "record.json"
            write_json(data={"z": 2, "a": 1}, output_path=path)
            self.assertEqual(
                path.read_text(encoding="utf-8"),
                '{\n  "a": 1,\n  "z": 2\n}\n',
            )

    def test_write_tsv_stream_returns_record_count(self) -> None:
        """Streaming TSV output should preserve blanks and count records."""
        with TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "stream.tsv"
            count = write_tsv_stream(
                records=({"a": index} for index in [1, 2]),
                output_path=path,
                fieldnames=["a", "missing"],
            )
            self.assertEqual(count, 2)
            self.assertEqual(
                path.read_text(encoding="utf-8"),
                "a\tmissing\n1\t\n2\t\n",
            )


class TestInputDecompression(unittest.TestCase):
    """Tests for optional input decompression helpers."""

    def test_open_text_reader_reads_python_gzip(self) -> None:
        """Python gzip decompression should preserve input text."""
        with TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "input.txt.gz"
            with gzip.open(path, "wt", encoding="utf-8") as handle:
                handle.write("a\nb\n")
            with open_text_reader(path, decompressor="python") as handle:
                text = handle.read()
        self.assertEqual(text, "a\nb\n")

    def test_auto_decompressor_falls_back_without_pigz(self) -> None:
        """Auto mode should use Python gzip when pigz is unavailable."""
        logger = Mock(spec=logging.Logger)
        with patch("kmersutra.io.shutil.which", return_value=None):
            resolved = resolve_decompressor(
                path="reads.fastq.gz",
                decompressor="auto",
                logger=logger,
            )
        self.assertEqual(resolved, "python")
        logger.info.assert_called_once()

    def test_explicit_pigz_requires_executable(self) -> None:
        """Explicit pigz mode should fail clearly when pigz is unavailable."""
        with patch("kmersutra.io.shutil.which", return_value=None):
            with self.assertRaises(RuntimeError):
                resolve_decompressor(path="reads.fastq.gz", decompressor="pigz")

    def test_plain_input_ignores_decompressor(self) -> None:
        """Plain input should resolve to a normal text reader."""
        self.assertEqual(
            resolve_decompressor(path="reads.fastq", decompressor="pigz"),
            "plain",
        )

    def test_pigz_decompression_matches_python_when_available(self) -> None:
        """pigz decompression should match Python gzip when pigz is installed."""
        import shutil

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

    def test_invalid_decompressor_is_rejected(self) -> None:
        """An unknown decompressor name should fail before opening input."""
        with self.assertRaisesRegex(ValueError, "decompressor must be"):
            resolve_decompressor(path="reads.fastq.gz", decompressor="xz")

    def test_auto_resolves_to_pigz_when_available(self) -> None:
        """Auto mode should select pigz for gzip input when available."""
        with patch("kmersutra.io.shutil.which", return_value="/usr/bin/pigz"):
            self.assertEqual(
                resolve_decompressor(
                    path="reads.fastq.gz",
                    decompressor="auto",
                ),
                "pigz",
            )

    def test_plain_reader_reports_path_to_logger(self) -> None:
        """Plain input should be read directly and logged."""
        logger = Mock(spec=logging.Logger)
        with TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "input.txt"
            path.write_text("plain\n", encoding="utf-8")
            with open_text_reader(path, logger=logger) as handle:
                self.assertEqual(handle.read(), "plain\n")
        logger.info.assert_called_once()

    def test_python_gzip_reader_reports_path_to_logger(self) -> None:
        """Python gzip selection should be visible through the logger."""
        logger = Mock(spec=logging.Logger)
        with TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "input.txt.gz"
            with gzip.open(path, "wt", encoding="utf-8") as handle:
                handle.write("compressed\n")
            with open_text_reader(
                path,
                decompressor="python",
                logger=logger,
            ) as handle:
                self.assertEqual(handle.read(), "compressed\n")
        logger.info.assert_called_once()

    def test_mocked_pigz_reader_returns_text(self) -> None:
        """The pigz path should wrap binary stdout as UTF-8 text."""
        process = Mock()
        process.stdout = io.BytesIO(b"one\ntwo\n")
        process.stderr = io.BytesIO(b"")
        process.wait.return_value = 0
        logger = Mock(spec=logging.Logger)
        with (
            patch("kmersutra.io.shutil.which", return_value="/usr/bin/pigz"),
            patch("kmersutra.io.subprocess.Popen", return_value=process),
        ):
            with open_text_reader(
                "reads.fastq.gz",
                decompressor="pigz",
                logger=logger,
            ) as handle:
                self.assertEqual(handle.read(), "one\ntwo\n")
        process.wait.assert_called_once()
        logger.info.assert_called_once()

    def test_mocked_pigz_failure_includes_stderr(self) -> None:
        """A failing pigz process should report its exit status and stderr."""
        process = Mock()
        process.stdout = io.BytesIO(b"")
        process.stderr = io.BytesIO(b"broken gzip")
        process.wait.return_value = 2
        with (
            patch("kmersutra.io.shutil.which", return_value="/usr/bin/pigz"),
            patch("kmersutra.io.subprocess.Popen", return_value=process),
        ):
            with self.assertRaisesRegex(RuntimeError, "broken gzip"):
                with open_text_reader(
                    "reads.fastq.gz",
                    decompressor="pigz",
                ):
                    pass

    def test_mocked_pigz_sigpipe_is_not_reported_as_failure(self) -> None:
        """A SIGPIPE after a consumer exits early should be tolerated."""
        process = Mock()
        process.stdout = io.BytesIO(b"first\nsecond\n")
        process.stderr = None
        process.wait.return_value = 141
        with (
            patch("kmersutra.io.shutil.which", return_value="/usr/bin/pigz"),
            patch("kmersutra.io.subprocess.Popen", return_value=process),
        ):
            with open_text_reader(
                "reads.fastq.gz",
                decompressor="pigz",
            ) as handle:
                self.assertEqual(handle.readline(), "first\n")

    def test_context_exception_is_not_masked_by_pigz_failure(self) -> None:
        """A consumer exception should remain the visible failure."""
        process = Mock()
        process.stdout = io.BytesIO(b"text\n")
        process.stderr = io.BytesIO(b"secondary error")
        process.wait.return_value = 2
        with (
            patch("kmersutra.io.shutil.which", return_value="/usr/bin/pigz"),
            patch("kmersutra.io.subprocess.Popen", return_value=process),
        ):
            with self.assertRaisesRegex(LookupError, "consumer failed"):
                with open_text_reader(
                    "reads.fastq.gz",
                    decompressor="pigz",
                ):
                    raise LookupError("consumer failed")

    def test_pigz_without_stdout_fails_clearly(self) -> None:
        """A malformed process handle should fail before yielding."""
        process = Mock()
        process.stdout = None
        process.stderr = io.BytesIO(b"")
        with (
            patch("kmersutra.io.shutil.which", return_value="/usr/bin/pigz"),
            patch("kmersutra.io.subprocess.Popen", return_value=process),
        ):
            with self.assertRaisesRegex(RuntimeError, "stdout"):
                with open_text_reader(
                    "reads.fastq.gz",
                    decompressor="pigz",
                ):
                    pass
