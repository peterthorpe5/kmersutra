"""Tests for generic TSV, TSV.GZ and Parquet-style table I/O."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from kmersutra.parquet_modules import OptionalParquetDependencyError
from kmersutra.panel_parquet import pyarrow_available
from kmersutra.table_io import (
    infer_table_format,
    read_records_table,
    validate_required_columns,
    write_records_table,
)


class TestGenericTableIo(unittest.TestCase):
    """Test generic table I/O helpers."""

    def test_infer_table_format_handles_supported_suffixes(self) -> None:
        """Format inference should support TSV, TSV.GZ and Parquet paths."""
        self.assertEqual(infer_table_format(path="a.tsv"), "tsv")
        self.assertEqual(infer_table_format(path="a.tsv.gz"), "tsv.gz")
        self.assertEqual(infer_table_format(path="a.parquet"), "parquet")

    def test_infer_table_format_rejects_unsupported_suffix(self) -> None:
        """Unsupported table suffixes should fail clearly."""
        with self.assertRaises(ValueError):
            infer_table_format(path="a.csv")

    def test_tsv_gz_round_trip_preserves_records(self) -> None:
        """Compressed TSV output should round-trip through generic I/O."""
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "records.tsv.gz"
            records = [{"sample_id": "s1", "value": 1}, {"sample_id": "s2", "value": 2}]
            n_written = write_records_table(
                records=records,
                output_path=path,
                fieldnames=["sample_id", "value"],
            )
            observed = read_records_table(
                input_path=path,
                required_columns=["sample_id", "value"],
            )
        self.assertEqual(n_written, 2)
        self.assertEqual(
            observed,
            [{"sample_id": "s1", "value": "1"}, {"sample_id": "s2", "value": "2"}],
        )

    def test_required_column_validation_reports_missing_columns(self) -> None:
        """Required-column validation should identify missing fields."""
        with self.assertRaisesRegex(ValueError, "missing"):
            validate_required_columns(
                records=[{"present": "yes"}],
                required_columns=["missing"],
            )

    def test_parquet_missing_dependency_fails_clearly(self) -> None:
        """Parquet output should request the optional dependency when unavailable."""
        if pyarrow_available():
            self.skipTest("pyarrow is installed; missing-dependency path is inactive")
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "records.parquet"
            with self.assertRaises(OptionalParquetDependencyError):
                write_records_table(records=[{"a": 1}], output_path=path)


if __name__ == "__main__":
    unittest.main()
