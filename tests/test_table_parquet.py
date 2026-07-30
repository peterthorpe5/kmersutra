"""Tests for optional generic Parquet table helpers."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from kmersutra.parquet_modules import OptionalParquetDependencyError
from kmersutra.panel_parquet import pyarrow_available
from kmersutra.table_parquet import (
    infer_fieldnames,
    records_to_string_columns,
    read_records_parquet,
    write_records_parquet,
)


class TestTableParquetHelpers(unittest.TestCase):
    """Test optional generic Parquet helpers."""

    def test_infer_fieldnames_preserves_first_seen_order(self) -> None:
        """Fieldnames should be inferred in deterministic first-seen order."""
        observed = infer_fieldnames(records=[{"b": 1, "a": 2}, {"c": 3, "a": 4}])
        self.assertEqual(observed, ["b", "a", "c"])

    def test_records_to_string_columns_uses_explicit_fieldnames(self) -> None:
        """Record conversion should preserve explicit output columns."""
        columns, values = records_to_string_columns(
            records=[{"a": 1, "b": None}],
            fieldnames=["b", "a", "missing"],
        )
        self.assertEqual(columns, ["b", "a", "missing"])
        self.assertEqual(values, {"b": [""], "a": ["1"], "missing": [""]})

    def test_records_to_string_columns_infers_union(self) -> None:
        """Implicit fields should use the first-seen union."""
        columns, values = records_to_string_columns(
            records=[{"b": 1}, {"a": 2, "b": None}],
        )
        self.assertEqual(columns, ["b", "a"])
        self.assertEqual(values, {"b": ["1", ""], "a": ["", "2"]})

    def test_read_write_requires_pyarrow_when_missing(self) -> None:
        """Parquet I/O should fail clearly without the optional dependency."""
        if pyarrow_available():
            self.skipTest("pyarrow is installed; missing-dependency path is inactive")
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "table.parquet"
            with self.assertRaises(OptionalParquetDependencyError):
                write_records_parquet(records=[{"a": 1}], output_path=path)
            with self.assertRaises(OptionalParquetDependencyError):
                read_records_parquet(input_path=path)

    def test_read_write_round_trip_when_pyarrow_available(self) -> None:
        """Parquet I/O should preserve columns and stringified values."""
        if not pyarrow_available():
            self.skipTest("pyarrow is not installed")
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "nested" / "table.parquet"
            count = write_records_parquet(
                records=[{"a": 1, "b": None}, {"a": 2, "b": "x"}],
                output_path=path,
                fieldnames=["b", "a"],
            )
            observed = read_records_parquet(input_path=path)
        self.assertEqual(count, 2)
        self.assertEqual(
            observed,
            [{"b": "", "a": "1"}, {"b": "x", "a": "2"}],
        )

    def test_read_rejects_missing_and_empty_files(self) -> None:
        """Missing and empty Parquet inputs should fail clearly."""
        if not pyarrow_available():
            self.skipTest("pyarrow is not installed")
        with tempfile.TemporaryDirectory() as tmpdir:
            missing = Path(tmpdir) / "missing.parquet"
            with self.assertRaisesRegex(FileNotFoundError, "missing or empty"):
                read_records_parquet(input_path=missing)
            empty = Path(tmpdir) / "empty.parquet"
            empty.touch()
            with self.assertRaisesRegex(FileNotFoundError, "missing or empty"):
                read_records_parquet(input_path=empty)

    def test_write_empty_table_with_explicit_schema(self) -> None:
        """An empty record set should work when columns are explicit."""
        if not pyarrow_available():
            self.skipTest("pyarrow is not installed")
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "empty.parquet"
            count = write_records_parquet(
                records=[],
                output_path=path,
                fieldnames=["sample_id"],
            )
            observed = read_records_parquet(input_path=path)
        self.assertEqual(count, 0)
        self.assertEqual(observed, [])


if __name__ == "__main__":
    unittest.main()
