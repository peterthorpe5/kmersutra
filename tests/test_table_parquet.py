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


if __name__ == "__main__":
    unittest.main()
