"""Tests for optional Parquet diagnostic panel helpers."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from kmersutra.panel_parquet import (
    OptionalParquetDependencyError,
    derive_panel_parquet_path,
    is_parquet_path,
    pyarrow_available,
    read_panel_parquet,
    select_screen_panel_path,
    write_panel_parquet,
)


class TestPanelParquetHelpers(unittest.TestCase):
    """Test optional Parquet panel helper behaviour."""

    def test_is_parquet_path_identifies_suffix(self) -> None:
        """Parquet suffix detection should be explicit and case-insensitive."""
        self.assertTrue(is_parquet_path(path="panel.parquet"))
        self.assertTrue(is_parquet_path(path="panel.PARQUET"))
        self.assertFalse(is_parquet_path(path="panel.tsv.gz"))

    def test_derive_panel_parquet_path_from_tsv_gz(self) -> None:
        """Default companion path should replace TSV suffixes cleanly."""
        self.assertEqual(
            derive_panel_parquet_path(panel_path="species_kmer_panel.tsv.gz").name,
            "species_kmer_panel.parquet",
        )
        self.assertEqual(
            derive_panel_parquet_path(panel_path="species_kmer_panel.tsv").name,
            "species_kmer_panel.parquet",
        )

    def test_select_screen_panel_path_prefers_existing_parquet_for_auto(self) -> None:
        """Auto selection should use Parquet when the companion exists."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            tsv_path = tmp_path / "panel.tsv.gz"
            parquet_path = tmp_path / "panel.parquet"
            tsv_path.write_text("placeholder\n", encoding="utf-8")
            parquet_path.write_text("placeholder\n", encoding="utf-8")
            observed = select_screen_panel_path(
                tsv_path=tsv_path,
                parquet_path=parquet_path,
                panel_storage_format="auto",
            )
            self.assertEqual(observed, parquet_path)

    def test_select_screen_panel_path_uses_tsv_without_parquet(self) -> None:
        """Auto selection should remain TSV-compatible when Parquet is absent."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            tsv_path = tmp_path / "panel.tsv.gz"
            tsv_path.write_text("placeholder\n", encoding="utf-8")
            observed = select_screen_panel_path(
                tsv_path=tsv_path,
                parquet_path=tmp_path / "panel.parquet",
                panel_storage_format="auto",
            )
            self.assertEqual(observed, tsv_path)

    def test_select_screen_panel_path_rejects_missing_required_parquet(self) -> None:
        """A forced Parquet selection should fail when the file is absent."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            with self.assertRaises(ValueError):
                select_screen_panel_path(
                    tsv_path=tmp_path / "panel.tsv.gz",
                    parquet_path=tmp_path / "panel.parquet",
                    panel_storage_format="parquet",
                )

    def test_parquet_read_write_requires_optional_dependency_when_missing(self) -> None:
        """Parquet I/O should fail clearly when pyarrow is unavailable."""
        if pyarrow_available():
            self.skipTest("pyarrow is installed; missing-dependency path is not active")
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "panel.parquet"
            with self.assertRaises(OptionalParquetDependencyError):
                write_panel_parquet(records=[], output_path=path)
            path.write_text("not parquet", encoding="utf-8")
            with self.assertRaises(OptionalParquetDependencyError):
                read_panel_parquet(input_path=path)


if __name__ == "__main__":
    unittest.main()
