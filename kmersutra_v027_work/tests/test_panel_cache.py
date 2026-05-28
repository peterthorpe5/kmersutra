"""Tests for panel-index caching."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from kmersutra.build_panel import DiagnosticKmer
from kmersutra.io import write_tsv
from kmersutra.panel_cache import (
    get_default_panel_cache_path,
    load_panel_index_cache,
    load_panel_with_cache,
    write_panel_index_cache,
)


class TestPanelCache(unittest.TestCase):
    """Tests for cached KmerSutra panel indexes."""

    def _write_panel(self, path: Path) -> None:
        """Write a tiny diagnostic panel for cache tests."""
        write_tsv(
            records=[
                {
                    "kmer": "AAAAA",
                    "k": 5,
                    "panel_type": "species_unique",
                    "species_name": "Alpha",
                    "clade": "Demo",
                    "source_genomes": "g1",
                    "source_contigs": "c1",
                    "example_position": 0,
                    "evidence_taxid": "1",
                    "evidence_name": "Alpha",
                    "evidence_rank": "species",
                    "lineage_taxids": "1",
                    "source_taxids": "1",
                }
            ],
            output_path=path,
            fieldnames=[
                "kmer",
                "k",
                "panel_type",
                "species_name",
                "clade",
                "source_genomes",
                "source_contigs",
                "example_position",
                "evidence_taxid",
                "evidence_name",
                "evidence_rank",
                "lineage_taxids",
                "source_taxids",
            ],
        )

    def test_default_panel_cache_path_appends_suffix(self) -> None:
        """Default cache path should append .index.pkl to the panel path."""
        cache_path = get_default_panel_cache_path(panel_path="panel.tsv.gz")
        self.assertEqual(str(cache_path), "panel.tsv.gz.index.pkl")

    def test_write_and_load_panel_index_cache(self) -> None:
        """Panel cache should round-trip a plain diagnostic index."""
        with tempfile.TemporaryDirectory() as tmpdir:
            panel_path = Path(tmpdir) / "panel.tsv"
            cache_path = Path(tmpdir) / "panel.index.pkl"
            self._write_panel(panel_path)
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
            index = {5: {"AAAAA": [diagnostic]}}
            write_panel_index_cache(
                panel_index=index,
                panel_path=panel_path,
                cache_path=cache_path,
            )
            loaded = load_panel_index_cache(
                cache_path=cache_path,
                panel_path=panel_path,
            )
            self.assertIn(5, loaded)
            self.assertEqual(loaded[5]["AAAAA"][0].species_name, "Alpha")

    def test_load_panel_with_cache_creates_cache_when_requested(self) -> None:
        """Panel loader should create a cache when use_cache is requested."""
        with tempfile.TemporaryDirectory() as tmpdir:
            panel_path = Path(tmpdir) / "panel.tsv"
            cache_path = Path(tmpdir) / "panel.index.pkl"
            self._write_panel(panel_path)
            index, source = load_panel_with_cache(
                panel_path=panel_path,
                cache_path=cache_path,
                use_cache=True,
            )
            self.assertEqual(source, "tsv")
            self.assertTrue(cache_path.exists())
            self.assertEqual(index[5]["AAAAA"][0].species_name, "Alpha")

            cached_index, cached_source = load_panel_with_cache(
                panel_path=panel_path,
                cache_path=cache_path,
                use_cache=True,
            )
            self.assertEqual(cached_source, "cache")
            self.assertEqual(cached_index[5]["AAAAA"][0].species_name, "Alpha")

    def test_stale_cache_is_rejected(self) -> None:
        """Panel cache should be rejected if the panel file changes."""
        with tempfile.TemporaryDirectory() as tmpdir:
            panel_path = Path(tmpdir) / "panel.tsv"
            cache_path = Path(tmpdir) / "panel.index.pkl"
            self._write_panel(panel_path)
            index, _ = load_panel_with_cache(
                panel_path=panel_path,
                cache_path=cache_path,
                use_cache=True,
            )
            self.assertIn(5, index)
            panel_path.write_text(panel_path.read_text() + "\n", encoding="utf-8")
            with self.assertRaises(ValueError):
                load_panel_index_cache(
                    cache_path=cache_path,
                    panel_path=panel_path,
                )
