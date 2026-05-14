"""Tests for KmerSutra screening command output options."""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from kmersutra.cli.screen_reads_for_clade_kmers import main as screen_main
from kmersutra.io import write_tsv


PANEL_COLUMNS = [
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
]


class TestScreenCliOutputs(unittest.TestCase):
    """Tests for screening CLI speed and diagnostics outputs."""

    def _write_panel_and_fastq(self, tmpdir: str) -> tuple[Path, Path]:
        """Write a tiny panel and FASTQ for CLI tests."""
        base = Path(tmpdir)
        panel = base / "panel.tsv"
        fastq = base / "reads.fastq"
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
            output_path=panel,
            fieldnames=PANEL_COLUMNS,
        )
        fastq.write_text("@read1\nGGGAAAAATTT\n+\nFFFFFFFFFFF\n", encoding="utf-8")
        return panel, fastq

    def test_screen_cli_can_skip_read_level_hits(self) -> None:
        """CLI should skip large read-level hit output when requested."""
        with tempfile.TemporaryDirectory() as tmpdir:
            panel, fastq = self._write_panel_and_fastq(tmpdir)
            out_dir = Path(tmpdir) / "out"
            argv = [
                "kmersutra-screen",
                "--input",
                str(fastq),
                "--input_format",
                "fastq",
                "--panel",
                str(panel),
                "--sample_id",
                "sample1",
                "--out_dir",
                str(out_dir),
                "--threads",
                "1",
                "--chunk_size",
                "1",
                "--no_read_level_hits",
            ]
            with patch.object(sys, "argv", argv):
                screen_main()
            self.assertFalse((out_dir / "read_level_species_kmer_hits.tsv.gz").exists())
            self.assertTrue((out_dir / "read_level_species_kmer_hits.disabled.tsv").exists())
            self.assertTrue((out_dir / "species_detection_calls.tsv").exists())
            self.assertTrue((out_dir / "sample_taxonomic_kmer_evidence.tsv").exists())

    def test_screen_cli_writes_profile_and_cache(self) -> None:
        """CLI should write profile timings and a reusable panel cache."""
        with tempfile.TemporaryDirectory() as tmpdir:
            panel, fastq = self._write_panel_and_fastq(tmpdir)
            out_dir = Path(tmpdir) / "out"
            cache_path = Path(tmpdir) / "panel.index.pkl"
            argv = [
                "kmersutra-screen",
                "--input",
                str(fastq),
                "--input_format",
                "fastq",
                "--panel",
                str(panel),
                "--sample_id",
                "sample1",
                "--out_dir",
                str(out_dir),
                "--threads",
                "1",
                "--chunk_size",
                "1",
                "--panel_cache",
                str(cache_path),
                "--write_panel_cache",
                "--profile",
            ]
            with patch.object(sys, "argv", argv):
                screen_main()
            self.assertTrue(cache_path.exists())
            self.assertTrue((out_dir / "profile_timing.tsv").exists())
            profile_text = (out_dir / "profile_timing.tsv").read_text(encoding="utf-8")
            self.assertIn("load_panel", profile_text)
            self.assertIn("screen_records", profile_text)
