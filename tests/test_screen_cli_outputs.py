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
            self.assertTrue((out_dir / "sample_lineage_interpretation.tsv").exists())

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

class TestScreenCliConsolidationOutputs(TestScreenCliOutputs):
    """Tests for CLI consolidation and optional Parquet output flags."""

    def test_screen_cli_can_consolidate_background_candidate_call(self) -> None:
        """CLI should write raw and consolidated calls when consolidation is enabled."""
        with tempfile.TemporaryDirectory() as tmpdir:
            base = Path(tmpdir)
            panel = base / "panel.tsv"
            fastq = base / "reads.fastq"
            write_tsv(
                records=[
                    {
                        "kmer": "AAAAA",
                        "k": 5,
                        "panel_type": "species_unique",
                        "species_name": "Hammondia hammondi",
                        "clade": "Apicomplexa",
                        "source_genomes": "g1",
                        "source_contigs": "c1",
                        "example_position": 0,
                        "evidence_taxid": "1",
                        "evidence_name": "Hammondia hammondi",
                        "evidence_rank": "species",
                        "lineage_taxids": "1",
                        "source_taxids": "1",
                    }
                ],
                output_path=panel,
                fieldnames=PANEL_COLUMNS,
            )
            fastq.write_text("@read1\nGGGAAAAATTT\n+\nFFFFFFFFFFF\n", encoding="utf-8")
            out_dir = base / "out"
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
                "--min_unique_kmers",
                "1",
                "--min_positive_sequences",
                "1",
                "--min_k_values_positive",
                "1",
                "--min_best_k",
                "5",
                "--min_exact_hits",
                "1",
                "--min_confidence_score",
                "0",
                "--consolidate_species_calls",
                "--background_candidate_taxa",
                "Hammondia hammondi",
                "--no_read_level_hits",
            ]
            with patch.object(sys, "argv", argv):
                screen_main()
            raw_text = (out_dir / "species_detection_calls_raw.tsv").read_text(
                encoding="utf-8"
            )
            consolidated_text = (out_dir / "species_detection_calls.tsv").read_text(
                encoding="utf-8"
            )
            self.assertIn("present_high_confidence", raw_text)
            self.assertIn("background_candidate_signal", consolidated_text)
            self.assertIn("pre_consolidation_call", consolidated_text)

    def test_screen_cli_parquet_flag_warns_without_pyarrow(self) -> None:
        """Parquet output should not make screening fail without pyarrow."""
        try:
            import pyarrow  # type: ignore  # noqa: F401
        except ImportError:
            pyarrow_available = False
        else:
            pyarrow_available = True
        if pyarrow_available:
            self.skipTest("pyarrow installed; warning path not active")
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
                "--write_parquet_outputs",
                "--no_read_level_hits",
            ]
            with patch.object(sys, "argv", argv):
                screen_main()
            self.assertTrue((out_dir / "species_detection_calls.tsv").exists())

