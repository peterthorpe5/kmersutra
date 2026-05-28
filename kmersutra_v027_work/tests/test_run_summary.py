"""Tests for KmerSutra run-level summary reporting."""

import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from kmersutra.io import read_tsv, write_tsv
from kmersutra.run_summary import (
    build_run_summary_reports,
    build_species_long_summary,
    find_species_labels,
    species_from_call_column,
    summarise_by_spike,
    summarise_call_counts,
)


class TestRunSummary(unittest.TestCase):
    """Tests for run-level summary helpers."""

    def _summary_records(self) -> list[dict[str, str]]:
        """Return a small KmerSutra spike-in summary table."""
        return [
            {
                "replicate": "1",
                "spike_n": "0",
                "total_spiked_reads": "0",
                "kmersutra_runtime_seconds": "10",
                "kmersutra_out_dir": "out0",
                "kmersutra_calls_tsv": "calls0.tsv",
                "kmersutra_evidence_tsv": "evidence0.tsv",
                "kmersutra_read_hits_tsv_gz": "hits0.tsv.gz",
                "kmersutra_Plasmodium_vivax_call": "not_detected",
                "kmersutra_Plasmodium_vivax_unique_kmers": "0",
                "kmersutra_Plasmodium_vivax_positive_reads": "0",
                "kmersutra_Plasmodium_vivax_confidence": "0",
            },
            {
                "replicate": "1",
                "spike_n": "5000",
                "total_spiked_reads": "15000",
                "kmersutra_runtime_seconds": "20",
                "kmersutra_out_dir": "out1",
                "kmersutra_calls_tsv": "calls1.tsv",
                "kmersutra_evidence_tsv": "evidence1.tsv",
                "kmersutra_read_hits_tsv_gz": "hits1.tsv.gz",
                "kmersutra_Plasmodium_vivax_call": "present_in_mixed_sample",
                "kmersutra_Plasmodium_vivax_unique_kmers": "2283",
                "kmersutra_Plasmodium_vivax_positive_reads": "307",
                "kmersutra_Plasmodium_vivax_confidence": "0.95",
            },
        ]

    def test_species_from_call_column(self) -> None:
        """Species label should be recovered from a wide call column name."""
        label = species_from_call_column(column="kmersutra_Plasmodium_vivax_call")
        self.assertEqual(label, "Plasmodium vivax")

    def test_find_species_labels(self) -> None:
        """Wide summary parsing should find KmerSutra species columns."""
        self.assertEqual(find_species_labels(records=self._summary_records()), ["Plasmodium vivax"])

    def test_build_species_long_summary(self) -> None:
        """Wide summary should convert to one row per species and spike."""
        long_rows = build_species_long_summary(records=self._summary_records())
        self.assertEqual(len(long_rows), 2)
        self.assertEqual(long_rows[1]["n_unique_kmers"], 2283)

    def test_summarise_call_counts(self) -> None:
        """Call count summary should count labels by species."""
        long_rows = build_species_long_summary(records=self._summary_records())
        counts = summarise_call_counts(species_long_rows=long_rows)
        self.assertEqual(sum(int(row["n_records"]) for row in counts), 2)

    def test_summarise_by_spike(self) -> None:
        """Spike summary should calculate present-call rate by spike level."""
        long_rows = build_species_long_summary(records=self._summary_records())
        by_spike = summarise_by_spike(species_long_rows=long_rows)
        high = [row for row in by_spike if row["spike_n"] == 5000][0]
        self.assertEqual(high["present_call_rate"], 1.0)

    def test_build_run_summary_reports_writes_html_and_xlsx(self) -> None:
        """Run summary builder should write formatted HTML and Excel outputs."""
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            summary_tsv = root / "summary.tsv"
            out_xlsx = root / "summary.xlsx"
            out_html = root / "summary.html"
            write_tsv(records=self._summary_records(), output_path=summary_tsv)
            written = build_run_summary_reports(
                summary_tsv=summary_tsv,
                out_xlsx=out_xlsx,
                out_html=out_html,
            )
            self.assertTrue(written["xlsx"].exists())
            self.assertTrue(written["html"].exists())
            self.assertIn("present_in_mixed_sample", out_html.read_text(encoding="utf-8"))

    def test_summary_tsv_roundtrip(self) -> None:
        """Summary records should remain parseable after TSV roundtrip."""
        with TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "summary.tsv"
            write_tsv(records=self._summary_records(), output_path=path)
            parsed = read_tsv(input_path=path)
            self.assertEqual(parsed[0]["kmersutra_Plasmodium_vivax_call"], "not_detected")
