"""End-to-end tests for the KmerSutra workflow."""

import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from kmersutra.build_panel import build_panel
from kmersutra.config import GenomeConfig
from kmersutra.io import read_tsv, write_tsv
from kmersutra.reporting import write_html_report
from kmersutra.screen_reads import screen_file_for_species_kmers
from kmersutra.summarise_hits import summarise_sample_species_evidence, summarise_species_hits
from kmersutra.thresholds import call_species_presence


class TestEndToEnd(unittest.TestCase):
    """End-to-end workflow tests using deterministic toy data."""

    def test_panel_screen_summary_and_report_outputs(self) -> None:
        """Workflow should produce valid TSV and HTML outputs."""
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            alpha = root / "alpha.fna"
            beta = root / "beta.fna"
            reads = root / "reads.fastq"
            alpha.write_text(">a\nAAAAACCCCCGGGGG\n", encoding="utf-8")
            beta.write_text(">b\nTTTTTGGGGGCCCCC\n", encoding="utf-8")
            reads.write_text("@r1\nAAAAACCCCC\n+\n!!!!!!!!!!\n", encoding="utf-8")
            diagnostics, panel_summary = build_panel(
                genome_configs=[
                    GenomeConfig(genome_fasta=alpha, species_name="Alpha", role="target_species", clade="Demo"),
                    GenomeConfig(genome_fasta=beta, species_name="Beta", role="target_species", clade="Demo"),
                ],
                k_values=[5],
            )
            panel_path = root / "panel.tsv"
            write_tsv(
                records=[item.to_record() for item in diagnostics],
                output_path=panel_path,
                fieldnames=[
                    "kmer", "k", "panel_type", "species_name", "clade",
                    "source_genomes", "source_contigs", "example_position",
                ],
            )
            hits = screen_file_for_species_kmers(
                input_path=reads,
                panel_path=panel_path,
                sample_id="sample1",
                input_format="fastq",
            )
            hit_summary = summarise_species_hits(hits=hits)
            evidence = summarise_sample_species_evidence(species_summary=hit_summary)
            calls = call_species_presence(
                evidence_records=evidence,
                min_unique_kmers=1,
                min_positive_sequences=1,
            )
            calls_path = root / "species_detection_calls.tsv"
            html_path = root / "report.html"
            write_tsv(records=calls, output_path=calls_path)
            write_html_report(
                output_path=html_path,
                title="Test report",
                panel_summary=panel_summary,
                hit_summary=hit_summary,
                detection_calls=calls,
            )
            parsed_calls = read_tsv(input_path=calls_path)
            self.assertTrue(parsed_calls)
            self.assertTrue(html_path.exists())
            self.assertIn("KmerSutra", html_path.read_text(encoding="utf-8"))
