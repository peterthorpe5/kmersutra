"""Tests for panel building."""

import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from kmersutra.build_panel import build_panel
from kmersutra.config import GenomeConfig


class TestBuildPanel(unittest.TestCase):
    """Tests for diagnostic panel construction."""

    def test_build_panel_species_unique(self) -> None:
        """Panel builder should remove k-mers found in outgroups."""
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            alpha = root / "alpha.fna"
            beta = root / "beta.fna"
            outgroup = root / "outgroup.fna"
            alpha.write_text(">a\nAAAAACCCCC\n", encoding="utf-8")
            beta.write_text(">b\nGGGGGCCCCC\n", encoding="utf-8")
            outgroup.write_text(">o\nTTTTTCCCCC\n", encoding="utf-8")
            configs = [
                GenomeConfig(genome_fasta=alpha, species_name="Alpha", role="target_species", clade="Demo"),
                GenomeConfig(genome_fasta=beta, species_name="Beta", role="target_species", clade="Demo"),
                GenomeConfig(genome_fasta=outgroup, species_name="Out", role="outgroup", clade="Out"),
            ]
            diagnostics, summary, collection_summary = build_panel(
                genome_configs=configs,
                k_values=[5],
            )
        species = {item.species_name for item in diagnostics if item.panel_type == "species_unique"}
        self.assertIn("Alpha", species)
        self.assertIn("Beta", species)
        self.assertTrue(all(item.species_name != "Out" for item in diagnostics))
        self.assertTrue(summary)
        self.assertEqual(len(collection_summary), 3)

    def test_build_panel_threads_matches_single_worker(self) -> None:
        """Parallel panel building should match single-worker output."""
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            alpha = root / "alpha.fna"
            beta = root / "beta.fna"
            alpha.write_text(">a\nAAAAACCCCCGGGGG\n", encoding="utf-8")
            beta.write_text(">b\nTTTTTGGGGGCCCCC\n", encoding="utf-8")
            configs = [
                GenomeConfig(genome_fasta=alpha, species_name="Alpha", role="target_species", clade="Demo"),
                GenomeConfig(genome_fasta=beta, species_name="Beta", role="target_species", clade="Demo"),
            ]
            serial, _, _ = build_panel(genome_configs=configs, k_values=[5], threads=1)
            parallel, _, _ = build_panel(genome_configs=configs, k_values=[5], threads=2)
        serial_keys = {(item.k, item.kmer, item.species_name) for item in serial}
        parallel_keys = {(item.k, item.kmer, item.species_name) for item in parallel}
        self.assertEqual(serial_keys, parallel_keys)

    def test_build_panel_thinning_limit(self) -> None:
        """Panel thinning should limit retained k-mers per species and k."""
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            alpha = root / "alpha.fna"
            beta = root / "beta.fna"
            alpha.write_text(">a\nAAAAACCCCCGGGGGTTTTT\n", encoding="utf-8")
            beta.write_text(">b\nCCCCCGGGGGTTTTTAAAAA\n", encoding="utf-8")
            configs = [
                GenomeConfig(genome_fasta=alpha, species_name="Alpha", role="target_species", clade="Demo"),
                GenomeConfig(genome_fasta=beta, species_name="Beta", role="target_species", clade="Demo"),
            ]
            diagnostics, _, _ = build_panel(
                genome_configs=configs,
                k_values=[5],
                max_per_species_per_k=2,
            )
        counts = {}
        for item in diagnostics:
            key = (item.panel_type, item.species_name, item.k)
            counts[key] = counts.get(key, 0) + 1
        self.assertTrue(all(count <= 2 for count in counts.values()))

    def test_build_panel_clade_core_kmers(self) -> None:
        """Panel builder should optionally include clade-core k-mers."""
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            alpha = root / "alpha.fna"
            beta = root / "beta.fna"
            outgroup = root / "outgroup.fna"
            alpha.write_text(">a\nAAAAACCCCC\n", encoding="utf-8")
            beta.write_text(">b\nAAAAAGGGGG\n", encoding="utf-8")
            outgroup.write_text(">o\nTTTTTGGGGG\n", encoding="utf-8")
            configs = [
                GenomeConfig(genome_fasta=alpha, species_name="Alpha", role="target_species", clade="Demo"),
                GenomeConfig(genome_fasta=beta, species_name="Beta", role="target_species", clade="Demo"),
                GenomeConfig(genome_fasta=outgroup, species_name="Out", role="outgroup", clade="Other"),
            ]
            diagnostics, _, _ = build_panel(
                genome_configs=configs,
                k_values=[5],
                target_clade="Demo",
            )
        panel_types = {item.panel_type for item in diagnostics}
        self.assertIn("clade_core", panel_types)
