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
            diagnostics, summary = build_panel(genome_configs=configs, k_values=[5])
        species = {item.species_name for item in diagnostics if item.panel_type == "species_unique"}
        self.assertIn("Alpha", species)
        self.assertIn("Beta", species)
        self.assertTrue(all(item.species_name != "Out" for item in diagnostics))
        self.assertTrue(summary)
