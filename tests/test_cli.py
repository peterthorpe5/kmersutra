"""Tests for command-line module availability."""

from __future__ import annotations

import importlib
import unittest


class TestCliModules(unittest.TestCase):
    """Test command-line modules and wrappers are importable."""

    def test_build_panel_cli_imports(self):
        """Build-panel CLI module should expose a main function."""
        module = importlib.import_module("kmersutra.cli.build_clade_kmer_panel")
        self.assertTrue(callable(module.main))

    def test_screen_cli_imports(self):
        """Screening CLI module should expose a main function."""
        module = importlib.import_module("kmersutra.cli.screen_reads_for_clade_kmers")
        self.assertTrue(callable(module.main))

    def test_ncbi_cli_imports(self):
        """NCBI downloader CLI module should expose a main function."""
        module = importlib.import_module("kmersutra.cli.download_ncbi_taxon_genomes")
        self.assertTrue(callable(module.main))
