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

    def test_extract_features_cli_imports(self):
        """Feature-extraction CLI module should expose a main function."""
        module = importlib.import_module("kmersutra.cli.extract_ml_features")
        self.assertTrue(callable(module.main))

    def test_train_classifier_cli_imports(self):
        """Classifier-training CLI module should expose a main function."""
        module = importlib.import_module("kmersutra.cli.train_classifier")
        self.assertTrue(callable(module.main))

    def test_predict_classifier_cli_imports(self):
        """Classifier-prediction CLI module should expose a main function."""
        module = importlib.import_module("kmersutra.cli.predict_classifier")
        self.assertTrue(callable(module.main))
    def test_summarise_run_cli_imports(self):
        """Run-summary CLI module should expose a main function."""
        module = importlib.import_module("kmersutra.cli.summarise_spikein_run")
        self.assertTrue(callable(module.main))


    def test_download_taxonomy_cli_imports(self):
        """Taxonomy downloader CLI module should expose a main function."""
        module = importlib.import_module("kmersutra.cli.download_taxonomy")
        self.assertTrue(callable(module.main))

    def test_merge_panels_cli_imports(self):
        """Panel-merge CLI module should expose a main function."""
        module = importlib.import_module("kmersutra.cli.merge_panels")
        self.assertTrue(callable(module.main))

    def test_validate_panel_cli_imports(self):
        """Panel-validation CLI module should expose a main function."""
        module = importlib.import_module("kmersutra.cli.validate_panel")
        self.assertTrue(callable(module.main))
