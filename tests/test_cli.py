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
        """Legacy NCBI downloader CLI module should expose a main function."""
        module = importlib.import_module("kmersutra.cli.download_ncbi_taxon_genomes")
        self.assertTrue(callable(module.main))

    def test_download_genomes_cli_imports(self):
        """KmerSutra genome downloader CLI module should expose a main function."""
        module = importlib.import_module("kmersutra.cli.download_genomes_for_kmersutra")
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


    def test_build_call_training_cli_imports(self):
        """AI call-training CLI module should expose a main function."""
        module = importlib.import_module("kmersutra.cli.build_call_training_table")
        self.assertTrue(callable(module.main))

    def test_train_call_calibrator_cli_imports(self):
        """AI call-calibrator CLI module should expose a main function."""
        module = importlib.import_module("kmersutra.cli.train_call_calibrator")
        self.assertTrue(callable(module.main))

    def test_summarise_run_cli_imports(self):
        """Run-summary CLI module should expose a main function."""
        module = importlib.import_module("kmersutra.cli.summarise_spikein_run")
        self.assertTrue(callable(module.main))

    def test_summarise_spikeins_cli_imports(self):
        """Multi-run spike-in summary CLI module should expose a main function."""
        module = importlib.import_module("kmersutra.cli.summarise_kmersutra_spikeins")
        self.assertTrue(callable(module.main))


    def test_summarise_lca_cli_imports(self):
        """LCA summary CLI module should expose a main function."""
        module = importlib.import_module("kmersutra.cli.summarise_lca")
        self.assertTrue(callable(module.main))


    def test_download_taxonomy_cli_imports(self):
        """Taxonomy downloader CLI module should expose a main function."""
        module = importlib.import_module("kmersutra.cli.download_taxonomy")
        self.assertTrue(callable(module.main))

    def test_merge_panels_cli_imports(self):
        """Panel-merge CLI module should expose a main function."""
        module = importlib.import_module("kmersutra.cli.merge_panels")
        self.assertTrue(callable(module.main))

    def test_merge_modules_cli_imports(self):
        """Parquet module-merge CLI should expose a main function."""
        module = importlib.import_module("kmersutra.cli.merge_modules")
        self.assertTrue(callable(module.main))

    def test_validate_panel_cli_imports(self):
        """Panel-validation CLI module should expose a main function."""
        module = importlib.import_module("kmersutra.cli.validate_panel")
        self.assertTrue(callable(module.main))


class TestBuildPanelCliDefaults(unittest.TestCase):
    """Test build-panel CLI defaults that affect publication builds."""

    def test_marker_selection_defaults_to_independent_multik_genome_spread(self) -> None:
        """Independent multi-k marker selection should be the default build behaviour."""
        from kmersutra.cli.build_clade_kmer_panel import parse_args
        import sys
        from unittest.mock import patch

        argv = [
            "kmersutra-build-panel",
            "--genome_config",
            "config.tsv",
            "--out_dir",
            "out",
        ]
        with patch.object(sys, "argv", argv):
            args = parse_args()
        self.assertEqual(args.marker_selection, "independent_multik_genome_spread")

    def test_global_index_progress_interval_defaults_to_production_value(self) -> None:
        """Global index progress logging should be quiet by default."""
        from kmersutra.cli.build_clade_kmer_panel import parse_args
        import sys
        from unittest.mock import patch

        argv = [
            "kmersutra-build-panel",
            "--genome_config",
            "config.tsv",
            "--out_dir",
            "out",
        ]
        with patch.object(sys, "argv", argv):
            args = parse_args()
        self.assertEqual(args.global_index_progress_interval, 5_000_000)
