"""Tests for genome configuration parsing."""

import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from kmersutra.config import load_genome_config


class TestConfig(unittest.TestCase):
    """Tests for configuration parsing."""

    def test_load_genome_config(self) -> None:
        """Genome config parser should retain required metadata."""
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            config_path = root / "config.tsv"
            fasta_path = root / "a.fna"
            fasta_path.write_text(">a\nACGT\n", encoding="utf-8")
            config_path.write_text(
                "genome_fasta\tspecies_name\trole\tclade\n"
                f"{fasta_path}\tSpecies alpha\ttarget_species\tDemo\n",
                encoding="utf-8",
            )
            records = load_genome_config(config_path=config_path)
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0].species_name, "Species alpha")
        self.assertTrue(records[0].is_target)
