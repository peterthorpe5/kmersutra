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

    def test_load_genome_config_rejects_missing_required_column(self) -> None:
        """Genome config parser should reject missing required columns."""
        with TemporaryDirectory() as tmpdir:
            config_path = Path(tmpdir) / "config.tsv"
            config_path.write_text(
                "genome_fasta\trole\n"
                "a.fna\ttarget_species\n",
                encoding="utf-8",
            )
            with self.assertRaises(ValueError):
                load_genome_config(config_path=config_path)


    def test_load_genome_config_accepts_specialised_outgroup_roles(self) -> None:
        """Genome config parser should accept biologically useful outgroup roles."""
        with TemporaryDirectory() as tmpdir:
            config_path = Path(tmpdir) / "config.tsv"
            config_path.write_text(
                "genome_fasta\tspecies_name\trole\tclade\n"
                "a.fna\tSpecies alpha\ttarget_species\tDemo\n"
                "b.fna\tSpecies beta\tapicomplexan_outgroup\tApicomplexa\n"
                "c.fna\tSpecies gamma\tdistant_outgroup\tCiliophora\n"
                "d.fna\tSpecies delta\thost_or_background\tHost\n",
                encoding="utf-8",
            )
            records = load_genome_config(config_path=config_path)
        self.assertEqual(len(records), 4)
        self.assertEqual(records[1].role, "apicomplexan_outgroup")
        self.assertEqual(records[2].role, "distant_outgroup")
        self.assertEqual(records[3].role, "host_or_background")

    def test_load_genome_config_accepts_downloaded_role_when_target_present(self) -> None:
        """Genome config parser should accept generic downloaded non-target rows."""
        with TemporaryDirectory() as tmpdir:
            config_path = Path(tmpdir) / "config.tsv"
            config_path.write_text(
                "genome_fasta\tspecies_name\trole\n"
                "a.fna\tSpecies alpha\ttarget_species\n"
                "b.fna\tSpecies beta\tdownloaded\n",
                encoding="utf-8",
            )
            records = load_genome_config(config_path=config_path)
        self.assertEqual(len(records), 2)
        self.assertEqual(records[1].role, "downloaded")

    def test_load_genome_config_rejects_invalid_role(self) -> None:
        """Genome config parser should reject unsupported roles."""
        with TemporaryDirectory() as tmpdir:
            config_path = Path(tmpdir) / "config.tsv"
            config_path.write_text(
                "genome_fasta\tspecies_name\trole\n"
                "a.fna\tSpecies alpha\tbanana\n",
                encoding="utf-8",
            )
            with self.assertRaises(ValueError):
                load_genome_config(config_path=config_path)

    def test_load_genome_config_skips_excluded_records(self) -> None:
        """Genome config parser should skip records with role exclude."""
        with TemporaryDirectory() as tmpdir:
            config_path = Path(tmpdir) / "config.tsv"
            config_path.write_text(
                "genome_fasta\tspecies_name\trole\n"
                "a.fna\tSpecies alpha\ttarget_species\n"
                "b.fna\tSpecies beta\texclude\n",
                encoding="utf-8",
            )
            records = load_genome_config(config_path=config_path)
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0].species_name, "Species alpha")

    def test_load_genome_config_requires_target(self) -> None:
        """Genome config parser should require at least one target species."""
        with TemporaryDirectory() as tmpdir:
            config_path = Path(tmpdir) / "config.tsv"
            config_path.write_text(
                "genome_fasta\tspecies_name\trole\n"
                "a.fna\tSpecies alpha\toutgroup\n",
                encoding="utf-8",
            )
            with self.assertRaises(ValueError):
                load_genome_config(config_path=config_path)
