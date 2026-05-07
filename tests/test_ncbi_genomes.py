"""Tests for KmerSutra NCBI genome downloader helpers."""

from __future__ import annotations

import argparse
import tempfile
import unittest
from pathlib import Path

from kmersutra import ncbi_genomes


class TestNcbiGenomeDownloader(unittest.TestCase):
    """Test NCBI genome downloader helper functions."""

    def test_safe_name_removes_spaces_and_symbols(self):
        """Filesystem labels should be safe and deterministic."""
        observed = ncbi_genomes.safe_name("Plasmodium vivax / strain A")
        self.assertEqual(observed, "Plasmodium_vivax_strain_A")

    def test_parse_optional_int_handles_empty_values(self):
        """Optional integer parser should keep blank config fields as None."""
        self.assertIsNone(ncbi_genomes.parse_optional_int(""))
        self.assertEqual(ncbi_genomes.parse_optional_int("3"), 3)

    def test_ftp_url_to_file_url_uses_assembly_stem(self):
        """Downloader should construct the expected NCBI assembly file URL."""
        observed = ncbi_genomes.ftp_url_to_file_url(
            "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.29_GRCh38.p14",
            "genomic.fna.gz",
        )
        self.assertTrue(observed.startswith("https://"))
        self.assertTrue(observed.endswith("GCA_000001405.29_GRCh38.p14_genomic.fna.gz"))

    def test_build_taxon_plan_reads_tsv(self):
        """Taxid plan files should produce per-taxon download plans."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            plan_path = Path(tmp_dir) / "plan.tsv"
            plan_path.write_text(
                "taxid\trole\tclade\tgroup_label\tmax_assemblies\tbest_per_species\n"
                "5820\tnear_neighbour\tPlasmodium\tplasmo\t10\t1\n",
                encoding="utf-8",
            )
            args = argparse.Namespace(
                taxid_plan=str(plan_path),
                taxids=None,
                default_role="downloaded",
                default_clade="",
                max_assemblies_per_taxid=None,
                best_per_species=None,
            )
            plans = ncbi_genomes.build_taxon_plan(args)
        self.assertEqual(len(plans), 1)
        self.assertEqual(plans[0].taxid, "5820")
        self.assertEqual(plans[0].role, "near_neighbour")
        self.assertEqual(plans[0].max_assemblies, 10)
        self.assertEqual(plans[0].best_per_species, 1)

    def test_select_best_per_species_uses_quality_order(self):
        """Best-per-species selection should retain the highest-quality assembly."""
        low = ncbi_genomes.AssemblyRecord(
            query_taxid="5820",
            assembly_uid="1",
            assembly_accession="GCA_low",
            assembly_name="low",
            organism_name="Plasmodium test",
            species_name="Plasmodium test",
            species_taxid="1",
            taxid="1",
            strain_name="low",
            assembly_level="contig",
            refseq_category="",
            scaffold_n50=10,
            contig_n50=10,
            total_length=100,
            ftp_path_refseq="",
            ftp_path_genbank="ftp://example/low",
            selected_source="genbank",
            selected_ftp_path="ftp://example/low",
            role="near_neighbour",
            clade="Plasmodium",
            group_label="",
        )
        high = ncbi_genomes.AssemblyRecord(
            query_taxid="5820",
            assembly_uid="2",
            assembly_accession="GCA_high",
            assembly_name="high",
            organism_name="Plasmodium test",
            species_name="Plasmodium test",
            species_taxid="1",
            taxid="1",
            strain_name="high",
            assembly_level="chromosome",
            refseq_category="",
            scaffold_n50=1000,
            contig_n50=1000,
            total_length=10000,
            ftp_path_refseq="",
            ftp_path_genbank="ftp://example/high",
            selected_source="genbank",
            selected_ftp_path="ftp://example/high",
            role="near_neighbour",
            clade="Plasmodium",
            group_label="",
        )
        selected = ncbi_genomes.select_best_per_species([low, high], best_per_species=1)
        self.assertEqual([record.assembly_accession for record in selected], ["GCA_high"])

    def test_metadata_to_kmersutra_config_rows_keeps_downloaded_fastas(self):
        """KmerSutra config rows should include records with usable genome FASTA paths."""
        rows = ncbi_genomes.metadata_to_kmersutra_config_rows(
            [
                {
                    "genome_fasta": "/tmp/a.fna",
                    "species_name": "Plasmodium vivax",
                    "strain_name": "test",
                    "taxid": "5855",
                    "role": "target_species",
                    "clade": "Plasmodium",
                    "assembly_accession": "GCA_test",
                    "query_taxid": "5820",
                    "assembly_level": "chromosome",
                    "scaffold_n50": "1000",
                    "contig_n50": "500",
                },
                {"genome_fasta": "", "species_name": "missing"},
            ]
        )
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["species_name"], "Plasmodium vivax")

    def test_module_import_does_not_require_biopython_configuration(self):
        """Importing the downloader helpers should not require live Entrez access."""
        self.assertTrue(hasattr(ncbi_genomes, "parse_args"))


if __name__ == "__main__":
    unittest.main()
