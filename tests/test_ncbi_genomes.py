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


    def test_build_taxon_plan_reads_quality_filters(self):
        """Taxid plan files should parse per-taxon assembly quality filters."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            plan_path = Path(tmp_dir) / "plan.tsv"
            plan_path.write_text(
                "taxid\trole\tclade\tgroup_label\tmax_assemblies\tbest_per_species\t"
                "min_total_length\tmax_total_length\tmin_scaffold_n50\tmin_contig_n50\n"
                "5911\tdistant_outgroup\tCiliophora\tTetrahymena\t2\t\t10000000\t50000000\t500000\t100000\n",
                encoding="utf-8",
            )
            args = argparse.Namespace(
                taxid_plan=str(plan_path),
                taxids=None,
                default_role="downloaded",
                default_clade="",
                max_assemblies_per_taxid=None,
                best_per_species=None,
                min_total_length=None,
                max_total_length=None,
                min_scaffold_n50=None,
                min_contig_n50=None,
            )
            plans = ncbi_genomes.build_taxon_plan(args)
        self.assertEqual(plans[0].min_total_length, 10000000)
        self.assertEqual(plans[0].max_total_length, 50000000)
        self.assertEqual(plans[0].min_scaffold_n50, 500000)
        self.assertEqual(plans[0].min_contig_n50, 100000)

    def test_filter_records_applies_length_and_n50_filters(self):
        """Assembly filters should remove very small or fragmented records."""
        small = ncbi_genomes.AssemblyRecord(
            query_taxid="5911",
            assembly_uid="1",
            assembly_accession="GCA_small",
            assembly_name="small",
            organism_name="Tetrahymena thermophila",
            species_name="Tetrahymena thermophila",
            species_taxid="5911",
            taxid="5911",
            strain_name="small",
            assembly_level="scaffold",
            refseq_category="",
            scaffold_n50=1000,
            contig_n50=1000,
            total_length=900000,
            ftp_path_refseq="",
            ftp_path_genbank="ftp://example/small",
            selected_source="genbank",
            selected_ftp_path="ftp://example/small",
            role="distant_outgroup",
            clade="Ciliophora",
            group_label="Tetrahymena",
        )
        large = ncbi_genomes.AssemblyRecord(
            query_taxid="5911",
            assembly_uid="2",
            assembly_accession="GCA_large",
            assembly_name="large",
            organism_name="Tetrahymena thermophila",
            species_name="Tetrahymena thermophila",
            species_taxid="5911",
            taxid="5911",
            strain_name="large",
            assembly_level="chromosome",
            refseq_category="",
            scaffold_n50=20000000,
            contig_n50=20000000,
            total_length=40000000,
            ftp_path_refseq="",
            ftp_path_genbank="ftp://example/large",
            selected_source="genbank",
            selected_ftp_path="ftp://example/large",
            role="distant_outgroup",
            clade="Ciliophora",
            group_label="Tetrahymena",
        )
        observed = ncbi_genomes.filter_records(
            records=[small, large],
            min_total_length=10000000,
            min_scaffold_n50=500000,
            min_contig_n50=500000,
        )
        self.assertEqual([record.assembly_accession for record in observed], ["GCA_large"])

    def test_filter_records_applies_max_total_length(self):
        """Assembly filters should optionally exclude unexpectedly large records."""
        record = ncbi_genomes.AssemblyRecord(
            query_taxid="1",
            assembly_uid="1",
            assembly_accession="GCA_large",
            assembly_name="large",
            organism_name="Species large",
            species_name="Species large",
            species_taxid="1",
            taxid="1",
            strain_name="",
            assembly_level="chromosome",
            refseq_category="",
            scaffold_n50=1000,
            contig_n50=1000,
            total_length=999999999,
            ftp_path_refseq="",
            ftp_path_genbank="ftp://example/large",
            selected_source="genbank",
            selected_ftp_path="ftp://example/large",
            role="outgroup",
            clade="Demo",
            group_label="",
        )
        observed = ncbi_genomes.filter_records(
            records=[record],
            max_total_length=1000000,
        )
        self.assertEqual(observed, [])

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


    def test_assembly_level_filter_is_case_insensitive(self):
        """Assembly-level filtering should accept NCBI case variants."""
        complete = ncbi_genomes.AssemblyRecord(
            query_taxid="1",
            assembly_uid="1",
            assembly_accession="GCA_complete",
            assembly_name="complete",
            organism_name="Example complete",
            species_name="Example complete",
            species_taxid="1",
            taxid="1",
            strain_name="",
            assembly_level="Complete Genome",
            refseq_category="",
            scaffold_n50=100,
            contig_n50=100,
            total_length=1000,
            ftp_path_refseq="ftp://example/complete",
            ftp_path_genbank="",
            selected_source="refseq",
            selected_ftp_path="ftp://example/complete",
            role="target",
            clade="demo",
            group_label="",
        )
        scaffold = ncbi_genomes.AssemblyRecord(
            query_taxid="1",
            assembly_uid="2",
            assembly_accession="GCA_scaffold",
            assembly_name="scaffold",
            organism_name="Example scaffold",
            species_name="Example scaffold",
            species_taxid="1",
            taxid="1",
            strain_name="",
            assembly_level="Scaffold",
            refseq_category="",
            scaffold_n50=100,
            contig_n50=100,
            total_length=1000,
            ftp_path_refseq="ftp://example/scaffold",
            ftp_path_genbank="",
            selected_source="refseq",
            selected_ftp_path="ftp://example/scaffold",
            role="target",
            clade="demo",
            group_label="",
        )
        observed = ncbi_genomes.filter_records(
            records=[complete, scaffold],
            assembly_levels=["complete genome", "scaffold"],
        )
        self.assertEqual(
            [record.assembly_accession for record in observed],
            ["GCA_complete", "GCA_scaffold"],
        )

    def test_candidate_filter_reason_reports_assembly_level(self):
        """Candidate audit should explain why a record was filtered out."""
        record = ncbi_genomes.AssemblyRecord(
            query_taxid="1",
            assembly_uid="1",
            assembly_accession="GCA_contig",
            assembly_name="contig",
            organism_name="Example contig",
            species_name="Example contig",
            species_taxid="1",
            taxid="1",
            strain_name="",
            assembly_level="Contig",
            refseq_category="",
            scaffold_n50=100,
            contig_n50=100,
            total_length=1000,
            ftp_path_refseq="ftp://example/contig",
            ftp_path_genbank="",
            selected_source="refseq",
            selected_ftp_path="ftp://example/contig",
            role="target",
            clade="demo",
            group_label="",
        )
        reason = ncbi_genomes.evaluate_record_filter_reason(
            record=record,
            assembly_levels=["complete genome"],
        )
        self.assertEqual(reason, "excluded_by_assembly_level")

    def test_prefer_refseq_falls_back_to_genbank_path(self):
        """Source preference should fall back to GenBank when RefSeq is absent."""
        source, path = ncbi_genomes.choose_ftp_path(
            summary={
                "FtpPath_RefSeq": "",
                "FtpPath_GenBank": "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/demo",
            },
            source="prefer_refseq",
        )
        self.assertEqual(source, "genbank")
        self.assertIn("GCA/demo", path)

    def test_entrez_read_retry_recovers_from_incomplete_read(self):
        """Entrez read retry should rerun the whole request after IncompleteRead."""
        import http.client

        class FakeEntrez:
            def __init__(self):
                self.calls = 0

            def read(self, handle, validate=False):
                self.calls += 1
                if self.calls == 1:
                    raise http.client.IncompleteRead(b"partial")
                return {"ok": "yes"}

        fake_entrez = FakeEntrez()
        old_entrez = ncbi_genomes.Entrez
        ncbi_genomes.Entrez = fake_entrez
        try:
            result = ncbi_genomes.entrez_read_retry(
                lambda: object(),
                retries=2,
                sleep_seconds=0,
                context="mock esummary batch",
            )
        finally:
            ncbi_genomes.Entrez = old_entrez
        self.assertEqual(result, {"ok": "yes"})
        self.assertEqual(fake_entrez.calls, 2)

    def test_zero_retained_case_writes_candidate_audit(self):
        """A zero-retained taxid should still produce a useful audit table."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_dir = Path(tmp_dir) / "out"
            old_configure = ncbi_genomes.configure_entrez
            old_search = ncbi_genomes.search_assembly_uids
            old_fetch = ncbi_genomes.fetch_assembly_summaries
            try:
                ncbi_genomes.configure_entrez = lambda email, api_key=None: None
                ncbi_genomes.search_assembly_uids = lambda **kwargs: ["1"]
                ncbi_genomes.fetch_assembly_summaries = lambda **kwargs: [
                    {
                        "AssemblyAccession": "GCA_missing_path",
                        "AssemblyName": "missing_path",
                        "Organism": "Example species",
                        "SpeciesName": "Example species",
                        "SpeciesTaxid": "123",
                        "Taxid": "123",
                        "AssemblyLevel": "Complete Genome",
                        "FtpPath_RefSeq": "",
                        "FtpPath_GenBank": "",
                    }
                ]
                status = ncbi_genomes.main(
                    [
                        "--taxids",
                        "123",
                        "--email",
                        "test@example.com",
                        "--out_dir",
                        str(out_dir),
                        "--metadata_only",
                    ]
                )
            finally:
                ncbi_genomes.configure_entrez = old_configure
                ncbi_genomes.search_assembly_uids = old_search
                ncbi_genomes.fetch_assembly_summaries = old_fetch
            self.assertEqual(status, 0)
            audit_path = out_dir / "ncbi_assembly_candidate_audit.tsv"
            config_path = out_dir / "kmersutra_genome_config.tsv"
            audit_text = audit_path.read_text(encoding="utf-8")
            self.assertIn("excluded_missing_download_path", audit_text)
            self.assertEqual(len(config_path.read_text(encoding="utf-8").splitlines()), 1)

    def test_custom_genome_config_can_be_loaded_without_entrez(self):
        """A local genome config should allow custom genomes without NCBI access."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)
            fasta = tmp_path / "custom.fna"
            fasta.write_text(">custom\nACGT\n", encoding="utf-8")
            config = tmp_path / "custom_config.tsv"
            config.write_text(
                "genome_fasta\tspecies_name\tstrain_name\ttaxid\trole\tclade\t"
                "assembly_accession\tquery_taxid\tassembly_level\tscaffold_n50\tcontig_n50\n"
                f"{fasta}\tCustom species\tstrain1\t999\tcustom\tCustom\tCUSTOM1\t999\tcustom\t\t\n",
                encoding="utf-8",
            )
            out_dir = tmp_path / "out"
            status = ncbi_genomes.main(
                [
                    "--custom_genome_config",
                    str(config),
                    "--out_dir",
                    str(out_dir),
                ]
            )
            self.assertEqual(status, 0)
            config_text = (out_dir / "kmersutra_genome_config.tsv").read_text(
                encoding="utf-8"
            )
            self.assertIn("Custom species", config_text)


if __name__ == "__main__":
    unittest.main()
