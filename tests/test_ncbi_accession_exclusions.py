"""Tests for truth-accession exclusions in NCBI genome preparation."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from kmersutra.ncbi_genomes import (
    AssemblyRecord,
    exclude_records_by_accession,
    fasta_excluded_accessions,
    load_excluded_accessions,
)


def assembly(accession: str) -> AssemblyRecord:
    """Build a minimal assembly record."""
    return AssemblyRecord(
        query_taxid="1",
        assembly_uid="uid",
        assembly_accession=accession,
        assembly_name="assembly",
        organism_name="Species alpha",
        species_name="Species alpha",
        species_taxid="1",
        taxid="1",
        strain_name="strain",
        assembly_level="Complete Genome",
        refseq_category="reference genome",
        scaffold_n50=100,
        contig_n50=100,
        total_length=100,
        ftp_path_refseq="",
        ftp_path_genbank="",
        selected_source="refseq",
        selected_ftp_path="",
        role="target_species",
        clade="alpha",
        group_label="alpha",
    )


class TestNcbiAccessionExclusions(unittest.TestCase):
    """Test manifest parsing and candidate filtering."""

    def test_load_exclusions_supports_semicolon_values(self) -> None:
        """Multiple truth accessions in one cell should be loaded."""
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "truth.tsv"
            path.write_text(
                "truth_sequence_accessions\ttruth_assembly_accessions\n"
                "NC_000001.2;NZ_000002.1\tGCF_000003.9\n",
                encoding="utf-8",
            )
            accessions = load_excluded_accessions(path)
        self.assertEqual(
            accessions,
            {"NC_000001", "NZ_000002", "GCF_000003"},
        )

    def test_candidate_assembly_is_excluded_without_version(self) -> None:
        """Assembly-version differences should not permit leakage."""
        audit = [
            {
                "assembly_accession": "GCF_000003.1",
                "filter_status": "retained",
                "filter_reason": "retained",
            }
        ]
        retained = exclude_records_by_accession(
            records=[assembly("GCF_000003.1")],
            audit_rows=audit,
            excluded_accessions={"GCF_000003"},
        )
        self.assertEqual(retained, [])
        self.assertEqual(
            audit[0]["filter_reason"],
            "excluded_by_accession_table",
        )

    def test_fasta_headers_are_checked_for_sequence_exclusions(self) -> None:
        """Sequence accessions should be caught after FASTA download."""
        with tempfile.TemporaryDirectory() as temporary:
            fasta = Path(temporary) / "genome.fna"
            fasta.write_text(">NC_000001.7 chromosome\nACGT\n", encoding="utf-8")
            matches = fasta_excluded_accessions(
                fasta_path=fasta,
                excluded_accessions={"NC_000001"},
            )
        self.assertEqual(matches, {"NC_000001"})


if __name__ == "__main__":
    unittest.main()
