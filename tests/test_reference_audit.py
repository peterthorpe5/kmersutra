"""Tests for reference-panel leakage auditing."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from kmersutra.reference_audit import (
    audit_reference_panel,
    fasta_header_accessions,
    write_reference_panel_audit,
)
from kmersutra.table_io import read_records_table


class TestReferenceAudit(unittest.TestCase):
    """Test exact assembly, sequence and checksum leakage checks."""

    def write_truth(self, root: Path, *, truth_fasta: str = "") -> Path:
        """Write a one-species truth manifest."""
        path = root / "truth.tsv"
        path.write_text(
            "benchmark_id\torganism_id\tspecies_name\texpected\t"
            "expected_abundance_fraction\tabundance_tier\t"
            "truth_assembly_accessions\ttruth_sequence_accessions\ttruth_fasta\n"
            f"example\tone\tSpecies alpha\ttrue\t1.0\thigh\t"
            f"GCF_000001.2\tNC_000001.4\t{truth_fasta}\n",
            encoding="utf-8",
        )
        return path

    def write_config(
        self,
        root: Path,
        *,
        fasta: Path,
        assembly_accession: str,
    ) -> Path:
        """Write a one-row genome configuration."""
        path = root / "genome_config.tsv"
        path.write_text(
            "genome_fasta\tspecies_name\tstrain_name\ttaxid\t"
            "assembly_accession\trole\tclade\tsource\n"
            f"{fasta}\tSpecies alpha\talternative\t1\t"
            f"{assembly_accession}\ttarget_species\talpha\ttest\n",
            encoding="utf-8",
        )
        return path

    def test_fasta_header_accessions_normalises_versions(self) -> None:
        """Header accessions should be returned without version suffixes."""
        with tempfile.TemporaryDirectory() as temporary:
            fasta = Path(temporary) / "panel.fna"
            fasta.write_text(">NC_000001.4 chromosome\nACGT\n", encoding="utf-8")
            accessions = fasta_header_accessions(fasta)
        self.assertIn("NC_000001", accessions)

    def test_sequence_accession_leakage_is_detected(self) -> None:
        """A truth sequence accession in a panel FASTA should fail the audit."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fasta = root / "panel.fna"
            fasta.write_text(">NC_000001.4 chromosome\nACGT\n", encoding="utf-8")
            rows = audit_reference_panel(
                genome_config=self.write_config(
                    root,
                    fasta=fasta,
                    assembly_accession="GCF_999999.1",
                ),
                truth_manifest=self.write_truth(root),
            )
        self.assertEqual(rows[0]["audit_status"], "leakage")
        self.assertEqual(rows[0]["sequence_accession_leakage"], "True")

    def test_same_species_alternative_is_allowed(self) -> None:
        """A held-out same-species assembly should be retained as legitimate."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fasta = root / "panel.fna"
            fasta.write_text(">NZ_OTHER.1 chromosome\nACGT\n", encoding="utf-8")
            rows = audit_reference_panel(
                genome_config=self.write_config(
                    root,
                    fasta=fasta,
                    assembly_accession="GCF_999999.1",
                ),
                truth_manifest=self.write_truth(root),
            )
        self.assertEqual(rows[0]["audit_status"], "same_species_reference")
        self.assertEqual(rows[0]["same_species_reference"], "True")

    def test_assembly_accession_leakage_is_detected_without_version(self) -> None:
        """Assembly versions should not evade truth exclusion."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fasta = root / "panel.fna"
            fasta.write_text(">NZ_OTHER.1 chromosome\nACGT\n", encoding="utf-8")
            rows = audit_reference_panel(
                genome_config=self.write_config(
                    root,
                    fasta=fasta,
                    assembly_accession="GCF_000001.99",
                ),
                truth_manifest=self.write_truth(root),
            )
        self.assertEqual(rows[0]["assembly_accession_leakage"], "True")

    def test_checksum_leakage_is_detected(self) -> None:
        """An exact truth FASTA copy should be detected by SHA-256."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            truth_fasta = root / "truth.fna"
            truth_fasta.write_text(">unlabelled\nACGTACGT\n", encoding="utf-8")
            panel_fasta = root / "panel.fna"
            panel_fasta.write_bytes(truth_fasta.read_bytes())
            rows = audit_reference_panel(
                genome_config=self.write_config(
                    root,
                    fasta=panel_fasta,
                    assembly_accession="GCF_999999.1",
                ),
                truth_manifest=self.write_truth(
                    root,
                    truth_fasta=str(truth_fasta),
                ),
            )
        self.assertEqual(rows[0]["checksum_leakage"], "True")

    def test_fail_on_leakage_writes_audit_before_raising(self) -> None:
        """The audit table should survive a deliberate leakage failure."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fasta = root / "panel.fna"
            fasta.write_text(">NC_000001.4\nACGT\n", encoding="utf-8")
            output = root / "audit.tsv"
            with self.assertRaisesRegex(ValueError, "leakage detected"):
                write_reference_panel_audit(
                    genome_config=self.write_config(
                        root,
                        fasta=fasta,
                        assembly_accession="GCF_999999.1",
                    ),
                    truth_manifest=self.write_truth(root),
                    output_table=output,
                    fail_on_leakage=True,
                )
            rows = read_records_table(input_path=output)
        self.assertEqual(rows[0]["audit_status"], "leakage")


if __name__ == "__main__":
    unittest.main()
