"""Tests for leakage-controlled ATCC reference preparation."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from kmersutra.atcc_reference import (
    create_taxid_plan,
    evaluate_reference_gate,
    finalise_reference_config,
)
from kmersutra.reference_audit import write_reference_panel_audit
from kmersutra.table_io import read_records_table
from kmersutra.taxonomy import TaxonomyDatabase


class TestAtccReference(unittest.TestCase):
    """Validate ATCC plan construction, finalisation and gating."""

    def setUp(self) -> None:
        """Create a minimal taxonomy, truth manifest and reference FASTAs."""
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary_directory.name)
        self.taxonomy_dir = self.root / "taxonomy"
        self.taxonomy_dir.mkdir()
        (self.taxonomy_dir / "nodes.dmp").write_text(
            "1\t|\t1\t|\tno rank\t|\n"
            "2\t|\t1\t|\tsuperkingdom\t|\n"
            "10\t|\t2\t|\tgenus\t|\n"
            "11\t|\t10\t|\tspecies\t|\n"
            "12\t|\t10\t|\tspecies\t|\n"
            "20\t|\t2\t|\tgenus\t|\n"
            "21\t|\t20\t|\tspecies\t|\n"
            "22\t|\t20\t|\tspecies\t|\n",
            encoding="utf-8",
        )
        (self.taxonomy_dir / "names.dmp").write_text(
            "1\t|\troot\t|\t\t|\tscientific name\t|\n"
            "2\t|\tBacteria\t|\t\t|\tscientific name\t|\n"
            "10\t|\tAlpha\t|\t\t|\tscientific name\t|\n"
            "11\t|\tAlpha one\t|\t\t|\tscientific name\t|\n"
            "12\t|\tAlpha neighbour\t|\t\t|\tscientific name\t|\n"
            "20\t|\tBeta\t|\t\t|\tscientific name\t|\n"
            "21\t|\tBeta two\t|\t\t|\tscientific name\t|\n"
            "22\t|\tBeta neighbour\t|\t\t|\tscientific name\t|\n",
            encoding="utf-8",
        )
        (self.taxonomy_dir / "merged.dmp").write_text("", encoding="utf-8")
        (self.taxonomy_dir / "delnodes.dmp").write_text("", encoding="utf-8")
        self.taxonomy = TaxonomyDatabase.from_taxdump(
            taxonomy_dir=self.taxonomy_dir
        )
        self.truth_manifest = self.root / "truth.tsv"
        self.truth_manifest.write_text(
            "benchmark_id\torganism_id\tspecies_name\taccepted_species_names\t"
            "expected\texpected_abundance_fraction\tabundance_tier\tncbi_taxid\t"
            "truth_strain\ttruth_assembly_accessions\ttruth_sequence_accessions\t"
            "truth_fasta\ttruth_fasta_sha256\tsource\n"
            "demo\ta\tAlpha one\t\ttrue\t0.5\thigh\t11\ttruth_a\t\tNC_123456.1\t\t\tdemo\n"
            "demo\tb\tBeta two\t\ttrue\t0.5\thigh\t21\ttruth_b\t\tNC_234567.1\t\t\tdemo\n",
            encoding="utf-8",
        )

    def tearDown(self) -> None:
        """Remove temporary test data."""
        self.temporary_directory.cleanup()

    def _fasta(self, name: str, header: str) -> Path:
        """Create a small FASTA file.

        Args:
            name: Output filename.
            header: FASTA record identifier.

        Returns:
            Created FASTA path.
        """
        path = self.root / name
        path.write_text(f">{header}\nACGTACGTACGT\n", encoding="utf-8")
        return path

    def test_plan_contains_targets_and_unique_genera(self) -> None:
        """The plan should contain two targets but only two genus queries."""
        output = self.root / "plan.tsv"
        rows = create_taxid_plan(
            truth_manifest=self.truth_manifest,
            taxonomy=self.taxonomy,
            output_path=output,
            target_assemblies=4,
            near_neighbour_assemblies_per_genus=7,
        )
        target_rows = [row for row in rows if row["role"] == "target_species"]
        neighbour_rows = [row for row in rows if row["role"] == "near_neighbour"]
        self.assertEqual({row["taxid"] for row in target_rows}, {"11", "21"})
        self.assertEqual({row["taxid"] for row in neighbour_rows}, {"10", "20"})
        self.assertTrue(all(row["max_assemblies"] == 4 for row in target_rows))
        self.assertTrue(all(row["max_assemblies"] == 7 for row in neighbour_rows))

    def test_finalise_filters_target_neighbours_and_background_duplicates(self) -> None:
        """Only explicit held-out downloads should represent ATCC targets."""
        target_a = self._fasta("target_a.fna", "ALT_A")
        target_b = self._fasta("target_b.fna", "ALT_B")
        neighbour = self._fasta("neighbour.fna", "NEIGHBOUR")
        background = self._fasta("background.fna", "BACKGROUND")
        background_target = self._fasta("background_target.fna", "OLD_TARGET")
        downloaded = self.root / "downloaded.tsv"
        downloaded.write_text(
            "genome_fasta\tspecies_name\tstrain_name\ttaxid\tassembly_accession\t"
            "role\tclade\n"
            f"{target_a}\tAlpha one\talt_a\t11\tGCF_000001.1\t"
            "target_species\tATCC\n"
            f"{target_b}\tBeta two\talt_b\t21\tGCF_000002.1\t"
            "target_species\tATCC\n"
            f"{target_a}\tAlpha one\tduplicate\t11\tGCF_000001.1\t"
            "near_neighbour\tAlpha\n"
            f"{neighbour}\tAlpha neighbour\tnear\t12\tGCF_000003.1\t"
            "near_neighbour\tAlpha\n",
            encoding="utf-8",
        )
        background_config = self.root / "background.tsv"
        background_config.write_text(
            "genome_fasta\tspecies_name\tstrain_name\ttaxid\tassembly_accession\t"
            "role\tclade\n"
            f"{background}\tBackground organism\tbg\t999\tGCF_000004.1\t"
            "outgroup\tBackground\n"
            f"{background_target}\tAlpha one\told\t11\tGCF_000005.1\t"
            "target_species\tLegacy\n",
            encoding="utf-8",
        )
        output_config = self.root / "final.tsv"
        coverage = self.root / "coverage.tsv"
        rows, coverage_rows = finalise_reference_config(
            downloaded_config=downloaded,
            background_config=background_config,
            truth_manifest=self.truth_manifest,
            taxonomy=self.taxonomy,
            output_config=output_config,
            coverage_table=coverage,
        )
        self.assertEqual(len(rows), 4)
        self.assertEqual(
            {row["assembly_accession"] for row in rows},
            {"GCF_000001.1", "GCF_000002.1", "GCF_000003.1", "GCF_000004.1"},
        )
        self.assertEqual(
            {row["reference_count"] for row in coverage_rows},
            {1},
        )

        audit = self.root / "audit.tsv"
        write_reference_panel_audit(
            genome_config=output_config,
            truth_manifest=self.truth_manifest,
            output_table=audit,
            fail_on_leakage=True,
        )
        gate_coverage = self.root / "gate_coverage.tsv"
        gate_summary = self.root / "gate_summary.tsv"
        self.assertTrue(
            evaluate_reference_gate(
                audit_table=audit,
                truth_manifest=self.truth_manifest,
                coverage_table=gate_coverage,
                summary_table=gate_summary,
            )
        )
        metrics = {
            row["metric"]: row["value"]
            for row in read_records_table(input_path=gate_summary)
        }
        self.assertEqual(metrics["represented_species"], "2")
        self.assertEqual(metrics["gate_status"], "PASS")

    def test_gate_blocks_truth_sequence_leakage(self) -> None:
        """A truth accession in a FASTA header must block the gate."""
        leaking_fasta = self._fasta("leaking.fna", "NC_123456.1")
        safe_fasta = self._fasta("safe.fna", "ALT_B")
        config = self.root / "leaking_config.tsv"
        config.write_text(
            "genome_fasta\tspecies_name\tassembly_accession\trole\n"
            f"{leaking_fasta}\tAlpha one\tGCF_100001.1\ttarget_species\n"
            f"{safe_fasta}\tBeta two\tGCF_100002.1\ttarget_species\n",
            encoding="utf-8",
        )
        audit = self.root / "leaking_audit.tsv"
        write_reference_panel_audit(
            genome_config=config,
            truth_manifest=self.truth_manifest,
            output_table=audit,
        )
        coverage = self.root / "leaking_coverage.tsv"
        summary = self.root / "leaking_summary.tsv"
        self.assertFalse(
            evaluate_reference_gate(
                audit_table=audit,
                truth_manifest=self.truth_manifest,
                coverage_table=coverage,
                summary_table=summary,
            )
        )
        metrics = {
            row["metric"]: row["value"]
            for row in read_records_table(input_path=summary)
        }
        self.assertEqual(metrics["leakage_rows"], "1")
        self.assertEqual(metrics["gate_status"], "BLOCKED")


if __name__ == "__main__":
    unittest.main()
