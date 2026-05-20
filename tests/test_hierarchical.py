"""Tests for hierarchical KmerSutra cascade screening."""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from kmersutra.cli.screen_reads_for_clade_kmers import main as screen_main
from kmersutra.hierarchical import (
    gate_summary_passes,
    load_module_manifest,
    order_modules_by_parentage,
    screen_file_hierarchical,
    summarise_gate_hits,
)
from kmersutra.io import read_tsv, write_tsv


PANEL_COLUMNS = [
    "kmer",
    "k",
    "panel_type",
    "species_name",
    "clade",
    "source_genomes",
    "source_contigs",
    "example_position",
    "evidence_taxid",
    "evidence_name",
    "evidence_rank",
    "lineage_taxids",
    "source_taxids",
]

MANIFEST_COLUMNS = [
    "module_id",
    "module_name",
    "rank",
    "parent_module_id",
    "gate_panel_path",
    "module_panel_path",
    "min_gate_unique_kmers",
    "min_gate_positive_sequences",
    "min_gate_k_values",
    "min_gate_best_k",
]


class TestHierarchicalScreening(unittest.TestCase):
    """Tests for hierarchical module activation and screening."""

    def _write_panel(
        self,
        *,
        path: Path,
        kmer: str,
        species_name: str,
        evidence_name: str,
        evidence_rank: str,
    ) -> None:
        """Write one minimal KmerSutra panel."""
        write_tsv(
            records=[
                {
                    "kmer": kmer,
                    "k": len(kmer),
                    "panel_type": "species_unique" if species_name else "clade_core",
                    "species_name": species_name,
                    "clade": "Apicomplexa",
                    "source_genomes": "g1",
                    "source_contigs": "c1",
                    "example_position": 0,
                    "evidence_taxid": "1",
                    "evidence_name": evidence_name,
                    "evidence_rank": evidence_rank,
                    "lineage_taxids": "1",
                    "source_taxids": "1",
                }
            ],
            output_path=path,
            fieldnames=PANEL_COLUMNS,
        )

    def _write_basic_inputs(self, tmpdir: str) -> tuple[Path, Path, Path, Path]:
        """Write gate panel, module panel, manifest and FASTQ."""
        base = Path(tmpdir)
        gate_panel = base / "gate.tsv"
        module_panel = base / "module.tsv"
        manifest = base / "modules.tsv"
        fastq = base / "reads.fastq"
        self._write_panel(
            path=gate_panel,
            kmer="AAAAA",
            species_name="",
            evidence_name="Apicomplexa",
            evidence_rank="phylum",
        )
        self._write_panel(
            path=module_panel,
            kmer="CCCCC",
            species_name="Plasmodium vivax",
            evidence_name="Plasmodium vivax",
            evidence_rank="species",
        )
        write_tsv(
            records=[
                {
                    "module_id": "apicomplexa",
                    "module_name": "Apicomplexa",
                    "rank": "phylum",
                    "parent_module_id": "",
                    "gate_panel_path": str(gate_panel),
                    "module_panel_path": str(module_panel),
                    "min_gate_unique_kmers": 1,
                    "min_gate_positive_sequences": 1,
                    "min_gate_k_values": 1,
                    "min_gate_best_k": 5,
                }
            ],
            output_path=manifest,
            fieldnames=MANIFEST_COLUMNS,
        )
        fastq.write_text("@read1\nGGGAAAAATTTCCCCC\n+\nFFFFFFFFFFFFFFFF\n", encoding="utf-8")
        return gate_panel, module_panel, manifest, fastq

    def test_load_module_manifest_resolves_modules(self) -> None:
        """The module manifest should resolve gate and module panels."""
        with tempfile.TemporaryDirectory() as tmpdir:
            _, _, manifest, _ = self._write_basic_inputs(tmpdir)
            modules = load_module_manifest(manifest_path=manifest)
            self.assertEqual(len(modules), 1)
            self.assertEqual(modules[0].module_id, "apicomplexa")
            self.assertTrue(Path(modules[0].gate_panel_path).is_file())
            self.assertTrue(Path(modules[0].module_panel_path).is_file())

    def test_gate_summary_thresholds(self) -> None:
        """Gate summaries should pass only when all thresholds are met."""
        with tempfile.TemporaryDirectory() as tmpdir:
            _, _, manifest, fastq = self._write_basic_inputs(tmpdir)
            result = screen_file_hierarchical(
                input_path=fastq,
                module_manifest_path=manifest,
                sample_id="sample1",
                input_format="fastq",
            )
            modules = load_module_manifest(manifest_path=manifest)
            gate_hits = [hit for hit in result.hits if hit.evidence_rank == "phylum"]
            summary = summarise_gate_hits(hits=gate_hits)
            self.assertTrue(gate_summary_passes(module=modules[0], summary=summary))

    def test_hierarchical_screen_activates_module_after_gate_hit(self) -> None:
        """A passing gate should activate the detailed module panel."""
        with tempfile.TemporaryDirectory() as tmpdir:
            _, _, manifest, fastq = self._write_basic_inputs(tmpdir)
            result = screen_file_hierarchical(
                input_path=fastq,
                module_manifest_path=manifest,
                sample_id="sample1",
                input_format="fastq",
            )
            species_hits = [hit for hit in result.hits if hit.species_name == "Plasmodium vivax"]
            self.assertEqual(len(species_hits), 1)
            self.assertTrue(any(row["activated"] == "True" for row in result.activation_records))

    def test_hierarchical_screen_skips_module_without_gate_hit(self) -> None:
        """A module panel should not be screened when its gate is absent."""
        with tempfile.TemporaryDirectory() as tmpdir:
            _, _, manifest, _ = self._write_basic_inputs(tmpdir)
            fastq = Path(tmpdir) / "no_gate.fastq"
            fastq.write_text("@read1\nGGGCCCCCGGGCCCCC\n+\nFFFFFFFFFFFFFFFF\n", encoding="utf-8")
            result = screen_file_hierarchical(
                input_path=fastq,
                module_manifest_path=manifest,
                sample_id="sample1",
                input_format="fastq",
            )
            self.assertFalse([hit for hit in result.hits if hit.species_name == "Plasmodium vivax"])
            self.assertEqual(result.activation_records[0]["activation_reason"], "gate_below_threshold")

    def test_hierarchical_fail_open_screens_weak_gate_modules(self) -> None:
        """Fail-open mode should screen modules with weak gate evidence."""
        with tempfile.TemporaryDirectory() as tmpdir:
            _, _, manifest, fastq = self._write_basic_inputs(tmpdir)
            records = read_tsv(input_path=manifest)
            records[0]["min_gate_unique_kmers"] = "2"
            write_tsv(records=records, output_path=manifest, fieldnames=MANIFEST_COLUMNS)
            result = screen_file_hierarchical(
                input_path=fastq,
                module_manifest_path=manifest,
                sample_id="sample1",
                input_format="fastq",
                hierarchical_fail_open=True,
            )
            self.assertTrue([hit for hit in result.hits if hit.species_name == "Plasmodium vivax"])
            self.assertTrue(
                any(row["activation_reason"] == "fail_open_weak_gate_signal" for row in result.activation_records)
            )


    def test_parent_child_order_is_independent_of_manifest_order(self) -> None:
        """Hierarchical screening should process parents before children."""
        with tempfile.TemporaryDirectory() as tmpdir:
            base = Path(tmpdir)
            root_gate = base / "root_gate.tsv"
            root_panel = base / "root_panel.tsv"
            child_gate = base / "child_gate.tsv"
            child_panel = base / "child_panel.tsv"
            manifest = base / "modules.tsv"
            fastq = base / "reads.fastq"
            self._write_panel(
                path=root_gate,
                kmer="AAAAA",
                species_name="",
                evidence_name="Apicomplexa",
                evidence_rank="phylum",
            )
            self._write_panel(
                path=root_panel,
                kmer="GGGGG",
                species_name="",
                evidence_name="Apicomplexa",
                evidence_rank="phylum",
            )
            self._write_panel(
                path=child_gate,
                kmer="CCCCC",
                species_name="",
                evidence_name="Plasmodium",
                evidence_rank="genus",
            )
            self._write_panel(
                path=child_panel,
                kmer="ATATA",
                species_name="Plasmodium vivax",
                evidence_name="Plasmodium vivax",
                evidence_rank="species",
            )
            write_tsv(
                records=[
                    {
                        "module_id": "plasmodium",
                        "module_name": "Plasmodium",
                        "rank": "genus",
                        "parent_module_id": "apicomplexa",
                        "gate_panel_path": str(child_gate),
                        "module_panel_path": str(child_panel),
                        "min_gate_unique_kmers": 1,
                        "min_gate_positive_sequences": 1,
                        "min_gate_k_values": 1,
                        "min_gate_best_k": 5,
                    },
                    {
                        "module_id": "apicomplexa",
                        "module_name": "Apicomplexa",
                        "rank": "phylum",
                        "parent_module_id": "",
                        "gate_panel_path": str(root_gate),
                        "module_panel_path": str(root_panel),
                        "min_gate_unique_kmers": 1,
                        "min_gate_positive_sequences": 1,
                        "min_gate_k_values": 1,
                        "min_gate_best_k": 5,
                    },
                ],
                output_path=manifest,
                fieldnames=MANIFEST_COLUMNS,
            )
            fastq.write_text(
                "@read1\nGGGAAAAACCCCCTATATGGG\n+\nFFFFFFFFFFFFFFFFFFFFF\n",
                encoding="utf-8",
            )
            modules = load_module_manifest(manifest_path=manifest)
            ordered = order_modules_by_parentage(modules=modules)
            self.assertEqual([module.module_id for module in ordered], ["apicomplexa", "plasmodium"])
            result = screen_file_hierarchical(
                input_path=fastq,
                module_manifest_path=manifest,
                sample_id="sample1",
                input_format="fastq",
            )
            self.assertTrue([hit for hit in result.hits if hit.species_name == "Plasmodium vivax"])

    def test_screen_cli_hierarchical_writes_activation_report(self) -> None:
        """The CLI should write module activation records in hierarchical mode."""
        with tempfile.TemporaryDirectory() as tmpdir:
            _, _, manifest, fastq = self._write_basic_inputs(tmpdir)
            out_dir = Path(tmpdir) / "out"
            argv = [
                "kmersutra-screen",
                "--input",
                str(fastq),
                "--input_format",
                "fastq",
                "--module_manifest",
                str(manifest),
                "--sample_id",
                "sample1",
                "--out_dir",
                str(out_dir),
                "--threads",
                "1",
                "--chunk_size",
                "1",
                "--no_read_level_hits",
            ]
            with patch.object(sys, "argv", argv):
                screen_main()
            activation = out_dir / "module_activation.tsv"
            self.assertTrue(activation.exists())
            text = activation.read_text(encoding="utf-8")
            self.assertIn("apicomplexa", text)
            calls = (out_dir / "species_detection_calls.tsv").read_text(encoding="utf-8")
            self.assertIn("Plasmodium vivax", calls)

    def test_screen_cli_flat_mode_requires_panel(self) -> None:
        """Flat mode should fail clearly without a panel."""
        with tempfile.TemporaryDirectory() as tmpdir:
            fastq = Path(tmpdir) / "reads.fastq"
            fastq.write_text("@read1\nAAAAA\n+\nFFFFF\n", encoding="utf-8")
            argv = [
                "kmersutra-screen",
                "--screen_mode",
                "flat",
                "--input",
                str(fastq),
                "--input_format",
                "fastq",
                "--sample_id",
                "sample1",
                "--out_dir",
                str(Path(tmpdir) / "out"),
            ]
            with patch.object(sys, "argv", argv):
                with self.assertRaises(ValueError):
                    screen_main()


if __name__ == "__main__":
    unittest.main()
