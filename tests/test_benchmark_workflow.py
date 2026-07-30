"""Tests for the restartable benchmark controller."""

from __future__ import annotations

import gzip
import json
import logging
import tempfile
import unittest
from pathlib import Path

from kmersutra.benchmark_workflow import (
    LOCKED_ATCC_BENCHMARK_ID,
    STAGES,
    BenchmarkWorkflow,
    StageResult,
    WorkflowOptions,
    filter_panel_by_k,
    load_benchmark_config,
    selected_stages,
    validate_locked_atcc_config,
)
from kmersutra.io import write_tsv


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


def locked_config() -> dict[str, object]:
    """Return the minimal locked settings used by validation tests."""
    return {
        "benchmark_id": LOCKED_ATCC_BENCHMARK_ID,
        "allow_primary_threshold_tuning": False,
        "dataset": {
            "sra_accession": "SRR9328980",
            "input_format": "fastq",
        },
        "screen": {
            "screen_preset": "exact",
            "max_mismatches": 0,
            "call_preset": "lineage_aware",
            "same_genus_reportable_min_fraction": 0.05,
            "k_values": [51, 77, 101, 151],
        },
        "depth": {
            "fractions": [0.01],
            "seeds": [1001],
        },
        "ai": {
            "enabled": False,
            "novelty_scale": 2.9,
        },
        "reference": {},
    }


class TestBenchmarkWorkflowHelpers(unittest.TestCase):
    """Test configuration, range and panel-filter helpers."""

    def test_locked_atcc_configuration_passes(self) -> None:
        """The registered settings should pass validation."""
        validate_locked_atcc_config(config=locked_config())

    def test_locked_atcc_rejects_threshold_tuning(self) -> None:
        """A primary tuning declaration should invalidate the lock."""
        config = locked_config()
        config["allow_primary_threshold_tuning"] = True
        with self.assertRaisesRegex(ValueError, "threshold tuning"):
            validate_locked_atcc_config(config=config)

    def test_locked_atcc_rejects_changed_k_values(self) -> None:
        """Changing the k ladder should invalidate the registered design."""
        config = locked_config()
        config["screen"]["k_values"] = [51, 77, 101]  # type: ignore[index]
        with self.assertRaisesRegex(ValueError, "k values"):
            validate_locked_atcc_config(config=config)

    def test_selected_stages_resolves_inclusive_range(self) -> None:
        """Stage ranges should include both named endpoints."""
        stages = selected_stages(
            start_at="03_screen_full",
            stop_after="05_screen_single_k",
        )
        self.assertEqual(
            stages,
            [
                "03_screen_full",
                "04_screen_depths",
                "05_screen_single_k",
            ],
        )

    def test_selected_stages_rejects_reverse_range(self) -> None:
        """A reversed stage range should fail before work begins."""
        with self.assertRaisesRegex(ValueError, "occurs after"):
            selected_stages(
                start_at="08_summarise",
                stop_after="03_screen_full",
            )

    def test_load_config_expands_repo_and_config_placeholders(self) -> None:
        """Configuration placeholders should be resolved deterministically."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config_path = root / "config.json"
            config_path.write_text(
                json.dumps(
                    {
                        "from_config": "${CONFIG_DIR}/file.tsv",
                        "from_repo": "${REPO_ROOT}/README.md",
                    }
                ),
                encoding="utf-8",
            )
            config, digest = load_benchmark_config(config_path=config_path)
        self.assertEqual(config["from_config"], str(root / "file.tsv"))
        self.assertTrue(str(config["from_repo"]).endswith("/README.md"))
        self.assertEqual(len(digest), 64)

    def test_filter_panel_by_k_supports_gzip(self) -> None:
        """Single-k filtering should retain only the requested panel rows."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "panel.tsv.gz"
            with gzip.open(source, "wt", encoding="utf-8") as handle:
                handle.write("k\tkmer\tspecies_name\n")
                handle.write("51\tAAAA\tSpecies alpha\n")
                handle.write("77\tCCCC\tSpecies alpha\n")
            output = root / "panel_k51.tsv.gz"
            retained = filter_panel_by_k(
                panel_path=source,
                output_path=output,
                k_value=51,
            )
            with gzip.open(output, "rt", encoding="utf-8") as handle:
                text = handle.read()
        self.assertEqual(retained, 1)
        self.assertIn("51\tAAAA", text)
        self.assertNotIn("77\tCCCC", text)

    def test_filter_panel_rejects_absent_k(self) -> None:
        """An empty single-k ablation should fail rather than appear complete."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "panel.tsv"
            source.write_text("k\tkmer\n51\tAAAA\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "no rows"):
                filter_panel_by_k(
                    panel_path=source,
                    output_path=root / "out.tsv",
                    k_value=77,
                )


class TestBenchmarkStageState(unittest.TestCase):
    """Test completion-state validation independently of scientific stages."""

    def build_workflow(self, root: Path) -> BenchmarkWorkflow:
        """Create a workflow with a temporary run root."""
        options = WorkflowOptions(
            output_root=root,
            run_name="run",
            threads=1,
            resume=True,
            start_at=None,
            stop_after=None,
            force_stages=frozenset(),
            dry_run=False,
        )
        return BenchmarkWorkflow(
            config=locked_config(),
            config_digest="a" * 64,
            options=options,
            logger=logging.getLogger("test"),
        )

    def test_stage_completion_requires_non_empty_output(self) -> None:
        """A success token must not validate an empty output file."""
        with tempfile.TemporaryDirectory() as temporary:
            workflow = self.build_workflow(Path(temporary))
            workflow.prepare_run_root()
            output = workflow.run_root / "empty.tsv"
            output.touch()
            workflow.write_stage_state(
                stage="00_preflight",
                status="success",
                started_at="2026-07-30T00:00:00+00:00",
                outputs=[output],
            )
            self.assertFalse(workflow.stage_is_complete("00_preflight"))

    def test_stage_completion_requires_current_digest(self) -> None:
        """A successful output from another configuration must not be skipped."""
        with tempfile.TemporaryDirectory() as temporary:
            workflow = self.build_workflow(Path(temporary))
            workflow.prepare_run_root()
            output = workflow.run_root / "output.tsv"
            output.write_text("header\n", encoding="utf-8")
            workflow.write_stage_state(
                stage="00_preflight",
                status="success",
                started_at="2026-07-30T00:00:00+00:00",
                outputs=[output],
            )
            self.assertTrue(workflow.stage_is_complete("00_preflight"))
            workflow.config_digest = "b" * 64
            self.assertFalse(workflow.stage_is_complete("00_preflight"))

    def test_stage_result_requires_declared_output(self) -> None:
        """The stage result dataclass should retain declared output paths."""
        path = Path("/tmp/example.tsv")
        result = StageResult(outputs=(path,), detail="test")
        self.assertEqual(result.outputs, (path,))


class TestBenchmarkWorkflowIntegration(unittest.TestCase):
    """Exercise the complete controller with a tiny independent benchmark."""

    def test_complete_workflow_is_restartable_and_writes_provenance(self) -> None:
        """All stages should complete and a resume should preserve outputs."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            reference_fasta = root / "alpha_alternative_reference.fna"
            reference_fasta.write_text(
                ">ALTERNATIVE_REFERENCE\n" + ("ACGT" * 100) + "\n",
                encoding="utf-8",
            )
            genome_config = root / "genomes.tsv"
            write_tsv(
                records=[
                    {
                        "genome_fasta": str(reference_fasta),
                        "species_name": "Alpha species",
                        "role": "target_species",
                        "assembly_accession": "GCF_999999999.1",
                    }
                ],
                output_path=genome_config,
                fieldnames=[
                    "genome_fasta",
                    "species_name",
                    "role",
                    "assembly_accession",
                ],
            )
            truth = root / "truth.tsv"
            write_tsv(
                records=[
                    {
                        "benchmark_id": "synthetic_mock_v1",
                        "organism_id": "alpha",
                        "species_name": "Alpha species",
                        "expected": "true",
                        "expected_abundance_fraction": "1.0",
                        "abundance_tier": "synthetic",
                    }
                ],
                output_path=truth,
                fieldnames=[
                    "benchmark_id",
                    "organism_id",
                    "species_name",
                    "expected",
                    "expected_abundance_fraction",
                    "abundance_tier",
                ],
            )
            panel = root / "panel.tsv"
            write_tsv(
                records=[
                    {
                        "kmer": "AAAAA",
                        "k": 5,
                        "panel_type": "species_unique",
                        "species_name": "Alpha species",
                        "clade": "Synthetic",
                        "source_genomes": "alternative",
                        "source_contigs": "contig1",
                        "example_position": 0,
                        "evidence_taxid": "1",
                        "evidence_name": "Alpha species",
                        "evidence_rank": "species",
                        "lineage_taxids": "1",
                        "source_taxids": "1",
                    },
                    {
                        "kmer": "AAAAAAA",
                        "k": 7,
                        "panel_type": "species_unique",
                        "species_name": "Alpha species",
                        "clade": "Synthetic",
                        "source_genomes": "alternative",
                        "source_contigs": "contig1",
                        "example_position": 20,
                        "evidence_taxid": "1",
                        "evidence_name": "Alpha species",
                        "evidence_rank": "species",
                        "lineage_taxids": "1",
                        "source_taxids": "1",
                    },
                ],
                output_path=panel,
                fieldnames=PANEL_COLUMNS,
            )
            reads = root / "reads.fastq"
            reads.write_text(
                "".join(
                    f"@read{index}\nCCCAAAAAAACCC\n+\nFFFFFFFFFFFFF\n"
                    for index in range(1, 9)
                ),
                encoding="utf-8",
            )
            config = {
                "benchmark_id": "synthetic_mock_v1",
                "sample_id": "synthetic",
                "allow_primary_threshold_tuning": False,
                "dataset": {
                    "input_reads": str(reads),
                    "input_format": "fastq",
                },
                "reference": {
                    "truth_manifest": str(truth),
                    "flat_panel": str(panel),
                    "panel_genome_config": str(genome_config),
                    "module_manifest": "",
                },
                "screen": {
                    "screen_preset": "exact",
                    "max_mismatches": 0,
                    "call_preset": "lineage_aware",
                    "same_genus_reportable_min_fraction": 0.05,
                    "k_values": [5, 7],
                    "chunk_size": 2,
                    "decompressor": "python",
                    "write_parquet_outputs": False,
                },
                "depth": {"fractions": [0.5], "seeds": [7]},
                "ai": {"enabled": False, "novelty_scale": 2.9},
            }
            options = WorkflowOptions(
                output_root=root / "runs",
                run_name="complete",
                threads=1,
                resume=True,
                start_at=None,
                stop_after=None,
                force_stages=frozenset(),
                dry_run=False,
            )
            workflow = BenchmarkWorkflow(
                config=config,
                config_digest="c" * 64,
                options=options,
                logger=logging.getLogger("workflow_integration"),
            )
            workflow.run()
            run_status = (
                workflow.stage_dir("09_provenance") / "run_status.tsv"
            ).read_text(encoding="utf-8")
            self.assertIn("09_provenance\tsuccess", run_status)
            for stage in STAGES:
                self.assertTrue(workflow.stage_is_complete(stage), stage)

            calls = (
                workflow.stage_dir("03_screen_full")
                / "synthetic"
                / "species_detection_calls.tsv"
            )
            before = calls.stat().st_mtime_ns
            workflow.run()
            self.assertEqual(calls.stat().st_mtime_ns, before)


if __name__ == "__main__":
    unittest.main()
