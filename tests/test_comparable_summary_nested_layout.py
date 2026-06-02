"""Tests for robust comparable summary sample-directory resolution."""

from __future__ import annotations

import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd


MODULE_PATH = Path(__file__).resolve().parents[1] / "summarise_kmersutra_comparable_benchmark.py"
SPEC = importlib.util.spec_from_file_location("kmersutra_comparable_summary", MODULE_PATH)
SUMMARY = importlib.util.module_from_spec(SPEC)
assert SPEC is not None and SPEC.loader is not None
sys.modules["kmersutra_comparable_summary"] = SUMMARY
SPEC.loader.exec_module(SUMMARY)


class TestComparableSummaryNestedLayout(unittest.TestCase):
    """Test nested-layout summary handling for historical manifest IDs."""

    def setUp(self) -> None:
        """Create a temporary comparable benchmark layout."""
        self.tmpdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tmpdir.name)
        self.out_root = self.root / "run"
        self.out_root.mkdir()
        self.manifest = self.out_root / "kmersutra_comparable_manifest.tsv"
        self.panel2 = self.root / "panel2.tsv"
        self.panel3 = self.root / "panel3.tsv"
        self.panel2.write_text(
            "fasta\ttarget_label\n"
            "pf.fa\tPlasmodium falciparum\n"
            "pv.fa\tPlasmodium vivax\n",
            encoding="utf-8",
        )
        self.panel3.write_text(
            "fasta\ttarget_label\n"
            "pf.fa\tPlasmodium falciparum\n"
            "pv.fa\tPlasmodium vivax\n"
            "pk.fa\tPlasmodium knowlesi\n",
            encoding="utf-8",
        )

    def tearDown(self) -> None:
        """Clean up temporary files."""
        self.tmpdir.cleanup()

    def write_manifest(self, rows: list[dict[str, str]]) -> None:
        """Write a manifest table.

        Parameters
        ----------
        rows : list[dict[str, str]]
            Manifest rows to write.
        """
        pd.DataFrame(rows).to_csv(self.manifest, sep="\t", index=False)

    def write_sample(
        self,
        *,
        family: str,
        output_name: str,
        call_label: str = "present_high_confidence",
    ) -> Path:
        """Write one minimal sample output directory.

        Parameters
        ----------
        family : str
            Benchmark family.
        output_name : str
            Output directory name.
        call_label : str
            Detection call label.

        Returns
        -------
        pathlib.Path
            Sample output directory.
        """
        sample_dir = self.out_root / "samples" / family / output_name
        sample_dir.mkdir(parents=True)
        calls = pd.DataFrame(
            [
                {
                    "species_name": "Plasmodium vivax",
                    "call": call_label,
                    "n_unique_kmers": "20",
                    "n_positive_sequences": "3",
                    "confidence_score": "0.95",
                    "conflict_ratio": "0",
                    "n_k_values_positive": "1",
                    "best_k": "151",
                    "evidence_rank": "species",
                }
            ]
        )
        calls.to_csv(sample_dir / "species_detection_calls.tsv", sep="\t", index=False)
        calls.to_csv(sample_dir / "sample_species_kmer_evidence.tsv", sep="\t", index=False)
        pd.DataFrame(
            [{"runtime_seconds": "1.25", "exit_status": "0"}]
        ).to_csv(sample_dir / "screen_task_timing.tsv", sep="\t", index=False)
        return sample_dir

    def test_sample_output_dir_matches_v014_manifest_to_v015_nested_output(self) -> None:
        """A v014 manifest ID should resolve to a v015 nested output directory."""
        row = {
            "sample_id": "kmersutra_v014_single_genome_spikein_single_read_kraken_confidence_0.10_20260427_165016_mix_rep1_n1",
            "input_fastq": "/fake/spikein_single_read_kraken_confidence_0.10_20260427_165016/mix_rep1_n1/mixed.fastq.gz",
            "benchmark_family": "single_genome",
            "panel": "panel1",
            "replicate": "1",
            "spike_n": "1",
            "source_run_dir": "/fake/spikein_single_read_kraken_confidence_0.10_20260427_165016",
            "source_relative_dir": "mix_rep1_n1",
        }
        expected_dir = self.write_sample(
            family="single_genome",
            output_name=(
                "kmersutra_v015_single_genome_spikein_single_read_kraken_"
                "confidence_0.10_20260427_165016_mix_rep1_n1"
            ),
        )
        observed = SUMMARY.sample_output_dir(out_root=self.out_root, row=row)
        self.assertEqual(observed, expected_dir)

    def test_sample_output_dir_uses_exact_suffix_not_prefix(self) -> None:
        """The n1 suffix must not accidentally match n10 or n100 outputs."""
        row = {
            "sample_id": "kmersutra_v014_single_genome_spikein_single_read_kraken_confidence_0.10_20260427_165016_mix_rep1_n1",
            "input_fastq": "/fake/spikein_single_read_kraken_confidence_0.10_20260427_165016/mix_rep1_n1/mixed.fastq.gz",
            "benchmark_family": "single_genome",
            "panel": "panel1",
            "replicate": "1",
            "spike_n": "1",
            "source_run_dir": "/fake/spikein_single_read_kraken_confidence_0.10_20260427_165016",
            "source_relative_dir": "mix_rep1_n1",
        }
        wrong_dir = self.write_sample(
            family="single_genome",
            output_name=(
                "kmersutra_v015_single_genome_spikein_single_read_kraken_"
                "confidence_0.10_20260427_165016_mix_rep1_n10"
            ),
        )
        correct_dir = self.write_sample(
            family="single_genome",
            output_name=(
                "kmersutra_v015_single_genome_spikein_single_read_kraken_"
                "confidence_0.10_20260427_165016_mix_rep1_n1"
            ),
        )
        observed = SUMMARY.sample_output_dir(out_root=self.out_root, row=row)
        self.assertEqual(observed, correct_dir)
        self.assertNotEqual(observed, wrong_dir)

    def test_run_summary_reads_nested_outputs_and_writes_long_tables(self) -> None:
        """The full summary should read nested calls and evidence outputs."""
        self.write_manifest(
            [
                {
                    "sample_id": "kmersutra_v014_single_genome_spikein_single_read_kraken_confidence_0.10_20260427_165016_mix_rep1_n1",
                    "input_fastq": "/fake/spikein_single_read_kraken_confidence_0.10_20260427_165016/mix_rep1_n1/mixed.fastq.gz",
                    "benchmark_family": "single_genome",
                    "panel": "panel1",
                    "replicate": "1",
                    "spike_n": "1",
                    "source_run_dir": "/fake/spikein_single_read_kraken_confidence_0.10_20260427_165016",
                    "source_relative_dir": "mix_rep1_n1",
                },
                {
                    "sample_id": "missing_sample",
                    "input_fastq": "/fake/spikein_single_read_kraken_confidence_0.10_20260427_165016/mix_rep2_n1/mixed.fastq.gz",
                    "benchmark_family": "single_genome",
                    "panel": "panel1",
                    "replicate": "2",
                    "spike_n": "1",
                    "source_run_dir": "/fake/spikein_single_read_kraken_confidence_0.10_20260427_165016",
                    "source_relative_dir": "mix_rep2_n1",
                },
            ]
        )
        self.write_sample(
            family="single_genome",
            output_name=(
                "kmersutra_v015_single_genome_spikein_single_read_kraken_"
                "confidence_0.10_20260427_165016_mix_rep1_n1"
            ),
        )
        out_dir = self.root / "summary"
        paths = SUMMARY.run_summary(
            out_root=self.out_root,
            manifest_path=self.manifest,
            out_dir=out_dir,
            panel1_targets=["Plasmodium vivax"],
            panel2_tsv=self.panel2,
            panel3_tsv=self.panel3,
            positive_calls={"present_high_confidence"},
            summary_name="unit_test_summary",
            allow_partial=True,
            strict=False,
            background_candidate_calls=set(),
            background_candidate_taxa=set(),
            demote_expected_genus_neighbours=True,
        )
        calls = pd.read_csv(paths.calls_long, sep="\t")
        status = pd.read_csv(paths.sample_status, sep="\t")
        self.assertEqual(len(calls), 1)
        self.assertEqual(int((status["screen_status"] == "ok").sum()), 1)
        self.assertEqual(int((status["screen_status"] != "ok").sum()), 1)

    def test_run_summary_raises_when_no_detection_calls_are_read(self) -> None:
        """A zero-call summary should fail rather than producing a false report."""
        self.write_manifest(
            [
                {
                    "sample_id": "missing_sample",
                    "input_fastq": "/fake/run/mix_rep1_n1/mixed.fastq.gz",
                    "benchmark_family": "single_genome",
                    "panel": "panel1",
                    "replicate": "1",
                    "spike_n": "1",
                    "source_run_dir": "/fake/run",
                    "source_relative_dir": "mix_rep1_n1",
                }
            ]
        )
        with self.assertRaises(RuntimeError):
            SUMMARY.run_summary(
                out_root=self.out_root,
                manifest_path=self.manifest,
                out_dir=self.root / "summary_empty",
                panel1_targets=["Plasmodium vivax"],
                panel2_tsv=self.panel2,
                panel3_tsv=self.panel3,
                positive_calls={"present_high_confidence"},
                summary_name="unit_test_summary_empty",
                allow_partial=True,
                strict=False,
                background_candidate_calls=set(),
                background_candidate_taxa=set(),
                demote_expected_genus_neighbours=True,
            )


if __name__ == "__main__":
    unittest.main()
