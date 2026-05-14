"""Unit tests for the KmerSutra comparable benchmark summariser."""

from __future__ import annotations

import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

SCRIPT_PATH = Path(__file__).resolve().parents[1] / "summarise_kmersutra_comparable_benchmark.py"
spec = importlib.util.spec_from_file_location("kmersutra_summary", SCRIPT_PATH)
kmersutra_summary = importlib.util.module_from_spec(spec)
sys.modules["kmersutra_summary"] = kmersutra_summary
assert spec.loader is not None
spec.loader.exec_module(kmersutra_summary)


class TestKmerSutraComparableSummary(unittest.TestCase):
    """Test the KmerSutra comparable summary workflow."""

    def setUp(self) -> None:
        """Create a temporary synthetic comparable benchmark run."""
        self.tmpdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tmpdir.name)
        self.out_root = self.root / "runs_kmersutra_v014_global_comparable_test"
        self.samples_dir = self.out_root / "samples"
        self.out_root.mkdir(parents=True)
        self.panel2_tsv = self.root / "pathogen_panel_2.tsv"
        self.panel3_tsv = self.root / "pathogen_panel_3.tsv"
        self.panel2_tsv.write_text(
            "fasta\ttarget_label\n"
            "pf.fa\tPlasmodium falciparum\n"
            "pv.fa\tPlasmodium vivax\n",
            encoding="utf-8",
        )
        self.panel3_tsv.write_text(
            "fasta\ttarget_label\n"
            "pf.fa\tPlasmodium falciparum\n"
            "pv.fa\tPlasmodium vivax\n"
            "pk.fa\tPlasmodium knowlesi\n",
            encoding="utf-8",
        )
        self.manifest = self.out_root / "kmersutra_v014_comparable_manifest.tsv"
        self._write_manifest()
        self._write_sample(
            family="single_genome",
            sample_id="single_rep1_n0",
            calls=[
                {
                    "species_name": "Plasmodium vivax",
                    "call": "absent",
                    "n_unique_kmers": "0",
                    "n_positive_sequences": "0",
                    "confidence_score": "0",
                    "conflict_ratio": "0",
                    "n_k_values_positive": "0",
                    "best_k": "",
                    "evidence_rank": "species",
                }
            ],
            runtime_seconds=4,
        )
        self._write_sample(
            family="single_genome",
            sample_id="single_rep1_n10",
            calls=[
                {
                    "species_name": "Plasmodium vivax",
                    "call": "present_high_confidence",
                    "n_unique_kmers": "12",
                    "n_positive_sequences": "2",
                    "confidence_score": "0.9",
                    "conflict_ratio": "0",
                    "n_k_values_positive": "2",
                    "best_k": "101",
                    "evidence_rank": "species",
                },
                {
                    "species_name": "Plasmodium falciparum",
                    "call": "present_low_confidence",
                    "n_unique_kmers": "4",
                    "n_positive_sequences": "1",
                    "confidence_score": "0.4",
                    "conflict_ratio": "0.2",
                    "n_k_values_positive": "1",
                    "best_k": "77",
                    "evidence_rank": "species",
                },
            ],
            runtime_seconds=5,
        )
        self._write_sample(
            family="two_genome",
            sample_id="panel2_rep1_n10",
            calls=[
                {
                    "species_name": "Plasmodium vivax",
                    "call": "present_high_confidence",
                    "n_unique_kmers": "20",
                    "n_positive_sequences": "3",
                    "confidence_score": "0.95",
                    "conflict_ratio": "0",
                    "n_k_values_positive": "2",
                    "best_k": "101",
                    "evidence_rank": "species",
                },
                {
                    "species_name": "Plasmodium falciparum",
                    "call": "present_high_confidence",
                    "n_unique_kmers": "22",
                    "n_positive_sequences": "3",
                    "confidence_score": "0.96",
                    "conflict_ratio": "0",
                    "n_k_values_positive": "2",
                    "best_k": "101",
                    "evidence_rank": "species",
                },
            ],
            runtime_seconds=6,
        )
        self._write_sample(
            family="shuffled_negative",
            sample_id="shuffle_rep1_n10",
            calls=[
                {
                    "species_name": "Plasmodium vivax",
                    "call": "present_high_confidence",
                    "n_unique_kmers": "10",
                    "n_positive_sequences": "1",
                    "confidence_score": "0.8",
                    "conflict_ratio": "0",
                    "n_k_values_positive": "2",
                    "best_k": "101",
                    "evidence_rank": "species",
                }
            ],
            runtime_seconds=7,
        )
        # Deliberately leave the three-genome sample missing to test partial summaries.

    def tearDown(self) -> None:
        """Clean up temporary files."""
        self.tmpdir.cleanup()

    def _write_manifest(self) -> None:
        """Write a synthetic comparable benchmark manifest."""
        rows = [
            [
                "single_rep1_n0",
                "/fake/single/mix_rep1_n0/mixed.fastq.gz",
                "single_genome",
                "panel1",
                "1",
                "0",
                "/fake/single",
                "mix_rep1_n0",
            ],
            [
                "single_rep1_n10",
                "/fake/single/mix_rep1_n10/mixed.fastq.gz",
                "single_genome",
                "panel1",
                "1",
                "10",
                "/fake/single",
                "mix_rep1_n10",
            ],
            [
                "panel2_rep1_n10",
                "/fake/panel2/mix_rep1_n10/mixed.fastq.gz",
                "two_genome",
                "panel2",
                "1",
                "10",
                "/fake/panel2",
                "mix_rep1_n10",
            ],
            [
                "shuffle_rep1_n10",
                "/fake/shuffle/mix_rep1_n10/mixed.fastq.gz",
                "shuffled_negative",
                "shuffled",
                "1",
                "10",
                "/fake/shuffle",
                "mix_rep1_n10",
            ],
            [
                "panel3_missing_rep1_n10",
                "/fake/panel3/mix_rep1_n10/mixed.fastq.gz",
                "three_genome",
                "panel3",
                "1",
                "10",
                "/fake/panel3",
                "mix_rep1_n10",
            ],
        ]
        dataframe = pd.DataFrame(
            rows,
            columns=[
                "sample_id",
                "input_fastq",
                "benchmark_family",
                "panel",
                "replicate",
                "spike_n",
                "source_run_dir",
                "source_relative_dir",
            ],
        )
        dataframe.to_csv(self.manifest, sep="\t", index=False)

    def _write_sample(
        self,
        *,
        family: str,
        sample_id: str,
        calls: list[dict[str, str]],
        runtime_seconds: int,
    ) -> None:
        """Write synthetic KmerSutra sample outputs."""
        sample_dir = self.samples_dir / family / sample_id
        sample_dir.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(calls).to_csv(
            sample_dir / "species_detection_calls.tsv",
            sep="\t",
            index=False,
        )
        pd.DataFrame(
            [
                {
                    "species_name": call["species_name"],
                    "n_hits": call.get("n_unique_kmers", "0"),
                    "n_unique_kmers": call.get("n_unique_kmers", "0"),
                    "n_positive_sequences": call.get("n_positive_sequences", "0"),
                }
                for call in calls
            ]
        ).to_csv(sample_dir / "sample_species_kmer_evidence.tsv", sep="\t", index=False)
        pd.DataFrame(
            [
                {
                    "sample_id": sample_id,
                    "runtime_seconds": runtime_seconds,
                    "exit_status": 0,
                }
            ]
        ).to_csv(sample_dir / "screen_timing.tsv", sep="\t", index=False)

    def test_load_manifest_requires_expected_columns(self) -> None:
        """Manifest loading should keep all required columns."""
        manifest = kmersutra_summary.load_manifest(path=self.manifest)
        self.assertEqual(len(manifest), 5)
        self.assertIn("sample_id", manifest.columns)

    def test_read_panel_targets_uses_named_label_column(self) -> None:
        """Panel target loading should use the target_label column."""
        labels = kmersutra_summary.read_panel_targets(path=self.panel2_tsv)
        self.assertEqual(
            labels,
            ["Plasmodium falciparum", "Plasmodium vivax"],
        )

    def test_shuffled_rows_are_negative_controls(self) -> None:
        """Shuffled rows should be negative even when spike_n is greater than zero."""
        manifest = kmersutra_summary.load_manifest(path=self.manifest)
        expected_map = kmersutra_summary.build_expected_target_map(
            panel1_targets=["Plasmodium vivax"],
            panel2_tsv=self.panel2_tsv,
            panel3_tsv=self.panel3_tsv,
        )
        row = manifest.loc[manifest["benchmark_family"] == "shuffled_negative"].iloc[0]
        metadata = kmersutra_summary.metadata_from_manifest_row(
            row=row.to_dict(),
            out_root=self.out_root,
            expected_targets=kmersutra_summary.expected_targets_for_row(
                row=row.to_dict(),
                expected_targets_by_panel=expected_map,
            ),
        )
        self.assertTrue(metadata["is_negative"])
        self.assertTrue(metadata["is_shuffled_control"])

    def test_sample_summary_counts_expected_and_off_target_calls(self) -> None:
        """Sample summary should count expected and off-target positive species."""
        paths = self._run_summary()
        sample_summary = pd.read_csv(paths.sample_summary, sep="\t")
        row = sample_summary.loc[sample_summary["sample_id"] == "single_rep1_n10"].iloc[0]
        self.assertEqual(int(row["n_expected_species_detected"]), 1)
        self.assertEqual(int(row["n_off_target_species"]), 1)
        self.assertEqual(int(row["all_expected_detected"]), 1)
        self.assertEqual(int(row["clean_expected_positive"]), 0)

    def test_missing_samples_are_retained_in_status_table(self) -> None:
        """Partial summaries should retain missing samples as missing_calls."""
        paths = self._run_summary()
        status = pd.read_csv(paths.sample_status, sep="\t")
        row = status.loc[status["sample_id"] == "panel3_missing_rep1_n10"].iloc[0]
        self.assertEqual(row["screen_status"], "missing_calls")

    def test_target_performance_has_lod_and_counts(self) -> None:
        """Tracked-target performance should include counts and LOD columns."""
        paths = self._run_summary()
        performance = pd.read_csv(paths.performance_by_target, sep="\t")
        single = performance.loc[
            (performance["benchmark_family"] == "single_genome")
            & (performance["target_label"] == "Plasmodium vivax")
        ].iloc[0]
        self.assertEqual(int(single["tp"]), 1)
        self.assertEqual(int(single["tn"]), 1)
        self.assertEqual(float(single["lod50_spike_n"]), 10.0)

    def test_real_world_summary_reports_off_target_rate(self) -> None:
        """Real-world summary should include positive off-target rates."""
        paths = self._run_summary()
        real_world = pd.read_csv(paths.real_world_summary, sep="\t")
        single = real_world.loc[real_world["benchmark_family"] == "single_genome"].iloc[0]
        self.assertGreater(float(single["positive_off_target_rate"]), 0.0)

    def test_off_target_summary_includes_non_expected_species(self) -> None:
        """Off-target summary should include the synthetic falciparum off-target."""
        paths = self._run_summary()
        off_target = pd.read_csv(paths.off_target_summary, sep="\t")
        self.assertIn("Plasmodium falciparum", set(off_target["report_label"]))

    def test_excel_and_html_outputs_are_written(self) -> None:
        """The summary workflow should write Excel and HTML reports."""
        paths = self._run_summary()
        self.assertTrue(paths.workbook.exists())
        self.assertGreater(paths.workbook.stat().st_size, 0)
        self.assertTrue(paths.html_report.exists())
        self.assertIn(
            "KmerSutra comparable spike-in summary",
            paths.html_report.read_text(encoding="utf-8"),
        )

    def test_strict_mode_fails_with_missing_calls(self) -> None:
        """Strict mode should fail when a manifest sample has no call table."""
        with self.assertRaises(RuntimeError):
            self._run_summary(strict=True, allow_partial=False)

    def _run_summary(
        self,
        *,
        strict: bool = False,
        allow_partial: bool = True,
    ) -> kmersutra_summary.SummaryPaths:
        """Run the summary workflow for the synthetic fixture."""
        summary_dir = self.out_root / "summary"
        paths = kmersutra_summary.build_summary_paths(
            out_dir=summary_dir,
            summary_name="kmersutra_comparable_summary",
        )
        kmersutra_summary.configure_logging(log_path=paths.log_path, verbose=False)
        return kmersutra_summary.run_summary(
            out_root=self.out_root,
            manifest_path=self.manifest,
            out_dir=summary_dir,
            panel1_targets=["Plasmodium vivax"],
            panel2_tsv=self.panel2_tsv,
            panel3_tsv=self.panel3_tsv,
            positive_calls=set(kmersutra_summary.DEFAULT_POSITIVE_CALLS),
            summary_name="kmersutra_comparable_summary",
            allow_partial=allow_partial,
            strict=strict,
        )


if __name__ == "__main__":
    unittest.main()
