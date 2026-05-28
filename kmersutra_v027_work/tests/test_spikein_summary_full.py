"""Tests for full KmerSutra spike-in summary workflow."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import pandas as pd

from kmersutra import spikein_summary


class TestSpikeinSummaryWorkflow(unittest.TestCase):
    """Test spike-in summary helpers and output generation."""

    def test_classify_call_state_groups_expected_calls(self):
        """Call classifier should identify positive, ambiguous, and negative states."""
        self.assertEqual(
            spikein_summary.classify_call_state("present_in_mixed_sample"),
            "positive",
        )
        self.assertEqual(
            spikein_summary.classify_call_state("ambiguous_conflicting_signal"),
            "ambiguous",
        )
        self.assertEqual(spikein_summary.classify_call_state("not_detected"), "not_detected")

    def test_parse_expected_spike_levels(self):
        """Expected spike levels should be parsed from space-separated text."""
        self.assertEqual(
            spikein_summary.parse_expected_spike_levels(raw_value="0 1 5"),
            [0.0, 1.0, 5.0],
        )

    def test_build_species_long_from_wide_extracts_shell_columns(self):
        """Wide shell summaries should become one row per species and sample."""
        dataframe = pd.DataFrame(
            [
                {
                    "run_name": "run_a",
                    "run_dir": "/tmp/run_a",
                    "replicate": 1,
                    "spike_n": 5,
                    "total_spiked_reads": 15,
                    "kmersutra_Plasmodium_vivax_call": "present_low_confidence",
                    "kmersutra_Plasmodium_vivax_unique_kmers": 4,
                    "kmersutra_Plasmodium_vivax_positive_reads": 1,
                    "kmersutra_Plasmodium_vivax_confidence": 0.5,
                }
            ]
        )
        observed = spikein_summary.build_species_long_from_wide(summary_df=dataframe)
        self.assertEqual(observed.shape[0], 1)
        self.assertEqual(observed.iloc[0]["species_name"], "Plasmodium vivax")
        self.assertEqual(observed.iloc[0]["n_unique_kmers"], 4)

    def test_run_summary_workflow_writes_tsv_excel_and_html(self):
        """Summary workflow should write combined TSV, Excel, and HTML outputs."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            root = Path(tmp_dir)
            run_dir = root / "spikein_multi_kmersutra_test"
            ks_dir = run_dir / "mix_rep1_n0" / "kmersutra"
            ks_dir.mkdir(parents=True)
            calls_path = ks_dir / "species_detection_calls.tsv"
            calls_path.write_text(
                "sample_id\tspecies_name\tclade\tn_hits\tn_unique_kmers\t"
                "n_positive_sequences\tn_k_values_positive\tbest_k\tn_exact_hits\t"
                "n_fuzzy_hits\tconflicting_unique_kmers\tconflict_ratio\t"
                "confidence_score\tcall\n"
                "rep1_n0\tPlasmodium vivax\tPlasmodium\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\tnot_detected\n",
                encoding="utf-8",
            )
            summary_path = run_dir / "spikein_multi_kmersutra_summary.tsv"
            summary_path.write_text(
                "replicate\tspike_n\ttotal_spiked_reads\tkmersutra_out_dir\t"
                "kmersutra_calls_tsv\tkmersutra_Plasmodium_vivax_call\t"
                "kmersutra_Plasmodium_vivax_unique_kmers\t"
                "kmersutra_Plasmodium_vivax_positive_reads\t"
                "kmersutra_Plasmodium_vivax_confidence\n"
                f"1\t0\t0\t{ks_dir}\t{calls_path}\tnot_detected\t0\t0\t0\n",
                encoding="utf-8",
            )
            out_dir = root / "summary"
            spikein_summary.run_summary_workflow(
                input_dirs=[str(root)],
                out_dir=out_dir,
                summary_name="spikein_multi_kmersutra_summary.tsv",
                run_glob="spikein_multi_kmersutra*",
                expected_replicates=1,
                expected_spike_levels=[0.0],
                max_html_rows=20,
            )
            self.assertTrue((out_dir / "combined_run_summary.tsv").exists())
            self.assertTrue((out_dir / "kmersutra_spikein_overall_summary.xlsx").exists())
            self.assertTrue((out_dir / "kmersutra_spikein_overall_summary.html").exists())


if __name__ == "__main__":
    unittest.main()
