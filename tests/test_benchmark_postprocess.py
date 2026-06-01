"""Tests for benchmark post-processing helpers."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from kmersutra.benchmark_postprocess import (
    BenchmarkPostprocessError,
    count_records_by_column,
    resolve_existing_table,
    run_benchmark_postprocess,
    summarise_lca_scopes,
    summarise_numeric_features,
    write_html_report,
)
from kmersutra.table_io import read_records_table, write_records_table


CALL_ROWS = [
    {
        "sample_id": "sample_1",
        "species_name": "Plasmodium vivax",
        "benchmark_family": "single_genome",
        "panel": "panel1",
        "replicate": "1",
        "spike_n": "500",
        "expected_targets": "Plasmodium vivax",
        "is_positive_expected": "True",
        "is_positive_call": "True",
        "n_hits": "100",
        "n_unique_kmers": "50",
        "n_positive_sequences": "20",
        "n_k_values_positive": "3",
        "best_k": "151",
        "n_exact_hits": "100",
        "conflict_ratio": "0.0",
        "confidence_score": "0.95",
    },
    {
        "sample_id": "sample_1",
        "species_name": "Plasmodium simium",
        "benchmark_family": "single_genome",
        "panel": "panel1",
        "replicate": "1",
        "spike_n": "500",
        "expected_targets": "Plasmodium vivax",
        "is_positive_neighbour_lineage": "True",
        "is_positive_call": "True",
        "n_hits": "70",
        "n_unique_kmers": "35",
        "n_positive_sequences": "10",
        "n_k_values_positive": "2",
        "best_k": "101",
        "n_exact_hits": "70",
        "conflict_ratio": "0.1",
        "confidence_score": "0.75",
    },
    {
        "sample_id": "sample_2",
        "species_name": "Hammondia hammondi",
        "benchmark_family": "shuffled_negative",
        "panel": "shuffled",
        "replicate": "1",
        "spike_n": "0",
        "expected_targets": "Plasmodium vivax",
        "is_background_candidate_signal": "True",
        "is_positive_call": "True",
        "n_hits": "80",
        "n_unique_kmers": "28",
        "n_positive_sequences": "73",
        "n_k_values_positive": "1",
        "best_k": "51",
        "n_exact_hits": "80",
        "conflict_ratio": "0.0",
        "confidence_score": "0.7",
    },
]


class TestBenchmarkPostprocess(unittest.TestCase):
    """Test v0.35 benchmark post-processing."""

    def _write_taxonomy(self, *, root: Path) -> Path:
        """Write a tiny taxonomy dump.

        Parameters
        ----------
        root : pathlib.Path
            Temporary root directory.

        Returns
        -------
        pathlib.Path
            Taxonomy directory.
        """
        taxonomy_dir = root / "taxonomy"
        taxonomy_dir.mkdir()
        (taxonomy_dir / "nodes.dmp").write_text(
            "1\t|\t1\t|\tno rank\t|\n"
            "2\t|\t1\t|\tsuperkingdom\t|\n"
            "10\t|\t2\t|\tphylum\t|\n"
            "20\t|\t10\t|\tclass\t|\n"
            "100\t|\t20\t|\tgenus\t|\n"
            "101\t|\t100\t|\tspecies\t|\n"
            "102\t|\t100\t|\tspecies\t|\n"
            "200\t|\t20\t|\tgenus\t|\n"
            "201\t|\t200\t|\tspecies\t|\n"
        )
        (taxonomy_dir / "names.dmp").write_text(
            "1\t|\troot\t|\t\t|\tscientific name\t|\n"
            "2\t|\tEukaryota\t|\t\t|\tscientific name\t|\n"
            "10\t|\tApicomplexa\t|\t\t|\tscientific name\t|\n"
            "20\t|\tAconoidasida\t|\t\t|\tscientific name\t|\n"
            "100\t|\tPlasmodium\t|\t\t|\tscientific name\t|\n"
            "101\t|\tPlasmodium vivax\t|\t\t|\tscientific name\t|\n"
            "102\t|\tPlasmodium simium\t|\t\t|\tscientific name\t|\n"
            "200\t|\tHammondia\t|\t\t|\tscientific name\t|\n"
            "201\t|\tHammondia hammondi\t|\t\t|\tscientific name\t|\n"
        )
        (taxonomy_dir / "merged.dmp").write_text("")
        (taxonomy_dir / "delnodes.dmp").write_text("")
        return taxonomy_dir

    def test_resolve_existing_table_uses_standard_summary_name(self) -> None:
        """Standard long-call table names should be auto-detected."""
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            calls_table = root / "kmersutra_detection_calls_long.tsv.gz"
            write_records_table(records=CALL_ROWS, output_path=calls_table)
            resolved = resolve_existing_table(
                explicit_path=None,
                search_dir=root,
                candidate_names=["kmersutra_detection_calls_long.tsv.gz"],
                description="calls table",
            )
        self.assertEqual(resolved.name, "kmersutra_detection_calls_long.tsv.gz")

    def test_resolve_existing_table_fails_for_missing_path(self) -> None:
        """Missing explicit paths should fail with a clear exception."""
        with self.assertRaises(BenchmarkPostprocessError):
            resolve_existing_table(
                explicit_path="/definitely/not/here.tsv",
                search_dir=None,
                candidate_names=[],
                description="missing table",
            )

    def test_summary_helpers_count_lca_and_ai_features(self) -> None:
        """Report helper functions should produce compact summaries."""
        lca_rows = [
            {
                "sample_id": "s1",
                "lca_scope": "dominant_lineage",
                "lca_rank": "genus",
                "lca_interpretation": "genus_resolved",
                "lca_name": "Plasmodium",
                "total_unique_kmers": "85",
                "total_positive_sequences": "30",
                "max_best_k": "151",
            },
            {
                "sample_id": "s2",
                "lca_scope": "background_candidate",
                "lca_rank": "species",
                "lca_interpretation": "species_resolved",
                "lca_name": "Hammondia hammondi",
                "total_unique_kmers": "28",
                "total_positive_sequences": "73",
                "max_best_k": "51",
            },
        ]
        ai_rows = [
            {"ml_report_label": "expected_target", "lca_dominant_lineage_rank_depth": "6"},
            {"ml_report_label": "expected_target", "lca_dominant_lineage_rank_depth": "7"},
            {"ml_report_label": "background", "lca_dominant_lineage_rank_depth": "0"},
        ]
        self.assertEqual(count_records_by_column(records=ai_rows, column="ml_report_label")[0]["n_records"], 2)
        self.assertEqual(len(summarise_lca_scopes(lca_records=lca_rows)), 2)
        feature_summary = summarise_numeric_features(records=ai_rows)
        self.assertEqual(feature_summary[0]["feature"], "lca_dominant_lineage_rank_depth")
        self.assertEqual(feature_summary[0]["n_non_zero"], 2)

    def test_write_html_report(self) -> None:
        """HTML reports should be written without requiring pandas."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output = Path(tmpdir) / "report.html"
            write_html_report(
                tables={"example": [{"a": "1", "b": "two"}]},
                output_path=output,
                metadata={"source": "unit-test"},
            )
            text = output.read_text()
        self.assertIn("KmerSutra benchmark post-processing report", text)
        self.assertIn("unit-test", text)

    def test_run_benchmark_postprocess_writes_lca_and_ai_tables(self) -> None:
        """The combined workflow should write LCA and AI training outputs."""
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            summary_dir = root / "summary"
            summary_dir.mkdir()
            calls_table = summary_dir / "kmersutra_detection_calls_long.tsv.gz"
            taxon_map = root / "kmersutra_genome_config.tsv"
            output_dir = root / "postprocess"
            taxonomy_dir = self._write_taxonomy(root=root)
            write_records_table(records=CALL_ROWS, output_path=calls_table)
            write_records_table(
                records=[
                    {"species_name": "Plasmodium vivax", "taxid": "101"},
                    {"species_name": "Plasmodium simium", "taxid": "102"},
                    {"species_name": "Hammondia hammondi", "taxid": "201"},
                ],
                output_path=taxon_map,
            )
            outputs = run_benchmark_postprocess(
                summary_dir=summary_dir,
                out_dir=output_dir,
                taxonomy_dir=taxonomy_dir,
                taxon_map_table=taxon_map,
                write_excel=False,
                write_html=True,
            )
            self.assertTrue(outputs["html_report"].exists())
            lca_records = read_records_table(input_path=outputs["lca_table"])
            ai_records = read_records_table(input_path=outputs["ai_training_table"])
            self.assertTrue(lca_records)
            self.assertTrue(ai_records)
            self.assertIn("lca_dominant_lineage_rank_depth", ai_records[0])


if __name__ == "__main__":
    unittest.main()
