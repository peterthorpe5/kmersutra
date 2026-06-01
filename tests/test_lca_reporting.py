"""Tests for conservative LCA reporting."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from kmersutra.lca_reporting import (
    build_taxon_name_map,
    summarise_lca_by_sample,
    summarise_lca_table,
)
from kmersutra.table_io import read_records_table, write_records_table
from kmersutra.taxonomy import TaxonomyDatabase, TaxonomyNode


class TestLcaReporting(unittest.TestCase):
    """Test LCA reporting from KmerSutra evidence rows."""

    def _taxonomy(self) -> TaxonomyDatabase:
        """Return a small taxonomy database for LCA tests.

        Returns
        -------
        TaxonomyDatabase
            Minimal taxonomy database.
        """
        nodes = {
            "1": TaxonomyNode("1", "1", "no rank", "root"),
            "2": TaxonomyNode("2", "1", "superkingdom", "Eukaryota"),
            "10": TaxonomyNode("10", "2", "phylum", "Apicomplexa"),
            "20": TaxonomyNode("20", "10", "class", "Aconoidasida"),
            "100": TaxonomyNode("100", "20", "genus", "Plasmodium"),
            "101": TaxonomyNode("101", "100", "species", "Plasmodium vivax"),
            "102": TaxonomyNode("102", "100", "species", "Plasmodium simium"),
            "200": TaxonomyNode("200", "20", "genus", "Hammondia"),
            "201": TaxonomyNode("201", "200", "species", "Hammondia hammondi"),
            "300": TaxonomyNode("300", "20", "genus", "Toxoplasma"),
            "301": TaxonomyNode("301", "300", "species", "Toxoplasma gondii"),
        }
        return TaxonomyDatabase(nodes=nodes)

    def test_name_map_detects_conflicting_taxids(self) -> None:
        """Name mapping should fail when one name has conflicting taxids."""
        records = [
            {"species_name": "Plasmodium vivax", "taxid": "101"},
            {"species_name": "Plasmodium vivax", "taxid": "999"},
        ]
        with self.assertRaises(ValueError):
            build_taxon_name_map(taxon_map_records=records)

    def test_dominant_lca_avoids_weak_distant_background_collapse(self) -> None:
        """Dominant-lineage LCA should not collapse because of weak background."""
        taxonomy = self._taxonomy()
        name_to_taxid = build_taxon_name_map(
            taxon_map_records=[
                {"species_name": "Plasmodium vivax", "taxid": "101"},
                {"species_name": "Plasmodium simium", "taxid": "102"},
                {"species_name": "Hammondia hammondi", "taxid": "201"},
            ]
        )
        rows = [
            {
                "sample_id": "sample_a",
                "species_name": "Plasmodium vivax",
                "n_unique_kmers": "400",
                "n_positive_sequences": "80",
                "n_k_values_positive": "4",
                "best_k": "151",
                "is_positive_call": "true",
            },
            {
                "sample_id": "sample_a",
                "species_name": "Plasmodium simium",
                "n_unique_kmers": "140",
                "n_positive_sequences": "20",
                "n_k_values_positive": "3",
                "best_k": "101",
                "is_positive_neighbour_lineage": "true",
            },
            {
                "sample_id": "sample_a",
                "species_name": "Hammondia hammondi",
                "n_unique_kmers": "28",
                "n_positive_sequences": "73",
                "n_k_values_positive": "1",
                "best_k": "51",
                "is_background_candidate_signal": "true",
            },
        ]
        records = summarise_lca_by_sample(
            evidence_records=rows,
            taxonomy=taxonomy,
            name_to_taxid=name_to_taxid,
            scopes=["dominant_lineage", "all_supported_evidence"],
        )
        by_scope = {record["lca_scope"]: record for record in records}
        self.assertEqual(by_scope["dominant_lineage"]["lca_name"], "Plasmodium")
        self.assertEqual(by_scope["dominant_lineage"]["lca_rank"], "genus")
        self.assertEqual(by_scope["all_supported_evidence"]["lca_name"], "Aconoidasida")
        self.assertEqual(by_scope["all_supported_evidence"]["lca_rank"], "class")

    def test_background_candidate_scope_reports_background_lca(self) -> None:
        """Background-candidate scope should use only background rows."""
        taxonomy = self._taxonomy()
        rows = [
            {
                "sample_id": "sample_b",
                "species_taxid": "201",
                "species_name": "Hammondia hammondi",
                "n_unique_kmers": "28",
                "n_positive_sequences": "73",
                "n_k_values_positive": "1",
                "best_k": "51",
                "is_background_candidate_signal": "true",
            },
            {
                "sample_id": "sample_b",
                "species_taxid": "101",
                "species_name": "Plasmodium vivax",
                "n_unique_kmers": "0",
                "n_positive_sequences": "0",
                "n_k_values_positive": "0",
                "best_k": "0",
                "is_positive_call": "false",
            },
        ]
        records = summarise_lca_by_sample(
            evidence_records=rows,
            taxonomy=taxonomy,
            scopes=["background_candidate"],
        )
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0]["lca_name"], "Hammondia hammondi")
        self.assertEqual(records[0]["lca_interpretation"], "species_resolved")

    def test_sample_without_supported_rows_is_retained(self) -> None:
        """A sample with no supported rows should still be reported."""
        records = summarise_lca_by_sample(
            evidence_records=[
                {
                    "sample_id": "sample_c",
                    "species_taxid": "101",
                    "species_name": "Plasmodium vivax",
                    "n_unique_kmers": "0",
                    "n_positive_sequences": "0",
                    "best_k": "0",
                }
            ],
            taxonomy=self._taxonomy(),
            scopes=["dominant_lineage"],
        )
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0]["sample_id"], "sample_c")
        self.assertEqual(records[0]["lca_interpretation"], "no_supported_taxonomic_placement")

    def test_table_lca_round_trip_with_taxdump_and_mapping(self) -> None:
        """LCA table writer should work with taxdump and name mapping inputs."""
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
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
            )
            (taxonomy_dir / "names.dmp").write_text(
                "1\t|\troot\t|\t\t|\tscientific name\t|\n"
                "2\t|\tEukaryota\t|\t\t|\tscientific name\t|\n"
                "10\t|\tApicomplexa\t|\t\t|\tscientific name\t|\n"
                "20\t|\tAconoidasida\t|\t\t|\tscientific name\t|\n"
                "100\t|\tPlasmodium\t|\t\t|\tscientific name\t|\n"
                "101\t|\tPlasmodium vivax\t|\t\t|\tscientific name\t|\n"
                "102\t|\tPlasmodium simium\t|\t\t|\tscientific name\t|\n"
            )
            (taxonomy_dir / "merged.dmp").write_text("")
            (taxonomy_dir / "delnodes.dmp").write_text("")
            evidence_table = root / "evidence.tsv"
            taxon_map = root / "taxon_map.tsv"
            output_table = root / "lca.tsv"
            write_records_table(
                records=[
                    {
                        "sample_id": "sample_d",
                        "species_name": "Plasmodium vivax",
                        "n_unique_kmers": "50",
                        "n_positive_sequences": "10",
                        "n_k_values_positive": "2",
                        "best_k": "101",
                    },
                    {
                        "sample_id": "sample_d",
                        "species_name": "Plasmodium simium",
                        "n_unique_kmers": "40",
                        "n_positive_sequences": "8",
                        "n_k_values_positive": "2",
                        "best_k": "101",
                    },
                ],
                output_path=evidence_table,
            )
            write_records_table(
                records=[
                    {"species_name": "Plasmodium vivax", "taxid": "101"},
                    {"species_name": "Plasmodium simium", "taxid": "102"},
                ],
                output_path=taxon_map,
            )
            summarise_lca_table(
                evidence_table=evidence_table,
                output_table=output_table,
                taxonomy_dir=taxonomy_dir,
                taxon_map_table=taxon_map,
                scopes=["dominant_lineage"],
            )
            output = read_records_table(input_path=output_table)
            self.assertEqual(len(output), 1)
            self.assertEqual(output[0]["lca_name"], "Plasmodium")
            self.assertEqual(output[0]["lca_rank"], "genus")


if __name__ == "__main__":
    unittest.main()
