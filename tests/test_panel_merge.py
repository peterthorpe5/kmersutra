"""Tests for KmerSutra panel merging and validation."""

from __future__ import annotations

import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from kmersutra.io import read_tsv, write_tsv
from kmersutra.panel_merge import (
    PANEL_REQUIRED_COLUMNS,
    load_panel_records,
    merge_panel_files,
    merge_panel_records,
    validate_panel_file,
    validate_panel_records,
)
from kmersutra.taxonomy import TaxonomyDatabase


def write_tiny_taxdump(root: Path) -> None:
    """Write a small NCBI-like taxonomy dump for panel merge tests."""
    root.mkdir(parents=True, exist_ok=True)
    (root / "nodes.dmp").write_text(
        "1\t|\t1\t|\tno rank\t|\n"
        "2\t|\t1\t|\tsuperkingdom\t|\n"
        "10\t|\t2\t|\tgenus\t|\n"
        "11\t|\t10\t|\tspecies\t|\n"
        "12\t|\t10\t|\tspecies\t|\n"
        "20\t|\t2\t|\tgenus\t|\n"
        "21\t|\t20\t|\tspecies\t|\n"
        "200\t|\t1\t|\tsuperkingdom\t|\n"
        "201\t|\t200\t|\tspecies\t|\n",
        encoding="utf-8",
    )
    (root / "names.dmp").write_text(
        "1\t|\troot\t|\t\t|\tscientific name\t|\n"
        "2\t|\tPathogenia\t|\t\t|\tscientific name\t|\n"
        "10\t|\tAlphaGenus\t|\t\t|\tscientific name\t|\n"
        "11\t|\tAlpha species one\t|\t\t|\tscientific name\t|\n"
        "12\t|\tAlpha species two\t|\t\t|\tscientific name\t|\n"
        "20\t|\tBetaGenus\t|\t\t|\tscientific name\t|\n"
        "21\t|\tBeta species one\t|\t\t|\tscientific name\t|\n"
        "200\t|\tHostia\t|\t\t|\tscientific name\t|\n"
        "201\t|\tHost species one\t|\t\t|\tscientific name\t|\n",
        encoding="utf-8",
    )
    (root / "merged.dmp").write_text("", encoding="utf-8")
    (root / "delnodes.dmp").write_text("", encoding="utf-8")


def panel_record(
    *,
    kmer: str,
    species_name: str,
    clade: str,
    taxid: str,
    panel_type: str = "species_unique",
    evidence_rank: str = "species",
    evidence_name: str | None = None,
) -> dict[str, object]:
    """Create a complete tiny panel record for tests."""
    return {
        "kmer": kmer,
        "k": len(kmer),
        "panel_type": panel_type,
        "species_name": species_name,
        "clade": clade,
        "source_genomes": f"genome_{taxid}",
        "source_contigs": f"contig_{taxid}",
        "example_position": 0,
        "evidence_taxid": taxid,
        "evidence_name": evidence_name or species_name,
        "evidence_rank": evidence_rank,
        "lineage_taxids": taxid,
        "source_taxids": taxid,
    }


class TestPanelMerge(unittest.TestCase):
    """Tests for merging independent KmerSutra panel modules."""

    def test_load_panel_records_adds_source_panel(self) -> None:
        """Panel loading should record which source panel each row came from."""
        with TemporaryDirectory() as tmpdir:
            panel_path = Path(tmpdir) / "panel.tsv"
            write_tsv(
                records=[
                    panel_record(
                        kmer="AAAAA",
                        species_name="Alpha species one",
                        clade="AlphaGenus",
                        taxid="11",
                    )
                ],
                output_path=panel_path,
                fieldnames=PANEL_REQUIRED_COLUMNS,
            )
            records = load_panel_records(panel_paths=[panel_path])
        self.assertEqual(len(records), 1)
        self.assertIn("source_panels", records[0])

    def test_merge_without_taxonomy_keeps_single_species_kmer(self) -> None:
        """Label-only merge should retain k-mers unique to one species."""
        result = merge_panel_records(
            records=[
                panel_record(
                    kmer="AAAAA",
                    species_name="Alpha species one",
                    clade="AlphaGenus",
                    taxid="11",
                )
            ]
        )
        self.assertEqual(len(result.master_records), 1)
        self.assertEqual(result.master_records[0]["panel_type"], "species_unique")

    def test_merge_without_taxonomy_downgrades_same_clade_kmer(self) -> None:
        """Label-only merge should downgrade shared same-clade k-mers."""
        result = merge_panel_records(
            records=[
                panel_record(
                    kmer="AAAAA",
                    species_name="Alpha species one",
                    clade="AlphaGenus",
                    taxid="11",
                ),
                panel_record(
                    kmer="AAAAA",
                    species_name="Alpha species two",
                    clade="AlphaGenus",
                    taxid="12",
                ),
            ]
        )
        self.assertEqual(len(result.master_records), 1)
        self.assertEqual(result.master_records[0]["panel_type"], "clade_core")
        self.assertEqual(len(result.downgraded_records), 1)

    def test_taxonomy_merge_downgrades_species_to_genus(self) -> None:
        """Taxonomy-aware merge should assign shared species k-mers to genus."""
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            write_tiny_taxdump(root)
            taxonomy = TaxonomyDatabase.from_taxdump(taxonomy_dir=root)
        result = merge_panel_records(
            records=[
                panel_record(
                    kmer="CCCCC",
                    species_name="Alpha species one",
                    clade="AlphaGenus",
                    taxid="11",
                ),
                panel_record(
                    kmer="CCCCC",
                    species_name="Alpha species two",
                    clade="AlphaGenus",
                    taxid="12",
                ),
            ],
            taxonomy_db=taxonomy,
            preferred_ranks=["species", "genus", "superkingdom"],
        )
        self.assertEqual(len(result.master_records), 1)
        self.assertEqual(result.master_records[0]["evidence_rank"], "genus")
        self.assertEqual(result.master_records[0]["evidence_name"], "AlphaGenus")
        self.assertEqual(len(result.downgraded_records), 1)

    def test_taxonomy_merge_removes_root_only_kmer(self) -> None:
        """Taxonomy-aware merge should remove k-mers with no useful evidence rank."""
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            write_tiny_taxdump(root)
            taxonomy = TaxonomyDatabase.from_taxdump(taxonomy_dir=root)
        result = merge_panel_records(
            records=[
                panel_record(
                    kmer="GGGGG",
                    species_name="Alpha species one",
                    clade="AlphaGenus",
                    taxid="11",
                ),
                panel_record(
                    kmer="GGGGG",
                    species_name="Host species one",
                    clade="Hostia",
                    taxid="201",
                ),
            ],
            taxonomy_db=taxonomy,
            preferred_ranks=["species", "genus", "superkingdom"],
        )
        self.assertEqual(len(result.master_records), 0)
        self.assertEqual(len(result.removed_records), 2)

    def test_merge_panel_files_writes_expected_outputs(self) -> None:
        """File-level merge should write master, removed, downgraded and summary files."""
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            panel_a = root / "panel_a.tsv.gz"
            panel_b = root / "panel_b.tsv.gz"
            out_dir = root / "merged"
            write_tsv(
                records=[
                    panel_record(
                        kmer="AAAAA",
                        species_name="Alpha species one",
                        clade="AlphaGenus",
                        taxid="11",
                    )
                ],
                output_path=panel_a,
                fieldnames=PANEL_REQUIRED_COLUMNS,
            )
            write_tsv(
                records=[
                    panel_record(
                        kmer="TTTTT",
                        species_name="Beta species one",
                        clade="BetaGenus",
                        taxid="21",
                    )
                ],
                output_path=panel_b,
                fieldnames=PANEL_REQUIRED_COLUMNS,
            )
            merge_panel_files(panel_paths=[panel_a, panel_b], out_dir=out_dir)
            self.assertTrue((out_dir / "master_kmer_panel.tsv.gz").exists())
            self.assertTrue((out_dir / "master_validation_summary.tsv").exists())
            self.assertTrue((out_dir / "taxonomic_level_summary.tsv").exists())
            self.assertEqual(len(read_tsv(input_path=out_dir / "master_kmer_panel.tsv.gz")), 2)

    def test_validate_panel_records_detects_conflicting_duplicates(self) -> None:
        """Validation should report duplicate k-mers with conflicting evidence."""
        summary, issues, _levels = validate_panel_records(
            records=[
                panel_record(
                    kmer="AAAAA",
                    species_name="Alpha species one",
                    clade="AlphaGenus",
                    taxid="11",
                ),
                panel_record(
                    kmer="AAAAA",
                    species_name="Beta species one",
                    clade="BetaGenus",
                    taxid="21",
                ),
            ]
        )
        issue_types = {str(issue["issue_type"]) for issue in issues}
        self.assertIn("conflicting_duplicate_kmer_key", issue_types)
        self.assertIn({"metric": "n_conflicting_duplicate_kmer_keys", "value": 1}, summary)

    def test_validate_panel_file_writes_validation_tables(self) -> None:
        """File-level validation should write summary, issue and level tables."""
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            panel = root / "panel.tsv"
            out_dir = root / "validation"
            write_tsv(
                records=[
                    panel_record(
                        kmer="AAAAA",
                        species_name="Alpha species one",
                        clade="AlphaGenus",
                        taxid="11",
                    )
                ],
                output_path=panel,
                fieldnames=PANEL_REQUIRED_COLUMNS,
            )
            validate_panel_file(panel_path=panel, out_dir=out_dir)
            self.assertTrue((out_dir / "panel_validation_summary.tsv").exists())
            self.assertTrue((out_dir / "panel_validation_issues.tsv").exists())
            self.assertTrue((out_dir / "taxonomic_level_summary.tsv").exists())


if __name__ == "__main__":
    unittest.main()
