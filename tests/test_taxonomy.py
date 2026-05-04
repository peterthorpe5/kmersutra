"""Tests for NCBI taxonomy support."""

from __future__ import annotations

import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from kmersutra.taxonomy import TaxonomyDatabase, ensure_taxdump_files


def write_tiny_taxdump(root: Path) -> None:
    """Write a tiny NCBI-like taxdump for tests."""
    root.mkdir(parents=True, exist_ok=True)
    (root / "nodes.dmp").write_text(
        "1\t|\t1\t|\tno rank\t|\n"
        "2\t|\t1\t|\tsuperkingdom\t|\n"
        "10\t|\t2\t|\tgenus\t|\n"
        "11\t|\t10\t|\tspecies\t|\n"
        "12\t|\t10\t|\tspecies\t|\n"
        "20\t|\t2\t|\tgenus\t|\n"
        "21\t|\t20\t|\tspecies\t|\n",
        encoding="utf-8",
    )
    (root / "names.dmp").write_text(
        "1\t|\troot\t|\t\t|\tscientific name\t|\n"
        "2\t|\tPathogenia\t|\t\t|\tscientific name\t|\n"
        "10\t|\tAlphaGenus\t|\t\t|\tscientific name\t|\n"
        "11\t|\tAlpha species one\t|\t\t|\tscientific name\t|\n"
        "12\t|\tAlpha species two\t|\t\t|\tscientific name\t|\n"
        "20\t|\tBetaGenus\t|\t\t|\tscientific name\t|\n"
        "21\t|\tBeta species one\t|\t\t|\tscientific name\t|\n",
        encoding="utf-8",
    )
    (root / "merged.dmp").write_text("111\t|\t11\t|\n", encoding="utf-8")
    (root / "delnodes.dmp").write_text("999\t|\n", encoding="utf-8")


class TestTaxonomy(unittest.TestCase):
    """Tests for taxonomy parsing and lineage operations."""

    def test_ensure_taxdump_files_reports_missing_files(self) -> None:
        """Missing taxdump files should raise a clear error."""
        with TemporaryDirectory() as tmpdir:
            with self.assertRaises(FileNotFoundError):
                ensure_taxdump_files(taxonomy_dir=tmpdir)

    def test_load_taxonomy_and_lineage(self) -> None:
        """Taxonomy parser should load names, ranks and lineages."""
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            write_tiny_taxdump(root)
            taxonomy = TaxonomyDatabase.from_taxdump(taxonomy_dir=root)
        self.assertEqual(taxonomy.get_name("11"), "Alpha species one")
        self.assertEqual(taxonomy.get_rank("10"), "genus")
        self.assertEqual(taxonomy.get_lineage("11"), ["1", "2", "10", "11"])

    def test_merged_and_deleted_taxids_are_handled(self) -> None:
        """Taxonomy parser should map merged taxids and reject deleted taxids."""
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            write_tiny_taxdump(root)
            taxonomy = TaxonomyDatabase.from_taxdump(taxonomy_dir=root)
        self.assertEqual(taxonomy.normalise_taxid("111"), "11")
        self.assertEqual(taxonomy.normalise_taxid("999"), "")

    def test_lowest_common_ancestor(self) -> None:
        """Lowest common ancestor should identify shared genus-level support."""
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            write_tiny_taxdump(root)
            taxonomy = TaxonomyDatabase.from_taxdump(taxonomy_dir=root)
        self.assertEqual(taxonomy.lowest_common_ancestor({"11", "12"}), "10")
        self.assertEqual(taxonomy.get_rank("10"), "genus")

    def test_descendant_check(self) -> None:
        """Descendant checks should find taxids inside an ancestor lineage."""
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            write_tiny_taxdump(root)
            taxonomy = TaxonomyDatabase.from_taxdump(taxonomy_dir=root)
        self.assertTrue(taxonomy.is_descendant(taxid="11", ancestor_taxid="10"))
        self.assertFalse(taxonomy.is_descendant(taxid="21", ancestor_taxid="10"))

    def test_best_named_ancestor_prefers_core_ranks(self) -> None:
        """Best named ancestor should report a useful evidence rank."""
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            write_tiny_taxdump(root)
            taxonomy = TaxonomyDatabase.from_taxdump(taxonomy_dir=root)
        node = taxonomy.best_named_ancestor(taxids={"11", "12"})
        self.assertIsNotNone(node)
        assert node is not None
        self.assertEqual(node.rank, "genus")
        self.assertEqual(node.name, "AlphaGenus")


if __name__ == "__main__":
    unittest.main()
