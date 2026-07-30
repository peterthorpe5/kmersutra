"""Tests for NCBI taxonomy support."""

from __future__ import annotations

import io
import logging
import unittest
import zipfile
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import Mock, patch

from kmersutra.taxonomy import (
    REQUIRED_TAXDUMP_FILES,
    TaxonomyDatabase,
    TaxonomyNode,
    download_taxdump,
    ensure_taxdump_files,
)


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

    def test_existing_taxdump_is_reused_and_logged(self) -> None:
        """A complete local taxdump should not trigger a download."""
        logger = Mock(spec=logging.Logger)
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            write_tiny_taxdump(root)
            observed = download_taxdump(
                taxonomy_dir=root,
                logger=logger,
            )
        self.assertEqual(observed, root)
        logger.info.assert_called_once()

    def test_download_taxdump_extracts_required_files(self) -> None:
        """A valid downloaded archive should extract the four required files."""
        archive = io.BytesIO()
        with zipfile.ZipFile(archive, "w") as zip_handle:
            zip_handle.writestr("nodes.dmp", "1\t|\t1\t|\tno rank\t|\n")
            zip_handle.writestr(
                "names.dmp",
                "1\t|\troot\t|\t\t|\tscientific name\t|\n",
            )
            zip_handle.writestr("merged.dmp", "")
            zip_handle.writestr("delnodes.dmp", "")
        logger = Mock(spec=logging.Logger)
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            with patch(
                "kmersutra.taxonomy.urllib.request.urlopen",
                return_value=io.BytesIO(archive.getvalue()),
            ):
                observed = download_taxdump(
                    taxonomy_dir=root,
                    url="https://example.invalid/taxdmp.zip",
                    logger=logger,
                )
            self.assertEqual(observed, root)
            self.assertEqual(
                {path.name for path in root.iterdir()},
                set(REQUIRED_TAXDUMP_FILES),
            )
        self.assertEqual(logger.info.call_count, 2)

    def test_download_taxdump_preserves_existing_file_without_overwrite(self) -> None:
        """A partial download should not replace an existing required file."""
        archive = io.BytesIO()
        with zipfile.ZipFile(archive, "w") as zip_handle:
            for filename in REQUIRED_TAXDUMP_FILES:
                zip_handle.writestr(filename, f"downloaded {filename}\n")
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            root.mkdir(exist_ok=True)
            existing = root / "nodes.dmp"
            existing.write_text("local nodes\n", encoding="utf-8")
            with patch(
                "kmersutra.taxonomy.urllib.request.urlopen",
                return_value=io.BytesIO(archive.getvalue()),
            ):
                download_taxdump(taxonomy_dir=root)
            self.assertEqual(
                existing.read_text(encoding="utf-8"),
                "local nodes\n",
            )

    def test_download_taxdump_rejects_incomplete_archive(self) -> None:
        """A corrupt logical taxdump should identify missing members."""
        archive = io.BytesIO()
        with zipfile.ZipFile(archive, "w") as zip_handle:
            zip_handle.writestr("nodes.dmp", "nodes\n")
        with TemporaryDirectory() as tmpdir:
            with (
                patch(
                    "kmersutra.taxonomy.urllib.request.urlopen",
                    return_value=io.BytesIO(archive.getvalue()),
                ),
                self.assertRaisesRegex(ValueError, "missing required files"),
            ):
                download_taxdump(taxonomy_dir=tmpdir)

    def test_ensure_taxdump_can_download_missing_files(self) -> None:
        """The ensure helper should delegate when downloading is enabled."""
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)

            def populate(**kwargs: object) -> Path:
                write_tiny_taxdump(Path(str(kwargs["taxonomy_dir"])))
                return Path(str(kwargs["taxonomy_dir"]))

            logger = Mock(spec=logging.Logger)
            with patch(
                "kmersutra.taxonomy.download_taxdump",
                side_effect=populate,
            ) as mocked:
                paths = ensure_taxdump_files(
                    taxonomy_dir=root,
                    download_if_missing=True,
                    logger=logger,
                )
            self.assertEqual(set(paths), set(REQUIRED_TAXDUMP_FILES))
            mocked.assert_called_once()
            logger.info.assert_called_once()

    def test_load_taxonomy_and_lineage(self) -> None:
        """Taxonomy parser should load names, ranks and lineages."""
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            write_tiny_taxdump(root)
            taxonomy = TaxonomyDatabase.from_taxdump(taxonomy_dir=root)
        self.assertEqual(taxonomy.get_name("11"), "Alpha species one")
        self.assertEqual(taxonomy.get_rank("10"), "genus")
        self.assertEqual(taxonomy.get_lineage("11"), ["1", "2", "10", "11"])

    def test_load_taxonomy_reports_counts_to_logger(self) -> None:
        """Taxdump loading should report parsed object counts."""
        logger = Mock(spec=logging.Logger)
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            write_tiny_taxdump(root)
            TaxonomyDatabase.from_taxdump(
                taxonomy_dir=root,
                logger=logger,
            )
        logger.info.assert_called_once()

    def test_taxdump_readers_ignore_malformed_lines(self) -> None:
        """Short dump lines should be ignored without inventing records."""
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            names = root / "names.dmp"
            nodes = root / "nodes.dmp"
            merged = root / "merged.dmp"
            deleted = root / "delnodes.dmp"
            names.write_text(
                "short\t|\n"
                "1\t|\tcommon root\t|\t\t|\tsynonym\t|\n"
                "1\t|\troot\t|\t\t|\tscientific name\t|\n",
                encoding="utf-8",
            )
            nodes.write_text(
                "short\t|\n1\t|\t1\t|\tno rank\t|\n",
                encoding="utf-8",
            )
            merged.write_text(
                "short\t|\n2\t|\t1\t|\n",
                encoding="utf-8",
            )
            deleted.write_text("3\t|\n", encoding="utf-8")
            parsed_names = TaxonomyDatabase._read_scientific_names(names)
            parsed_nodes = TaxonomyDatabase._read_nodes(nodes, parsed_names)
            parsed_merged = TaxonomyDatabase._read_merged(merged)
            parsed_deleted = TaxonomyDatabase._read_deleted(deleted)
        self.assertEqual(parsed_names, {"1": "root"})
        self.assertEqual(parsed_nodes["1"].name, "root")
        self.assertEqual(parsed_merged, {"2": "1"})
        self.assertEqual(parsed_deleted, {"3"})

    def test_merged_and_deleted_taxids_are_handled(self) -> None:
        """Taxonomy parser should map merged taxids and reject deleted taxids."""
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            write_tiny_taxdump(root)
            taxonomy = TaxonomyDatabase.from_taxdump(taxonomy_dir=root)
        self.assertEqual(taxonomy.normalise_taxid("111"), "11")
        self.assertEqual(taxonomy.normalise_taxid("999"), "")
        self.assertEqual(taxonomy.normalise_taxid(None), "")
        self.assertEqual(taxonomy.normalise_taxid("  "), "")
        self.assertIsNone(taxonomy.get_node("missing"))
        self.assertEqual(taxonomy.get_name("missing"), "")
        self.assertEqual(taxonomy.get_rank("missing"), "")
        self.assertEqual(taxonomy.get_lineage("missing"), [])

    def test_lowest_common_ancestor(self) -> None:
        """Lowest common ancestor should identify shared genus-level support."""
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            write_tiny_taxdump(root)
            taxonomy = TaxonomyDatabase.from_taxdump(taxonomy_dir=root)
        self.assertEqual(taxonomy.lowest_common_ancestor({"11", "12"}), "10")
        self.assertEqual(taxonomy.get_rank("10"), "genus")
        self.assertEqual(taxonomy.lowest_common_ancestor({"missing"}), "")
        self.assertEqual(taxonomy.lowest_common_ancestor({"11", "21"}), "2")

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

    def test_best_named_ancestor_handles_missing_and_custom_ranks(self) -> None:
        """Named-ancestor selection should support missing and custom ranks."""
        taxonomy = TaxonomyDatabase(
            nodes={
                "1": TaxonomyNode("1", "1", "no rank", "root"),
                "2": TaxonomyNode("2", "1", "strain", "isolate"),
            }
        )
        self.assertIsNone(taxonomy.best_named_ancestor(taxids={"missing"}))
        fallback = taxonomy.best_named_ancestor(
            taxids={"2"},
            preferred_ranks=["species"],
        )
        custom = taxonomy.best_named_ancestor(
            taxids={"2"},
            preferred_ranks=["strain"],
        )
        self.assertEqual(fallback, taxonomy.nodes["2"])
        self.assertEqual(custom, taxonomy.nodes["2"])

    def test_lineage_cycle_is_safely_bounded(self) -> None:
        """Malformed cyclic taxonomies should not loop forever."""
        taxonomy = TaxonomyDatabase(
            nodes={
                "1": TaxonomyNode("1", "2", "no rank"),
                "2": TaxonomyNode("2", "1", "species"),
            }
        )
        self.assertEqual(taxonomy.get_lineage("1"), ["2", "1"])


if __name__ == "__main__":
    unittest.main()
