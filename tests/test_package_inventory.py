"""Tests for deterministic package inventories."""

from __future__ import annotations

import csv
import tempfile
import unittest
from pathlib import Path

from kmersutra.package_inventory import (
    build_parser,
    create_package_inventory,
    sha256_file,
)


class TestPackageInventory(unittest.TestCase):
    """Test inventory content, exclusions and validation."""

    def test_inventory_is_sorted_and_excludes_generated_results(self) -> None:
        """Generated test outputs and caches should not enter inventories."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "z.txt").write_text("z\n", encoding="utf-8")
            (root / "a.txt").write_text("a\n", encoding="utf-8")
            (root / "test_results").mkdir()
            (root / "test_results" / "run.log").write_text(
                "generated\n",
                encoding="utf-8",
            )
            (root / "__pycache__").mkdir()
            (root / "__pycache__" / "module.pyc").write_bytes(b"cache")
            output = root / "outputs" / "inventory.tsv"

            count = create_package_inventory(
                repo_root=root,
                output_path=output,
            )

            with output.open("r", encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(count, 2)
        self.assertEqual(
            [row["relative_path"] for row in rows],
            ["a.txt", "z.txt"],
        )
        self.assertTrue(all(len(row["sha256"]) == 64 for row in rows))

    def test_inventory_rejects_missing_repository(self) -> None:
        """A missing repository root should fail clearly."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            with self.assertRaisesRegex(NotADirectoryError, "not a directory"):
                create_package_inventory(
                    repo_root=root / "missing",
                    output_path=root / "inventory.tsv",
                )

    def test_sha256_file_is_deterministic(self) -> None:
        """Repeated file hashing should return the same known digest."""
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "value.txt"
            path.write_bytes(b"abc")
            observed = sha256_file(path=path, chunk_size=1)
        self.assertEqual(
            observed,
            "ba7816bf8f01cfea414140de5dae2223"
            "b00361a396177a9cb410ff61f20015ad",
        )

    def test_parser_requires_named_paths(self) -> None:
        """The command should expose named input and output paths."""
        args = build_parser().parse_args(
            [
                "--repo-root",
                "/repository",
                "--output",
                "/results/inventory.tsv",
            ]
        )
        self.assertEqual(args.repo_root, "/repository")
        self.assertEqual(args.output, "/results/inventory.tsv")


if __name__ == "__main__":
    unittest.main()
