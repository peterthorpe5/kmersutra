"""Tests for deterministic benchmark depth subsets."""

from __future__ import annotations

import gzip
import tempfile
import unittest
from pathlib import Path

from kmersutra.depth_subsets import (
    create_depth_subsets,
    deterministic_unit_interval,
    fraction_label,
    validate_fractions,
)


class TestDepthSubsets(unittest.TestCase):
    """Test deterministic and nested record selection."""

    def write_fastq(self, root: Path, *, n_records: int = 100) -> Path:
        """Write a small valid FASTQ."""
        path = root / "reads.fastq"
        with path.open("w", encoding="utf-8") as handle:
            for index in range(n_records):
                handle.write(
                    f"@read_{index} description\n"
                    "ACGTACGT\n"
                    "+\n"
                    "IIIIIIII\n"
                )
        return path

    def count_fastq(self, path: Path) -> int:
        """Count records in a gzip FASTQ."""
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            lines = handle.readlines()
        self.assertEqual(len(lines) % 4, 0)
        return len(lines) // 4

    def read_identifiers(self, path: Path) -> set[str]:
        """Read FASTQ identifiers from a gzip subset."""
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            lines = handle.readlines()
        return {lines[index].split()[0] for index in range(0, len(lines), 4)}

    def test_fraction_validation(self) -> None:
        """Fractions should be sorted, unique and bounded."""
        self.assertEqual(validate_fractions([0.25, 0.01, 0.25]), [0.01, 0.25])
        with self.assertRaises(ValueError):
            validate_fractions([0.0])

    def test_deterministic_draw_is_repeatable_and_bounded(self) -> None:
        """The same record identity and seed should yield the same draw."""
        first = deterministic_unit_interval(
            identifier="read_1",
            record_index=1,
            seed=1001,
        )
        second = deterministic_unit_interval(
            identifier="read_1",
            record_index=1,
            seed=1001,
        )
        self.assertEqual(first, second)
        self.assertGreaterEqual(first, 0.0)
        self.assertLess(first, 1.0)

    def test_depth_subsets_are_nested_for_one_seed(self) -> None:
        """A lower fraction should be a subset of a higher fraction."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            manifest = create_depth_subsets(
                input_path=self.write_fastq(root, n_records=500),
                input_format="fastq",
                output_dir=root / "out",
                fractions=[0.1, 0.25],
                seeds=[1001],
                sample_prefix="sample",
            )
            paths = {
                row["fraction"]: Path(str(row["reads_path"])) for row in manifest
            }
            lower = self.read_identifiers(paths["0.10000000"])
            higher = self.read_identifiers(paths["0.25000000"])
        self.assertTrue(lower.issubset(higher))

    def test_repeat_runs_are_byte_identical(self) -> None:
        """Repeated runs should produce identical compressed FASTQ bytes."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = self.write_fastq(root)
            first = create_depth_subsets(
                input_path=source,
                input_format="fastq",
                output_dir=root / "one",
                fractions=[0.25],
                seeds=[7],
                sample_prefix="sample",
            )
            second = create_depth_subsets(
                input_path=source,
                input_format="fastq",
                output_dir=root / "two",
                fractions=[0.25],
                seeds=[7],
                sample_prefix="sample",
            )
            first_bytes = gzip.open(first[0]["reads_path"], "rb").read()
            second_bytes = gzip.open(second[0]["reads_path"], "rb").read()
        self.assertEqual(first_bytes, second_bytes)

    def test_fastq_quality_is_preserved(self) -> None:
        """Subsetting must retain original qualities and descriptions."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            manifest = create_depth_subsets(
                input_path=self.write_fastq(root, n_records=2),
                input_format="fastq",
                output_dir=root / "out",
                fractions=[1.0],
                seeds=[1],
                sample_prefix="sample",
            )
            with gzip.open(
                manifest[0]["reads_path"],
                "rt",
                encoding="utf-8",
            ) as handle:
                text = handle.read()
            self.assertEqual(
                self.count_fastq(Path(str(manifest[0]["reads_path"]))),
                2,
            )
        self.assertIn("@read_0 description", text)
        self.assertIn("IIIIIIII", text)

    def test_fraction_label_is_stable(self) -> None:
        """Fraction labels should not depend on locale."""
        self.assertEqual(fraction_label(0.01), "p010000")
        self.assertEqual(fraction_label(0.25), "p250000")


if __name__ == "__main__":
    unittest.main()
