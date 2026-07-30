"""Tests for deterministic benchmark depth subsets."""

from __future__ import annotations

import gzip
import io
import logging
import tempfile
import unittest
from pathlib import Path
from unittest.mock import Mock

from kmersutra.depth_subsets import (
    RawSequenceRecord,
    _write_record,
    create_depth_subsets,
    deterministic_unit_interval,
    fraction_label,
    iter_fasta_raw,
    iter_fastq_raw,
    validate_fractions,
    validate_seeds,
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
        with self.assertRaises(ValueError):
            validate_fractions([1.01])
        with self.assertRaisesRegex(ValueError, "At least one fraction"):
            validate_fractions([])

    def test_seed_validation_retains_first_seen_order(self) -> None:
        """Duplicate seeds should be removed without sorting."""
        self.assertEqual(validate_seeds([3, 1, 3]), [3, 1])
        with self.assertRaisesRegex(ValueError, "At least one seed"):
            validate_seeds([])

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
        with self.assertRaisesRegex(ValueError, "zero or greater"):
            deterministic_unit_interval(
                identifier="read_1",
                record_index=-1,
                seed=1001,
            )

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

    def test_fastq_iterator_rejects_truncated_record(self) -> None:
        """A partial four-line FASTQ record should be rejected."""
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "reads.fastq"
            path.write_text("@read\nACGT\n+\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "Truncated FASTQ"):
                list(iter_fastq_raw(input_path=path, decompressor="python"))

    def test_fastq_iterator_rejects_malformed_markers(self) -> None:
        """FASTQ headers and separator lines should be validated."""
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "reads.fastq"
            path.write_text("read\nACGT\n-\nIIII\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "Malformed FASTQ"):
                list(iter_fastq_raw(input_path=path, decompressor="python"))

    def test_fastq_iterator_rejects_length_mismatch(self) -> None:
        """Sequence and quality lengths must match."""
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "reads.fastq"
            path.write_text("@read\nACGT\n+\nIII\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "lengths differ"):
                list(iter_fastq_raw(input_path=path, decompressor="python"))

    def test_fastq_iterator_rejects_blank_identifier(self) -> None:
        """A FASTQ record must have a non-empty identifier."""
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "reads.fastq"
            path.write_text("@   \nACGT\n+\nIIII\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "Blank FASTQ identifier"):
                list(iter_fastq_raw(input_path=path, decompressor="python"))

    def test_fasta_iterator_preserves_description_and_sequence(self) -> None:
        """Multiline FASTA records should be joined and described."""
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "reads.fasta"
            path.write_text(
                "\n>read_1 description\nACGT\nTGCA\n>read_2\nNN\n",
                encoding="utf-8",
            )
            records = list(
                iter_fasta_raw(input_path=path, decompressor="python")
            )
        self.assertEqual(
            records,
            [
                RawSequenceRecord(
                    identifier="read_1",
                    description="read_1 description",
                    sequence="ACGTTGCA",
                    quality=None,
                ),
                RawSequenceRecord(
                    identifier="read_2",
                    description="read_2",
                    sequence="NN",
                    quality=None,
                ),
            ],
        )

    def test_fasta_iterator_rejects_blank_header(self) -> None:
        """Blank FASTA descriptions should fail clearly."""
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "reads.fasta"
            path.write_text(">   \nACGT\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "Blank FASTA header"):
                list(iter_fasta_raw(input_path=path, decompressor="python"))

    def test_fasta_iterator_rejects_sequence_before_header(self) -> None:
        """Sequence text before the first FASTA header is invalid."""
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "reads.fasta"
            path.write_text("ACGT\n>read\nACGT\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "before first header"):
                list(iter_fasta_raw(input_path=path, decompressor="python"))

    def test_fastq_writer_requires_quality(self) -> None:
        """FASTQ output cannot silently invent missing qualities."""
        handle = io.StringIO()
        record = RawSequenceRecord(
            identifier="read",
            description="read",
            sequence="ACGT",
            quality=None,
        )
        with self.assertRaisesRegex(ValueError, "quality string"):
            _write_record(
                handle=handle,
                record=record,
                input_format="fastq",
            )

    def test_uncompressed_fasta_subsets_and_explicit_manifest(self) -> None:
        """FASTA subsets should support plain text and named manifest paths."""
        logger = Mock(spec=logging.Logger)
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "reads.fasta"
            source.write_text(">one\nACGT\n>two detail\nTGCA\n", encoding="utf-8")
            manifest_path = root / "reports" / "manifest.tsv"
            manifest = create_depth_subsets(
                input_path=source,
                input_format="fasta",
                output_dir=root / "out",
                fractions=[1.0],
                seeds=[7],
                sample_prefix="sample",
                compress=False,
                manifest_path=manifest_path,
                progress_interval=1,
                logger=logger,
            )
            subset_path = Path(str(manifest[0]["reads_path"]))
            observed = subset_path.read_text(encoding="utf-8")
            self.assertTrue(manifest_path.is_file())
        self.assertEqual(observed, ">one\nACGT\n>two detail\nTGCA\n")
        self.assertEqual(manifest[0]["n_input_records"], 2)
        self.assertEqual(manifest[0]["n_selected_records"], 2)
        self.assertGreaterEqual(logger.info.call_count, 3)

    def test_create_depth_subsets_validates_inputs(self) -> None:
        """Missing inputs and invalid run settings should fail early."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            missing = root / "missing.fastq"
            with self.assertRaisesRegex(FileNotFoundError, "missing or empty"):
                create_depth_subsets(
                    input_path=missing,
                    input_format="fastq",
                    output_dir=root / "out",
                    fractions=[1.0],
                    seeds=[1],
                    sample_prefix="sample",
                )
            source = self.write_fastq(root, n_records=1)
            with self.assertRaisesRegex(ValueError, "input_format"):
                create_depth_subsets(
                    input_path=source,
                    input_format="bam",
                    output_dir=root / "out",
                    fractions=[1.0],
                    seeds=[1],
                    sample_prefix="sample",
                )
            with self.assertRaisesRegex(ValueError, "progress_interval"):
                create_depth_subsets(
                    input_path=source,
                    input_format="fastq",
                    output_dir=root / "out",
                    fractions=[1.0],
                    seeds=[1],
                    sample_prefix="sample",
                    progress_interval=0,
                )

    def test_no_record_input_removes_temporary_outputs(self) -> None:
        """A record-free input should leave no partial subset files."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "reads.fasta"
            source.write_text("\n\n", encoding="utf-8")
            output = root / "out"
            with self.assertRaisesRegex(ValueError, "contains no records"):
                create_depth_subsets(
                    input_path=source,
                    input_format="fasta",
                    output_dir=output,
                    fractions=[1.0],
                    seeds=[1],
                    sample_prefix="sample",
                    compress=False,
                )
            self.assertEqual(list(output.glob(".*.tmp")), [])


if __name__ == "__main__":
    unittest.main()
