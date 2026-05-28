"""Tests for optional Parquet module helpers."""

from __future__ import annotations

import sqlite3
import tempfile
import unittest
from pathlib import Path

from kmersutra.global_candidate_evidence import initialise_global_candidate_database
from kmersutra.parquet_modules import (
    GLOBAL_KMER_COLUMNS,
    OptionalParquetDependencyError,
    deduplicate_global_kmer_sqlite,
    export_sqlite_table_to_parquet,
    import_global_kmer_parquets_to_sqlite,
    join_semicolon,
    merge_global_kmer_record_group,
    merge_global_kmer_records,
    normalise_global_kmer_record,
    require_pyarrow,
    resolve_global_kmer_parquet_paths,
    split_semicolon,
    upsert_global_kmer_record,
)
from kmersutra.target_evidence import _connect


class TestParquetModuleHelpers(unittest.TestCase):
    """Test pure helpers for module-aware k-mer source tables."""

    def test_split_semicolon_removes_empty_values(self) -> None:
        """Semicolon splitting should trim values and remove blanks."""
        self.assertEqual(split_semicolon(" a ;b;; a ; "), {"a", "b"})

    def test_join_semicolon_is_sorted_and_unique(self) -> None:
        """Semicolon joining should be deterministic."""
        self.assertEqual(join_semicolon(["b", "a", "b", ""]), "a;b")

    def test_normalise_global_kmer_record_fills_schema(self) -> None:
        """Global-kmer records should be normalised to the source schema."""
        record = normalise_global_kmer_record(
            {
                "k": "77",
                "kmer": "ACGT",
                "species_names": "Species B;Species A;Species A",
                "example_position": "12",
            }
        )
        self.assertEqual(set(record), set(GLOBAL_KMER_COLUMNS))
        self.assertEqual(record["k"], 77)
        self.assertEqual(record["species_names"], "Species A;Species B")
        self.assertEqual(record["example_position"], 12)

    def test_merge_global_kmer_record_group_unions_source_metadata(self) -> None:
        """Duplicate k-mer records should be merged by source metadata."""
        merged = merge_global_kmer_record_group(
            [
                {
                    "k": 77,
                    "kmer": "AAAA",
                    "species_names": "Species A",
                    "genome_ids": "G1",
                    "contig_ids": "c1",
                    "taxids": "1",
                    "clades": "CladeA",
                    "roles": "near_neighbour",
                    "example_position": 50,
                },
                {
                    "k": 77,
                    "kmer": "AAAA",
                    "species_names": "Species B",
                    "genome_ids": "G2",
                    "contig_ids": "c2",
                    "taxids": "2",
                    "clades": "CladeA",
                    "roles": "outgroup",
                    "example_position": 20,
                },
            ]
        )
        self.assertEqual(merged["species_names"], "Species A;Species B")
        self.assertEqual(merged["genome_ids"], "G1;G2")
        self.assertEqual(merged["taxids"], "1;2")
        self.assertEqual(merged["roles"], "near_neighbour;outgroup")
        self.assertEqual(merged["example_position"], 20)

    def test_merge_global_kmer_record_group_rejects_mixed_keys(self) -> None:
        """A record group should not silently merge different k-mers."""
        with self.assertRaises(ValueError):
            merge_global_kmer_record_group(
                [
                    {"k": 77, "kmer": "AAAA"},
                    {"k": 77, "kmer": "CCCC"},
                ]
            )

    def test_merge_global_kmer_records_groups_duplicate_keys(self) -> None:
        """In-memory merging should group by k and kmer."""
        records = merge_global_kmer_records(
            [
                {"k": 77, "kmer": "AAAA", "species_names": "A"},
                {"k": 101, "kmer": "AAAA", "species_names": "LongA"},
                {"k": 77, "kmer": "AAAA", "species_names": "B"},
            ]
        )
        self.assertEqual(len(records), 2)
        self.assertEqual(records[0]["species_names"], "A;B")
        self.assertEqual(records[1]["species_names"], "LongA")

    def test_resolve_global_kmer_parquet_paths_requires_input(self) -> None:
        """Module path resolution should require at least one input."""
        with self.assertRaises(ValueError):
            resolve_global_kmer_parquet_paths(module_dirs=[], global_kmer_parquets=[])

    def test_resolve_global_kmer_parquet_paths_checks_missing_files(self) -> None:
        """Missing module parquet files should fail loudly."""
        with tempfile.TemporaryDirectory() as tmpdir:
            with self.assertRaises(FileNotFoundError):
                resolve_global_kmer_parquet_paths(module_dirs=[tmpdir])

    def test_upsert_global_kmer_record_merges_in_sqlite(self) -> None:
        """SQLite upsert should merge duplicate module source metadata."""
        with tempfile.TemporaryDirectory() as tmpdir:
            sqlite_path = Path(tmpdir) / "module.sqlite"
            initialise_global_candidate_database(sqlite_path=sqlite_path)
            connection = _connect(sqlite_path=sqlite_path)
            try:
                upsert_global_kmer_record(
                    connection=connection,
                    record={
                        "k": 77,
                        "kmer": "AAAA",
                        "species_names": "A",
                        "genome_ids": "G1",
                        "taxids": "1",
                        "roles": "near_neighbour",
                        "example_position": 100,
                    },
                )
                upsert_global_kmer_record(
                    connection=connection,
                    record={
                        "k": 77,
                        "kmer": "AAAA",
                        "species_names": "B",
                        "genome_ids": "G2",
                        "taxids": "2",
                        "roles": "outgroup",
                        "example_position": 10,
                    },
                )
                connection.commit()
            finally:
                connection.close()
            deduplicate_global_kmer_sqlite(sqlite_path=sqlite_path)
            connection = sqlite3.connect(sqlite_path)
            try:
                row = connection.execute(
                    "SELECT species_names, genome_ids, taxids, roles, example_position FROM global_kmers"
                ).fetchone()
            finally:
                connection.close()
            self.assertEqual(row[0], "A;B")
            self.assertEqual(row[1], "G1;G2")
            self.assertEqual(row[2], "1;2")
            self.assertEqual(row[3], "near_neighbour;outgroup")
            self.assertEqual(row[4], 10)


class TestParquetOptionalDependency(unittest.TestCase):
    """Test behaviour when optional pyarrow is unavailable or available."""

    def test_require_pyarrow_returns_modules_or_clear_error(self) -> None:
        """The optional dependency helper should either import or explain."""
        try:
            pa, pq = require_pyarrow()
        except OptionalParquetDependencyError as exc:
            self.assertIn("pyarrow", str(exc))
        else:
            self.assertTrue(hasattr(pa, "Table"))
            self.assertTrue(hasattr(pq, "ParquetFile"))

    def test_export_table_without_pyarrow_gives_clear_error(self) -> None:
        """Parquet export should fail clearly when pyarrow is unavailable."""
        with tempfile.TemporaryDirectory() as tmpdir:
            sqlite_path = Path(tmpdir) / "module.sqlite"
            initialise_global_candidate_database(sqlite_path=sqlite_path)
            out_path = Path(tmpdir) / "global_kmers.parquet"
            try:
                require_pyarrow()
            except OptionalParquetDependencyError:
                with self.assertRaises(OptionalParquetDependencyError):
                    export_sqlite_table_to_parquet(
                        sqlite_path=sqlite_path,
                        table_name="global_kmers",
                        output_path=out_path,
                    )
            else:
                n_records = export_sqlite_table_to_parquet(
                    sqlite_path=sqlite_path,
                    table_name="global_kmers",
                    output_path=out_path,
                )
                self.assertEqual(n_records, 0)
                self.assertTrue(out_path.exists())

    def test_import_without_pyarrow_gives_clear_error(self) -> None:
        """Parquet import should fail clearly when pyarrow is unavailable."""
        try:
            require_pyarrow()
        except OptionalParquetDependencyError:
            with tempfile.TemporaryDirectory() as tmpdir:
                with self.assertRaises(OptionalParquetDependencyError):
                    import_global_kmer_parquets_to_sqlite(
                        parquet_paths=[Path(tmpdir) / "missing.parquet"],
                        sqlite_path=Path(tmpdir) / "merged.sqlite",
                    )
        else:
            self.skipTest("pyarrow is installed; missing-file behaviour is exercised elsewhere")


if __name__ == "__main__":
    unittest.main()
