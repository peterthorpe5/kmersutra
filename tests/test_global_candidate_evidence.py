"""Tests for the scalable global all-candidate evidence builder."""

from __future__ import annotations

import sqlite3
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from kmersutra.config import GenomeConfig
from kmersutra.global_candidate_evidence import (
    build_global_candidate_evidence_sqlite,
    collect_global_kmer_sources_sqlite,
    iter_retained_global_candidate_diagnostics,
)
from kmersutra.taxonomy import TaxonomyDatabase
from tests.test_target_evidence import write_test_taxdump


class TestGlobalCandidateEvidenceBuild(unittest.TestCase):
    """Tests for global query-agnostic evidence building."""

    def test_global_builder_indexes_each_genome_once(self) -> None:
        """Global collection should make one collection record per genome."""
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            alpha = root / "alpha.fna"
            beta = root / "beta.fna"
            alpha.write_text(">alpha\nAAAAACCCCC\n", encoding="utf-8")
            beta.write_text(">beta\nTTTTTGGGGG\n", encoding="utf-8")
            configs = [
                GenomeConfig(
                    genome_fasta=alpha,
                    species_name="Alpha target",
                    taxid="11",
                    role="target_species",
                    clade="AlphaGenus",
                ),
                GenomeConfig(
                    genome_fasta=beta,
                    species_name="Beta outgroup",
                    taxid="21",
                    role="outgroup",
                    clade="BetaGenus",
                ),
            ]
            sqlite_path = root / "global.sqlite"
            collect_global_kmer_sources_sqlite(
                genome_configs=configs,
                k_values=[5],
                sqlite_path=sqlite_path,
                batch_size=2,
            )
            connection = sqlite3.connect(sqlite_path)
            try:
                n_events = connection.execute(
                    "SELECT COUNT(*) FROM build_events WHERE stage='collect_global_sources'"
                ).fetchone()[0]
                n_keys = connection.execute("SELECT COUNT(*) FROM global_kmers").fetchone()[0]
            finally:
                connection.close()
        self.assertEqual(n_events, 2)
        self.assertGreater(n_keys, 0)

    def test_global_builder_retains_species_genus_and_outgroup_evidence(self) -> None:
        """Global mode should retain evidence for all reportable taxa."""
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            taxdump = root / "taxdump"
            write_test_taxdump(taxdump)
            taxonomy = TaxonomyDatabase.from_taxdump(taxonomy_dir=taxdump)

            target = root / "target.fna"
            neighbour = root / "neighbour.fna"
            outgroup = root / "outgroup.fna"
            target.write_text(">target\nAAAAACCCCCGGGGG\n", encoding="utf-8")
            neighbour.write_text(">neighbour\nAAAAATTTTTCCCCC\n", encoding="utf-8")
            outgroup.write_text(">outgroup\nACGTACGTACGT\n", encoding="utf-8")
            configs = [
                GenomeConfig(
                    genome_fasta=target,
                    species_name="Alpha target",
                    taxid="11",
                    role="target_species",
                    clade="AlphaGenus",
                ),
                GenomeConfig(
                    genome_fasta=neighbour,
                    species_name="Alpha neighbour",
                    taxid="12",
                    role="near_neighbour",
                    clade="AlphaGenus",
                ),
                GenomeConfig(
                    genome_fasta=outgroup,
                    species_name="Beta outgroup",
                    taxid="21",
                    role="outgroup",
                    clade="BetaGenus",
                ),
            ]
            result = build_global_candidate_evidence_sqlite(
                genome_configs=configs,
                k_values=[5],
                sqlite_path=root / "global.sqlite",
                taxonomy_db=taxonomy,
                preferred_ranks=["species", "genus", "superkingdom"],
                batch_size=3,
                max_per_evidence_per_k=20,
            )
            diagnostics = list(
                iter_retained_global_candidate_diagnostics(sqlite_path=result.sqlite_path)
            )
            species_names = {item.species_name for item in diagnostics if item.species_name}
            evidence_names = {item.evidence_name for item in diagnostics}

        self.assertIn("Alpha target", species_names)
        self.assertIn("Alpha neighbour", species_names)
        self.assertIn("Beta outgroup", species_names)
        self.assertIn("AlphaGenus", evidence_names)
        self.assertGreater(len(result.panel_summary), 0)
        self.assertTrue(
            any(row["summary_name"] == "global_distinct_kmer_keys" for row in result.build_summary)
        )

    def test_global_builder_applies_evidence_limit(self) -> None:
        """Per-evidence caps should be enforced during global assignment."""
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            taxdump = root / "taxdump"
            write_test_taxdump(taxdump)
            taxonomy = TaxonomyDatabase.from_taxdump(taxonomy_dir=taxdump)
            target = root / "target.fna"
            target.write_text(">target\nAAAAACCCCCGGGGGTTTTT\n", encoding="utf-8")
            configs = [
                GenomeConfig(
                    genome_fasta=target,
                    species_name="Alpha target",
                    taxid="11",
                    role="target_species",
                    clade="AlphaGenus",
                ),
            ]
            result = build_global_candidate_evidence_sqlite(
                genome_configs=configs,
                k_values=[5],
                sqlite_path=root / "global.sqlite",
                taxonomy_db=taxonomy,
                preferred_ranks=["species", "genus"],
                batch_size=2,
                max_per_evidence_per_k=2,
            )
            connection = sqlite3.connect(result.sqlite_path)
            try:
                max_count = connection.execute(
                    """
                    SELECT MAX(n_records) FROM (
                        SELECT evidence_taxid, evidence_rank, k, COUNT(*) AS n_records
                        FROM retained_kmers
                        GROUP BY evidence_taxid, evidence_rank, k
                    )
                    """
                ).fetchone()[0]
            finally:
                connection.close()
        self.assertLessEqual(max_count, 2)


if __name__ == "__main__":
    unittest.main()

class TestGlobalCandidateEvidenceSourceRows(unittest.TestCase):
    """Tests for the faster source-row global index mode."""

    def test_source_row_mode_materialises_global_kmers(self) -> None:
        """Source-row mode should build source rows and aggregated keys."""
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            genome = root / "alpha.fna"
            genome.write_text(">alpha\nAAAAACCCCCAAAAA\n", encoding="utf-8")
            configs = [
                GenomeConfig(
                    genome_fasta=genome,
                    species_name="Alpha target",
                    taxid="11",
                    role="target_species",
                    clade="AlphaGenus",
                )
            ]
            sqlite_path = root / "global.sqlite"
            summary = collect_global_kmer_sources_sqlite(
                genome_configs=configs,
                k_values=[5],
                sqlite_path=sqlite_path,
                batch_size=2,
                source_index_mode="source_rows",
                progress_interval=3,
            )
            connection = sqlite3.connect(sqlite_path)
            try:
                n_source_rows = connection.execute(
                    "SELECT COUNT(*) FROM global_kmer_sources"
                ).fetchone()[0]
                n_global_rows = connection.execute(
                    "SELECT COUNT(*) FROM global_kmers"
                ).fetchone()[0]
                materialise_events = connection.execute(
                    """
                    SELECT COUNT(*) FROM build_events
                    WHERE stage='materialise_global_sources'
                    """
                ).fetchone()[0]
            finally:
                connection.close()

        self.assertGreater(n_source_rows, 0)
        self.assertGreater(n_global_rows, 0)
        self.assertGreaterEqual(n_source_rows, n_global_rows)
        self.assertEqual(materialise_events, 1)
        self.assertTrue(
            any(row["stage"] == "materialise_global_sources" for row in summary)
        )

    def test_source_row_and_aggregated_modes_retain_same_evidence(self) -> None:
        """The faster source-row mode should preserve retained evidence calls."""
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            taxdump = root / "taxdump"
            write_test_taxdump(taxdump)
            taxonomy = TaxonomyDatabase.from_taxdump(taxonomy_dir=taxdump)
            target = root / "target.fna"
            neighbour = root / "neighbour.fna"
            target.write_text(">target\nAAAAACCCCCGGGGG\n", encoding="utf-8")
            neighbour.write_text(">neighbour\nAAAAATTTTTCCCCC\n", encoding="utf-8")
            configs = [
                GenomeConfig(
                    genome_fasta=target,
                    species_name="Alpha target",
                    taxid="11",
                    role="target_species",
                    clade="AlphaGenus",
                ),
                GenomeConfig(
                    genome_fasta=neighbour,
                    species_name="Alpha neighbour",
                    taxid="12",
                    role="near_neighbour",
                    clade="AlphaGenus",
                ),
            ]
            source_result = build_global_candidate_evidence_sqlite(
                genome_configs=configs,
                k_values=[5],
                sqlite_path=root / "source.sqlite",
                taxonomy_db=taxonomy,
                preferred_ranks=["species", "genus"],
                batch_size=3,
                max_per_evidence_per_k=50,
                source_index_mode="source_rows",
            )
            aggregated_result = build_global_candidate_evidence_sqlite(
                genome_configs=configs,
                k_values=[5],
                sqlite_path=root / "aggregated.sqlite",
                taxonomy_db=taxonomy,
                preferred_ranks=["species", "genus"],
                batch_size=3,
                max_per_evidence_per_k=50,
                source_index_mode="aggregated",
            )
            source_evidence = {
                (item.kmer, item.k, item.evidence_name, item.evidence_rank)
                for item in iter_retained_global_candidate_diagnostics(
                    sqlite_path=source_result.sqlite_path
                )
            }
            aggregated_evidence = {
                (item.kmer, item.k, item.evidence_name, item.evidence_rank)
                for item in iter_retained_global_candidate_diagnostics(
                    sqlite_path=aggregated_result.sqlite_path
                )
            }

        self.assertEqual(source_evidence, aggregated_evidence)

    def test_global_builder_rejects_bad_source_index_mode(self) -> None:
        """Invalid source-index mode values should fail clearly."""
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            taxdump = root / "taxdump"
            write_test_taxdump(taxdump)
            taxonomy = TaxonomyDatabase.from_taxdump(taxonomy_dir=taxdump)
            genome = root / "target.fna"
            genome.write_text(">target\nAAAAACCCCC\n", encoding="utf-8")
            configs = [
                GenomeConfig(
                    genome_fasta=genome,
                    species_name="Alpha target",
                    taxid="11",
                    role="target_species",
                    clade="AlphaGenus",
                )
            ]
            with self.assertRaisesRegex(ValueError, "source_index_mode"):
                build_global_candidate_evidence_sqlite(
                    genome_configs=configs,
                    k_values=[5],
                    sqlite_path=root / "bad.sqlite",
                    taxonomy_db=taxonomy,
                    source_index_mode="not_a_mode",
                )

class TestGlobalCandidateEvidenceAssignmentOptimisation(unittest.TestCase):
    """Tests for cached, batched global evidence assignment."""

    def test_assignment_reuses_taxonomy_result_for_repeated_taxid_sets(self) -> None:
        """Repeated taxid sets should use the cached evidence assignment."""
        from kmersutra.global_candidate_evidence import (
            assign_global_candidate_evidence_sqlite,
            initialise_global_candidate_database,
        )

        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            taxdump = root / "taxdump"
            write_test_taxdump(taxdump)
            taxonomy = TaxonomyDatabase.from_taxdump(taxonomy_dir=taxdump)
            sqlite_path = root / "global.sqlite"
            initialise_global_candidate_database(sqlite_path=sqlite_path)
            connection = sqlite3.connect(sqlite_path)
            try:
                connection.executemany(
                    """
                    INSERT INTO global_kmers(
                        k, kmer, species_names, genome_ids, contig_ids,
                        taxids, clades, roles, example_position
                    )
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
                    """,
                    [
                        (5, "AAAAA", "Alpha target", "G1", "c1", "11", "AlphaGenus", "target_species", 1),
                        (5, "AAAAC", "Alpha target", "G1", "c1", "11", "AlphaGenus", "target_species", 2),
                        (5, "AAACC", "Alpha target", "G1", "c1", "11", "AlphaGenus", "target_species", 3),
                    ],
                )
                connection.commit()
            finally:
                connection.close()

            original = taxonomy.best_named_ancestor
            call_count = {"n": 0}

            def counted_best_named_ancestor(*, taxids, preferred_ranks=None):
                call_count["n"] += 1
                return original(taxids=taxids, preferred_ranks=preferred_ranks)

            taxonomy.best_named_ancestor = counted_best_named_ancestor  # type: ignore[method-assign]
            summary = assign_global_candidate_evidence_sqlite(
                sqlite_path=sqlite_path,
                taxonomy_db=taxonomy,
                preferred_ranks=["species", "genus"],
                batch_size=2,
                max_per_evidence_per_k=10,
            )
            retained = {
                row["metric"]: row["value"]
                for row in summary
                if row["stage"] == "assign_global_evidence"
            }

        self.assertEqual(call_count["n"], 1)
        self.assertEqual(retained["diagnostics_retained"], 3)

