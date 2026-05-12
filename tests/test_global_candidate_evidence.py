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
