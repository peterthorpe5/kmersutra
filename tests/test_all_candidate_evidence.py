"""Tests for low-memory all-candidate evidence building."""

from __future__ import annotations

import sqlite3
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from kmersutra.all_candidate_evidence import (
    build_all_candidate_evidence_sqlite,
    iter_retained_all_candidate_diagnostics,
    select_candidate_genomes,
)
from kmersutra.config import GenomeConfig
from kmersutra.taxonomy import TaxonomyDatabase
from tests.test_target_evidence import write_test_taxdump


class TestAllCandidateEvidenceBuild(unittest.TestCase):
    """Tests for query-agnostic all-candidate evidence building."""

    def test_select_candidate_genomes_excludes_background_by_default(self) -> None:
        """Default candidate selection should exclude background roles."""
        configs = [
            GenomeConfig(
                genome_fasta=Path("a.fna"),
                species_name="Alpha target",
                role="target_species",
            ),
            GenomeConfig(
                genome_fasta=Path("b.fna"),
                species_name="Host background",
                role="host_or_background",
            ),
        ]
        selected = select_candidate_genomes(genome_configs=configs)
        self.assertEqual([item.species_name for item in selected], ["Alpha target"])

    def test_all_candidate_build_retains_multiple_species(self) -> None:
        """All-candidate mode should retain reportable evidence for each species."""
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            taxdump = root / "taxdump"
            write_test_taxdump(taxdump)
            taxonomy = TaxonomyDatabase.from_taxdump(taxonomy_dir=taxdump)

            target = root / "target.fna"
            neighbour = root / "neighbour.fna"
            outgroup = root / "outgroup.fna"
            target.write_text(">target\nAAAAACCCCCGGGGG\n", encoding="utf-8")
            neighbour.write_text(">neighbour\nTTTTTCCCCCAAAAA\n", encoding="utf-8")
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

            result = build_all_candidate_evidence_sqlite(
                genome_configs=configs,
                k_values=[5],
                retained_sqlite_path=root / "retained.sqlite",
                work_sqlite_path=root / "work.sqlite",
                taxonomy_db=taxonomy,
                preferred_ranks=["species", "genus", "superkingdom"],
                batch_size=3,
                max_per_evidence_per_k=10,
            )
            diagnostics = list(
                iter_retained_all_candidate_diagnostics(
                    sqlite_path=result.retained_sqlite_path,
                )
            )
            species_names = {item.species_name for item in diagnostics if item.species_name}
            evidence_names = {item.evidence_name for item in diagnostics}

        self.assertIn("Alpha target", species_names)
        self.assertIn("Alpha neighbour", species_names)
        self.assertIn("Beta outgroup", species_names)
        self.assertIn("AlphaGenus", evidence_names)
        self.assertGreater(len(result.panel_summary), 0)
        self.assertTrue(
            any(row["summary_name"] == "candidate_species_groups" for row in result.build_summary)
        )

    def test_all_candidate_build_applies_global_limit(self) -> None:
        """Per-evidence limits should be applied across candidate rounds."""
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            taxdump = root / "taxdump"
            write_test_taxdump(taxdump)
            taxonomy = TaxonomyDatabase.from_taxdump(taxonomy_dir=taxdump)
            target = root / "target.fna"
            neighbour = root / "neighbour.fna"
            target.write_text(">target\nAAAAACCCCCGGGGGTTTTT\n", encoding="utf-8")
            neighbour.write_text(">neighbour\nTTTTTGGGGGCCCCCAAAAA\n", encoding="utf-8")
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
            result = build_all_candidate_evidence_sqlite(
                genome_configs=configs,
                k_values=[5],
                retained_sqlite_path=root / "retained.sqlite",
                work_sqlite_path=root / "work.sqlite",
                taxonomy_db=taxonomy,
                preferred_ranks=["species", "genus"],
                batch_size=2,
                max_per_evidence_per_k=2,
            )
            connection = sqlite3.connect(result.retained_sqlite_path)
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
