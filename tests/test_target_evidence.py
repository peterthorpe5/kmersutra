"""Tests for SQLite-backed target-evidence panel building."""

from __future__ import annotations

import sqlite3
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from kmersutra.config import GenomeConfig
from kmersutra.target_evidence import (
    build_target_evidence_sqlite,
    iter_target_evidence_diagnostics,
    summarise_diagnostics_stream,
)
from kmersutra.taxonomy import TaxonomyDatabase


def write_test_taxdump(root: Path) -> None:
    """Write a minimal taxonomy dump for target-evidence tests.

    Parameters
    ----------
    root : pathlib.Path
        Directory where taxdump files are written.
    """
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
        "11\t|\tAlpha target\t|\t\t|\tscientific name\t|\n"
        "12\t|\tAlpha neighbour\t|\t\t|\tscientific name\t|\n"
        "20\t|\tBetaGenus\t|\t\t|\tscientific name\t|\n"
        "21\t|\tBeta outgroup\t|\t\t|\tscientific name\t|\n",
        encoding="utf-8",
    )
    (root / "merged.dmp").write_text("", encoding="utf-8")
    (root / "delnodes.dmp").write_text("", encoding="utf-8")


class TestTargetEvidenceBuild(unittest.TestCase):
    """Tests for low-memory target-evidence building."""

    def test_sqlite_build_marks_non_target_overlap(self) -> None:
        """Non-target genomes should mark matching target candidates only."""
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            target = root / "target.fna"
            neighbour = root / "neighbour.fna"
            outgroup = root / "outgroup.fna"
            target.write_text(">target\nAAAAACCCCCGGGGG\n", encoding="utf-8")
            neighbour.write_text(">neighbour\nAAAAATTTTT\n", encoding="utf-8")
            outgroup.write_text(">outgroup\nACGTGACGTG\n", encoding="utf-8")
            sqlite_path = root / "candidates.sqlite"
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
            result = build_target_evidence_sqlite(
                genome_configs=configs,
                k_values=[5],
                sqlite_path=sqlite_path,
                batch_size=3,
            )
            connection = sqlite3.connect(sqlite_path)
            try:
                n_candidates = connection.execute(
                    "SELECT COUNT(*) FROM target_kmers"
                ).fetchone()[0]
                n_shared = connection.execute(
                    "SELECT COUNT(*) FROM target_kmers WHERE non_target_taxids != ''"
                ).fetchone()[0]
            finally:
                connection.close()
        self.assertGreater(n_candidates, 0)
        self.assertGreater(n_shared, 0)
        self.assertEqual(result.sqlite_path, sqlite_path)
        self.assertEqual(len(result.collection_summary), 3)

    def test_sqlite_build_yields_species_and_genus_evidence(self) -> None:
        """Target-only and shared target-neighbour k-mers should be ranked."""
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            taxdump = root / "taxdump"
            write_test_taxdump(taxdump)
            taxonomy = TaxonomyDatabase.from_taxdump(taxonomy_dir=taxdump)
            target = root / "target.fna"
            neighbour = root / "neighbour.fna"
            outgroup = root / "outgroup.fna"
            target.write_text(">target\nAAAAACCCCCGGGGG\n", encoding="utf-8")
            neighbour.write_text(">neighbour\nAAAAATTTTT\n", encoding="utf-8")
            outgroup.write_text(">outgroup\nACGTGACGTG\n", encoding="utf-8")
            sqlite_path = root / "candidates.sqlite"
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
            build_target_evidence_sqlite(
                genome_configs=configs,
                k_values=[5],
                sqlite_path=sqlite_path,
                batch_size=2,
            )
            diagnostics = list(
                iter_target_evidence_diagnostics(
                    sqlite_path=sqlite_path,
                    taxonomy_db=taxonomy,
                    target_taxid="10",
                    preferred_ranks=["species", "genus"],
                    chunk_size=2,
                )
            )
            ranks = {item.evidence_rank for item in diagnostics}
        self.assertIn("species", ranks)
        self.assertIn("genus", ranks)
        self.assertTrue(all(item.evidence_taxid in {"10", "11"} for item in diagnostics))

    def test_summarise_diagnostics_stream_applies_limit(self) -> None:
        """Diagnostic stream summarisation should enforce retention limits."""
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            target = root / "target.fna"
            target.write_text(">target\nAAAAACCCCCGGGGGTTTTT\n", encoding="utf-8")
            sqlite_path = root / "candidates.sqlite"
            configs = [
                GenomeConfig(
                    genome_fasta=target,
                    species_name="Alpha target",
                    taxid="11",
                    role="target_species",
                    clade="AlphaGenus",
                )
            ]
            build_target_evidence_sqlite(
                genome_configs=configs,
                k_values=[5],
                sqlite_path=sqlite_path,
                batch_size=2,
            )
            diagnostics, summary = summarise_diagnostics_stream(
                diagnostics=iter_target_evidence_diagnostics(sqlite_path=sqlite_path),
                max_per_evidence_per_k=2,
            )
        self.assertEqual(len(diagnostics), 2)
        self.assertEqual(summary[0]["diagnostic_kmers"], 2)


if __name__ == "__main__":
    unittest.main()
