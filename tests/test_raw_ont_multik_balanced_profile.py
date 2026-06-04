"""Tests for raw ONT multi-k balanced candidate sampling."""

from __future__ import annotations

import sqlite3
import tempfile
import unittest
from pathlib import Path

from kmersutra.cli.build_clade_kmer_panel import (
    resolve_candidate_k_order,
    resolve_max_per_genome_bin_by_k,
    resolve_min_cross_k_marker_distance,
)
from kmersutra.config import GenomeConfig
from kmersutra.global_candidate_evidence import (
    collect_candidate_universe_sqlite,
    parse_max_per_genome_bin_by_k_spec,
    resolve_candidate_bin_quota,
)


class TestRawOntMultikBalancedProfile(unittest.TestCase):
    """Test the raw ONT multi-k balanced build controls."""

    def test_parse_max_per_genome_bin_by_k_spec(self) -> None:
        """Per-k quota strings should parse to integer mappings."""
        observed = parse_max_per_genome_bin_by_k_spec(
            spec="51:3,77:3,101:2,151:2"
        )
        self.assertEqual(observed, {51: 3, 77: 3, 101: 2, 151: 2})

    def test_parse_max_per_genome_bin_by_k_rejects_bad_values(self) -> None:
        """Malformed quota strings should fail with useful errors."""
        with self.assertRaisesRegex(ValueError, "k:quota"):
            parse_max_per_genome_bin_by_k_spec(spec="51=3")
        with self.assertRaisesRegex(ValueError, "positive"):
            parse_max_per_genome_bin_by_k_spec(spec="51:0")

    def test_profile_defaults_are_multik_balanced(self) -> None:
        """The raw ONT multi-k profile should use per-k quotas by default."""
        self.assertEqual(
            resolve_candidate_k_order(
                marker_profile="raw_ont_multik_balanced",
                candidate_k_order="auto",
            ),
            "input",
        )
        self.assertEqual(
            resolve_min_cross_k_marker_distance(
                marker_profile="raw_ont_multik_balanced",
                min_cross_k_marker_distance=None,
            ),
            0,
        )
        self.assertEqual(
            resolve_max_per_genome_bin_by_k(
                marker_profile="raw_ont_multik_balanced",
                max_per_genome_bin_by_k="",
            ),
            {51: 3, 77: 3, 101: 2, 151: 2},
        )

    def test_explicit_quota_overrides_profile_default(self) -> None:
        """User-supplied quotas should override profile defaults."""
        self.assertEqual(
            resolve_max_per_genome_bin_by_k(
                marker_profile="raw_ont_multik_balanced",
                max_per_genome_bin_by_k="51:4,77:3,101:2,151:1",
            ),
            {51: 4, 77: 3, 101: 2, 151: 1},
        )

    def test_resolve_candidate_bin_quota_falls_back(self) -> None:
        """K values absent from the mapping should use the global fallback."""
        self.assertEqual(
            resolve_candidate_bin_quota(
                k=77,
                max_per_genome_bin=10,
                max_per_genome_bin_by_k={51: 4},
            ),
            10,
        )
        self.assertEqual(
            resolve_candidate_bin_quota(
                k=51,
                max_per_genome_bin=10,
                max_per_genome_bin_by_k={51: 4},
            ),
            4,
        )

    def test_candidate_universe_retains_every_k_with_per_k_quotas(self) -> None:
        """Short k values should not starve longer k values under per-k quotas."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)
            fasta_path = tmp_path / "genome.fna"
            sqlite_path = tmp_path / "candidates.sqlite"
            fasta_path.write_text(
                ">contig1\n" + "ACGT" * 80 + "\n",
                encoding="utf-8",
            )
            config = GenomeConfig(
                genome_fasta=fasta_path,
                species_name="Example species",
                strain_name="strain1",
                taxid="12345",
                role="target_species",
                clade="example",
                assembly_accession="ASM1",
            )

            collect_candidate_universe_sqlite(
                genome_configs=[config],
                k_values=[5, 7, 9, 11],
                sqlite_path=sqlite_path,
                batch_size=100,
                genome_bin_size=50,
                max_per_genome_bin=1,
                max_per_genome_bin_by_k={5: 2, 7: 2, 9: 2, 11: 2},
                min_cross_k_marker_distance=0,
                assembly_aware_binning=False,
                progress_interval=1000000,
                candidate_k_order="input",
            )

            with sqlite3.connect(sqlite_path) as connection:
                rows = connection.execute(
                    """
                    SELECT k, COUNT(*)
                    FROM candidate_kmers
                    GROUP BY k
                    ORDER BY k
                    """
                ).fetchall()

            counts = dict(rows)
            self.assertEqual(set(counts), {5, 7, 9, 11})
            for count in counts.values():
                self.assertGreater(count, 0)


if __name__ == "__main__":
    unittest.main()
