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


class TestRawOntMultikLocusBalancedProfile(unittest.TestCase):
    """Test the raw ONT multi-k locus-balanced build controls."""

    def test_locus_balanced_profile_defaults(self) -> None:
        """The locus-balanced profile should avoid nested cross-k markers."""
        self.assertEqual(
            resolve_candidate_k_order(
                marker_profile="raw_ont_multik_locus_balanced",
                candidate_k_order="auto",
            ),
            "input",
        )
        self.assertEqual(
            resolve_min_cross_k_marker_distance(
                marker_profile="raw_ont_multik_locus_balanced",
                min_cross_k_marker_distance=None,
            ),
            250,
        )
        self.assertEqual(
            resolve_max_per_genome_bin_by_k(
                marker_profile="raw_ont_multik_locus_balanced",
                max_per_genome_bin_by_k="",
            ),
            {51: 3, 77: 3, 101: 2, 151: 2},
        )

    def test_candidate_universe_locus_separation_keeps_all_k_values(self) -> None:
        """Per-k quotas plus small locus spacing should retain all k values."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)
            fasta_path = tmp_path / "genome.fna"
            sqlite_path = tmp_path / "candidates.sqlite"
            fasta_path.write_text(
                ">contig1\n" + "ACGT" * 300 + "\n",
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
                genome_bin_size=100,
                max_per_genome_bin=1,
                max_per_genome_bin_by_k={5: 2, 7: 2, 9: 2, 11: 2},
                min_cross_k_marker_distance=1,
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
                positions = connection.execute(
                    """
                    SELECT first_contig_id, first_position, k
                    FROM candidate_kmers
                    ORDER BY first_position, k
                    """
                ).fetchall()

            counts = dict(rows)
            self.assertEqual(set(counts), {5, 7, 9, 11})
            for count in counts.values():
                self.assertGreater(count, 0)

            for index, left in enumerate(positions):
                for right in positions[index + 1:]:
                    left_contig, left_pos, left_k = left
                    right_contig, right_pos, right_k = right
                    if left_contig != right_contig or int(left_k) == int(right_k):
                        continue
                    left_start = int(left_pos)
                    left_end = left_start + int(left_k)
                    right_start = int(right_pos)
                    right_end = right_start + int(right_k)
                    gap = max(0, max(left_start, right_start) - min(left_end, right_end))
                    self.assertGreaterEqual(gap, 1)


class TestRawOntLodBalancedProfile(unittest.TestCase):
    """Test the raw ONT LOD-balanced exact-screening build controls."""

    def test_lod_balanced_profile_defaults(self) -> None:
        """The LOD-balanced profile should favour k=77/k=101 exact evidence."""
        self.assertEqual(
            resolve_candidate_k_order(
                marker_profile="raw_ont_lod_balanced",
                candidate_k_order="auto",
            ),
            "input",
        )
        self.assertEqual(
            resolve_min_cross_k_marker_distance(
                marker_profile="raw_ont_lod_balanced",
                min_cross_k_marker_distance=None,
            ),
            150,
        )
        self.assertEqual(
            resolve_max_per_genome_bin_by_k(
                marker_profile="raw_ont_lod_balanced",
                max_per_genome_bin_by_k="",
            ),
            {51: 4, 77: 6, 101: 4, 151: 1},
        )

    def test_lod_balanced_profile_explicit_quota_override(self) -> None:
        """Explicit quotas should override LOD-balanced defaults."""
        self.assertEqual(
            resolve_max_per_genome_bin_by_k(
                marker_profile="raw_ont_lod_balanced",
                max_per_genome_bin_by_k="51:2,77:2,101:2,151:2",
            ),
            {51: 2, 77: 2, 101: 2, 151: 2},
        )

    def test_lod_balanced_candidate_universe_retains_mid_k_values(self) -> None:
        """LOD-balanced quotas should retain comparatively dense mid-k evidence."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)
            fasta_path = tmp_path / "genome.fna"
            sqlite_path = tmp_path / "candidates.sqlite"
            fasta_path.write_text(
                ">contig1\n" + "ACGT" * 400 + "\n",
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
                genome_bin_size=100,
                max_per_genome_bin=1,
                max_per_genome_bin_by_k={5: 4, 7: 6, 9: 4, 11: 1},
                min_cross_k_marker_distance=1,
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
            self.assertIn(7, counts)
            self.assertIn(9, counts)
            self.assertGreater(counts[7], 0)
            self.assertGreater(counts[9], 0)
