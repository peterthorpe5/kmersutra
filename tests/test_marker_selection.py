"""Tests for genome-aware marker selection."""

from __future__ import annotations

import unittest

from kmersutra.build_panel import DiagnosticKmer
from kmersutra.marker_selection import (
    MarkerSelectionConfig,
    diagnostic_retention_key,
    first_semicolon_value,
    genome_bin_key,
    select_genome_spread_markers,
)


def make_marker(
    *,
    kmer: str,
    position: int,
    contig: str = "contig1",
    species_name: str = "Alpha target",
    evidence_taxid: str = "11",
    evidence_name: str = "Alpha target",
    evidence_rank: str = "species",
    k: int = 5,
) -> DiagnosticKmer:
    """Create a DiagnosticKmer for marker-selection tests.

    Parameters
    ----------
    kmer : str
        K-mer sequence.
    position : int
        Example source position.
    contig : str, optional
        Source contig ID.
    species_name : str, optional
        Species name.
    evidence_taxid : str, optional
        Evidence taxid.
    evidence_name : str, optional
        Evidence taxon name.
    evidence_rank : str, optional
        Evidence rank.
    k : int, optional
        K-mer length.

    Returns
    -------
    DiagnosticKmer
        Diagnostic k-mer test record.
    """
    panel_type = "species_unique" if evidence_rank == "species" else f"{evidence_rank}_core"
    return DiagnosticKmer(
        kmer=kmer,
        k=k,
        panel_type=panel_type,
        species_name=species_name if evidence_rank == "species" else "",
        clade="AlphaGenus",
        source_genomes="GCA_TEST_1",
        source_contigs=contig,
        example_position=position,
        evidence_taxid=evidence_taxid,
        evidence_name=evidence_name,
        evidence_rank=evidence_rank,
        lineage_taxids="1;10;11",
        source_taxids=evidence_taxid,
    )


class TestMarkerSelection(unittest.TestCase):
    """Tests for deterministic genome-spread marker selection."""

    def test_first_semicolon_value_returns_first_non_empty_value(self) -> None:
        """The helper should parse semicolon-separated metadata."""
        self.assertEqual(first_semicolon_value("alpha;beta"), "alpha")
        self.assertEqual(first_semicolon_value(";beta"), "beta")
        self.assertEqual(first_semicolon_value(""), "")

    def test_genome_bin_key_uses_position_bins(self) -> None:
        """Genome bins should be based on source genome, contig and position."""
        marker = make_marker(kmer="AAAAA", position=2500, contig="chr1")
        self.assertEqual(
            genome_bin_key(item=marker, genome_bin_size=1000),
            ("GCA_TEST_1", "chr1", 2),
        )

    def test_diagnostic_retention_key_groups_by_evidence_bucket(self) -> None:
        """Species and genus evidence should have distinct retention buckets."""
        species_marker = make_marker(kmer="AAAAA", position=0)
        genus_marker = make_marker(
            kmer="CCCCC",
            position=0,
            species_name="",
            evidence_taxid="10",
            evidence_name="AlphaGenus",
            evidence_rank="genus",
        )
        self.assertNotEqual(
            diagnostic_retention_key(species_marker),
            diagnostic_retention_key(genus_marker),
        )

    def test_config_rejects_invalid_values(self) -> None:
        """MarkerSelectionConfig should reject unsafe parameters."""
        with self.assertRaises(ValueError):
            MarkerSelectionConfig(strategy="unknown").validate()
        with self.assertRaises(ValueError):
            MarkerSelectionConfig(max_per_bucket=0).validate()
        with self.assertRaises(ValueError):
            MarkerSelectionConfig(genome_bin_size=0).validate()
        with self.assertRaises(ValueError):
            MarkerSelectionConfig(max_per_genome_bin=0).validate()

    def test_genome_spread_limits_dense_same_bin_markers(self) -> None:
        """Selection should not retain many adjacent markers from one bin."""
        markers = [
            make_marker(kmer=f"AAAA{i:02d}", position=i)
            for i in range(10)
        ]
        selected = list(
            select_genome_spread_markers(
                diagnostic_kmers=markers,
                config=MarkerSelectionConfig(
                    strategy="genome_spread",
                    max_per_bucket=10,
                    genome_bin_size=1000,
                    max_per_genome_bin=2,
                ),
            )
        )
        self.assertEqual(len(selected), 2)
        self.assertTrue(all(item.example_position < 10 for item in selected))

    def test_genome_spread_selects_across_multiple_bins(self) -> None:
        """Selection should retain markers from several genomic bins."""
        markers = []
        for bin_index in range(5):
            for offset in range(4):
                position = (bin_index * 1000) + offset
                markers.append(make_marker(kmer=f"A{bin_index}{offset}CCC", position=position))
        selected = list(
            select_genome_spread_markers(
                diagnostic_kmers=markers,
                config=MarkerSelectionConfig(
                    strategy="genome_spread",
                    max_per_bucket=5,
                    genome_bin_size=1000,
                    max_per_genome_bin=1,
                ),
            )
        )
        selected_bins = {
            genome_bin_key(item=item, genome_bin_size=1000)[2]
            for item in selected
        }
        self.assertEqual(len(selected), 5)
        self.assertEqual(selected_bins, {0, 1, 2, 3, 4})

    def test_genome_spread_keeps_buckets_separate(self) -> None:
        """Caps should be applied independently to each evidence bucket."""
        species_markers = [
            make_marker(kmer=f"AAAA{i}", position=i * 1000)
            for i in range(3)
        ]
        genus_markers = [
            make_marker(
                kmer=f"CCCC{i}",
                position=i * 1000,
                species_name="",
                evidence_taxid="10",
                evidence_name="AlphaGenus",
                evidence_rank="genus",
            )
            for i in range(3)
        ]
        selected = list(
            select_genome_spread_markers(
                diagnostic_kmers=species_markers + genus_markers,
                config=MarkerSelectionConfig(
                    strategy="genome_spread",
                    max_per_bucket=2,
                    genome_bin_size=1000,
                    max_per_genome_bin=1,
                ),
            )
        )
        species_selected = [item for item in selected if item.evidence_rank == "species"]
        genus_selected = [item for item in selected if item.evidence_rank == "genus"]
        self.assertEqual(len(species_selected), 2)
        self.assertEqual(len(genus_selected), 2)

    def test_genome_spread_is_deterministic(self) -> None:
        """Repeated selection over the same markers should be reproducible."""
        markers = [
            make_marker(kmer=f"AAAA{i:02d}", position=i * 100)
            for i in range(20)
        ]
        config = MarkerSelectionConfig(
            strategy="genome_spread",
            max_per_bucket=5,
            genome_bin_size=500,
            max_per_genome_bin=2,
        )
        first = [item.kmer for item in select_genome_spread_markers(diagnostic_kmers=markers, config=config)]
        second = [item.kmer for item in select_genome_spread_markers(diagnostic_kmers=markers, config=config)]
        self.assertEqual(first, second)


if __name__ == "__main__":
    unittest.main()

class TestMarkerSelectionScalability(unittest.TestCase):
    """Tests for scalable genome-spread selection behaviour."""

    def test_large_bucket_respects_bucket_and_bin_limits(self) -> None:
        """Large buckets should be thinned without violating caps."""
        markers = [
            make_marker(kmer=f"A{i:05d}", position=i * 10)
            for i in range(1000)
        ]
        config = MarkerSelectionConfig(
            strategy="genome_spread",
            max_per_bucket=50,
            genome_bin_size=100,
            max_per_genome_bin=2,
        )
        selected = list(
            select_genome_spread_markers(
                diagnostic_kmers=markers,
                config=config,
            )
        )
        bin_counts = {}
        for item in selected:
            key = genome_bin_key(item=item, genome_bin_size=100)
            bin_counts[key] = bin_counts.get(key, 0) + 1

        self.assertLessEqual(len(selected), 50)
        self.assertTrue(bin_counts)
        self.assertLessEqual(max(bin_counts.values()), 2)



class TestIndependentMultiKMarkerSelection(unittest.TestCase):
    """Tests for v0.31 independent multi-k marker sampling."""

    def test_shifted_bins_differ_between_k_values(self) -> None:
        """Different k values should use different bin phases by default."""
        from kmersutra.marker_selection import calculate_k_bin_offset

        offsets = [
            calculate_k_bin_offset(k=k, k_values=[51, 77, 101, 151], genome_bin_size=10000)
            for k in [51, 77, 101, 151]
        ]
        self.assertEqual(len(set(offsets)), 4)
        self.assertEqual(offsets[0], 0)

    def test_independent_mode_reduces_cross_k_same_region_selection(self) -> None:
        """Nested cross-k markers from one region should be de-correlated."""
        markers = []
        for k in [51, 77, 101]:
            for offset in [0, 5, 10, 15]:
                markers.append(
                    make_marker(
                        kmer=("A" * k)[:-len(str(offset))] + str(offset),
                        position=1000 + offset,
                        k=k,
                    )
                )
        selected = list(
            select_genome_spread_markers(
                diagnostic_kmers=markers,
                config=MarkerSelectionConfig(
                    strategy="independent_multik_genome_spread",
                    max_per_bucket=10,
                    genome_bin_size=10000,
                    max_per_genome_bin=10,
                    min_cross_k_marker_distance=500,
                ),
            )
        )
        selected_k_values = {item.k for item in selected}
        self.assertGreaterEqual(len(selected), 1)
        self.assertEqual(len(selected_k_values), 1)

    def test_independent_mode_allows_distant_cross_k_regions(self) -> None:
        """Different k values should be retained when they come from distant regions."""
        markers = [
            make_marker(kmer="A" * 51, position=1000, k=51),
            make_marker(kmer="C" * 77, position=10000, k=77),
            make_marker(kmer="G" * 101, position=20000, k=101),
        ]
        selected = list(
            select_genome_spread_markers(
                diagnostic_kmers=markers,
                config=MarkerSelectionConfig(
                    strategy="independent_multik_genome_spread",
                    max_per_bucket=10,
                    genome_bin_size=10000,
                    max_per_genome_bin=10,
                    min_cross_k_marker_distance=500,
                ),
            )
        )
        self.assertEqual({item.k for item in selected}, {51, 77, 101})

    def test_independent_mode_is_default_strategy(self) -> None:
        """The v0.31 default should use independent multi-k sampling."""
        self.assertEqual(
            MarkerSelectionConfig().strategy,
            "independent_multik_genome_spread",
        )
