"""Tests for assembly-aware candidate-sampling helpers."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from kmersutra.assembly_stats import (
    build_contig_offsets,
    calculate_assembly_stats,
    calculate_n50,
    choose_assembly_aware_bin_plan,
    collect_fasta_contig_lengths,
)


class TestAssemblyStats(unittest.TestCase):
    """Test assembly statistics and adaptive bin planning."""

    def test_calculate_n50_uses_length_weighted_threshold(self) -> None:
        """N50 should be calculated from sorted cumulative contig lengths."""
        self.assertEqual(calculate_n50(lengths=[90, 80, 30]), 80)

    def test_small_viral_like_assembly_uses_larger_bin_size(self) -> None:
        """Small assemblies should not be over-binned with tiny requested bins."""
        stats = calculate_assembly_stats(
            contig_lengths={"viral_contig": 30000},
            assembly_id="viral",
            k_values=[51, 77, 101],
        )
        plan = choose_assembly_aware_bin_plan(
            stats=stats,
            requested_bin_size=1000,
            small_assembly_min_bin_size=10000,
            small_assembly_target_bins=25,
        )
        self.assertTrue(plan.is_small_assembly)
        self.assertEqual(plan.effective_bin_size, 10000)
        self.assertEqual(plan.estimated_global_bins, 3)

    def test_fragmented_assembly_uses_global_cumulative_bins(self) -> None:
        """Fragmented assemblies should use an effective global bin plan."""
        contigs = {f"contig_{index}": 1000 for index in range(600)}
        stats = calculate_assembly_stats(
            contig_lengths=contigs,
            assembly_id="fragmented",
            k_values=[51],
        )
        plan = choose_assembly_aware_bin_plan(
            stats=stats,
            requested_bin_size=1000,
            fragmented_contig_count=500,
            fragmented_max_global_bins=100,
        )
        self.assertTrue(plan.is_fragmented_assembly)
        self.assertGreaterEqual(plan.effective_bin_size, 6000)
        self.assertLessEqual(plan.estimated_global_bins, 100)

    def test_contig_offsets_skip_too_short_contigs(self) -> None:
        """Contigs shorter than the minimum k should not receive offsets."""
        offsets = build_contig_offsets(
            contig_lengths={"short": 20, "long": 200, "long2": 100},
            min_contig_length=51,
        )
        self.assertNotIn("short", offsets)
        self.assertEqual(offsets["long"], 0)
        self.assertEqual(offsets["long2"], 200)

    def test_collect_fasta_contig_lengths_reads_all_records(self) -> None:
        """FASTA length collection should report all contigs."""
        with tempfile.TemporaryDirectory() as tmpdir:
            fasta = Path(tmpdir) / "test.fna"
            fasta.write_text(">a\nAAAA\n>b\nCCCCCC\n", encoding="utf-8")
            lengths = collect_fasta_contig_lengths(fasta_path=fasta)
        self.assertEqual(lengths, {"a": 4, "b": 6})


if __name__ == "__main__":
    unittest.main()
