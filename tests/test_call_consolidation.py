"""Tests for post-threshold KmerSutra species-call consolidation."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from kmersutra.call_consolidation import (
    BACKGROUND_CANDIDATE_CALL,
    NEIGHBOUR_LINEAGE_CALL,
    consolidate_species_calls,
    extract_genus,
    merge_background_taxa,
    normalise_taxon_name,
)


class TestCallConsolidation(unittest.TestCase):
    """Test report-layer consolidation for species calls."""

    def test_normalise_taxon_name_handles_underscores_and_whitespace(self) -> None:
        """Taxon-name normalisation should support robust matching."""
        self.assertEqual(
            normalise_taxon_name(value="  Hammondia_hammondi  "),
            "hammondia hammondi",
        )

    def test_extract_genus_returns_first_token(self) -> None:
        """Genus extraction should return the normalised first word."""
        self.assertEqual(extract_genus(species_name="Plasmodium vivax"), "plasmodium")
        self.assertEqual(extract_genus(species_name=""), "")

    def test_merge_background_taxa_reads_direct_and_file_values(self) -> None:
        """Background taxa should merge direct CLI labels and text-file labels."""
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "background_taxa.txt"
            path.write_text("# comment\nHammondia hammondi\nToxoplasma_gondii\n", encoding="utf-8")
            observed = merge_background_taxa(
                background_candidate_taxa=["Neospora caninum"],
                background_candidate_file=path,
            )
        self.assertEqual(
            observed,
            {"hammondia hammondi", "toxoplasma gondii", "neospora caninum"},
        )

    def test_background_candidate_taxon_is_relabelled_not_dropped(self) -> None:
        """A plausible background taxon should be retained in a separate layer."""
        rows = [
            {
                "sample_id": "s1",
                "species_name": "Hammondia hammondi",
                "call": "present_high_confidence",
                "n_unique_kmers": 28,
                "n_positive_sequences": 73,
                "best_k": 151,
            }
        ]
        observed = consolidate_species_calls(
            species_calls=rows,
            background_candidate_taxa=["Hammondia hammondi"],
        )
        self.assertEqual(observed[0]["pre_consolidation_call"], "present_high_confidence")
        self.assertEqual(observed[0]["call"], BACKGROUND_CANDIDATE_CALL)
        self.assertEqual(observed[0]["report_layer"], "background_candidate")
        self.assertEqual(observed[0]["is_background_candidate"], 1)

    def test_dominated_same_genus_species_is_demoted_to_neighbour_lineage(self) -> None:
        """A weak same-genus species should be demoted when dominated by a primary."""
        rows = [
            {
                "sample_id": "s1",
                "species_name": "Plasmodium vivax",
                "call": "present_high_confidence",
                "n_unique_kmers": 200,
                "n_positive_sequences": 30,
                "best_k": 151,
            },
            {
                "sample_id": "s1",
                "species_name": "Plasmodium simium",
                "call": "present_high_confidence",
                "n_unique_kmers": 50,
                "n_positive_sequences": 8,
                "best_k": 101,
            },
        ]
        observed = consolidate_species_calls(
            species_calls=rows,
            dominant_species_min_margin=25,
            dominant_species_min_ratio=2.0,
        )
        calls = {row["species_name"]: row for row in observed}
        self.assertEqual(calls["Plasmodium vivax"]["call"], "present_high_confidence")
        self.assertEqual(calls["Plasmodium simium"]["call"], NEIGHBOUR_LINEAGE_CALL)
        self.assertEqual(
            calls["Plasmodium simium"]["consolidation_reason"],
            "dominated_same_genus_neighbour",
        )

    def test_co_dominant_same_genus_species_is_not_demoted(self) -> None:
        """Possible true mixed species should survive if support is co-dominant."""
        rows = [
            {
                "sample_id": "s1",
                "species_name": "Plasmodium alpha",
                "call": "present_in_mixed_sample",
                "n_unique_kmers": 100,
                "n_positive_sequences": 20,
                "best_k": 151,
            },
            {
                "sample_id": "s1",
                "species_name": "Plasmodium beta",
                "call": "present_in_mixed_sample",
                "n_unique_kmers": 80,
                "n_positive_sequences": 18,
                "best_k": 151,
            },
        ]
        observed = consolidate_species_calls(
            species_calls=rows,
            dominant_species_min_margin=25,
            dominant_species_min_ratio=2.0,
        )
        self.assertEqual(
            {row["call"] for row in observed},
            {"present_in_mixed_sample"},
        )

    def test_consolidation_can_be_disabled(self) -> None:
        """Same-genus demotion should be optional for reproducible old runs."""
        rows = [
            {
                "sample_id": "s1",
                "species_name": "Plasmodium vivax",
                "call": "present_high_confidence",
                "n_unique_kmers": 200,
            },
            {
                "sample_id": "s1",
                "species_name": "Plasmodium simium",
                "call": "present_high_confidence",
                "n_unique_kmers": 50,
            },
        ]
        observed = consolidate_species_calls(
            species_calls=rows,
            demote_same_genus_neighbours=False,
        )
        self.assertEqual(
            {row["call"] for row in observed},
            {"present_high_confidence"},
        )


if __name__ == "__main__":
    unittest.main()
