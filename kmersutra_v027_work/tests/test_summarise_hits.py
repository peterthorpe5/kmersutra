"""Tests for KmerSutra hit summarisation."""

import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from kmersutra.io import write_tsv
from kmersutra.screen_reads import KmerHit
from kmersutra.summarise_hits import (
    complete_sample_species_evidence,
    load_panel_species_metadata,
    summarise_sample_species_evidence,
    summarise_species_hits,
)


class TestSummariseHits(unittest.TestCase):
    """Tests for summarising diagnostic k-mer hits."""

    def test_summarise_species_hits_groups_by_species_and_k(self) -> None:
        """Hit summary should aggregate hits by species, panel type and k."""
        hits = [
            KmerHit("s1", "r1", "read", 5, 0, "AAAAA", "AAAAA", 0, "species_unique", "Alpha", "Demo"),
            KmerHit("s1", "r2", "read", 5, 0, "AAAAC", "AAAAC", 0, "species_unique", "Alpha", "Demo"),
        ]
        summary = summarise_species_hits(hits=hits)
        self.assertEqual(summary[0]["n_hits"], 2)
        self.assertEqual(summary[0]["n_positive_sequences"], 2)

    def test_summarise_sample_species_evidence_collapses_k_values(self) -> None:
        """Species evidence should collapse k-specific rows into one record."""
        summary = [
            {
                "sample_id": "s1", "panel_type": "species_unique", "label": "Alpha",
                "clade": "Demo", "k": 5, "n_hits": 2, "n_unique_kmers": 2,
                "n_positive_sequences": 1, "n_exact_hits": 2, "n_fuzzy_hits": 0,
            },
            {
                "sample_id": "s1", "panel_type": "species_unique", "label": "Alpha",
                "clade": "Demo", "k": 7, "n_hits": 3, "n_unique_kmers": 3,
                "n_positive_sequences": 2, "n_exact_hits": 3, "n_fuzzy_hits": 0,
            },
        ]
        evidence = summarise_sample_species_evidence(species_summary=summary)
        self.assertEqual(evidence[0]["n_hits"], 5)
        self.assertEqual(evidence[0]["n_k_values_positive"], 2)
        self.assertEqual(evidence[0]["best_k"], 7)

    def test_load_panel_species_metadata(self) -> None:
        """Panel metadata loader should recover expected species labels."""
        with TemporaryDirectory() as tmpdir:
            panel_path = Path(tmpdir) / "panel.tsv"
            write_tsv(
                records=[
                    {
                        "kmer": "AAAAA", "k": 5, "panel_type": "species_unique",
                        "species_name": "Alpha", "clade": "Demo", "source_genomes": "g1",
                        "source_contigs": "c1", "example_position": 0,
                    },
                    {
                        "kmer": "CCCCC", "k": 5, "panel_type": "clade_core",
                        "species_name": "", "clade": "Demo", "source_genomes": "g1;g2",
                        "source_contigs": "c1", "example_position": 0,
                    },
                ],
                output_path=panel_path,
            )
            species = load_panel_species_metadata(panel_path=panel_path)
            self.assertEqual(species, [{"species_name": "Alpha", "clade": "Demo"}])

    def test_complete_sample_species_evidence_adds_zero_rows(self) -> None:
        """Completed evidence should include not-observed expected species."""
        completed = complete_sample_species_evidence(
            evidence_records=[
                {
                    "sample_id": "s1", "species_name": "Alpha", "clade": "Demo",
                    "n_hits": 5, "n_unique_kmers": 5, "n_positive_sequences": 2,
                    "n_k_values_positive": 1, "best_k": 71, "n_exact_hits": 5,
                    "n_fuzzy_hits": 0,
                }
            ],
            expected_species=[
                {"species_name": "Alpha", "clade": "Demo"},
                {"species_name": "Beta", "clade": "Demo"},
            ],
            sample_id="s1",
        )
        beta = [row for row in completed if row["species_name"] == "Beta"][0]
        self.assertEqual(beta["n_unique_kmers"], 0)
        self.assertEqual(beta["best_k"], 0)


class TestTaxonomicEvidenceSummaries(unittest.TestCase):
    """Tests for v0.15 taxonomic evidence retention."""

    def test_summarise_taxonomic_hits_keeps_genus_evidence(self) -> None:
        """Genus-level hits should be summarised without forcing species."""
        from kmersutra.screen_reads import KmerHit
        from kmersutra.summarise_hits import summarise_taxonomic_hits

        hits = [
            KmerHit(
                sample_id="s1",
                sequence_id="read1",
                sequence_type="read",
                k=101,
                query_position=0,
                matched_kmer="A" * 101,
                query_kmer="A" * 101,
                mismatches=0,
                panel_type="genus_core",
                species_name="",
                clade="Plasmodium",
                evidence_taxid="5820",
                evidence_name="Plasmodium",
                evidence_rank="genus",
            )
        ]
        summary = summarise_taxonomic_hits(hits=hits)
        self.assertEqual(len(summary), 1)
        self.assertEqual(summary[0]["evidence_rank"], "genus")
        self.assertEqual(summary[0]["evidence_name"], "Plasmodium")
        self.assertEqual(summary[0]["n_unique_kmers"], 1)

    def test_summarise_sample_taxonomic_evidence_collapses_k_values(self) -> None:
        """Taxonomic evidence should collapse across k values."""
        from kmersutra.summarise_hits import summarise_sample_taxonomic_evidence

        taxonomic_summary = [
            {
                "sample_id": "s1",
                "evidence_rank": "genus",
                "evidence_name": "Plasmodium",
                "evidence_taxid": "5820",
                "panel_type": "genus_core",
                "clade": "Plasmodium",
                "k": 77,
                "n_hits": 5,
                "n_unique_kmers": 4,
                "n_positive_sequences": 2,
                "n_exact_hits": 5,
                "n_fuzzy_hits": 0,
            },
            {
                "sample_id": "s1",
                "evidence_rank": "genus",
                "evidence_name": "Plasmodium",
                "evidence_taxid": "5820",
                "panel_type": "genus_core",
                "clade": "Plasmodium",
                "k": 101,
                "n_hits": 7,
                "n_unique_kmers": 6,
                "n_positive_sequences": 3,
                "n_exact_hits": 7,
                "n_fuzzy_hits": 0,
            },
        ]
        evidence = summarise_sample_taxonomic_evidence(
            taxonomic_summary=taxonomic_summary,
        )
        self.assertEqual(len(evidence), 1)
        self.assertEqual(evidence[0]["n_hits"], 12)
        self.assertEqual(evidence[0]["n_unique_kmers"], 10)
        self.assertEqual(evidence[0]["n_k_values_positive"], 2)
        self.assertEqual(evidence[0]["best_k"], 101)
