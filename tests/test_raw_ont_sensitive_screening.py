"""Tests for raw-ONT sensitive screening and calling presets."""

import logging
import unittest
from argparse import Namespace
from unittest.mock import patch

from kmersutra.build_panel import DiagnosticKmer
from kmersutra.cli.screen_reads_for_clade_kmers import (
    resolve_screening_preset,
    select_fuzzy_rescue_species,
)
from kmersutra.fasta import SequenceRecord
from kmersutra.screen_reads import (
    build_one_mismatch_seed_index,
    deduplicate_hits,
    filter_panel_index_by_species,
    screen_records_for_species_kmers,
    screen_sequence_for_kmers,
)
from kmersutra.summarise_hits import (
    summarise_sample_species_evidence,
    summarise_species_hits,
)
from kmersutra.thresholds import apply_species_call_preset, call_species_presence


class TestRawOntSensitiveScreening(unittest.TestCase):
    """Tests for sensitivity-improving raw-read behaviour."""

    def test_screening_preset_sets_one_mismatch_long_k_defaults(self) -> None:
        """raw_ont_sensitive should enable one-mismatch long-k screening."""
        args = Namespace(
            screen_preset="raw_ont_sensitive",
            max_mismatches=None,
            fuzzy_min_k=None,
        )
        resolve_screening_preset(args=args, logger=logging.getLogger(__name__))

        self.assertEqual(args.max_mismatches, 1)
        self.assertEqual(args.fuzzy_min_k, 101)

    def test_screening_preset_preserves_explicit_user_values(self) -> None:
        """Explicit mismatch settings should override preset defaults."""
        args = Namespace(
            screen_preset="raw_ont_sensitive",
            max_mismatches=0,
            fuzzy_min_k=151,
        )
        resolve_screening_preset(args=args, logger=logging.getLogger(__name__))

        self.assertEqual(args.max_mismatches, 0)
        self.assertEqual(args.fuzzy_min_k, 151)

    def test_raw_ont_sensitive_call_preset_allows_fuzzy_long_k_evidence(self) -> None:
        """Fuzzy long-k evidence should be able to support a reportable call."""
        settings = apply_species_call_preset(preset_name="raw_ont_sensitive")
        self.assertEqual(settings["min_exact_hits"], 0)
        self.assertEqual(settings["min_total_hits"], 8)
        self.assertEqual(settings["min_best_k"], 77)

        evidence = [
            {
                "sample_id": "s1",
                "species_name": "Plasmodium example",
                "clade": "Plasmodium",
                "n_hits": 9,
                "n_unique_kmers": 9,
                "n_positive_sequences": 3,
                "n_k_values_positive": 1,
                "best_k": 101,
                "n_exact_hits": 0,
                "n_fuzzy_hits": 9,
            }
        ]
        calls = call_species_presence(evidence_records=evidence, **settings)
        self.assertEqual(calls[0]["call"], "present_high_confidence")

    def test_raw_ont_sensitive_call_preset_rejects_sparse_fuzzy_evidence(self) -> None:
        """Sparse fuzzy hits should remain below threshold."""
        settings = apply_species_call_preset(preset_name="raw_ont_sensitive")
        evidence = [
            {
                "sample_id": "s1",
                "species_name": "Plasmodium example",
                "clade": "Plasmodium",
                "n_hits": 3,
                "n_unique_kmers": 3,
                "n_positive_sequences": 2,
                "n_k_values_positive": 1,
                "best_k": 101,
                "n_exact_hits": 0,
                "n_fuzzy_hits": 3,
            }
        ]
        calls = call_species_presence(evidence_records=evidence, **settings)
        self.assertEqual(calls[0]["call"], "observed_below_threshold")

    def test_one_mismatch_long_k_hit_contributes_to_species_evidence(self) -> None:
        """One-mismatch k=101 hits should be counted as fuzzy evidence."""
        panel_kmer = "A" * 101
        query_kmer = "A" * 100 + "C"
        diagnostic = DiagnosticKmer(
            kmer=panel_kmer,
            k=101,
            panel_type="species_unique",
            species_name="Plasmodium example",
            clade="Plasmodium",
            source_genomes="g1",
            source_contigs="c1",
            example_position=0,
        )
        record = SequenceRecord(
            identifier="read1",
            description="read1",
            sequence=query_kmer,
        )
        hits = screen_sequence_for_kmers(
            sequence_record=record,
            panel_index={101: {panel_kmer: [diagnostic]}},
            sample_id="s1",
            sequence_type="read",
            max_mismatches=1,
            fuzzy_min_k=101,
        )
        hit_summary = summarise_species_hits(hits=hits)
        evidence = summarise_sample_species_evidence(species_summary=hit_summary)

        self.assertEqual(len(hits), 1)
        self.assertEqual(hits[0].mismatches, 1)
        self.assertEqual(evidence[0]["n_exact_hits"], 0)
        self.assertEqual(evidence[0]["n_fuzzy_hits"], 1)
        self.assertEqual(evidence[0]["best_k"], 101)


class TestRawOntSensitiveSeedAcceleration(unittest.TestCase):
    """Tests for accelerated one-mismatch raw-read screening."""

    def test_one_mismatch_seed_index_finds_long_k_hit(self) -> None:
        """Seeded one-mismatch screening should preserve fuzzy sensitivity."""
        panel_kmer = "A" * 101
        query_kmer = "A" * 100 + "C"
        diagnostic = DiagnosticKmer(
            kmer=panel_kmer,
            k=101,
            panel_type="species_unique",
            species_name="Plasmodium example",
            clade="Plasmodium",
            source_genomes="g1",
            source_contigs="c1",
            example_position=0,
        )
        panel_index = {101: {panel_kmer: [diagnostic]}}
        seed_index = build_one_mismatch_seed_index(
            panel_index=panel_index,
            fuzzy_min_k=101,
        )
        record = SequenceRecord(
            identifier="read1",
            description="read1",
            sequence=query_kmer,
        )

        hits = screen_sequence_for_kmers(
            sequence_record=record,
            panel_index=panel_index,
            sample_id="s1",
            sequence_type="read",
            max_mismatches=1,
            fuzzy_min_k=101,
            one_mismatch_seed_indices=seed_index,
        )

        self.assertEqual(len(hits), 1)
        self.assertEqual(hits[0].mismatches, 1)
        self.assertEqual(hits[0].species_name, "Plasmodium example")

    def test_one_mismatch_seed_index_rejects_two_mismatch_candidate(self) -> None:
        """Seeded one-mismatch screening should still verify Hamming distance."""
        panel_kmer = "A" * 101
        query_kmer = "A" * 99 + "CC"
        diagnostic = DiagnosticKmer(
            kmer=panel_kmer,
            k=101,
            panel_type="species_unique",
            species_name="Plasmodium example",
            clade="Plasmodium",
            source_genomes="g1",
            source_contigs="c1",
            example_position=0,
        )
        panel_index = {101: {panel_kmer: [diagnostic]}}
        seed_index = build_one_mismatch_seed_index(
            panel_index=panel_index,
            fuzzy_min_k=101,
        )
        record = SequenceRecord(
            identifier="read1",
            description="read1",
            sequence=query_kmer,
        )

        hits = screen_sequence_for_kmers(
            sequence_record=record,
            panel_index=panel_index,
            sample_id="s1",
            sequence_type="read",
            max_mismatches=1,
            fuzzy_min_k=101,
            one_mismatch_seed_indices=seed_index,
        )

        self.assertEqual(hits, [])

    def test_record_screening_uses_seed_path_for_one_mismatch(self) -> None:
        """Record-level one-mismatch screening should avoid neighbour generation."""
        panel_kmer = "A" * 101
        query_kmer = "A" * 100 + "C"
        diagnostic = DiagnosticKmer(
            kmer=panel_kmer,
            k=101,
            panel_type="species_unique",
            species_name="Plasmodium example",
            clade="Plasmodium",
            source_genomes="g1",
            source_contigs="c1",
            example_position=0,
        )
        records = [
            SequenceRecord(
                identifier="read1",
                description="read1",
                sequence=query_kmer,
            )
        ]

        with patch("kmersutra.screen_reads.iter_mismatch_neighbourhood") as mocked:
            mocked.side_effect = AssertionError(
                "Neighbour-generation fallback should not be used for max_mismatches=1"
            )
            hits = screen_records_for_species_kmers(
                records=records,
                panel_index={101: {panel_kmer: [diagnostic]}},
                sample_id="s1",
                sequence_type="read",
                max_mismatches=1,
                fuzzy_min_k=101,
                threads=1,
                chunk_size=1,
            )

        self.assertEqual(len(hits), 1)
        self.assertEqual(hits[0].mismatches, 1)

class TestRawOntSensitiveTwoStageRescue(unittest.TestCase):
    """Tests for exact-first fuzzy rescue candidate selection."""

    def test_select_fuzzy_rescue_species_limits_to_exact_evidence_candidates(self) -> None:
        """Fuzzy rescue should select supported candidates and ignore zero rows."""
        evidence = [
            {
                "species_name": "Alpha",
                "n_hits": 5,
                "n_unique_kmers": 5,
                "n_positive_sequences": 2,
                "best_k": 77,
            },
            {
                "species_name": "Beta",
                "n_hits": 0,
                "n_unique_kmers": 0,
                "n_positive_sequences": 0,
                "best_k": 0,
            },
            {
                "species_name": "Gamma",
                "n_hits": 9,
                "n_unique_kmers": 9,
                "n_positive_sequences": 3,
                "best_k": 51,
            },
        ]
        selected = select_fuzzy_rescue_species(
            evidence_records=evidence,
            max_species=1,
            min_unique_kmers=1,
            min_positive_sequences=1,
            logger=logging.getLogger(__name__),
        )
        self.assertEqual(selected, ["Gamma"])

    def test_filter_panel_index_by_species_keeps_requested_long_k_only(self) -> None:
        """Rescue panels should include only selected species and long k values."""
        alpha_short = DiagnosticKmer(
            kmer="A" * 77,
            k=77,
            panel_type="species_unique",
            species_name="Alpha",
            clade="Demo",
            source_genomes="g1",
            source_contigs="c1",
            example_position=0,
        )
        alpha_long = DiagnosticKmer(
            kmer="A" * 101,
            k=101,
            panel_type="species_unique",
            species_name="Alpha",
            clade="Demo",
            source_genomes="g1",
            source_contigs="c1",
            example_position=0,
        )
        beta_long = DiagnosticKmer(
            kmer="C" * 101,
            k=101,
            panel_type="species_unique",
            species_name="Beta",
            clade="Demo",
            source_genomes="g2",
            source_contigs="c1",
            example_position=0,
        )
        panel = {
            77: {"A" * 77: [alpha_short]},
            101: {"A" * 101: [alpha_long], "C" * 101: [beta_long]},
        }
        filtered = filter_panel_index_by_species(
            panel_index=panel,
            species_names={"Alpha"},
            min_k=101,
        )
        self.assertNotIn(77, filtered)
        self.assertIn(101, filtered)
        self.assertIn("A" * 101, filtered[101])
        self.assertNotIn("C" * 101, filtered[101])

    def test_deduplicate_hits_removes_repeated_exact_rescue_records(self) -> None:
        """Combined exact and rescue hits should not duplicate identical records."""
        diagnostic = DiagnosticKmer(
            kmer="A" * 101,
            k=101,
            panel_type="species_unique",
            species_name="Alpha",
            clade="Demo",
            source_genomes="g1",
            source_contigs="c1",
            example_position=0,
        )
        record = SequenceRecord(
            identifier="read1",
            description="read1",
            sequence="A" * 101,
        )
        hits = screen_sequence_for_kmers(
            sequence_record=record,
            panel_index={101: {"A" * 101: [diagnostic]}},
            sample_id="s1",
            sequence_type="read",
        )
        combined = deduplicate_hits(hits=[*hits, *hits])
        self.assertEqual(len(combined), 1)


if __name__ == "__main__":
    unittest.main()
