"""Tests for raw-ONT sensitive screening and calling presets."""

import logging
import unittest
from argparse import Namespace

from kmersutra.build_panel import DiagnosticKmer
from kmersutra.cli.screen_reads_for_clade_kmers import resolve_screening_preset
from kmersutra.fasta import SequenceRecord
from kmersutra.screen_reads import screen_sequence_for_kmers
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


if __name__ == "__main__":
    unittest.main()
