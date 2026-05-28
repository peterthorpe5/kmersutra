"""Tests for KmerSutra-ML feature extraction."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from kmersutra.features import (
    FEATURE_FIELDNAMES,
    hit_from_record,
    infer_numeric_feature_columns,
    load_hits_as_features,
    summarise_sequence_features,
    write_sequence_features,
)
from kmersutra.io import read_tsv, write_tsv
from kmersutra.screen_reads import KmerHit


class TestFeatures(unittest.TestCase):
    """Test sequence-level feature extraction from diagnostic hit rows."""

    def test_hit_from_record_converts_types(self):
        """Convert a read-level TSV row into a typed KmerHit object."""
        record = {
            "sample_id": "s1",
            "sequence_id": "read1",
            "sequence_type": "read",
            "k": "71",
            "query_position": "5",
            "matched_kmer": "AAAA",
            "query_kmer": "AAAT",
            "mismatches": "1",
            "panel_type": "species_unique",
            "species_name": "Species A",
            "clade": "Clade X",
        }
        hit = hit_from_record(record)
        self.assertEqual(hit.k, 71)
        self.assertEqual(hit.query_position, 5)
        self.assertEqual(hit.mismatches, 1)
        self.assertEqual(hit.species_name, "Species A")

    def test_summarise_sequence_features_best_species(self):
        """Summarise one read and identify the strongest species support."""
        hits = [
            KmerHit("s1", "read1", "read", 71, 0, "AAA", "AAA", 0, "species_unique", "Species A", "Clade X"),
            KmerHit("s1", "read1", "read", 101, 10, "CCC", "CCC", 0, "species_unique", "Species A", "Clade X"),
            KmerHit("s1", "read1", "read", 71, 20, "GGG", "GGA", 1, "species_unique", "Species B", "Clade X"),
        ]
        features = summarise_sequence_features(hits=hits)
        self.assertEqual(len(features), 1)
        row = features[0]
        self.assertEqual(row["best_species"], "Species A")
        self.assertEqual(row["best_species_hits"], 2)
        self.assertEqual(row["second_species"], "Species B")
        self.assertEqual(row["n_fuzzy_hits"], 1)
        self.assertAlmostEqual(float(row["conflict_ratio"]), 1 / 3, places=6)

    def test_summarise_sequence_features_unresolved_clade(self):
        """Mark clade-only evidence as unresolved at species level."""
        hits = [
            KmerHit("s1", "read1", "read", 51, 0, "AAA", "AAA", 0, "clade_core", "", "Plasmodium"),
        ]
        features = summarise_sequence_features(hits=hits)
        self.assertEqual(features[0]["best_clade"], "Plasmodium")
        self.assertEqual(features[0]["best_species"], "")
        self.assertEqual(features[0]["unresolved_clade_signal"], 1)

    def test_load_hits_as_features_from_gzip_tsv(self):
        """Load gzipped read-level hits and write feature output."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            hits_tsv = tmp_path / "hits.tsv.gz"
            out_tsv = tmp_path / "features.tsv"
            write_tsv(
                records=[
                    {
                        "sample_id": "s1",
                        "sequence_id": "read1",
                        "sequence_type": "read",
                        "k": 71,
                        "query_position": 0,
                        "matched_kmer": "AAA",
                        "query_kmer": "AAA",
                        "mismatches": 0,
                        "panel_type": "species_unique",
                        "species_name": "Species A",
                        "clade": "Clade X",
                    }
                ],
                output_path=hits_tsv,
                fieldnames=[
                    "sample_id",
                    "sequence_id",
                    "sequence_type",
                    "k",
                    "query_position",
                    "matched_kmer",
                    "query_kmer",
                    "mismatches",
                    "panel_type",
                    "species_name",
                    "clade",
                ],
            )
            features = load_hits_as_features(hits_tsv=hits_tsv)
            write_sequence_features(features=features, output_path=out_tsv)
            loaded = read_tsv(input_path=out_tsv)
            self.assertEqual(loaded[0]["best_species"], "Species A")
            self.assertEqual(set(loaded[0]).issuperset(FEATURE_FIELDNAMES), True)

    def test_infer_numeric_feature_columns_excludes_labels(self):
        """Infer numeric feature columns while excluding metadata and labels."""
        records = [
            {"sample_id": "s1", "true_species": "A", "n_total_hits": "3", "best_species": "A"},
            {"sample_id": "s2", "true_species": "B", "n_total_hits": "5", "best_species": "B"},
        ]
        columns = infer_numeric_feature_columns(records=records, label_column="true_species")
        self.assertEqual(columns, ["n_total_hits"])


if __name__ == "__main__":
    unittest.main()
