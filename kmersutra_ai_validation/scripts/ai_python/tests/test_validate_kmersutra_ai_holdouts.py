#!/usr/bin/env python3
"""Unit tests for internal holdout-validation helper functions."""

from __future__ import annotations

import unittest

from validate_kmersutra_ai_holdouts import (
    build_validation_splits,
    make_sample_hash_split,
    should_skip_split,
    summarise_metrics_by_split,
)


class TestInternalHoldoutHelpers(unittest.TestCase):
    """Tests for validation split and summary helpers."""

    def setUp(self) -> None:
        """Create small artificial records."""
        self.records = [
            {
                "sample_id": "s1",
                "benchmark_family": "single_genome",
                "panel": "panel1",
                "species_name": "Plasmodium vivax",
                "ml_report_label": "expected_target",
            },
            {
                "sample_id": "s2",
                "benchmark_family": "two_genome",
                "panel": "panel2",
                "species_name": "Plasmodium falciparum",
                "ml_report_label": "expected_target",
            },
            {
                "sample_id": "s3",
                "benchmark_family": "shuffled_negative",
                "panel": "shuffled",
                "species_name": "Toxoplasma gondii",
                "ml_report_label": "not_detected",
            },
            {
                "sample_id": "s4",
                "benchmark_family": "single_genome",
                "panel": "panel1",
                "species_name": "Plasmodium simium",
                "ml_report_label": "observed_below_threshold",
            },
        ]

    def test_make_sample_hash_split_keeps_all_records(self) -> None:
        """The sample-hash split should not drop records."""
        train, test = make_sample_hash_split(records=self.records, test_fraction=0.5)
        self.assertEqual(len(train) + len(test), len(self.records))

    def test_build_validation_splits_includes_expected_species_splits(self) -> None:
        """Expected-target species splits should be included."""
        splits = build_validation_splits(
            records=self.records,
            label_column="ml_report_label",
            sample_test_fraction=0.5,
        )
        names = {split[0] for split in splits}
        self.assertIn("sample_group_hash_20pct", names)
        self.assertIn(
            "leave_one_expected_target_species_Plasmodium_vivax",
            names,
        )

    def test_should_skip_split(self) -> None:
        """Split validation should report empty train/test cases."""
        self.assertEqual(
            should_skip_split(
                train_records=[],
                test_records=self.records,
                label_column="ml_report_label",
            ),
            "no_training_records",
        )
        self.assertEqual(
            should_skip_split(
                train_records=self.records,
                test_records=[],
                label_column="ml_report_label",
            ),
            "no_test_records",
        )

    def test_summarise_metrics_by_split(self) -> None:
        """Metric summary should extract overall and expected-target rows."""
        rows = [
            {
                "validation_name": "split1",
                "holdout_column": "panel",
                "holdout_value": "panel1",
                "n_train": 10,
                "n_test": 3,
                "label": "expected_target",
                "n": 2,
                "precision": 1.0,
                "recall": 0.5,
            },
            {
                "validation_name": "split1",
                "holdout_column": "panel",
                "holdout_value": "panel1",
                "n_train": 10,
                "n_test": 3,
                "label": "overall",
                "f1": 0.8,
            },
        ]
        summary = summarise_metrics_by_split(rows)
        self.assertEqual(len(summary), 1)
        self.assertEqual(summary[0]["overall_accuracy"], 0.8)
        self.assertEqual(summary[0]["expected_target_recall"], 0.5)


if __name__ == "__main__":
    unittest.main()
