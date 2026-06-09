#!/usr/bin/env python3
"""Unit tests for external Zymo validation helpers."""

from __future__ import annotations

import unittest

from apply_internal_model_to_zymo import validate_feature_compatibility


class FakeModel:
    """Small fake model for compatibility tests."""

    feature_columns = ["n_hits", "n_unique_kmers"]


class TestExternalZymoHelpers(unittest.TestCase):
    """Tests for external-validation helper behaviour."""

    def test_validate_feature_compatibility_passes(self) -> None:
        """Compatible feature records should pass."""
        validate_feature_compatibility(
            feature_records=[{"n_hits": 1, "n_unique_kmers": 2}],
            model=FakeModel(),
        )

    def test_validate_feature_compatibility_fails_on_missing(self) -> None:
        """Missing model features should raise ValueError."""
        with self.assertRaises(ValueError):
            validate_feature_compatibility(
                feature_records=[{"n_hits": 1}],
                model=FakeModel(),
            )

    def test_validate_feature_compatibility_fails_on_empty(self) -> None:
        """Empty external feature tables should raise ValueError."""
        with self.assertRaises(ValueError):
            validate_feature_compatibility(feature_records=[], model=FakeModel())


if __name__ == "__main__":
    unittest.main()
