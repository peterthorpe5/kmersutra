"""Tests for summary CLI dispatch helpers."""

from __future__ import annotations

import unittest

from kmersutra.cli.summarise_spikein_run import should_use_comparable_mode


class TestSummaryCliDispatch(unittest.TestCase):
    """Test backwards-compatible summary command dispatch."""

    def test_out_root_uses_comparable_mode(self) -> None:
        """Comparable arguments should dispatch to comparable summary mode."""
        self.assertTrue(should_use_comparable_mode(argv=["--out_root", "run"]))

    def test_legacy_summary_tsv_does_not_use_comparable_mode(self) -> None:
        """Legacy summary arguments should keep the old report formatter."""
        self.assertFalse(
            should_use_comparable_mode(
                argv=["--summary_tsv", "summary.tsv", "--out_xlsx", "out.xlsx"]
            )
        )


if __name__ == "__main__":
    unittest.main()
