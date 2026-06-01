"""Tests for production-scale progress logging throttling."""

from __future__ import annotations

import unittest

from kmersutra.global_candidate_evidence import _should_log_progress_update


class TestProgressLogging(unittest.TestCase):
    """Test progress logging decisions for large database builds."""

    def test_does_not_log_before_interval(self) -> None:
        """Progress logging should wait until the configured interval."""
        result = _should_log_progress_update(
            attempted=999,
            last_logged_at=0,
            current_count=10,
            last_logged_count=0,
            progress_interval=1000,
        )
        self.assertFalse(result)

    def test_logs_when_count_changes_after_interval(self) -> None:
        """Progress logging should emit updates when the informative count changes."""
        result = _should_log_progress_update(
            attempted=1000,
            last_logged_at=0,
            current_count=10,
            last_logged_count=9,
            progress_interval=1000,
        )
        self.assertTrue(result)

    def test_suppresses_unchanged_count_until_heartbeat(self) -> None:
        """Unchanged candidate or hit counts should not be logged every interval."""
        result = _should_log_progress_update(
            attempted=5000,
            last_logged_at=0,
            current_count=10,
            last_logged_count=10,
            progress_interval=1000,
            heartbeat_multiplier=10,
        )
        self.assertFalse(result)

    def test_logs_unchanged_count_at_heartbeat(self) -> None:
        """A heartbeat should still be emitted when counts stay unchanged for a long scan."""
        result = _should_log_progress_update(
            attempted=10000,
            last_logged_at=0,
            current_count=10,
            last_logged_count=10,
            progress_interval=1000,
            heartbeat_multiplier=10,
        )
        self.assertTrue(result)

    def test_rejects_non_positive_interval(self) -> None:
        """A non-positive progress interval should fail clearly."""
        with self.assertRaises(ValueError):
            _should_log_progress_update(
                attempted=10,
                last_logged_at=0,
                current_count=1,
                last_logged_count=0,
                progress_interval=0,
            )


if __name__ == "__main__":
    unittest.main()
