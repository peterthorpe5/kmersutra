"""Tests for dependency-free RAM monitoring."""

from __future__ import annotations

import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from kmersutra.resource_monitor import (
    ResourceMonitor,
    get_current_rss_bytes,
    get_peak_rss_bytes,
)


class TestResourceMonitor(unittest.TestCase):
    """Tests for process RAM monitoring helpers."""

    def test_rss_helpers_return_non_negative_integers(self) -> None:
        """RSS helpers should return non-negative integer values."""
        self.assertIsInstance(get_current_rss_bytes(), int)
        self.assertIsInstance(get_peak_rss_bytes(), int)
        self.assertGreaterEqual(get_current_rss_bytes(), 0)
        self.assertGreaterEqual(get_peak_rss_bytes(), 0)

    def test_resource_monitor_writes_header_and_rows(self) -> None:
        """ResourceMonitor should write a parseable TSV file."""
        with TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "ram_usage.tsv"
            monitor = ResourceMonitor(
                output_path=output_path,
                interval_seconds=0.01,
            )
            monitor.start()
            _ = [str(index) for index in range(100)]
            monitor.stop()
            lines = output_path.read_text(encoding="utf-8").strip().splitlines()
        self.assertGreaterEqual(len(lines), 2)
        self.assertEqual(
            lines[0].split("\t"),
            [
                "elapsed_seconds",
                "timestamp_epoch",
                "rss_bytes",
                "rss_mb",
                "peak_rss_bytes",
                "peak_rss_mb",
            ],
        )

    def test_resource_monitor_rejects_bad_interval(self) -> None:
        """ResourceMonitor should reject non-positive sampling intervals."""
        with TemporaryDirectory() as tmpdir:
            with self.assertRaises(ValueError):
                ResourceMonitor(
                    output_path=Path(tmpdir) / "ram_usage.tsv",
                    interval_seconds=0,
                )


if __name__ == "__main__":
    unittest.main()
