"""Tests for lightweight profiling helpers."""

from __future__ import annotations

import unittest

from kmersutra.profiling import WorkflowProfiler


class TestProfiling(unittest.TestCase):
    """Tests for workflow timing records."""

    def test_profiler_records_named_stage(self) -> None:
        """Profiler should record stage name, elapsed seconds and detail."""
        profiler = WorkflowProfiler()
        with profiler.time_stage(stage="example", detail="demo"):
            value = sum([1, 2, 3])
        self.assertEqual(value, 6)
        records = profiler.to_records()
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0]["stage"], "example")
        self.assertEqual(records[0]["detail"], "demo")
        self.assertGreaterEqual(float(records[0]["seconds"]), 0.0)
