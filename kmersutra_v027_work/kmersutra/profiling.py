"""Lightweight profiling helpers for KmerSutra workflows."""

from __future__ import annotations

import time
from collections.abc import Iterator
from contextlib import contextmanager
from dataclasses import dataclass


@dataclass(frozen=True)
class TimingRecord:
    """A single named timing measurement.

    Attributes
    ----------
    stage : str
        Name of the workflow stage.
    seconds : float
        Elapsed wall-clock time in seconds.
    detail : str
        Optional human-readable details about the stage.
    """

    stage: str
    seconds: float
    detail: str = ""

    def to_record(self) -> dict[str, object]:
        """Convert the timing record to a serialisable dictionary.

        Returns
        -------
        dict[str, object]
            Dictionary representation of the timing record.
        """
        return {
            "stage": self.stage,
            "seconds": f"{self.seconds:.6f}",
            "detail": self.detail,
        }


class WorkflowProfiler:
    """Collect named wall-clock timing records.

    The profiler is deliberately simple so it can be used in command-line
    workflows without adding extra dependencies.
    """

    def __init__(self) -> None:
        """Initialise an empty profiler."""
        self.records: list[TimingRecord] = []

    @contextmanager
    def time_stage(self, *, stage: str, detail: str = "") -> Iterator[None]:
        """Measure a block of code as a named workflow stage.

        Parameters
        ----------
        stage : str
            Name of the workflow stage.
        detail : str, optional
            Optional details to write alongside the timing.

        Yields
        ------
        None
            Context manager body.
        """
        start = time.perf_counter()
        try:
            yield
        finally:
            end = time.perf_counter()
            self.records.append(
                TimingRecord(stage=stage, seconds=end - start, detail=detail)
            )

    def to_records(self) -> list[dict[str, object]]:
        """Return profiler records as serialisable dictionaries.

        Returns
        -------
        list[dict[str, object]]
            Timing records suitable for TSV output.
        """
        return [record.to_record() for record in self.records]
