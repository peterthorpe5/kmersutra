"""Resource-monitoring utilities for KmerSutra workflows."""

from __future__ import annotations

import logging
import os
import platform
import resource
import threading
import time
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO


@dataclass(frozen=True)
class ResourceSnapshot:
    """A point-in-time resource-usage measurement.

    Attributes
    ----------
    elapsed_seconds : float
        Seconds elapsed since monitoring started.
    timestamp_epoch : float
        Unix epoch timestamp at which the snapshot was taken.
    rss_bytes : int
        Current resident set size in bytes when available.
    peak_rss_bytes : int
        Peak resident set size in bytes when available.
    """

    elapsed_seconds: float
    timestamp_epoch: float
    rss_bytes: int
    peak_rss_bytes: int

    def to_record(self) -> dict[str, object]:
        """Convert the resource snapshot to a serialisable record.

        Returns
        -------
        dict[str, object]
            Dictionary representation suitable for TSV output.
        """
        return {
            "elapsed_seconds": f"{self.elapsed_seconds:.3f}",
            "timestamp_epoch": f"{self.timestamp_epoch:.3f}",
            "rss_bytes": self.rss_bytes,
            "rss_mb": f"{self.rss_bytes / (1024 ** 2):.3f}",
            "peak_rss_bytes": self.peak_rss_bytes,
            "peak_rss_mb": f"{self.peak_rss_bytes / (1024 ** 2):.3f}",
        }


def get_current_rss_bytes() -> int:
    """Return the current process resident set size in bytes.

    Returns
    -------
    int
        Current resident set size in bytes. Returns 0 if the platform does not
        expose the value.
    """
    statm_path = Path("/proc/self/statm")
    if statm_path.exists():
        try:
            fields = statm_path.read_text(encoding="utf-8").split()
            resident_pages = int(fields[1])
            return resident_pages * os.sysconf("SC_PAGE_SIZE")
        except (IndexError, OSError, ValueError):
            return 0
    return 0


def get_peak_rss_bytes() -> int:
    """Return the peak resident set size in bytes.

    Returns
    -------
    int
        Peak resident set size in bytes, or 0 if unavailable.
    """
    usage = resource.getrusage(resource.RUSAGE_SELF)
    value = int(usage.ru_maxrss)
    system_name = platform.system().lower()
    if system_name == "darwin":
        return value
    return value * 1024


class ResourceMonitor:
    """Background RAM monitor for long-running KmerSutra commands.

    The monitor writes one TSV row per interval. It is dependency-free and is
    intended for cluster jobs where external tools such as ``/usr/bin/time`` may
    not be available.
    """

    def __init__(
        self,
        *,
        output_path: str | Path,
        interval_seconds: float = 60.0,
        logger: logging.Logger | None = None,
    ) -> None:
        """Initialise the resource monitor.

        Parameters
        ----------
        output_path : str or pathlib.Path
            TSV file to write.
        interval_seconds : float, optional
            Sampling interval in seconds.
        logger : logging.Logger | None, optional
            Optional logger for start and stop messages.
        """
        if interval_seconds <= 0:
            raise ValueError("interval_seconds must be positive")
        self.output_path = Path(output_path)
        self.interval_seconds = float(interval_seconds)
        self.logger = logger
        self._start_time = 0.0
        self._stop_event = threading.Event()
        self._thread: threading.Thread | None = None
        self._handle: TextIO | None = None

    def __enter__(self) -> "ResourceMonitor":
        """Start monitoring and return the monitor instance.

        Returns
        -------
        ResourceMonitor
            Started monitor.
        """
        self.start()
        return self

    def __exit__(self, exc_type, exc_value, traceback) -> None:  # noqa: ANN001
        """Stop monitoring when leaving a context manager."""
        self.stop()

    def start(self) -> None:
        """Start the background monitoring thread."""
        self.output_path.parent.mkdir(parents=True, exist_ok=True)
        self._handle = self.output_path.open("w", encoding="utf-8")
        self._handle.write(
            "elapsed_seconds\ttimestamp_epoch\trss_bytes\trss_mb\t"
            "peak_rss_bytes\tpeak_rss_mb\n"
        )
        self._handle.flush()
        self._start_time = time.time()
        self._write_snapshot()
        self._thread = threading.Thread(target=self._run, daemon=True)
        self._thread.start()
        if self.logger:
            self.logger.info("Started RAM monitor: %s", self.output_path)

    def stop(self) -> None:
        """Stop the background monitoring thread and close the TSV file."""
        self._stop_event.set()
        if self._thread is not None:
            self._thread.join(timeout=max(1.0, self.interval_seconds + 1.0))
        self._write_snapshot()
        if self._handle is not None:
            self._handle.close()
            self._handle = None
        if self.logger:
            self.logger.info("Stopped RAM monitor: %s", self.output_path)

    def _run(self) -> None:
        """Run the background sampling loop."""
        while not self._stop_event.wait(self.interval_seconds):
            self._write_snapshot()

    def _write_snapshot(self) -> None:
        """Write a single resource snapshot to the TSV file."""
        if self._handle is None:
            return
        now = time.time()
        snapshot = ResourceSnapshot(
            elapsed_seconds=now - self._start_time,
            timestamp_epoch=now,
            rss_bytes=get_current_rss_bytes(),
            peak_rss_bytes=get_peak_rss_bytes(),
        )
        record = snapshot.to_record()
        self._handle.write(
            f"{record['elapsed_seconds']}\t{record['timestamp_epoch']}\t"
            f"{record['rss_bytes']}\t{record['rss_mb']}\t"
            f"{record['peak_rss_bytes']}\t{record['peak_rss_mb']}\n"
        )
        self._handle.flush()
