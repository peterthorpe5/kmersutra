"""Command-line entry point for restartable KmerSutra benchmarks."""

from __future__ import annotations

from kmersutra.benchmark_workflow import run_benchmark


def main() -> int:
    """Run the benchmark controller."""
    return run_benchmark()


if __name__ == "__main__":
    raise SystemExit(main())
