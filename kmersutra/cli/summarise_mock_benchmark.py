"""Summarise a manifest-driven mock-community benchmark."""

from __future__ import annotations

import argparse
from pathlib import Path

from kmersutra.logging_utils import configure_logging
from kmersutra.mock_benchmark_summary import summarise_mock_benchmark


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments.

    Args:
        argv: Optional explicit argument list.

    Returns:
        Parsed arguments.
    """
    parser = argparse.ArgumentParser(
        description="Summarise KmerSutra mock-community benchmark tasks."
    )
    parser.add_argument("--task_manifest", required=True)
    parser.add_argument("--truth_manifest", required=True)
    parser.add_argument("--out_dir", required=True)
    parser.add_argument("--verbose", action="store_true")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    """Run mock-community benchmark summarisation.

    Args:
        argv: Optional explicit argument list.

    Returns:
        Process exit status.
    """
    args = parse_args(argv)
    output_dir = Path(args.out_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    logger = configure_logging(
        log_file=output_dir / "summarise_mock_benchmark.log",
        verbose=args.verbose,
    )
    paths = summarise_mock_benchmark(
        task_manifest=args.task_manifest,
        truth_manifest=args.truth_manifest,
        output_dir=output_dir,
        logger=logger,
    )
    logger.info("Wrote %d mock-community summary tables", len(paths))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
