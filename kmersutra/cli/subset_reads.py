"""Create deterministic benchmark depth subsets."""

from __future__ import annotations

import argparse
from pathlib import Path

from kmersutra.depth_subsets import create_depth_subsets
from kmersutra.logging_utils import configure_logging


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments.

    Args:
        argv: Optional explicit argument list.

    Returns:
        Parsed arguments.
    """
    parser = argparse.ArgumentParser(
        description="Create deterministic nested FASTA/FASTQ benchmark subsets."
    )
    parser.add_argument("--input", required=True)
    parser.add_argument(
        "--input_format",
        choices=["fastq", "fasta"],
        required=True,
    )
    parser.add_argument("--out_dir", required=True)
    parser.add_argument("--fractions", nargs="+", type=float, required=True)
    parser.add_argument("--seeds", nargs="+", type=int, required=True)
    parser.add_argument("--sample_prefix", required=True)
    parser.add_argument(
        "--compress",
        action=argparse.BooleanOptionalAction,
        default=True,
    )
    parser.add_argument(
        "--decompressor",
        choices=["python", "pigz", "auto"],
        default="auto",
    )
    parser.add_argument("--manifest", default=None)
    parser.add_argument("--progress_interval", type=int, default=1_000_000)
    parser.add_argument("--verbose", action="store_true")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    """Create depth subsets.

    Args:
        argv: Optional explicit argument list.

    Returns:
        Process exit status.
    """
    args = parse_args(argv)
    output_dir = Path(args.out_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    logger = configure_logging(
        log_file=output_dir / "subset_reads.log",
        verbose=args.verbose,
    )
    create_depth_subsets(
        input_path=args.input,
        input_format=args.input_format,
        output_dir=output_dir,
        fractions=args.fractions,
        seeds=args.seeds,
        sample_prefix=args.sample_prefix,
        compress=args.compress,
        decompressor=args.decompressor,
        manifest_path=args.manifest,
        progress_interval=args.progress_interval,
        logger=logger,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
