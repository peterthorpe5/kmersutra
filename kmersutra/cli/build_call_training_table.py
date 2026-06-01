"""Build an AI-ready KmerSutra call-training table."""

from __future__ import annotations

import argparse
from pathlib import Path

from kmersutra.ai_calibration import write_call_training_table_from_table
from kmersutra.logging_utils import configure_logging


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed arguments.
    """
    parser = argparse.ArgumentParser(
        description="Build an AI-ready call calibration table from KmerSutra calls."
    )
    parser.add_argument(
        "--calls_table",
        default=None,
        help="Input call table. Supported suffixes: .tsv, .tsv.gz and .parquet.",
    )
    parser.add_argument(
        "--calls_tsv",
        default=None,
        help=(
            "Backward-compatible alias for --calls_table. Despite the name, "
            ".tsv, .tsv.gz and .parquet are supported."
        ),
    )
    parser.add_argument(
        "--out_table",
        default=None,
        help="Output training table. Supported suffixes: .tsv, .tsv.gz and .parquet.",
    )
    parser.add_argument(
        "--out_tsv",
        default=None,
        help=(
            "Backward-compatible alias for --out_table. Despite the name, "
            ".tsv, .tsv.gz and .parquet are supported."
        ),
    )
    parser.add_argument(
        "--lca_table",
        default=None,
        help=(
            "Optional LCA summary table from kmersutra-summarise-lca. "
            "Supported suffixes: .tsv, .tsv.gz and .parquet."
        ),
    )
    parser.add_argument(
        "--lca_sample_column",
        default="sample_id",
        help="Sample identifier column used to join LCA features.",
    )
    parser.add_argument("--label_column", default="ml_report_label")
    parser.add_argument("--drop_not_detected", action="store_true")
    parser.add_argument("--max_not_detected", type=int, default=None)
    parser.add_argument("--verbose", action="store_true")
    return parser.parse_args()


def main() -> None:
    """Run call-training table construction."""
    args = parse_args()
    calls_table = args.calls_table or args.calls_tsv
    output_table = args.out_table or args.out_tsv
    if calls_table is None:
        raise ValueError("Either --calls_table or --calls_tsv is required")
    if output_table is None:
        raise ValueError("Either --out_table or --out_tsv is required")
    out_path = Path(output_table)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    logger = configure_logging(log_file=out_path.with_suffix(".log"), verbose=args.verbose)
    logger.info("Starting KmerSutra AI call-training table construction")
    logger.info("Call table: %s", calls_table)
    logger.info("Output table: %s", output_table)
    logger.info("LCA table: %s", args.lca_table or "not supplied")
    records = write_call_training_table_from_table(
        calls_table=calls_table,
        output_table=out_path,
        label_column=args.label_column,
        include_not_detected=not args.drop_not_detected,
        max_not_detected=args.max_not_detected,
        lca_table=args.lca_table,
        lca_sample_column=args.lca_sample_column,
        logger=logger,
    )
    logger.info("Wrote %d AI training records to %s", len(records), out_path)
    logger.info("Done")


if __name__ == "__main__":
    main()
