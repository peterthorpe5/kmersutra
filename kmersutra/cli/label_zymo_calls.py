"""Build clean Zymo AI validation labels from KmerSutra calls."""

from __future__ import annotations

import argparse
from pathlib import Path

from kmersutra.logging_utils import configure_logging
from kmersutra.zymo_truth import write_zymo_ai_feature_table


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed arguments.
    """
    parser = argparse.ArgumentParser(
        description="Label ZymoBIOMICS D6300 KmerSutra calls for AI validation."
    )
    parser.add_argument("--calls_table", required=True)
    parser.add_argument("--reference_label_map", default=None)
    parser.add_argument("--out_table", required=True)
    parser.add_argument("--out_category_counts", required=True)
    parser.add_argument("--out_coarse_label_counts", required=True)
    parser.add_argument("--verbose", action="store_true")
    return parser.parse_args()


def main() -> None:
    """Run Zymo truth labelling."""
    args = parse_args()
    out_path = Path(args.out_table)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    logger = configure_logging(
        log_file=out_path.with_suffix(".log"),
        verbose=args.verbose,
    )
    logger.info("Starting Zymo call truth labelling")
    logger.info("Calls table: %s", args.calls_table)
    logger.info("Reference-label map: %s", args.reference_label_map or "not supplied")
    records = write_zymo_ai_feature_table(
        calls_table=args.calls_table,
        output_table=args.out_table,
        reference_label_map=args.reference_label_map,
        category_counts_table=args.out_category_counts,
        coarse_label_counts_table=args.out_coarse_label_counts,
        logger=logger,
    )
    logger.info("Wrote %d Zymo AI feature records", len(records))
    logger.info("Done")


if __name__ == "__main__":
    main()
