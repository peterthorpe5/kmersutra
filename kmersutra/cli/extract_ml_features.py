"""Extract KmerSutra-ML feature tables from read-level hit tables."""

from __future__ import annotations

import argparse
from pathlib import Path

from kmersutra.features import FEATURE_FIELDNAMES, load_hits_as_features
from kmersutra.io import write_tsv
from kmersutra.logging_utils import configure_logging


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed arguments.
    """
    parser = argparse.ArgumentParser(
        description="Extract sequence-level KmerSutra-ML features from hit tables."
    )
    parser.add_argument("--hits_tsv", required=True)
    parser.add_argument("--out_tsv", required=True)
    parser.add_argument("--verbose", action="store_true")
    return parser.parse_args()


def main() -> None:
    """Run feature extraction."""
    args = parse_args()
    out_path = Path(args.out_tsv)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    logger = configure_logging(log_file=out_path.with_suffix(".log"), verbose=args.verbose)
    logger.info("Starting KmerSutra-ML feature extraction")
    logger.info("Input hits TSV: %s", args.hits_tsv)
    logger.info("Output feature TSV: %s", args.out_tsv)
    features = load_hits_as_features(hits_tsv=args.hits_tsv, logger=logger)
    write_tsv(records=features, output_path=out_path, fieldnames=FEATURE_FIELDNAMES)
    logger.info("Wrote %d feature records", len(features))
    logger.info("Done")


if __name__ == "__main__":
    main()
