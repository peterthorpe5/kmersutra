"""Predict labels with a trained KmerSutra classifier."""

from __future__ import annotations

import argparse
from pathlib import Path

from kmersutra.logging_utils import configure_logging
from kmersutra.ml import predict_tsv


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed arguments.
    """
    parser = argparse.ArgumentParser(
        description="Predict labels using a trained open-set KmerSutra-ML model."
    )
    parser.add_argument("--features_tsv", required=True)
    parser.add_argument("--model_json", required=True)
    parser.add_argument("--out_tsv", required=True)
    parser.add_argument("--novelty_scale", type=float, default=1.0)
    parser.add_argument("--verbose", action="store_true")
    return parser.parse_args()


def main() -> None:
    """Run classifier prediction."""
    args = parse_args()
    out_path = Path(args.out_tsv)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    logger = configure_logging(log_file=out_path.with_suffix(".log"), verbose=args.verbose)
    logger.info("Starting KmerSutra-ML prediction")
    logger.info("Feature table: %s", args.features_tsv)
    logger.info("Model JSON: %s", args.model_json)
    logger.info("Novelty scale: %.3f", args.novelty_scale)
    predictions = predict_tsv(
        features_tsv=args.features_tsv,
        model_json=args.model_json,
        output_tsv=args.out_tsv,
        novelty_scale=args.novelty_scale,
        logger=logger,
    )
    logger.info("Wrote %d prediction records", len(predictions))
    logger.info("Done")


if __name__ == "__main__":
    main()
