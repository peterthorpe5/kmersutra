"""Train a lightweight open-set KmerSutra classifier."""

from __future__ import annotations

import argparse
from pathlib import Path

from kmersutra.logging_utils import configure_logging
from kmersutra.ml import train_model_from_tsv


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed arguments.
    """
    parser = argparse.ArgumentParser(
        description="Train a centroid-based open-set KmerSutra-ML classifier."
    )
    parser.add_argument("--features_tsv", required=True)
    parser.add_argument("--label_column", required=True)
    parser.add_argument("--out_model_json", required=True)
    parser.add_argument("--out_summary_tsv", required=True)
    parser.add_argument("--distance_quantile", type=float, default=0.95)
    parser.add_argument("--unknown_label", default="unknown_or_unresolved")
    parser.add_argument("--verbose", action="store_true")
    return parser.parse_args()


def main() -> None:
    """Run classifier training."""
    args = parse_args()
    model_path = Path(args.out_model_json)
    model_path.parent.mkdir(parents=True, exist_ok=True)
    logger = configure_logging(log_file=model_path.with_suffix(".log"), verbose=args.verbose)
    logger.info("Starting KmerSutra-ML classifier training")
    logger.info("Feature table: %s", args.features_tsv)
    logger.info("Label column: %s", args.label_column)
    logger.info("Distance quantile: %.3f", args.distance_quantile)
    model = train_model_from_tsv(
        features_tsv=args.features_tsv,
        label_column=args.label_column,
        model_json=args.out_model_json,
        summary_tsv=args.out_summary_tsv,
        distance_quantile=args.distance_quantile,
        unknown_label=args.unknown_label,
        logger=logger,
    )
    logger.info("Model labels: %s", "; ".join(sorted(model.class_counts)))
    logger.info("Model written to %s", args.out_model_json)
    logger.info("Done")


if __name__ == "__main__":
    main()
