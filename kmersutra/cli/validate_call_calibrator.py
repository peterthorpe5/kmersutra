"""Run KmerSutra call-calibrator holdout validation."""

from __future__ import annotations

import argparse
from pathlib import Path

from kmersutra.call_validation import validate_call_calibrator_holdouts
from kmersutra.logging_utils import configure_logging


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed arguments.
    """
    parser = argparse.ArgumentParser(
        description="Run full holdout validation for a KmerSutra call calibrator."
    )
    parser.add_argument("--training_table", required=True)
    parser.add_argument("--out_dir", required=True)
    parser.add_argument(
        "--feature_profile",
        choices=["legacy", "safe_raw", "safe_transformed"],
        default="safe_transformed",
    )
    parser.add_argument("--feature_columns", nargs="+", default=None)
    parser.add_argument("--label_column", default="ml_report_label")
    parser.add_argument("--distance_quantile", type=float, default=0.95)
    parser.add_argument("--sample_test_fraction", type=float, default=0.20)
    parser.add_argument("--unknown_label", default="unknown_or_unresolved")
    parser.add_argument("--verbose", action="store_true")
    return parser.parse_args()


def main() -> None:
    """Run holdout validation."""
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    logger = configure_logging(
        log_file=out_dir / "validate_call_calibrator.log",
        verbose=args.verbose,
    )
    logger.info("Starting KmerSutra call-calibrator validation")
    logger.info("Training table: %s", args.training_table)
    logger.info("Output directory: %s", out_dir)
    logger.info("Feature profile: %s", args.feature_profile)
    manifest = validate_call_calibrator_holdouts(
        training_table=args.training_table,
        out_dir=out_dir,
        feature_profile=args.feature_profile,
        feature_columns=args.feature_columns,
        label_column=args.label_column,
        distance_quantile=args.distance_quantile,
        sample_test_fraction=args.sample_test_fraction,
        unknown_label=args.unknown_label,
        logger=logger,
    )
    logger.info("Validation splits completed: %s", manifest["n_completed_splits"])
    logger.info("Done")


if __name__ == "__main__":
    main()
