"""Train and evaluate an interpretable KmerSutra AI call calibrator."""

from __future__ import annotations

import argparse
from pathlib import Path

from kmersutra.ai_calibration import train_evaluate_call_calibrator
from kmersutra.logging_utils import configure_logging


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed arguments.
    """
    parser = argparse.ArgumentParser(
        description="Train an interpretable AI call calibrator for KmerSutra outputs."
    )
    parser.add_argument(
        "--training_table",
        default=None,
        help="Input training table. Supported suffixes: .tsv, .tsv.gz and .parquet.",
    )
    parser.add_argument(
        "--training_tsv",
        default=None,
        help=(
            "Backward-compatible alias for --training_table. Despite the name, "
            ".tsv, .tsv.gz and .parquet are supported."
        ),
    )
    parser.add_argument("--out_model_json", required=True)
    parser.add_argument(
        "--out_summary_table",
        default=None,
        help="Output summary table. Supported suffixes: .tsv, .tsv.gz and .parquet.",
    )
    parser.add_argument(
        "--out_summary_tsv",
        default=None,
        help=(
            "Backward-compatible alias for --out_summary_table. Despite the name, "
            ".tsv, .tsv.gz and .parquet are supported."
        ),
    )
    parser.add_argument(
        "--out_evaluation_table",
        default=None,
        help="Output evaluation table. Supported suffixes: .tsv, .tsv.gz and .parquet.",
    )
    parser.add_argument(
        "--out_evaluation_tsv",
        default=None,
        help=(
            "Backward-compatible alias for --out_evaluation_table. Despite the name, "
            ".tsv, .tsv.gz and .parquet are supported."
        ),
    )
    parser.add_argument("--label_column", default="ml_report_label")
    parser.add_argument("--feature_columns", nargs="+", default=None)
    parser.add_argument("--group_columns", nargs="+", default=["sample_id"])
    parser.add_argument("--test_fraction", type=float, default=0.2)
    parser.add_argument("--distance_quantile", type=float, default=0.95)
    parser.add_argument("--unknown_label", default="unknown_or_unresolved")
    parser.add_argument("--verbose", action="store_true")
    return parser.parse_args()


def main() -> None:
    """Run call-calibrator training and evaluation."""
    args = parse_args()
    training_table = args.training_table or args.training_tsv
    summary_table = args.out_summary_table or args.out_summary_tsv
    evaluation_table = args.out_evaluation_table or args.out_evaluation_tsv
    if training_table is None:
        raise ValueError("Either --training_table or --training_tsv is required")
    if summary_table is None:
        raise ValueError("Either --out_summary_table or --out_summary_tsv is required")
    if evaluation_table is None:
        raise ValueError("Either --out_evaluation_table or --out_evaluation_tsv is required")
    model_path = Path(args.out_model_json)
    model_path.parent.mkdir(parents=True, exist_ok=True)
    logger = configure_logging(log_file=model_path.with_suffix(".log"), verbose=args.verbose)
    logger.info("Starting KmerSutra AI call-calibrator training")
    logger.info("Training table: %s", training_table)
    model, predictions, metrics = train_evaluate_call_calibrator(
        training_table=training_table,
        model_json=args.out_model_json,
        summary_table=summary_table,
        evaluation_table=evaluation_table,
        label_column=args.label_column,
        feature_columns=args.feature_columns,
        test_fraction=args.test_fraction,
        group_columns=args.group_columns,
        distance_quantile=args.distance_quantile,
        unknown_label=args.unknown_label,
        logger=logger,
    )
    logger.info("Model labels: %s", "; ".join(sorted(model.class_counts)))
    logger.info("Held-out predictions: %d", len(predictions))
    logger.info("Evaluation rows: %d", len(metrics))
    logger.info("Done")


if __name__ == "__main__":
    main()
