"""Evaluate KmerSutra call-calibrator prediction tables."""

from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path

from kmersutra.ai_calibration import evaluate_predictions
from kmersutra.logging_utils import configure_logging
from kmersutra.table_io import read_records_table, write_records_table


def count_column(*, records: list[dict[str, object]], column: str) -> list[dict[str, object]]:
    """Count records by a column.

    Parameters
    ----------
    records : list of dict
        Input records.
    column : str
        Column name.

    Returns
    -------
    list[dict[str, object]]
        Count records.
    """
    counts = Counter(str(row.get(column, "")) for row in records)
    return [
        {column: value, "n_records": counts[value]}
        for value in sorted(counts)
    ]


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed arguments.
    """
    parser = argparse.ArgumentParser(
        description="Evaluate KmerSutra call-calibrator predictions."
    )
    parser.add_argument("--predictions_table", required=True)
    parser.add_argument("--out_metrics", required=True)
    parser.add_argument("--out_prediction_counts", required=True)
    parser.add_argument("--out_label_counts", required=True)
    parser.add_argument("--label_column", default="ml_report_label")
    parser.add_argument("--prediction_column", default="prediction")
    parser.add_argument(
        "--include_label",
        action="append",
        default=[],
        help="Optional truth label to include. May be supplied more than once.",
    )
    parser.add_argument(
        "--exclude_label",
        action="append",
        default=[],
        help="Optional truth label to exclude. May be supplied more than once.",
    )
    parser.add_argument("--verbose", action="store_true")
    return parser.parse_args()


def main() -> None:
    """Evaluate predictions."""
    args = parse_args()
    out_path = Path(args.out_metrics)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    logger = configure_logging(
        log_file=out_path.with_suffix(".log"),
        verbose=args.verbose,
    )
    logger.info("Reading prediction table: %s", args.predictions_table)
    records = read_records_table(input_path=args.predictions_table, logger=logger)
    include_labels = set(args.include_label or [])
    exclude_labels = set(args.exclude_label or [])
    if include_labels:
        records = [
            row for row in records
            if str(row.get(args.label_column, "")) in include_labels
        ]
        logger.info("Filtered to %d rows using include labels", len(records))
    if exclude_labels:
        records = [
            row for row in records
            if str(row.get(args.label_column, "")) not in exclude_labels
        ]
        logger.info("Filtered to %d rows after excluding labels", len(records))
    metrics = evaluate_predictions(
        predictions=records,
        label_column=args.label_column,
        prediction_column=args.prediction_column,
    )
    write_records_table(
        records=metrics,
        output_path=args.out_metrics,
        fieldnames=["label", "n", "tp", "fp", "fn", "precision", "recall", "f1"],
        logger=logger,
    )
    write_records_table(
        records=count_column(records=records, column=args.prediction_column),
        output_path=args.out_prediction_counts,
        fieldnames=[args.prediction_column, "n_records"],
        logger=logger,
    )
    write_records_table(
        records=count_column(records=records, column=args.label_column),
        output_path=args.out_label_counts,
        fieldnames=[args.label_column, "n_records"],
        logger=logger,
    )
    logger.info("Evaluation rows: %d", len(metrics))
    logger.info("Done")


if __name__ == "__main__":
    main()
