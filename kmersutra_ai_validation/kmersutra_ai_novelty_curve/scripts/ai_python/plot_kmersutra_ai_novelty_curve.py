#!/usr/bin/env python3
"""Summarise and plot KmerSutra AI novelty-scale sweep results.

This script is intended for the external Zymo AI-calibration sweep produced by
KmerSutra v0.50.1-style validation runs. It scans a directory containing one or
more run folders, reads each run's ``external_zymo_predictions.tsv.gz`` table,
extracts the novelty scale from the corresponding prediction log, and writes
summary tables plus publication/supplement-friendly operating-point plots.

The primary plot is an open-set novelty-scale operating curve rather than a
classical probabilistic ROC curve. This is because the KmerSutra AI layer is an
open-set evidence calibrator: rows are assigned to the nearest learned evidence
class and then accepted or rejected according to a novelty-distance threshold.
For completeness, the script can also draw precision-recall and ROC-style
operating-point plots for the strict expected-target label.

Outputs are tab-separated, with PNG and PDF figures when matplotlib is
available.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import logging
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

STRICT_CATEGORIES = {
    "expected_species",
    "expected_reference_label",
    "near_neighbour_evidence",
    "true_off_target_reportable",
    "not_detected",
}

EXPECTED_CATEGORIES = {
    "expected_species",
    "expected_reference_label",
}

NON_EXPECTED_STRICT_CATEGORIES = STRICT_CATEGORIES - EXPECTED_CATEGORIES

FLOAT_TOLERANCE = 1e-9


@dataclass(frozen=True)
class SweepSummary:
    """Summary metrics for one novelty-scale run.

    Attributes
    ----------
    run_name : str
        Name of the run directory.
    novelty_scale : float
        Novelty scale extracted from the run log or settings.
    prediction_table : pathlib.Path
        Input prediction table used for the run.
    n_all_rows : int
        Total rows in the prediction table.
    n_strict_rows : int
        Rows included in the strict operating-point summary.
    n_strict_expected : int
        Strict expected organism rows.
    n_strict_expected_tp : int
        Strict expected rows predicted as expected_target.
    n_strict_expected_fn : int
        Strict expected rows not predicted as expected_target.
    n_strict_expected_fp : int
        Strict non-expected rows predicted as expected_target.
    n_strict_not_detected : int
        Strict not_detected rows.
    n_strict_not_detected_tp : int
        Strict not_detected rows predicted as not_detected.
    n_strict_not_detected_fn : int
        Strict not_detected rows not predicted as not_detected.
    n_neighbour_rows : int
        Strict near-neighbour rows.
    n_neighbour_expected_overpromoted : int
        Strict near-neighbour rows predicted as expected_target.
    n_strict_unknown_or_unresolved : int
        Strict rows predicted as unknown_or_unresolved.
    n_strict_correct : int
        Strict rows correctly predicted using the expected/not-detected/
        observed-below-threshold coarse labels.
    strict_accuracy : float
        Accuracy over strict rows.
    expected_precision : float
        Precision for expected_target over strict rows.
    expected_recall : float
        Recall for expected_target over strict rows.
    expected_f1 : float
        F1 for expected_target over strict rows.
    expected_fpr : float
        False-positive rate for expected_target over strict non-expected rows.
    """

    run_name: str
    novelty_scale: float
    prediction_table: Path
    n_all_rows: int
    n_strict_rows: int
    n_strict_expected: int
    n_strict_expected_tp: int
    n_strict_expected_fn: int
    n_strict_expected_fp: int
    n_strict_not_detected: int
    n_strict_not_detected_tp: int
    n_strict_not_detected_fn: int
    n_neighbour_rows: int
    n_neighbour_expected_overpromoted: int
    n_strict_unknown_or_unresolved: int
    n_strict_correct: int
    strict_accuracy: float
    expected_precision: float
    expected_recall: float
    expected_f1: float
    expected_fpr: float

    def as_dict(self) -> dict[str, object]:
        """Return this summary as a TSV-friendly dictionary.

        Returns
        -------
        dict[str, object]
            Summary values.
        """
        return {
            "run_name": self.run_name,
            "novelty_scale": format_float(self.novelty_scale),
            "prediction_table": str(self.prediction_table),
            "n_all_rows": self.n_all_rows,
            "n_strict_rows": self.n_strict_rows,
            "n_strict_expected": self.n_strict_expected,
            "n_strict_expected_tp": self.n_strict_expected_tp,
            "n_strict_expected_fn": self.n_strict_expected_fn,
            "n_strict_expected_fp": self.n_strict_expected_fp,
            "strict_expected_precision": format_float(self.expected_precision),
            "strict_expected_recall": format_float(self.expected_recall),
            "strict_expected_f1": format_float(self.expected_f1),
            "strict_expected_fpr": format_float(self.expected_fpr),
            "n_strict_not_detected": self.n_strict_not_detected,
            "n_strict_not_detected_tp": self.n_strict_not_detected_tp,
            "n_strict_not_detected_fn": self.n_strict_not_detected_fn,
            "strict_not_detected_recall": format_float(
                safe_divide(
                    self.n_strict_not_detected_tp,
                    self.n_strict_not_detected,
                )
            ),
            "n_neighbour_rows": self.n_neighbour_rows,
            "n_neighbour_expected_overpromoted": (
                self.n_neighbour_expected_overpromoted
            ),
            "n_strict_unknown_or_unresolved": self.n_strict_unknown_or_unresolved,
            "n_strict_correct": self.n_strict_correct,
            "strict_accuracy": format_float(self.strict_accuracy),
        }


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Summarise KmerSutra AI external Zymo novelty-scale sweep results "
            "and draw operating-point curves."
        )
    )
    parser.add_argument(
        "--sweep_root",
        required=True,
        help=(
            "Directory containing KmerSutra AI external Zymo run folders, or "
            "a single run folder. The script searches recursively for "
            "external_zymo_predictions.tsv.gz."
        ),
    )
    parser.add_argument(
        "--out_dir",
        required=True,
        help="Output directory for TSV summaries and figures.",
    )
    parser.add_argument(
        "--run_glob",
        default="external_zymo_predictions.tsv.gz",
        help=(
            "Prediction table basename to search for. Default: "
            "external_zymo_predictions.tsv.gz."
        ),
    )
    parser.add_argument(
        "--exclude_run_regex",
        default="",
        help=(
            "Optional regular expression. Any prediction table whose path "
            "matches this expression is skipped. Useful for excluding known "
            "provenance-collision runs."
        ),
    )
    parser.add_argument(
        "--figure_prefix",
        default="kmersutra_ai_novelty_sweep",
        help="Prefix for output figure filenames.",
    )
    parser.add_argument(
        "--no_plots",
        action="store_true",
        help="Only write TSV summaries; do not draw figures.",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable verbose logging.",
    )
    return parser.parse_args()


def configure_logging(*, verbose: bool, out_dir: Path) -> logging.Logger:
    """Configure console and file logging.

    Parameters
    ----------
    verbose : bool
        Whether to emit debug logs to the console.
    out_dir : pathlib.Path
        Directory where the log file will be written.

    Returns
    -------
    logging.Logger
        Configured logger.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger("kmersutra_ai_novelty_sweep")
    logger.setLevel(logging.DEBUG)
    logger.handlers.clear()

    formatter = logging.Formatter(
        "%(asctime)s %(levelname)s [%(name)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.DEBUG if verbose else logging.INFO)
    stream_handler.setFormatter(formatter)
    logger.addHandler(stream_handler)

    file_handler = logging.FileHandler(out_dir / "plot_kmersutra_ai_novelty_curve.log")
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    return logger


def safe_divide(numerator: float, denominator: float) -> float:
    """Safely divide two values.

    Parameters
    ----------
    numerator : float
        Numerator.
    denominator : float
        Denominator.

    Returns
    -------
    float
        Division result, or 0.0 when the denominator is zero.
    """
    if abs(denominator) < FLOAT_TOLERANCE:
        return 0.0
    return numerator / denominator


def format_float(value: float) -> str:
    """Format a float for TSV output.

    Parameters
    ----------
    value : float
        Numeric value.

    Returns
    -------
    str
        Compact numeric string.
    """
    if math.isnan(value):
        return "nan"
    return f"{value:.6g}"


def open_text(path: Path):
    """Open a plain text or gzip-compressed text file.

    Parameters
    ----------
    path : pathlib.Path
        Input path.

    Returns
    -------
    file-like
        Open text handle.
    """
    if path.suffix == ".gz":
        return gzip.open(path, "rt", newline="")
    return path.open("r", newline="")


def read_tsv_records(path: Path) -> list[dict[str, str]]:
    """Read a TSV or TSV.GZ file into dictionaries.

    Parameters
    ----------
    path : pathlib.Path
        Input table.

    Returns
    -------
    list[dict[str, str]]
        Table records.

    Raises
    ------
    ValueError
        If the table has no header.
    """
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"No header found in {path}")
        return [dict(row) for row in reader]


def write_tsv(path: Path, rows: Iterable[dict[str, object]], fieldnames: list[str]) -> None:
    """Write records to a TSV file.

    Parameters
    ----------
    path : pathlib.Path
        Output path.
    rows : iterable of dict
        Output records.
    fieldnames : list[str]
        Column order.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=fieldnames,
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def find_prediction_tables(
    *,
    sweep_root: Path,
    basename: str,
    exclude_run_regex: str,
    logger: logging.Logger,
) -> list[Path]:
    """Find prediction tables under a sweep directory.

    Parameters
    ----------
    sweep_root : pathlib.Path
        Root directory or single run directory.
    basename : str
        Prediction table basename.
    exclude_run_regex : str
        Optional exclusion regular expression.
    logger : logging.Logger
        Logger.

    Returns
    -------
    list[pathlib.Path]
        Sorted prediction table paths.
    """
    if not sweep_root.exists():
        raise FileNotFoundError(f"Sweep root does not exist: {sweep_root}")

    paths = sorted(sweep_root.rglob(basename))
    if exclude_run_regex:
        pattern = re.compile(exclude_run_regex)
        kept = []
        for path in paths:
            if pattern.search(str(path)):
                logger.info("Excluding run path by regex: %s", path)
                continue
            kept.append(path)
        paths = kept

    if not paths:
        raise FileNotFoundError(
            f"No {basename} files found below {sweep_root}"
        )

    return paths


def get_run_name(prediction_table: Path) -> str:
    """Infer the run directory name for a prediction table.

    Parameters
    ----------
    prediction_table : pathlib.Path
        Path to external_zymo_predictions.tsv.gz.

    Returns
    -------
    str
        Run directory name.
    """
    output_dir = prediction_table.parent
    return output_dir.parent.name if output_dir.name == "outputs" else output_dir.name


def extract_novelty_scale(prediction_table: Path) -> float:
    """Extract the novelty scale for a run.

    The function first reads ``external_zymo_predictions.tsv.log`` and searches
    for a line such as ``Novelty scale: 3.000``. If that is unavailable, it
    falls back to ``run_submission_settings.tsv`` and looks for a
    ``novelty_scale`` setting.

    Parameters
    ----------
    prediction_table : pathlib.Path
        Prediction table path.

    Returns
    -------
    float
        Novelty scale.

    Raises
    ------
    ValueError
        If the novelty scale cannot be found.
    """
    output_dir = prediction_table.parent
    log_path = output_dir / "external_zymo_predictions.tsv.log"
    if log_path.exists():
        text = log_path.read_text(errors="replace")
        match = re.search(r"Novelty scale:\s*([0-9.]+)", text)
        if match:
            return float(match.group(1))

    settings_path = output_dir / "run_submission_settings.tsv"
    if settings_path.exists():
        with settings_path.open("r", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                if row.get("setting", "").strip().lower() == "novelty_scale":
                    return float(row.get("value", ""))

    raise ValueError(
        "Could not extract novelty scale for "
        f"{prediction_table}. Expected external_zymo_predictions.tsv.log "
        "or run_submission_settings.tsv."
    )


def strict_expected_truth(row: dict[str, str]) -> bool:
    """Return whether a row is a strict expected organism row.

    Parameters
    ----------
    row : dict[str, str]
        Prediction table row.

    Returns
    -------
    bool
        True if the row is strict expected.
    """
    return row.get("zymo_truth_category", "") in EXPECTED_CATEGORIES


def strict_non_expected_truth(row: dict[str, str]) -> bool:
    """Return whether a row is a strict non-expected row.

    Parameters
    ----------
    row : dict[str, str]
        Prediction table row.

    Returns
    -------
    bool
        True if the row is strict but not expected.
    """
    return row.get("zymo_truth_category", "") in NON_EXPECTED_STRICT_CATEGORIES


def row_is_correct_strict(row: dict[str, str]) -> bool:
    """Return whether a strict-row prediction is correct.

    Parameters
    ----------
    row : dict[str, str]
        Prediction row.

    Returns
    -------
    bool
        True for correct strict predictions.
    """
    category = row.get("zymo_truth_category", "")
    prediction = row.get("prediction", "")

    if category in EXPECTED_CATEGORIES:
        return prediction == "expected_target"
    if category == "not_detected":
        return prediction == "not_detected"
    if category in {"near_neighbour_evidence", "true_off_target_reportable"}:
        return prediction == "observed_below_threshold"
    return False


def summarise_prediction_table(prediction_table: Path) -> SweepSummary:
    """Summarise one prediction table.

    Parameters
    ----------
    prediction_table : pathlib.Path
        Path to ``external_zymo_predictions.tsv.gz``.

    Returns
    -------
    SweepSummary
        Summary metrics.
    """
    records = read_tsv_records(prediction_table)
    novelty_scale = extract_novelty_scale(prediction_table)
    run_name = get_run_name(prediction_table)

    strict_records = [
        row
        for row in records
        if row.get("zymo_truth_category", "") in STRICT_CATEGORIES
    ]
    expected_records = [row for row in strict_records if strict_expected_truth(row)]
    non_expected_strict_records = [
        row for row in strict_records if strict_non_expected_truth(row)
    ]
    not_detected_records = [
        row
        for row in strict_records
        if row.get("zymo_truth_category", "") == "not_detected"
    ]
    neighbour_records = [
        row
        for row in strict_records
        if row.get("zymo_truth_category", "") == "near_neighbour_evidence"
    ]

    expected_tp = sum(
        1 for row in expected_records
        if row.get("prediction", "") == "expected_target"
    )
    expected_fn = len(expected_records) - expected_tp
    expected_fp = sum(
        1 for row in non_expected_strict_records
        if row.get("prediction", "") == "expected_target"
    )

    not_detected_tp = sum(
        1 for row in not_detected_records
        if row.get("prediction", "") == "not_detected"
    )
    not_detected_fn = len(not_detected_records) - not_detected_tp

    neighbour_overpromoted = sum(
        1 for row in neighbour_records
        if row.get("prediction", "") == "expected_target"
    )
    unresolved = sum(
        1 for row in strict_records
        if row.get("prediction", "") == "unknown_or_unresolved"
    )
    strict_correct = sum(1 for row in strict_records if row_is_correct_strict(row))

    precision = safe_divide(expected_tp, expected_tp + expected_fp)
    recall = safe_divide(expected_tp, len(expected_records))
    f1 = safe_divide(2 * precision * recall, precision + recall)
    fpr = safe_divide(expected_fp, len(non_expected_strict_records))

    return SweepSummary(
        run_name=run_name,
        novelty_scale=novelty_scale,
        prediction_table=prediction_table,
        n_all_rows=len(records),
        n_strict_rows=len(strict_records),
        n_strict_expected=len(expected_records),
        n_strict_expected_tp=expected_tp,
        n_strict_expected_fn=expected_fn,
        n_strict_expected_fp=expected_fp,
        n_strict_not_detected=len(not_detected_records),
        n_strict_not_detected_tp=not_detected_tp,
        n_strict_not_detected_fn=not_detected_fn,
        n_neighbour_rows=len(neighbour_records),
        n_neighbour_expected_overpromoted=neighbour_overpromoted,
        n_strict_unknown_or_unresolved=unresolved,
        n_strict_correct=strict_correct,
        strict_accuracy=safe_divide(strict_correct, len(strict_records)),
        expected_precision=precision,
        expected_recall=recall,
        expected_f1=f1,
        expected_fpr=fpr,
    )


def deduplicate_by_scale(
    *,
    summaries: list[SweepSummary],
    logger: logging.Logger,
) -> list[SweepSummary]:
    """Deduplicate summaries with the same novelty scale.

    When multiple runs have the same novelty scale, the function retains the
    lexicographically last run name. This is useful after exploratory reruns.

    Parameters
    ----------
    summaries : list[SweepSummary]
        Input summaries.
    logger : logging.Logger
        Logger.

    Returns
    -------
    list[SweepSummary]
        Deduplicated summaries sorted by novelty scale.
    """
    by_scale: dict[float, SweepSummary] = {}
    for summary in sorted(summaries, key=lambda item: item.run_name):
        previous = by_scale.get(summary.novelty_scale)
        if previous is not None:
            logger.warning(
                "Duplicate novelty scale %.3f: replacing %s with %s",
                summary.novelty_scale,
                previous.run_name,
                summary.run_name,
            )
        by_scale[summary.novelty_scale] = summary
    return [by_scale[scale] for scale in sorted(by_scale)]


def draw_figures(
    *,
    summaries: list[SweepSummary],
    out_dir: Path,
    prefix: str,
    logger: logging.Logger,
) -> None:
    """Draw novelty-scale, precision-recall and ROC-style plots.

    Parameters
    ----------
    summaries : list[SweepSummary]
        Sorted sweep summaries.
    out_dir : pathlib.Path
        Output directory.
    prefix : str
        Output filename prefix.
    logger : logging.Logger
        Logger.
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        logger.warning("matplotlib is not available; skipping figures")
        return

    scales = [item.novelty_scale for item in summaries]
    expected_recovery = [item.expected_recall for item in summaries]
    not_detected_recovery = [
        safe_divide(item.n_strict_not_detected_tp, item.n_strict_not_detected)
        for item in summaries
    ]
    strict_accuracy = [item.strict_accuracy for item in summaries]
    unknown_rows = [item.n_strict_unknown_or_unresolved for item in summaries]
    near_neighbour_overpromotions = [
        item.n_neighbour_expected_overpromoted for item in summaries
    ]

    fig, primary_axis = plt.subplots(figsize=(8, 5))
    primary_axis.plot(scales, expected_recovery, marker="o", label="Strict expected recovery")
    primary_axis.plot(scales, not_detected_recovery, marker="o", label="Strict not-detected recovery")
    primary_axis.plot(scales, strict_accuracy, marker="o", label="Strict overall accuracy")
    primary_axis.set_xlabel("Novelty scale")
    primary_axis.set_ylabel("Fraction")
    primary_axis.set_ylim(-0.02, 1.05)
    primary_axis.grid(True, alpha=0.3)

    secondary_axis = primary_axis.twinx()
    secondary_axis.plot(
        scales,
        unknown_rows,
        marker="s",
        linestyle="--",
        label="Strict rows unresolved",
    )
    secondary_axis.plot(
        scales,
        near_neighbour_overpromotions,
        marker="s",
        linestyle="--",
        label="Near-neighbour over-promotions",
    )
    secondary_axis.set_ylabel("Row count")

    handles_1, labels_1 = primary_axis.get_legend_handles_labels()
    handles_2, labels_2 = secondary_axis.get_legend_handles_labels()
    primary_axis.legend(
        handles_1 + handles_2,
        labels_1 + labels_2,
        loc="lower right",
        frameon=True,
    )
    fig.tight_layout()
    for suffix in ["png", "pdf"]:
        fig.savefig(out_dir / f"{prefix}_operating_curve.{suffix}", dpi=300)
    plt.close(fig)

    fig, axis = plt.subplots(figsize=(6, 5))
    precision = [item.expected_precision for item in summaries]
    recall = [item.expected_recall for item in summaries]
    axis.plot(recall, precision, marker="o")
    for item in summaries:
        axis.annotate(
            format_float(item.novelty_scale),
            (item.expected_recall, item.expected_precision),
            textcoords="offset points",
            xytext=(5, 5),
        )
    axis.set_xlabel("Strict expected-target recall")
    axis.set_ylabel("Strict expected-target precision")
    axis.set_xlim(-0.02, 1.05)
    axis.set_ylim(-0.02, 1.05)
    axis.grid(True, alpha=0.3)
    fig.tight_layout()
    for suffix in ["png", "pdf"]:
        fig.savefig(out_dir / f"{prefix}_expected_pr_points.{suffix}", dpi=300)
    plt.close(fig)

    fig, axis = plt.subplots(figsize=(6, 5))
    fpr = [item.expected_fpr for item in summaries]
    tpr = [item.expected_recall for item in summaries]
    axis.plot(fpr, tpr, marker="o")
    for item in summaries:
        axis.annotate(
            format_float(item.novelty_scale),
            (item.expected_fpr, item.expected_recall),
            textcoords="offset points",
            xytext=(5, 5),
        )
    axis.set_xlabel("Strict expected-target false-positive rate")
    axis.set_ylabel("Strict expected-target true-positive rate")
    axis.set_xlim(-0.02, max(0.05, max(fpr) + 0.02))
    axis.set_ylim(-0.02, 1.05)
    axis.grid(True, alpha=0.3)
    fig.tight_layout()
    for suffix in ["png", "pdf"]:
        fig.savefig(out_dir / f"{prefix}_expected_roc_points.{suffix}", dpi=300)
    plt.close(fig)


def main() -> None:
    """Run the command-line workflow."""
    args = parse_args()
    out_dir = Path(args.out_dir).resolve()
    logger = configure_logging(verbose=args.verbose, out_dir=out_dir)

    sweep_root = Path(args.sweep_root).resolve()
    logger.info("Scanning sweep root: %s", sweep_root)

    prediction_tables = find_prediction_tables(
        sweep_root=sweep_root,
        basename=args.run_glob,
        exclude_run_regex=args.exclude_run_regex,
        logger=logger,
    )
    logger.info("Prediction tables found: %d", len(prediction_tables))

    summaries = []
    for table in prediction_tables:
        logger.info("Summarising %s", table)
        summaries.append(summarise_prediction_table(table))

    summaries = deduplicate_by_scale(summaries=summaries, logger=logger)
    logger.info("Unique novelty scales retained: %d", len(summaries))

    rows = [summary.as_dict() for summary in summaries]
    fieldnames = list(rows[0].keys())
    write_tsv(out_dir / "ai_novelty_sweep_summary.tsv", rows, fieldnames)

    pr_rows = [
        {
            "novelty_scale": format_float(summary.novelty_scale),
            "strict_expected_precision": format_float(summary.expected_precision),
            "strict_expected_recall": format_float(summary.expected_recall),
            "strict_expected_f1": format_float(summary.expected_f1),
            "strict_expected_fpr": format_float(summary.expected_fpr),
            "n_strict_expected_tp": summary.n_strict_expected_tp,
            "n_strict_expected_fp": summary.n_strict_expected_fp,
            "n_strict_expected_fn": summary.n_strict_expected_fn,
        }
        for summary in summaries
    ]
    write_tsv(
        out_dir / "ai_expected_target_pr_roc_operating_points.tsv",
        pr_rows,
        list(pr_rows[0].keys()),
    )

    if not args.no_plots:
        draw_figures(
            summaries=summaries,
            out_dir=out_dir,
            prefix=args.figure_prefix,
            logger=logger,
        )

    logger.info("Wrote outputs to %s", out_dir)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:  # pragma: no cover - command-line safety net
        logging.getLogger("kmersutra_ai_novelty_sweep").exception(
            "Fatal error: %s", exc
        )
        sys.exit(1)
