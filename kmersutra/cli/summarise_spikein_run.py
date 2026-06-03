"""Build KmerSutra run-level summaries.

This entry point preserves the older ``--summary_tsv`` report-formatting mode
and now also dispatches comparable benchmark arguments such as ``--out_root`` to
``kmersutra.comparable_benchmark_summary``. The dual behaviour keeps existing
workflows working while allowing users to run comparable summaries through the
installed ``kmersutra-summarise-run`` command.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from kmersutra.comparable_benchmark_summary import main as comparable_main
from kmersutra.logging_utils import configure_logging
from kmersutra.run_summary import build_run_summary_reports


COMPARABLE_MODE_FLAGS = {
    "--out_root",
    "--manifest",
    "--out_dir",
    "--panel1_targets",
    "--panel2_tsv",
    "--panel3_tsv",
    "--background_candidate_taxa",
    "--demote_expected_genus_neighbours",
}


def parse_args() -> argparse.Namespace:
    """Parse arguments for the legacy summary formatter.

    Returns
    -------
    argparse.Namespace
        Parsed legacy summary arguments.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Build Excel and HTML summaries from a KmerSutra spike-in summary "
            "TSV. For comparable benchmark runs, pass --out_root and related "
            "arguments; this command will dispatch to the comparable summary "
            "workflow."
        )
    )
    parser.add_argument("--summary_tsv", required=True)
    parser.add_argument("--out_xlsx", required=True)
    parser.add_argument("--out_html", required=True)
    parser.add_argument("--verbose", action="store_true")
    return parser.parse_args()


def should_use_comparable_mode(*, argv: list[str]) -> bool:
    """Return whether arguments request comparable benchmark summary mode.

    Parameters
    ----------
    argv : list[str]
        Command-line arguments excluding the executable name.

    Returns
    -------
    bool
        True if comparable mode should be used.
    """
    return any(argument in COMPARABLE_MODE_FLAGS for argument in argv)


def run_legacy_summary() -> None:
    """Run the legacy summary-TSV-to-Excel/HTML workflow."""
    args = parse_args()
    out_xlsx = Path(args.out_xlsx)
    out_xlsx.parent.mkdir(parents=True, exist_ok=True)
    logger = configure_logging(log_file=out_xlsx.with_suffix(".log"), verbose=args.verbose)
    logger.info("Starting KmerSutra run summary")
    logger.info("Summary TSV: %s", args.summary_tsv)
    logger.info("Output Excel: %s", args.out_xlsx)
    logger.info("Output HTML: %s", args.out_html)
    build_run_summary_reports(
        summary_tsv=args.summary_tsv,
        out_xlsx=args.out_xlsx,
        out_html=args.out_html,
        logger=logger,
    )
    logger.info("Done")


def main() -> None:
    """Run the appropriate KmerSutra summary workflow."""
    argv = sys.argv[1:]
    if should_use_comparable_mode(argv=argv):
        comparable_main()
        return
    run_legacy_summary()


if __name__ == "__main__":  # pragma: no cover
    main()
