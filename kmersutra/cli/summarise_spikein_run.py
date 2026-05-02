"""Build run-level KmerSutra spike-in summary reports."""

from __future__ import annotations

import argparse
from pathlib import Path

from kmersutra.logging_utils import configure_logging
from kmersutra.run_summary import build_run_summary_reports


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed arguments.
    """
    parser = argparse.ArgumentParser(
        description="Build Excel and HTML summaries from a KmerSutra spike-in summary TSV."
    )
    parser.add_argument("--summary_tsv", required=True)
    parser.add_argument("--out_xlsx", required=True)
    parser.add_argument("--out_html", required=True)
    parser.add_argument("--verbose", action="store_true")
    return parser.parse_args()


def main() -> None:
    """Run the KmerSutra spike-in summary workflow."""
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


if __name__ == "__main__":
    main()
