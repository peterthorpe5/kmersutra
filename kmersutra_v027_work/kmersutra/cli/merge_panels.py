"""Command-line interface for merging KmerSutra panels."""

from __future__ import annotations

import argparse

from kmersutra.logging_utils import configure_logging
from kmersutra.panel_merge import merge_panel_files
from kmersutra.taxonomy import CORE_RANK_ORDER, TaxonomyDatabase


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Merge independently built KmerSutra panels into a validated master panel."
    )
    parser.add_argument("--panels", nargs="+", required=True, help="Input panel TSV/TSV.GZ files.")
    parser.add_argument("--out_dir", required=True, help="Output directory for the master panel.")
    parser.add_argument("--taxonomy_dir", default=None, help="Optional NCBI taxonomy dump directory.")
    parser.add_argument(
        "--download_taxonomy_if_missing",
        action="store_true",
        help="Download NCBI taxdump files if they are missing from --taxonomy_dir.",
    )
    parser.add_argument(
        "--evidence_ranks",
        nargs="+",
        default=CORE_RANK_ORDER,
        help="Taxonomic ranks that may be retained as evidence levels.",
    )
    parser.add_argument("--verbose", action="store_true", help="Enable detailed logging.")
    return parser.parse_args()


def main() -> None:
    """Run the panel merge command."""
    args = parse_args()
    logger = configure_logging(verbose=args.verbose)
    taxonomy_db = None
    if args.taxonomy_dir:
        logger.info("Loading taxonomy from %s", args.taxonomy_dir)
        taxonomy_db = TaxonomyDatabase.from_taxdump(
            taxonomy_dir=args.taxonomy_dir,
            download_if_missing=args.download_taxonomy_if_missing,
            logger=logger,
        )
    merge_panel_files(
        panel_paths=args.panels,
        out_dir=args.out_dir,
        taxonomy_db=taxonomy_db,
        preferred_ranks=args.evidence_ranks,
        logger=logger,
    )
    logger.info("Panel merge complete: %s", args.out_dir)


if __name__ == "__main__":
    main()
