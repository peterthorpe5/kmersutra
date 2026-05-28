"""Download NCBI taxonomy files for KmerSutra."""

from __future__ import annotations

import argparse
from pathlib import Path

from kmersutra.logging_utils import configure_logging
from kmersutra.taxonomy import NCBI_TAXONOMY_URL, download_taxdump


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed arguments.
    """
    parser = argparse.ArgumentParser(description="Download NCBI taxdmp.zip for KmerSutra.")
    parser.add_argument("--taxonomy_dir", required=True)
    parser.add_argument("--url", default=NCBI_TAXONOMY_URL)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--verbose", action="store_true")
    return parser.parse_args()


def main() -> None:
    """Download and extract NCBI taxonomy files."""
    args = parse_args()
    taxonomy_dir = Path(args.taxonomy_dir)
    taxonomy_dir.mkdir(parents=True, exist_ok=True)
    logger = configure_logging(log_file=taxonomy_dir / "download_taxonomy.log", verbose=args.verbose)
    logger.info("Starting NCBI taxonomy download")
    logger.info("Taxonomy directory: %s", taxonomy_dir)
    logger.info("URL: %s", args.url)
    download_taxdump(
        taxonomy_dir=taxonomy_dir,
        url=args.url,
        overwrite=args.overwrite,
        logger=logger,
    )
    logger.info("Done")


if __name__ == "__main__":
    main()
