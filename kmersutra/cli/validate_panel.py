"""Command-line interface for validating a KmerSutra panel."""

from __future__ import annotations

import argparse

from kmersutra.logging_utils import configure_logging
from kmersutra.panel_merge import validate_panel_file


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(description="Validate a KmerSutra panel file.")
    parser.add_argument("--panel", required=True, help="Panel TSV/TSV.GZ file to validate.")
    parser.add_argument("--out_dir", required=True, help="Output directory for validation tables.")
    parser.add_argument("--verbose", action="store_true", help="Enable detailed logging.")
    return parser.parse_args()


def main() -> None:
    """Run the panel validation command."""
    args = parse_args()
    logger = configure_logging(verbose=args.verbose)
    validate_panel_file(panel_path=args.panel, out_dir=args.out_dir, logger=logger)
    logger.info("Panel validation complete: %s", args.out_dir)


if __name__ == "__main__":
    main()
