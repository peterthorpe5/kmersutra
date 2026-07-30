"""Audit a KmerSutra reference panel for exact benchmark-truth leakage."""

from __future__ import annotations

import argparse
from pathlib import Path

from kmersutra.logging_utils import configure_logging
from kmersutra.reference_audit import write_reference_panel_audit


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments.

    Args:
        argv: Optional explicit argument list.

    Returns:
        Parsed arguments.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Audit a KmerSutra genome configuration for exact mock-community "
            "truth accessions and FASTA identity."
        )
    )
    parser.add_argument("--genome_config", required=True)
    parser.add_argument("--truth_manifest", required=True)
    parser.add_argument("--out_table", required=True)
    parser.add_argument("--allow_missing_fastas", action="store_true")
    parser.add_argument("--fail_on_leakage", action="store_true")
    parser.add_argument("--verbose", action="store_true")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    """Run the reference-panel leakage audit.

    Args:
        argv: Optional explicit argument list.

    Returns:
        Process exit status.
    """
    args = parse_args(argv)
    output_path = Path(args.out_table).expanduser().resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    logger = configure_logging(
        log_file=output_path.with_suffix(".log"),
        verbose=args.verbose,
    )
    rows = write_reference_panel_audit(
        genome_config=args.genome_config,
        truth_manifest=args.truth_manifest,
        output_table=output_path,
        allow_missing_fastas=args.allow_missing_fastas,
        fail_on_leakage=args.fail_on_leakage,
        logger=logger,
    )
    logger.info("Wrote %d reference-audit rows", len(rows))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
