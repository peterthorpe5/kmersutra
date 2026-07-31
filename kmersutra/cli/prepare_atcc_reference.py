"""Command-line interface for ATCC reference-panel preparation."""

from __future__ import annotations

import argparse
from pathlib import Path

from kmersutra.atcc_reference import (
    create_taxid_plan,
    evaluate_reference_gate,
    finalise_reference_config,
)
from kmersutra.logging_utils import configure_logging
from kmersutra.taxonomy import TaxonomyDatabase


def build_parser() -> argparse.ArgumentParser:
    """Build the command-line parser.

    Returns:
        Configured argument parser.
    """
    parser = argparse.ArgumentParser(
        description="Prepare a leakage-controlled ATCC MSA-1003 reference panel."
    )
    parser.add_argument("--verbose", action="store_true")
    subparsers = parser.add_subparsers(dest="command", required=True)

    plan = subparsers.add_parser("plan", help="Create the NCBI taxid download plan.")
    plan.add_argument("--truth-manifest", required=True)
    plan.add_argument("--taxonomy-dir", required=True)
    plan.add_argument("--output", required=True)
    plan.add_argument("--target-assemblies", type=int, default=5)
    plan.add_argument("--near-neighbour-assemblies-per-genus", type=int, default=25)

    finalise = subparsers.add_parser(
        "finalise",
        help="Combine ATCC downloads with the established background config.",
    )
    finalise.add_argument("--downloaded-config", required=True)
    finalise.add_argument("--background-config", required=True)
    finalise.add_argument("--truth-manifest", required=True)
    finalise.add_argument("--taxonomy-dir", required=True)
    finalise.add_argument("--output-config", required=True)
    finalise.add_argument("--coverage-table", required=True)
    finalise.add_argument("--minimum-target-references", type=int, default=1)

    gate = subparsers.add_parser("gate", help="Evaluate the completed audit.")
    gate.add_argument("--audit-table", required=True)
    gate.add_argument("--truth-manifest", required=True)
    gate.add_argument("--coverage-table", required=True)
    gate.add_argument("--summary-table", required=True)
    gate.add_argument("--minimum-target-references", type=int, default=1)
    return parser


def main(argv: list[str] | None = None) -> int:
    """Run ATCC reference preparation.

    Args:
        argv: Optional explicit command-line arguments.

    Returns:
        Process exit status. A blocked gate returns status 3.
    """
    args = build_parser().parse_args(argv)
    output_hint = Path(
        getattr(args, "output", "")
        or getattr(args, "output_config", "")
        or getattr(args, "summary_table", "")
    ).expanduser()
    log_path = output_hint.parent / "atcc_reference_preparation.log"
    logger = configure_logging(log_file=log_path, verbose=args.verbose)

    if args.command == "plan":
        taxonomy = TaxonomyDatabase.from_taxdump(
            taxonomy_dir=args.taxonomy_dir,
            logger=logger,
        )
        create_taxid_plan(
            truth_manifest=args.truth_manifest,
            taxonomy=taxonomy,
            output_path=args.output,
            target_assemblies=args.target_assemblies,
            near_neighbour_assemblies_per_genus=(
                args.near_neighbour_assemblies_per_genus
            ),
            logger=logger,
        )
        return 0
    if args.command == "finalise":
        taxonomy = TaxonomyDatabase.from_taxdump(
            taxonomy_dir=args.taxonomy_dir,
            logger=logger,
        )
        finalise_reference_config(
            downloaded_config=args.downloaded_config,
            background_config=args.background_config,
            truth_manifest=args.truth_manifest,
            taxonomy=taxonomy,
            output_config=args.output_config,
            coverage_table=args.coverage_table,
            minimum_target_references=args.minimum_target_references,
            logger=logger,
        )
        return 0
    passed = evaluate_reference_gate(
        audit_table=args.audit_table,
        truth_manifest=args.truth_manifest,
        coverage_table=args.coverage_table,
        summary_table=args.summary_table,
        minimum_target_references=args.minimum_target_references,
        logger=logger,
    )
    return 0 if passed else 3


if __name__ == "__main__":
    raise SystemExit(main())
