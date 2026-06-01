"""Command-line interface for KmerSutra LCA reporting."""

from __future__ import annotations

import argparse
from pathlib import Path

from kmersutra.lca_reporting import summarise_lca_table
from kmersutra.logging_utils import configure_logging


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Summarise KmerSutra evidence at lowest-common-ancestor levels. "
            "Supported table suffixes are .tsv, .tsv.gz and .parquet."
        )
    )
    parser.add_argument(
        "--evidence_table",
        required=True,
        help="Input KmerSutra detection/evidence table.",
    )
    parser.add_argument(
        "--out_table",
        required=True,
        help="Output LCA report table.",
    )
    parser.add_argument(
        "--taxonomy_dir",
        required=True,
        help="Directory containing NCBI nodes.dmp, names.dmp, merged.dmp and delnodes.dmp.",
    )
    parser.add_argument(
        "--taxon_map_table",
        default=None,
        help=(
            "Optional name-to-taxid mapping table. Required when the evidence "
            "table does not already contain taxid columns."
        ),
    )
    parser.add_argument(
        "--sample_column",
        default="sample_id",
        help="Sample identifier column in the input evidence table.",
    )
    parser.add_argument(
        "--taxid_column",
        default=None,
        help="Optional explicit taxid column in the evidence table.",
    )
    parser.add_argument(
        "--taxon_name_column",
        default=None,
        help="Optional explicit taxon-name column in the evidence table.",
    )
    parser.add_argument(
        "--taxon_map_name_column",
        default=None,
        help="Optional explicit taxon-name column in the taxon map table.",
    )
    parser.add_argument(
        "--taxon_map_taxid_column",
        default=None,
        help="Optional explicit taxid column in the taxon map table.",
    )
    parser.add_argument("--min_unique_kmers", type=int, default=1)
    parser.add_argument("--min_positive_sequences", type=int, default=1)
    parser.add_argument("--min_best_k", type=int, default=0)
    parser.add_argument(
        "--dominant_min_score_fraction",
        type=float,
        default=0.25,
        help="Minimum score fraction for non-top taxa in dominant-lineage LCA.",
    )
    parser.add_argument(
        "--scopes",
        nargs="+",
        default=[
            "dominant_lineage",
            "all_supported_evidence",
            "background_candidate",
            "neighbour_lineage",
            "reportable_positive",
        ],
        help="LCA scopes to write.",
    )
    parser.add_argument(
        "--download_taxonomy_if_missing",
        action="store_true",
        help="Download NCBI taxdump if required files are missing.",
    )
    parser.add_argument("--verbose", action="store_true")
    return parser.parse_args()


def main() -> None:
    """Run KmerSutra LCA reporting."""
    args = parse_args()
    output_path = Path(args.out_table)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    logger = configure_logging(
        log_file=output_path.with_suffix(".log"),
        verbose=args.verbose,
    )
    logger.info("Starting KmerSutra LCA reporting")
    logger.info("Evidence table: %s", args.evidence_table)
    logger.info("Output table: %s", args.out_table)
    logger.info("Taxonomy directory: %s", args.taxonomy_dir)
    records = summarise_lca_table(
        evidence_table=args.evidence_table,
        output_table=args.out_table,
        taxonomy_dir=args.taxonomy_dir,
        taxon_map_table=args.taxon_map_table,
        download_taxonomy_if_missing=args.download_taxonomy_if_missing,
        sample_column=args.sample_column,
        taxid_column=args.taxid_column,
        taxon_name_column=args.taxon_name_column,
        taxon_map_name_column=args.taxon_map_name_column,
        taxon_map_taxid_column=args.taxon_map_taxid_column,
        min_unique_kmers=args.min_unique_kmers,
        min_positive_sequences=args.min_positive_sequences,
        min_best_k=args.min_best_k,
        dominant_min_score_fraction=args.dominant_min_score_fraction,
        scopes=args.scopes,
        logger=logger,
    )
    logger.info("Wrote %d LCA report rows", len(records))
    logger.info("Done")


if __name__ == "__main__":
    main()
