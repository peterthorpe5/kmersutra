"""Command-line interface for merging Parquet KmerSutra modules."""

from __future__ import annotations

import argparse
from pathlib import Path

from kmersutra.cli.build_clade_kmer_panel import _write_panel_streaming
from kmersutra.global_candidate_evidence import (
    assign_global_candidate_evidence_sqlite,
    iter_retained_global_candidate_diagnostics,
    summarise_retained_global_candidate_evidence,
)
from kmersutra.io import write_json, write_tsv
from kmersutra.logging_utils import configure_logging
from kmersutra.parquet_modules import (
    import_global_kmer_parquets_to_sqlite,
    export_global_candidate_module,
    resolve_global_kmer_parquet_paths,
)
from kmersutra.reporting import write_html_report
from kmersutra.taxonomy import CORE_RANK_ORDER, TaxonomyDatabase


PANEL_FIELDNAMES = [
    "kmer",
    "k",
    "panel_type",
    "species_name",
    "clade",
    "source_genomes",
    "source_contigs",
    "example_position",
    "evidence_taxid",
    "evidence_name",
    "evidence_rank",
    "lineage_taxids",
    "source_taxids",
]


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Merge KmerSutra Parquet module source tables and globally "
            "revalidate k-mer evidence across all modules."
        )
    )
    parser.add_argument(
        "--module_dirs",
        nargs="*",
        default=[],
        help="Module directories containing global_kmers.parquet.",
    )
    parser.add_argument(
        "--global_kmer_parquets",
        nargs="*",
        default=[],
        help="Explicit global_kmers.parquet files to merge.",
    )
    parser.add_argument("--out_dir", required=True, help="Output directory.")
    parser.add_argument("--taxonomy_dir", required=True, help="NCBI taxonomy dump directory.")
    parser.add_argument("--download_taxonomy_if_missing", action="store_true")
    parser.add_argument(
        "--evidence_ranks",
        nargs="+",
        default=CORE_RANK_ORDER,
        help="Taxonomic ranks that may be retained as evidence levels.",
    )
    parser.add_argument(
        "--target_taxid",
        default="",
        help="Optional subtree restriction. Leave empty for query-agnostic panels.",
    )
    parser.add_argument(
        "--candidate_roles",
        nargs="*",
        default=None,
        help="Optional whitelist of reportable source roles.",
    )
    parser.add_argument(
        "--excluded_candidate_roles",
        nargs="*",
        default=None,
        help="Optional roles excluded from reportable candidate status.",
    )
    parser.add_argument(
        "--sqlite_path",
        default="",
        help="Optional merged SQLite path. Defaults to merged_global_candidate_evidence.sqlite in --out_dir.",
    )
    parser.add_argument("--sqlite_batch_size", type=int, default=50000)
    parser.add_argument("--max_per_species_per_k", type=int, default=None)
    parser.add_argument(
        "--marker_selection",
        choices=["first_seen", "genome_spread", "independent_multik_genome_spread"],
        default="independent_multik_genome_spread",
        help="Marker-selection strategy for capped merged panels.",
    )
    parser.add_argument("--genome_bin_size", type=int, default=10000)
    parser.add_argument("--max_per_genome_bin", type=int, default=10)
    parser.add_argument("--min_cross_k_marker_distance", type=int, default=5000)
    parser.add_argument(
        "--write_merged_module_parquet",
        action="store_true",
        help="Export the merged source index and retained evidence as a new module Parquet directory.",
    )
    parser.add_argument(
        "--merged_module_name",
        default="merged_module",
        help="Module name recorded if --write_merged_module_parquet is set.",
    )
    parser.add_argument("--verbose", action="store_true")
    return parser.parse_args()


def main() -> None:
    """Run Parquet module merging and global revalidation."""
    args = parse_args()
    logger = configure_logging(verbose=args.verbose)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    module_parquets = resolve_global_kmer_parquet_paths(
        module_dirs=args.module_dirs,
        global_kmer_parquets=args.global_kmer_parquets,
    )
    logger.info("Merging %d global-kmer module table(s)", len(module_parquets))

    taxonomy_db = TaxonomyDatabase.from_taxdump(
        taxonomy_dir=args.taxonomy_dir,
        download_if_missing=args.download_taxonomy_if_missing,
        logger=logger,
    )

    sqlite_path = Path(args.sqlite_path) if args.sqlite_path else out_dir / "merged_global_candidate_evidence.sqlite"
    import_summary = import_global_kmer_parquets_to_sqlite(
        parquet_paths=module_parquets,
        sqlite_path=sqlite_path,
        batch_size=args.sqlite_batch_size,
        logger=logger,
    )
    write_tsv(
        records=import_summary,
        output_path=out_dir / "module_import_summary.tsv",
        fieldnames=["module_parquet", "imported_records"],
    )

    assignment_summary = assign_global_candidate_evidence_sqlite(
        sqlite_path=sqlite_path,
        taxonomy_db=taxonomy_db,
        preferred_ranks=args.evidence_ranks,
        target_taxid=args.target_taxid,
        candidate_roles=args.candidate_roles,
        excluded_candidate_roles=args.excluded_candidate_roles,
        batch_size=args.sqlite_batch_size,
        max_per_evidence_per_k=(
            None
            if args.marker_selection != "first_seen"
            else args.max_per_species_per_k
        ),
        logger=logger,
    )
    write_tsv(
        records=assignment_summary,
        output_path=out_dir / "global_revalidation_summary.tsv",
        fieldnames=["stage", "metric", "value"],
    )

    panel_path = out_dir / "species_kmer_panel.tsv.gz"
    n_kmers, panel_summary = _write_panel_streaming(
        diagnostic_kmers=iter_retained_global_candidate_diagnostics(sqlite_path=sqlite_path),
        panel_path=panel_path,
        max_per_species_per_k=(
            args.max_per_species_per_k
            if args.marker_selection != "first_seen"
            else None
        ),
        logger=logger,
        marker_selection=args.marker_selection,
        genome_bin_size=args.genome_bin_size,
        max_per_genome_bin=args.max_per_genome_bin,
        min_cross_k_marker_distance=args.min_cross_k_marker_distance,
    )
    write_tsv(
        records=panel_summary,
        output_path=out_dir / "kmer_uniqueness_summary.tsv",
        fieldnames=[
            "panel_type",
            "species_name",
            "clade",
            "evidence_taxid",
            "evidence_rank",
            "k",
            "diagnostic_kmers",
        ],
    )

    retained_summary = summarise_retained_global_candidate_evidence(sqlite_path=sqlite_path)
    write_tsv(
        records=retained_summary,
        output_path=out_dir / "panel_evidence_rank_summary.tsv",
        fieldnames=["evidence_rank", "evidence_name", "species_name", "k", "n_kmers"],
    )

    if args.write_merged_module_parquet:
        export_global_candidate_module(
            sqlite_path=sqlite_path,
            module_dir=out_dir / "merged_module_parquet",
            module_name=args.merged_module_name,
            metadata={
                "source_module_count": len(module_parquets),
                "marker_selection": args.marker_selection,
                "genome_bin_size": args.genome_bin_size,
                "max_per_genome_bin": args.max_per_genome_bin,
                "min_cross_k_marker_distance": args.min_cross_k_marker_distance,
                "max_per_species_per_k": args.max_per_species_per_k,
                "evidence_ranks": ";".join(args.evidence_ranks),
            },
            batch_size=args.sqlite_batch_size,
            logger=logger,
        )

    write_json(
        data={
            "module_parquets": [str(path) for path in module_parquets],
            "sqlite_path": str(sqlite_path),
            "panel_path": str(panel_path),
            "n_diagnostic_kmers": n_kmers,
            "marker_selection": args.marker_selection,
            "genome_bin_size": args.genome_bin_size,
            "max_per_genome_bin": args.max_per_genome_bin,
                "min_cross_k_marker_distance": args.min_cross_k_marker_distance,
            "max_per_species_per_k": args.max_per_species_per_k,
            "evidence_ranks": args.evidence_ranks,
        },
        output_path=out_dir / "module_merge_metadata.json",
    )
    write_html_report(
        output_path=out_dir / "species_detection_report.html",
        title="KmerSutra module merge and global revalidation report",
        panel_summary=panel_summary,
    )
    logger.info("Merged panel written to %s", panel_path)
    logger.info("Retained %d diagnostic k-mers", n_kmers)


if __name__ == "__main__":
    main()
