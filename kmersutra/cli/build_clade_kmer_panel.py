"""Build a KmerSutra clade/species diagnostic k-mer panel."""

from __future__ import annotations

import argparse
from pathlib import Path

from kmersutra.build_panel import build_panel
from kmersutra.config import load_genome_config
from kmersutra.io import write_json, write_tsv
from kmersutra.logging_utils import configure_logging
from kmersutra.profiling import WorkflowProfiler
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
        Parsed arguments.
    """
    parser = argparse.ArgumentParser(
        description="Build an outgroup-aware KmerSutra diagnostic k-mer panel."
    )
    parser.add_argument("--genome_config", required=True)
    parser.add_argument("--out_dir", required=True)
    parser.add_argument("--k_values", nargs="+", type=int, default=[31, 51, 71, 101])
    parser.add_argument("--target_clade", default="")
    parser.add_argument("--target_taxid", default="")
    parser.add_argument("--taxonomy_dir", default="")
    parser.add_argument("--download_taxonomy_if_missing", action="store_true")
    parser.add_argument(
        "--taxonomy_url",
        default="https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip",
    )
    parser.add_argument(
        "--evidence_ranks",
        nargs="+",
        default=CORE_RANK_ORDER,
        help="Taxonomic ranks retained as evidence levels when taxonomy is used.",
    )
    parser.add_argument("--max_per_species_per_k", type=int, default=None)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument(
        "--compact_build",
        dest="compact_build",
        action="store_true",
        default=True,
        help="Use compact set-based k-mer grouping. This is the default and is recommended for larger databases.",
    )
    parser.add_argument(
        "--legacy_observation_build",
        dest="compact_build",
        action="store_false",
        help="Use the older occurrence-level builder for debugging/regression checks.",
    )
    parser.add_argument(
        "--profile",
        action="store_true",
        help="Write build_profile_timing.tsv with wall-clock timings for build stages.",
    )
    parser.add_argument("--verbose", action="store_true")
    return parser.parse_args()


def main() -> None:
    """Run the panel builder."""
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    logger = configure_logging(log_file=out_dir / "build_panel.log", verbose=args.verbose)
    logger.info("Starting KmerSutra panel build")
    logger.info("Genome config: %s", args.genome_config)
    logger.info("k values: %s", args.k_values)
    logger.info("Worker processes: %d", args.threads)
    logger.info("Compact build: %s", args.compact_build)
    logger.info("Build profiling: %s", args.profile)
    profiler = WorkflowProfiler()

    taxonomy_db = None
    if args.taxonomy_dir:
        logger.info("Loading NCBI taxonomy from %s", args.taxonomy_dir)
        taxonomy_db = TaxonomyDatabase.from_taxdump(
            taxonomy_dir=args.taxonomy_dir,
            download_if_missing=args.download_taxonomy_if_missing,
            url=args.taxonomy_url,
            logger=logger,
        )
        logger.info("Taxonomic evidence ranks: %s", ", ".join(args.evidence_ranks))

    with profiler.time_stage(stage="load_genome_config", detail=str(args.genome_config)):
        genome_configs = load_genome_config(config_path=args.genome_config)
    logger.info("Loaded %d genome records", len(genome_configs))

    with profiler.time_stage(
        stage="build_panel",
        detail=f"compact_build={args.compact_build};k_values={','.join(map(str, args.k_values))}",
    ):
        diagnostic_kmers, summary, collection_summary = build_panel(
            genome_configs=genome_configs,
            k_values=args.k_values,
            target_clade=args.target_clade,
            max_per_species_per_k=args.max_per_species_per_k,
            threads=args.threads,
            taxonomy_db=taxonomy_db,
            target_taxid=args.target_taxid,
            preferred_ranks=args.evidence_ranks,
            compact_build=args.compact_build,
            logger=logger,
        )
    logger.info("Retained %d diagnostic k-mers", len(diagnostic_kmers))

    panel_path = out_dir / "species_kmer_panel.tsv.gz"
    summary_path = out_dir / "kmer_uniqueness_summary.tsv"
    collection_summary_path = out_dir / "kmer_collection_summary.tsv"
    metadata_path = out_dir / "species_kmer_panel_metadata.json"
    html_path = out_dir / "species_detection_report.html"

    with profiler.time_stage(stage="write_panel", detail=str(panel_path)):
        write_tsv(
            records=[item.to_record() for item in diagnostic_kmers],
            output_path=panel_path,
            fieldnames=PANEL_FIELDNAMES,
        )
    with profiler.time_stage(stage="write_summaries", detail=str(summary_path)):
        write_tsv(
            records=summary,
            output_path=summary_path,
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
        collection_fieldnames = sorted({key for record in collection_summary for key in record})
        write_tsv(
            records=collection_summary,
            output_path=collection_summary_path,
            fieldnames=collection_fieldnames,
        )
    write_json(
        data={
            "genome_config": str(args.genome_config),
            "k_values": args.k_values,
            "target_clade": args.target_clade,
            "target_taxid": args.target_taxid,
            "taxonomy_dir": args.taxonomy_dir,
            "download_taxonomy_if_missing": args.download_taxonomy_if_missing,
            "evidence_ranks": args.evidence_ranks,
            "max_per_species_per_k": args.max_per_species_per_k,
            "threads": args.threads,
            "compact_build": args.compact_build,
            "profile": args.profile,
            "n_genomes": len(genome_configs),
            "n_diagnostic_kmers": len(diagnostic_kmers),
        },
        output_path=metadata_path,
    )
    if args.profile:
        write_tsv(
            records=profiler.to_records(),
            output_path=out_dir / "build_profile_timing.tsv",
            fieldnames=["stage", "seconds", "detail"],
        )
    write_html_report(
        output_path=html_path,
        title="KmerSutra panel build report",
        panel_summary=summary,
    )
    logger.info("Panel written to %s", panel_path)
    logger.info("Done")


if __name__ == "__main__":
    main()
