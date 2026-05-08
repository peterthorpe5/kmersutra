"""Build a KmerSutra clade/species diagnostic k-mer panel."""

from __future__ import annotations

import argparse
from pathlib import Path

from kmersutra.build_panel import DiagnosticKmer, build_panel
from kmersutra.config import load_genome_config
from kmersutra.io import write_json, write_tsv
from kmersutra.logging_utils import configure_logging
from kmersutra.profiling import WorkflowProfiler
from kmersutra.reporting import write_html_report
from kmersutra.resource_monitor import ResourceMonitor
from kmersutra.target_evidence import (
    build_target_evidence_sqlite,
    iter_target_evidence_diagnostics,
)
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
    parser.add_argument(
        "--target_evidence_only",
        action="store_true",
        help=(
            "Use the SQLite-backed target-evidence builder. This low-memory "
            "mode stores candidate k-mers from genomes marked target_species "
            "on disk and streams all other genomes only as filters or "
            "downgrade evidence. It is recommended for larger near-neighbour "
            "and outgroup panels when the aim is to test named targets."
        ),
    )
    parser.add_argument(
        "--sqlite_path",
        default="",
        help=(
            "Optional path for the SQLite candidate database used by "
            "--target_evidence_only. Defaults to target_evidence_candidates.sqlite "
            "inside --out_dir."
        ),
    )
    parser.add_argument(
        "--sqlite_batch_size",
        type=int,
        default=50000,
        help="Commit interval for the SQLite-backed target-evidence builder.",
    )
    parser.add_argument(
        "--ram_log_interval_seconds",
        type=float,
        default=60.0,
        help=(
            "Interval for writing RAM measurements to ram_usage.tsv. Use 0 to "
            "disable the background RAM monitor."
        ),
    )
    parser.add_argument(
        "--ram_log_path",
        default="",
        help="Optional RAM log path. Defaults to ram_usage.tsv inside --out_dir.",
    )
    parser.add_argument("--verbose", action="store_true")
    return parser.parse_args()



def _diagnostic_retention_key(item: DiagnosticKmer) -> tuple[str, str, str, str, int]:
    """Return the throttling key for one diagnostic k-mer.

    Parameters
    ----------
    item : DiagnosticKmer
        Diagnostic k-mer record.

    Returns
    -------
    tuple[str, str, str, str, int]
        Retention key based on evidence type, taxon and k value.
    """
    return (
        item.panel_type,
        item.evidence_taxid,
        item.species_name,
        item.clade,
        item.k,
    )


def _summary_key(item: DiagnosticKmer) -> tuple[str, str, str, str, str, int]:
    """Return the panel-summary key for one diagnostic k-mer.

    Parameters
    ----------
    item : DiagnosticKmer
        Diagnostic k-mer record.

    Returns
    -------
    tuple[str, str, str, str, str, int]
        Summary key.
    """
    return (
        item.panel_type,
        item.species_name,
        item.clade,
        item.evidence_taxid,
        item.evidence_rank,
        item.k,
    )


def _write_panel_streaming(
    *,
    diagnostic_kmers,
    panel_path: Path,
    max_per_species_per_k: int | None,
    logger,
) -> tuple[int, list[dict[str, object]]]:
    """Write diagnostic k-mers to a panel file while summarising counts.

    Parameters
    ----------
    diagnostic_kmers : iterable
        Iterable of DiagnosticKmer records.
    panel_path : pathlib.Path
        Output TSV or TSV.GZ path.
    max_per_species_per_k : int | None
        Optional maximum number of retained records per evidence bucket and k.
    logger : logging.Logger
        Logger for progress messages.

    Returns
    -------
    tuple[int, list[dict[str, object]]]
        Number of retained diagnostic k-mers and panel summary rows.
    """
    from collections import defaultdict

    from kmersutra.io import open_text

    if max_per_species_per_k is not None and max_per_species_per_k <= 0:
        raise ValueError("max_per_species_per_k must be positive")

    retained_by_key = defaultdict(int)
    summary_counts = defaultdict(int)
    n_written = 0

    with open_text(panel_path, "wt") as handle:
        handle.write("\t".join(PANEL_FIELDNAMES) + "\n")
        for item in diagnostic_kmers:
            retention_key = _diagnostic_retention_key(item)
            if max_per_species_per_k is not None:
                if retained_by_key[retention_key] >= max_per_species_per_k:
                    continue
            retained_by_key[retention_key] += 1
            summary_counts[_summary_key(item)] += 1
            record = item.to_record()
            handle.write(
                "\t".join(str(record.get(column, "")) for column in PANEL_FIELDNAMES)
                + "\n"
            )
            n_written += 1
            if n_written % 100000 == 0:
                logger.info("Wrote %d diagnostic k-mers", n_written)

    summary = [
        {
            "panel_type": panel_type,
            "species_name": species_name,
            "clade": clade,
            "evidence_taxid": evidence_taxid,
            "evidence_rank": evidence_rank,
            "k": k,
            "diagnostic_kmers": count,
        }
        for (panel_type, species_name, clade, evidence_taxid, evidence_rank, k), count
        in sorted(summary_counts.items())
    ]
    return n_written, summary


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
    logger.info("Target-evidence-only build: %s", args.target_evidence_only)
    logger.info("Build profiling: %s", args.profile)
    profiler = WorkflowProfiler()

    ram_log_path = Path(args.ram_log_path) if args.ram_log_path else out_dir / "ram_usage.tsv"
    monitor = None
    if args.ram_log_interval_seconds > 0:
        monitor = ResourceMonitor(
            output_path=ram_log_path,
            interval_seconds=args.ram_log_interval_seconds,
            logger=logger,
        )
        monitor.start()

    try:
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

        panel_path = out_dir / "species_kmer_panel.tsv.gz"
        summary_path = out_dir / "kmer_uniqueness_summary.tsv"
        collection_summary_path = out_dir / "kmer_collection_summary.tsv"
        metadata_path = out_dir / "species_kmer_panel_metadata.json"
        html_path = out_dir / "species_detection_report.html"
        target_evidence_summary_path = out_dir / "target_evidence_build_summary.tsv"

        if args.target_evidence_only:
            sqlite_path = (
                Path(args.sqlite_path)
                if args.sqlite_path
                else out_dir / "target_evidence_candidates.sqlite"
            )
            with profiler.time_stage(
                stage="build_target_evidence_sqlite",
                detail=f"sqlite_path={sqlite_path};k_values={','.join(map(str, args.k_values))}",
            ):
                target_result = build_target_evidence_sqlite(
                    genome_configs=genome_configs,
                    k_values=args.k_values,
                    sqlite_path=sqlite_path,
                    batch_size=args.sqlite_batch_size,
                    logger=logger,
                )
            with profiler.time_stage(stage="write_panel", detail=str(panel_path)):
                diagnostic_iter = iter_target_evidence_diagnostics(
                    sqlite_path=target_result.sqlite_path,
                    taxonomy_db=taxonomy_db,
                    target_taxid=args.target_taxid,
                    preferred_ranks=args.evidence_ranks,
                )
                n_diagnostic_kmers, summary = _write_panel_streaming(
                    diagnostic_kmers=diagnostic_iter,
                    panel_path=panel_path,
                    max_per_species_per_k=args.max_per_species_per_k,
                    logger=logger,
                )
            collection_summary = target_result.collection_summary
            write_tsv(
                records=target_result.build_summary,
                output_path=target_evidence_summary_path,
                fieldnames=["summary_name", "summary_value"],
            )
        else:
            with profiler.time_stage(
                stage="build_panel",
                detail=(
                    f"compact_build={args.compact_build};"
                    f"k_values={','.join(map(str, args.k_values))}"
                ),
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
            n_diagnostic_kmers = len(diagnostic_kmers)
            with profiler.time_stage(stage="write_panel", detail=str(panel_path)):
                write_tsv(
                    records=[item.to_record() for item in diagnostic_kmers],
                    output_path=panel_path,
                    fieldnames=PANEL_FIELDNAMES,
                )

        logger.info("Retained %d diagnostic k-mers", n_diagnostic_kmers)

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
                "target_evidence_only": args.target_evidence_only,
                "sqlite_path": args.sqlite_path,
                "sqlite_batch_size": args.sqlite_batch_size,
                "profile": args.profile,
                "ram_log_path": str(ram_log_path) if args.ram_log_interval_seconds > 0 else "",
                "ram_log_interval_seconds": args.ram_log_interval_seconds,
                "n_genomes": len(genome_configs),
                "n_diagnostic_kmers": n_diagnostic_kmers,
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
    finally:
        if monitor is not None:
            monitor.stop()


if __name__ == "__main__":
    main()
