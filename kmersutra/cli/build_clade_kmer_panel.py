"""Build a KmerSutra clade/species diagnostic k-mer panel."""

from __future__ import annotations

import argparse
from pathlib import Path

from kmersutra.build_panel import build_panel
from kmersutra.config import load_genome_config
from kmersutra.io import write_json, write_tsv
from kmersutra.logging_utils import configure_logging
from kmersutra.reporting import write_html_report


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
    parser.add_argument("--max_per_species_per_k", type=int, default=None)
    parser.add_argument("--threads", type=int, default=1)
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

    genome_configs = load_genome_config(config_path=args.genome_config)
    logger.info("Loaded %d genome records", len(genome_configs))

    diagnostic_kmers, summary, collection_summary = build_panel(
        genome_configs=genome_configs,
        k_values=args.k_values,
        target_clade=args.target_clade,
        max_per_species_per_k=args.max_per_species_per_k,
        threads=args.threads,
        logger=logger,
    )
    logger.info("Retained %d diagnostic k-mers", len(diagnostic_kmers))

    panel_path = out_dir / "species_kmer_panel.tsv.gz"
    summary_path = out_dir / "kmer_uniqueness_summary.tsv"
    collection_summary_path = out_dir / "kmer_collection_summary.tsv"
    metadata_path = out_dir / "species_kmer_panel_metadata.json"
    html_path = out_dir / "species_detection_report.html"

    write_tsv(
        records=[item.to_record() for item in diagnostic_kmers],
        output_path=panel_path,
        fieldnames=[
            "kmer",
            "k",
            "panel_type",
            "species_name",
            "clade",
            "source_genomes",
            "source_contigs",
            "example_position",
        ],
    )
    write_tsv(
        records=summary,
        output_path=summary_path,
        fieldnames=["panel_type", "species_name", "clade", "k", "diagnostic_kmers"],
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
            "max_per_species_per_k": args.max_per_species_per_k,
            "threads": args.threads,
            "n_genomes": len(genome_configs),
            "n_diagnostic_kmers": len(diagnostic_kmers),
        },
        output_path=metadata_path,
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
