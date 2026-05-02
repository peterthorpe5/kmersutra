"""Screen reads or assemblies against a KmerSutra panel."""

from __future__ import annotations

import argparse
from pathlib import Path

from kmersutra.features import FEATURE_FIELDNAMES, summarise_sequence_features
from kmersutra.io import write_tsv
from kmersutra.logging_utils import configure_logging
from kmersutra.reporting import write_html_report
from kmersutra.screen_reads import screen_file_for_species_kmers
from kmersutra.summarise_hits import (
    SPECIES_EVIDENCE_FIELDNAMES,
    complete_sample_species_evidence,
    load_panel_species_metadata,
    summarise_sample_species_evidence,
    summarise_species_hits,
)
from kmersutra.thresholds import call_species_presence


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed arguments.
    """
    parser = argparse.ArgumentParser(
        description="Screen FASTQ/FASTA sequences against a KmerSutra k-mer panel."
    )
    parser.add_argument("--input", required=True)
    parser.add_argument("--panel", required=True)
    parser.add_argument("--sample_id", required=True)
    parser.add_argument("--input_format", choices=["fastq", "fasta"], required=True)
    parser.add_argument("--out_dir", required=True)
    parser.add_argument("--max_mismatches", type=int, default=0)
    parser.add_argument("--fuzzy_min_k", type=int, default=71)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--chunk_size", type=int, default=1000)
    parser.add_argument("--min_unique_kmers", type=int, default=3)
    parser.add_argument("--min_positive_sequences", type=int, default=2)
    parser.add_argument("--min_k_values_positive", type=int, default=1)
    parser.add_argument("--max_conflict_ratio", type=float, default=0.10)
    parser.add_argument(
        "--disallow_mixed_species",
        action="store_true",
        help=(
            "Treat multiple species passing evidence thresholds as conflicting "
            "rather than as a true mixed-species sample."
        ),
    )
    parser.add_argument("--verbose", action="store_true")
    return parser.parse_args()


def main() -> None:
    """Run the sequence screening workflow."""
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    logger = configure_logging(log_file=out_dir / "screen_reads.log", verbose=args.verbose)
    logger.info("Starting KmerSutra screening")
    logger.info("Input: %s", args.input)
    logger.info("Panel: %s", args.panel)
    logger.info("Worker processes: %d", args.threads)
    logger.info("Chunk size: %d", args.chunk_size)
    logger.info("Maximum mismatches: %d", args.max_mismatches)
    logger.info("Fuzzy minimum k: %d", args.fuzzy_min_k)
    logger.info("Mixed-species calls allowed: %s", not args.disallow_mixed_species)

    expected_species = load_panel_species_metadata(panel_path=args.panel)
    logger.info("Loaded %d expected species labels from panel", len(expected_species))

    hits = screen_file_for_species_kmers(
        input_path=args.input,
        panel_path=args.panel,
        sample_id=args.sample_id,
        input_format=args.input_format,
        max_mismatches=args.max_mismatches,
        fuzzy_min_k=args.fuzzy_min_k,
        threads=args.threads,
        chunk_size=args.chunk_size,
        logger=logger,
    )
    logger.info("Detected %d diagnostic k-mer hits", len(hits))

    sequence_features = summarise_sequence_features(hits=hits, logger=logger)
    hit_summary = summarise_species_hits(hits=hits)
    observed_evidence = summarise_sample_species_evidence(species_summary=hit_summary)
    evidence = complete_sample_species_evidence(
        evidence_records=observed_evidence,
        expected_species=expected_species,
        sample_id=args.sample_id,
    )
    logger.info("Built %d completed species evidence rows", len(evidence))
    detection_calls = call_species_presence(
        evidence_records=evidence,
        min_unique_kmers=args.min_unique_kmers,
        min_positive_sequences=args.min_positive_sequences,
        min_k_values_positive=args.min_k_values_positive,
        max_conflict_ratio=args.max_conflict_ratio,
        allow_mixed_species=not args.disallow_mixed_species,
    )
    logger.info("Built %d species detection-call rows", len(detection_calls))

    write_tsv(
        records=[hit.to_record() for hit in hits],
        output_path=out_dir / "read_level_species_kmer_hits.tsv.gz",
        fieldnames=[
            "sample_id",
            "sequence_id",
            "sequence_type",
            "k",
            "query_position",
            "matched_kmer",
            "query_kmer",
            "mismatches",
            "panel_type",
            "species_name",
            "clade",
        ],
    )
    write_tsv(
        records=sequence_features,
        output_path=out_dir / "sequence_ml_features.tsv",
        fieldnames=FEATURE_FIELDNAMES,
    )
    write_tsv(
        records=hit_summary,
        output_path=out_dir / "sample_species_kmer_hits.tsv",
        fieldnames=[
            "sample_id",
            "panel_type",
            "label",
            "clade",
            "k",
            "n_hits",
            "n_unique_kmers",
            "n_positive_sequences",
            "n_exact_hits",
            "n_fuzzy_hits",
            "min_mismatches",
            "max_mismatches",
        ],
    )
    write_tsv(
        records=evidence,
        output_path=out_dir / "sample_species_kmer_evidence.tsv",
        fieldnames=SPECIES_EVIDENCE_FIELDNAMES,
    )
    write_tsv(
        records=detection_calls,
        output_path=out_dir / "species_detection_calls.tsv",
        fieldnames=[
            "sample_id",
            "species_name",
            "clade",
            "n_hits",
            "n_unique_kmers",
            "n_positive_sequences",
            "n_k_values_positive",
            "best_k",
            "n_exact_hits",
            "n_fuzzy_hits",
            "conflicting_unique_kmers",
            "conflict_ratio",
            "confidence_score",
            "call",
        ],
    )
    write_html_report(
        output_path=out_dir / "species_detection_report.html",
        title="KmerSutra species detection report",
        hit_summary=hit_summary,
        detection_calls=detection_calls,
    )
    logger.info("Done")


if __name__ == "__main__":
    main()
