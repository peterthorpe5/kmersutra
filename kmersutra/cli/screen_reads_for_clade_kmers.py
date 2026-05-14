"""Screen reads or assemblies against a KmerSutra panel."""

from __future__ import annotations

import argparse
from contextlib import nullcontext
from pathlib import Path

from kmersutra.features import FEATURE_FIELDNAMES, summarise_sequence_features
from kmersutra.io import write_tsv
from kmersutra.profiling import WorkflowProfiler
from kmersutra.logging_utils import configure_logging
from kmersutra.reporting import write_html_report
from kmersutra.screen_reads import screen_file_for_species_kmers
from kmersutra.summarise_hits import (
    SPECIES_EVIDENCE_FIELDNAMES,
    TAXONOMIC_EVIDENCE_FIELDNAMES,
    complete_sample_species_evidence,
    load_panel_species_metadata,
    summarise_sample_species_evidence,
    summarise_sample_taxonomic_evidence,
    summarise_species_hits,
    summarise_taxonomic_hits,
)
from kmersutra.thresholds import apply_species_call_preset, call_species_presence


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
    parser.add_argument("--chunk_size", type=int, default=5000)
    parser.add_argument(
        "--max_pending_chunks",
        type=int,
        default=None,
        help="Maximum queued chunks during threaded screening. Default: twice --threads.",
    )
    parser.add_argument(
        "--panel_cache",
        default=None,
        help="Optional pickled panel-index cache path.",
    )
    parser.add_argument(
        "--use_panel_cache",
        action="store_true",
        help="Use a current panel-index cache if available; otherwise build one.",
    )
    parser.add_argument(
        "--write_panel_cache",
        action="store_true",
        help="Write a pickled panel-index cache after loading the TSV panel.",
    )
    parser.add_argument(
        "--no_read_level_hits",
        action="store_true",
        help="Do not write read_level_species_kmer_hits.tsv.gz. This speeds up large screens.",
    )
    parser.add_argument(
        "--profile",
        action="store_true",
        help="Write profile_timing.tsv with wall-clock timings for major stages.",
    )
    parser.add_argument(
        "--call_preset",
        choices=["legacy", "conservative", "strict"],
        default="legacy",
        help=(
            "Species-call threshold preset. The legacy preset preserves older "
            "behaviour. Conservative and strict presets make species calls "
            "harder to earn and label weak evidence as observed_below_threshold."
        ),
    )
    parser.add_argument("--min_unique_kmers", type=int, default=None)
    parser.add_argument("--min_positive_sequences", type=int, default=None)
    parser.add_argument("--min_k_values_positive", type=int, default=None)
    parser.add_argument("--max_conflict_ratio", type=float, default=None)
    parser.add_argument(
        "--min_best_k",
        type=int,
        default=None,
        help="Minimum longest k-mer length required for a reportable species call.",
    )
    parser.add_argument(
        "--min_exact_hits",
        type=int,
        default=None,
        help="Minimum exact-hit count required for a reportable species call.",
    )
    parser.add_argument(
        "--min_confidence_score",
        type=float,
        default=None,
        help="Minimum heuristic confidence score for a reportable species call.",
    )
    parser.add_argument(
        "--min_unique_kmer_margin",
        type=int,
        default=0,
        help="Optional unique-k-mer margin over the second-best species.",
    )
    parser.add_argument(
        "--min_unique_kmer_ratio",
        type=float,
        default=0.0,
        help="Optional focal-to-second-best unique-k-mer ratio requirement.",
    )
    parser.add_argument(
        "--low_evidence_call",
        choices=["present_low_confidence", "observed_below_threshold"],
        default=None,
        help="Call label for species evidence observed below reportable thresholds.",
    )
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
    logger.info("Maximum pending chunks: %s", args.max_pending_chunks or "auto")
    logger.info("Panel cache: %s", args.panel_cache or "default")
    logger.info("Use panel cache: %s", args.use_panel_cache)
    logger.info("Write panel cache: %s", args.write_panel_cache)
    logger.info("Write read-level hits: %s", not args.no_read_level_hits)
    logger.info("Profiling enabled: %s", args.profile)
    logger.info("Maximum mismatches: %d", args.max_mismatches)
    logger.info("Fuzzy minimum k: %d", args.fuzzy_min_k)
    logger.info("Mixed-species calls allowed: %s", not args.disallow_mixed_species)
    call_settings = apply_species_call_preset(
        preset_name=args.call_preset,
        min_unique_kmers=args.min_unique_kmers,
        min_positive_sequences=args.min_positive_sequences,
        min_k_values_positive=args.min_k_values_positive,
        max_conflict_ratio=args.max_conflict_ratio,
        min_best_k=args.min_best_k,
        min_exact_hits=args.min_exact_hits,
        min_confidence_score=args.min_confidence_score,
        low_evidence_call=args.low_evidence_call,
    )
    logger.info("Species-call preset: %s", args.call_preset)
    for setting_name, setting_value in sorted(call_settings.items()):
        logger.info("Species-call setting %s: %s", setting_name, setting_value)
    logger.info("Minimum unique-kmer margin: %d", args.min_unique_kmer_margin)
    logger.info("Minimum unique-kmer ratio: %.4f", args.min_unique_kmer_ratio)

    profiler = WorkflowProfiler() if args.profile else None
    screen_profile_records: list[dict[str, object]] = []

    with (profiler.time_stage(stage="load_expected_species", detail="panel_metadata") if profiler else nullcontext()):
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
        max_pending_chunks=args.max_pending_chunks,
        panel_cache_path=args.panel_cache,
        use_panel_cache=args.use_panel_cache,
        write_panel_cache=args.write_panel_cache,
        profile_records=screen_profile_records if profiler else None,
        logger=logger,
    )
    logger.info("Detected %d diagnostic k-mer hits", len(hits))

    with (profiler.time_stage(stage="summarise_sequence_features", detail="ml_features") if profiler else nullcontext()):
        sequence_features = summarise_sequence_features(hits=hits, logger=logger)
    with (profiler.time_stage(stage="summarise_hits", detail="sample_species") if profiler else nullcontext()):
        hit_summary = summarise_species_hits(hits=hits)
        observed_evidence = summarise_sample_species_evidence(species_summary=hit_summary)
        evidence = complete_sample_species_evidence(
            evidence_records=observed_evidence,
            expected_species=expected_species,
            sample_id=args.sample_id,
        )
        taxonomic_hit_summary = summarise_taxonomic_hits(hits=hits)
        taxonomic_evidence = summarise_sample_taxonomic_evidence(
            taxonomic_summary=taxonomic_hit_summary,
        )
    logger.info("Built %d completed species evidence rows", len(evidence))
    logger.info("Built %d taxonomic evidence rows", len(taxonomic_evidence))
    with (profiler.time_stage(stage="call_species_presence", detail="thresholds") if profiler else nullcontext()):
        detection_calls = call_species_presence(
            evidence_records=evidence,
            min_unique_kmers=int(call_settings["min_unique_kmers"]),
            min_positive_sequences=int(call_settings["min_positive_sequences"]),
            min_k_values_positive=int(call_settings["min_k_values_positive"]),
            max_conflict_ratio=float(call_settings["max_conflict_ratio"]),
            allow_mixed_species=not args.disallow_mixed_species,
            min_best_k=int(call_settings["min_best_k"]),
            min_exact_hits=int(call_settings["min_exact_hits"]),
            min_confidence_score=float(call_settings["min_confidence_score"]),
            min_unique_kmer_margin=args.min_unique_kmer_margin,
            min_unique_kmer_ratio=args.min_unique_kmer_ratio,
            low_evidence_call=str(call_settings["low_evidence_call"]),
        )
    logger.info("Built %d species detection-call rows", len(detection_calls))

    read_level_fieldnames = [
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
        "evidence_taxid",
        "evidence_name",
        "evidence_rank",
    ]
    if args.no_read_level_hits:
        logger.info("Skipping read-level hit output by user request")
        write_tsv(
            records=[{
                "sample_id": args.sample_id,
                "note": "read-level hit output disabled; use without --no_read_level_hits to write this file",
            }],
            output_path=out_dir / "read_level_species_kmer_hits.disabled.tsv",
            fieldnames=["sample_id", "note"],
        )
    else:
        with (profiler.time_stage(stage="write_read_level_hits", detail=f"hits={len(hits)}") if profiler else nullcontext()):
            write_tsv(
                records=[hit.to_record() for hit in hits],
                output_path=out_dir / "read_level_species_kmer_hits.tsv.gz",
                fieldnames=read_level_fieldnames,
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
        records=taxonomic_evidence,
        output_path=out_dir / "sample_taxonomic_kmer_evidence.tsv",
        fieldnames=TAXONOMIC_EVIDENCE_FIELDNAMES,
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
    with (profiler.time_stage(stage="write_html_report", detail="species_detection_report") if profiler else nullcontext()):
        write_html_report(
            output_path=out_dir / "species_detection_report.html",
            title="KmerSutra species detection report",
            hit_summary=hit_summary,
            detection_calls=detection_calls,
        )
    if profiler:
        write_tsv(
            records=[*screen_profile_records, *profiler.to_records()],
            output_path=out_dir / "profile_timing.tsv",
            fieldnames=["stage", "seconds", "detail"],
        )
        logger.info("Wrote profile timings: %s", out_dir / "profile_timing.tsv")
    logger.info("Done")


if __name__ == "__main__":
    main()
