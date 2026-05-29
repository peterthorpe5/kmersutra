"""Screen reads or assemblies against a KmerSutra panel."""

from __future__ import annotations

import argparse
from contextlib import nullcontext
from pathlib import Path

from kmersutra.call_consolidation import (
    CONSOLIDATED_CALL_FIELDNAMES,
    consolidate_species_calls,
    merge_background_taxa,
)
from kmersutra.features import FEATURE_FIELDNAMES, summarise_sequence_features
from kmersutra.io import write_tsv
from kmersutra.parquet_modules import OptionalParquetDependencyError
from kmersutra.table_parquet import write_records_parquet
from kmersutra.hierarchical import (
    MODULE_ACTIVATION_FIELDNAMES,
    load_module_manifest,
    screen_file_hierarchical,
)
from kmersutra.lineage_interpretation import (
    LINEAGE_INTERPRETATION_FIELDNAMES,
    interpret_lineage_evidence,
)
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
    parser.add_argument("--panel", default=None)
    parser.add_argument("--sample_id", required=True)
    parser.add_argument("--input_format", choices=["fastq", "fasta"], required=True)
    parser.add_argument(
        "--screen_mode",
        choices=["hierarchical", "flat"],
        default="hierarchical",
        help=(
            "Screening mode. Hierarchical is the default and activates "
            "taxonomic modules from broad gate panels when --module_manifest "
            "is supplied. If no module manifest is supplied, a single --panel "
            "is screened as a compatibility fallback."
        ),
    )
    parser.add_argument(
        "--module_manifest",
        default=None,
        help=(
            "Optional hierarchical module manifest TSV. Required for true "
            "cascade screening; otherwise --panel is used as a single-panel "
            "compatibility fallback."
        ),
    )
    parser.add_argument(
        "--hierarchical_fail_open",
        action="store_true",
        help=(
            "In hierarchical mode, activate detailed modules when weak gate "
            "evidence is observed but no gate passes strict activation "
            "thresholds. This retains sensitivity for unresolved or unsampled "
            "lineage signals."
        ),
    )
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
        choices=["legacy", "conservative", "lineage_aware", "strict"],
        default="legacy",
        help=(
            "Species-call threshold preset. The legacy preset preserves older "
            "behaviour. Conservative and strict presets make species calls "
            "harder to earn. The lineage_aware preset keeps weak neighbouring "
            "species evidence visible but demotes it from reportable species calls."
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
        "--min_mixed_species_fraction",
        type=float,
        default=None,
        help=(
            "Minimum fraction of strongest species-level support required for "
            "a passing species to be promoted to a reportable mixed-species "
            "call. Lower-support passing species are retained as "
            "neighbour_lineage_evidence."
        ),
    )
    parser.add_argument(
        "--low_evidence_call",
        choices=["present_low_confidence", "observed_below_threshold"],
        default=None,
        help="Call label for species evidence observed below reportable thresholds.",
    )
    parser.add_argument(
        "--min_taxonomic_unique_kmers",
        type=int,
        default=20,
        help="Minimum unique k-mers for unresolved lineage evidence.",
    )
    parser.add_argument(
        "--min_taxonomic_positive_sequences",
        type=int,
        default=5,
        help="Minimum independent sequences for unresolved lineage evidence.",
    )
    parser.add_argument(
        "--min_taxonomic_k_values",
        type=int,
        default=1,
        help="Minimum positive k values for unresolved lineage evidence.",
    )
    parser.add_argument(
        "--min_taxonomic_best_k",
        type=int,
        default=77,
        help="Minimum longest k value for unresolved lineage evidence.",
    )
    parser.add_argument(
        "--min_taxonomic_confidence_score",
        type=float,
        default=0.40,
        help="Minimum heuristic confidence score for unresolved lineage evidence.",
    )
    parser.add_argument(
        "--min_neighbour_species_for_novelty",
        type=int,
        default=2,
        help=(
            "Minimum number of weak neighbouring species supporting a possible "
            "novel or unsampled lineage interpretation."
        ),
    )

    parser.add_argument(
        "--consolidate_species_calls",
        action="store_true",
        help=(
            "After raw species thresholding, demote dominated same-genus "
            "species to neighbour_lineage_evidence and separate supplied "
            "empirical-background candidate taxa from ordinary species calls."
        ),
    )
    parser.add_argument(
        "--background_candidate_taxa",
        nargs="*",
        default=[],
        help=(
            "Species labels that should be retained as plausible empirical-"
            "background candidate signals rather than ordinary off-target "
            "species calls, for example 'Hammondia hammondi'."
        ),
    )
    parser.add_argument(
        "--background_candidate_file",
        default=None,
        help="Optional text file containing one background-candidate taxon per line.",
    )
    parser.add_argument(
        "--disable_same_genus_neighbour_demotion",
        action="store_true",
        help=(
            "With --consolidate_species_calls, keep co-reported same-genus "
            "species as reportable species rather than demoting dominated "
            "neighbours."
        ),
    )
    parser.add_argument(
        "--dominant_species_min_margin",
        type=int,
        default=25,
        help=(
            "Minimum unique-k-mer margin required to demote a same-genus "
            "neighbour when --consolidate_species_calls is enabled."
        ),
    )
    parser.add_argument(
        "--dominant_species_min_ratio",
        type=float,
        default=2.0,
        help=(
            "Minimum primary/candidate unique-k-mer ratio required to demote "
            "a same-genus neighbour when --consolidate_species_calls is enabled."
        ),
    )
    parser.add_argument(
        "--write_parquet_outputs",
        action="store_true",
        help=(
            "Also write key screening output tables as Parquet files. Requires "
            "the optional pyarrow dependency."
        ),
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


def _deduplicate_species_metadata(
    *,
    records: list[dict[str, str]],
) -> list[dict[str, str]]:
    """Deduplicate panel species metadata records.

    Parameters
    ----------
    records : list[dict[str, str]]
        Species metadata records.

    Returns
    -------
    list[dict[str, str]]
        Deduplicated records in first-seen order.
    """
    seen: set[tuple[str, str]] = set()
    output: list[dict[str, str]] = []
    for record in records:
        key = (record.get("species_name", ""), record.get("clade", ""))
        if key in seen:
            continue
        seen.add(key)
        output.append(record)
    return output


def load_expected_species_for_screening(
    *,
    panel_path: str | None,
    module_manifest_path: str | None,
) -> list[dict[str, str]]:
    """Load expected species metadata from a panel or module manifest.

    Parameters
    ----------
    panel_path : str or None
        Single flat panel path.
    module_manifest_path : str or None
        Hierarchical module manifest path.

    Returns
    -------
    list[dict[str, str]]
        Deduplicated species metadata records.
    """
    records: list[dict[str, str]] = []
    if module_manifest_path:
        modules = load_module_manifest(manifest_path=module_manifest_path)
        panel_paths = []
        for module in modules:
            if module.module_panel_path:
                panel_paths.append(module.module_panel_path)
        for module_panel in sorted(set(panel_paths)):
            records.extend(load_panel_species_metadata(panel_path=module_panel))
    elif panel_path:
        records.extend(load_panel_species_metadata(panel_path=panel_path))
    return _deduplicate_species_metadata(records=records)



def maybe_write_parquet_table(
    *,
    records: list[dict[str, object]],
    output_path: Path,
    fieldnames: list[str],
    enabled: bool,
    logger,
) -> None:
    """Write an optional Parquet companion table.

    Parameters
    ----------
    records : list of dict
        Records to write.
    output_path : pathlib.Path
        Parquet output path.
    fieldnames : list of str
        Stable output column order.
    enabled : bool
        Whether Parquet output was requested.
    logger : logging.Logger
        Logger used for diagnostics.
    """
    if not enabled:
        return
    try:
        n_written = write_records_parquet(
            records=records,
            output_path=output_path,
            fieldnames=fieldnames,
        )
    except OptionalParquetDependencyError as exc:
        logger.warning(
            "Could not write Parquet output %s because pyarrow is unavailable: %s",
            output_path,
            exc,
        )
        return
    logger.info("Wrote Parquet output %s (%d rows)", output_path, n_written)

def main() -> None:
    """Run the sequence screening workflow."""
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    logger = configure_logging(log_file=out_dir / "screen_reads.log", verbose=args.verbose)
    logger.info("Starting KmerSutra screening")
    logger.info("Input: %s", args.input)
    logger.info("Screen mode: %s", args.screen_mode)
    logger.info("Panel: %s", args.panel or "not supplied")
    logger.info("Module manifest: %s", args.module_manifest or "not supplied")
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
    if args.screen_mode == "flat" and not args.panel:
        raise ValueError("--panel is required when --screen_mode flat is used")
    if args.screen_mode == "hierarchical" and not args.module_manifest and not args.panel:
        raise ValueError(
            "--module_manifest or --panel is required for hierarchical screening"
        )
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
        min_mixed_species_fraction=args.min_mixed_species_fraction,
    )
    logger.info("Species-call preset: %s", args.call_preset)
    for setting_name, setting_value in sorted(call_settings.items()):
        logger.info("Species-call setting %s: %s", setting_name, setting_value)
    logger.info("Minimum unique-kmer margin: %d", args.min_unique_kmer_margin)
    logger.info("Minimum unique-kmer ratio: %.4f", args.min_unique_kmer_ratio)
    logger.info(
        "Minimum mixed-species support fraction: %.4f",
        float(call_settings["min_mixed_species_fraction"]),
    )

    logger.info("Minimum taxonomic unique k-mers: %d", args.min_taxonomic_unique_kmers)
    logger.info(
        "Minimum taxonomic positive sequences: %d",
        args.min_taxonomic_positive_sequences,
    )
    logger.info("Minimum taxonomic k values: %d", args.min_taxonomic_k_values)
    logger.info("Minimum taxonomic best k: %d", args.min_taxonomic_best_k)
    logger.info(
        "Minimum taxonomic confidence score: %.4f",
        args.min_taxonomic_confidence_score,
    )
    logger.info(
        "Minimum weak neighbour species for possible novelty: %d",
        args.min_neighbour_species_for_novelty,
    )
    logger.info("Consolidate species calls: %s", args.consolidate_species_calls)
    logger.info(
        "Same-genus neighbour demotion enabled: %s",
        not args.disable_same_genus_neighbour_demotion,
    )
    logger.info(
        "Dominant species demotion margin: %d; ratio: %.4f",
        args.dominant_species_min_margin,
        args.dominant_species_min_ratio,
    )
    logger.info("Write Parquet outputs: %s", args.write_parquet_outputs)

    profiler = WorkflowProfiler() if args.profile else None
    screen_profile_records: list[dict[str, object]] = []

    with (profiler.time_stage(stage="load_expected_species", detail="panel_metadata") if profiler else nullcontext()):
        expected_species = load_expected_species_for_screening(
            panel_path=args.panel,
            module_manifest_path=args.module_manifest,
        )
    logger.info("Loaded %d expected species labels from panel metadata", len(expected_species))

    module_activation_records = []
    if args.screen_mode == "hierarchical" and args.module_manifest:
        logger.info("Running hierarchical cascade screening")
        hierarchical_result = screen_file_hierarchical(
            input_path=args.input,
            module_manifest_path=args.module_manifest,
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
            hierarchical_fail_open=args.hierarchical_fail_open,
            profile_records=screen_profile_records if profiler else None,
            logger=logger,
        )
        hits = hierarchical_result.hits
        module_activation_records = hierarchical_result.activation_records
    else:
        if args.screen_mode == "hierarchical":
            logger.info(
                "No module manifest supplied; using single-panel compatibility fallback"
            )
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
            min_mixed_species_fraction=float(
                call_settings["min_mixed_species_fraction"]
            ),
        )
    logger.info("Built %d raw species detection-call rows", len(detection_calls))

    raw_detection_calls = [dict(record) for record in detection_calls]
    if args.consolidate_species_calls:
        background_taxa = merge_background_taxa(
            background_candidate_taxa=args.background_candidate_taxa,
            background_candidate_file=args.background_candidate_file,
        )
        logger.info(
            "Background candidate taxa supplied for consolidation: %s",
            "; ".join(sorted(background_taxa)) or "none",
        )
        detection_calls = consolidate_species_calls(
            species_calls=raw_detection_calls,
            background_candidate_taxa=background_taxa,
            demote_same_genus_neighbours=not args.disable_same_genus_neighbour_demotion,
            dominant_species_min_margin=args.dominant_species_min_margin,
            dominant_species_min_ratio=args.dominant_species_min_ratio,
            logger=logger,
        )
    logger.info("Built %d report-layer species detection-call rows", len(detection_calls))

    with (profiler.time_stage(stage="interpret_lineage_evidence", detail="unresolved_novelty") if profiler else nullcontext()):
        lineage_interpretation = interpret_lineage_evidence(
            species_calls=detection_calls,
            taxonomic_evidence=taxonomic_evidence,
            min_taxonomic_unique_kmers=args.min_taxonomic_unique_kmers,
            min_taxonomic_positive_sequences=args.min_taxonomic_positive_sequences,
            min_taxonomic_k_values=args.min_taxonomic_k_values,
            min_taxonomic_best_k=args.min_taxonomic_best_k,
            min_taxonomic_confidence_score=args.min_taxonomic_confidence_score,
            min_neighbour_species_for_novelty=args.min_neighbour_species_for_novelty,
        )
    logger.info(
        "Built %d sample lineage-interpretation rows",
        len(lineage_interpretation),
    )

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
    maybe_write_parquet_table(
        records=taxonomic_evidence,
        output_path=out_dir / "sample_taxonomic_kmer_evidence.parquet",
        fieldnames=TAXONOMIC_EVIDENCE_FIELDNAMES,
        enabled=args.write_parquet_outputs,
        logger=logger,
    )
    raw_call_fieldnames = [
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
        "reportable_conflicting_unique_kmers",
        "reportable_conflict_ratio",
        "mixed_species_support_fraction",
        "confidence_score",
        "signal_confidence_score",
        "call",
    ]
    if args.consolidate_species_calls:
        write_tsv(
            records=raw_detection_calls,
            output_path=out_dir / "species_detection_calls_raw.tsv",
            fieldnames=raw_call_fieldnames,
        )
        maybe_write_parquet_table(
            records=raw_detection_calls,
            output_path=out_dir / "species_detection_calls_raw.parquet",
            fieldnames=raw_call_fieldnames,
            enabled=args.write_parquet_outputs,
            logger=logger,
        )
        call_fieldnames = CONSOLIDATED_CALL_FIELDNAMES
    else:
        call_fieldnames = raw_call_fieldnames
    write_tsv(
        records=detection_calls,
        output_path=out_dir / "species_detection_calls.tsv",
        fieldnames=call_fieldnames,
    )
    maybe_write_parquet_table(
        records=detection_calls,
        output_path=out_dir / "species_detection_calls.parquet",
        fieldnames=call_fieldnames,
        enabled=args.write_parquet_outputs,
        logger=logger,
    )
    write_tsv(
        records=lineage_interpretation,
        output_path=out_dir / "sample_lineage_interpretation.tsv",
        fieldnames=LINEAGE_INTERPRETATION_FIELDNAMES,
    )
    maybe_write_parquet_table(
        records=lineage_interpretation,
        output_path=out_dir / "sample_lineage_interpretation.parquet",
        fieldnames=LINEAGE_INTERPRETATION_FIELDNAMES,
        enabled=args.write_parquet_outputs,
        logger=logger,
    )
    write_tsv(
        records=module_activation_records,
        output_path=out_dir / "module_activation.tsv",
        fieldnames=MODULE_ACTIVATION_FIELDNAMES,
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
