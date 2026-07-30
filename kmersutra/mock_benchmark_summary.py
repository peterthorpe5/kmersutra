"""Summaries for manifest-driven mock-community benchmarks."""

from __future__ import annotations

import logging
from collections import defaultdict
from pathlib import Path

from kmersutra.mock_community import (
    MockCommunityTruth,
    load_truth_manifest,
    normalise_species_key,
)
from kmersutra.table_io import read_records_table, write_records_table
from kmersutra.zymo_truth import infer_species_name, row_has_evidence, row_is_reportable

EXPECTED_RESULT_FIELDS = [
    "sample_id",
    "analysis_type",
    "fraction",
    "seed",
    "k_value",
    "benchmark_id",
    "truth_organism_id",
    "species_name",
    "expected_abundance_fraction",
    "abundance_tier",
    "raw_evidence_detected",
    "species_reportable",
    "n_matching_call_rows",
    "max_unique_kmers",
    "max_positive_sequences",
    "best_k",
]

OFF_TARGET_FIELDS = [
    "sample_id",
    "analysis_type",
    "fraction",
    "seed",
    "k_value",
    "species_name",
    "raw_evidence_detected",
    "species_reportable",
    "n_unique_kmers",
    "n_positive_sequences",
    "best_k",
    "call",
]

SAMPLE_SUMMARY_FIELDS = [
    "sample_id",
    "analysis_type",
    "fraction",
    "seed",
    "k_value",
    "n_expected_species",
    "n_expected_raw_detected",
    "n_expected_reportable",
    "expected_raw_recall",
    "expected_reportable_recall",
    "n_unexpected_raw_evidence_species",
    "n_reportable_off_target_species",
    "reportable_precision",
    "strictly_clean",
]

TIER_SUMMARY_FIELDS = [
    "analysis_type",
    "fraction",
    "k_value",
    "abundance_tier",
    "n_sample_species_opportunities",
    "n_raw_detected",
    "n_reportable",
    "raw_recall",
    "reportable_recall",
]


def safe_int(value: object) -> int:
    """Parse an integer-like table value.

    Args:
        value: Input value.

    Returns:
        Parsed integer, or zero if unavailable.
    """
    try:
        return int(float(str(value or "0")))
    except (TypeError, ValueError):
        return 0


def _task_metadata(task: dict[str, str]) -> dict[str, str]:
    """Return stable metadata fields from a task manifest row."""
    return {
        "sample_id": str(task.get("sample_id", "")),
        "analysis_type": str(task.get("analysis_type", "")),
        "fraction": str(task.get("fraction", "")),
        "seed": str(task.get("seed", "")),
        "k_value": str(task.get("k_value", "")),
    }


def _call_rows_by_species(
    call_rows: list[dict[str, str]],
) -> dict[str, list[dict[str, str]]]:
    """Group species-call rows by normalised species name."""
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in call_rows:
        key = normalise_species_key(infer_species_name(record=row))
        if key:
            grouped[key].append(row)
    return grouped


def _matching_rows(
    *,
    truth: MockCommunityTruth,
    grouped_calls: dict[str, list[dict[str, str]]],
) -> list[dict[str, str]]:
    """Collect call rows matching any accepted truth species name."""
    rows: list[dict[str, str]] = []
    for name in truth.accepted_species_names:
        rows.extend(grouped_calls.get(name, []))
    return rows


def summarise_one_sample(
    *,
    task: dict[str, str],
    truth_records: list[MockCommunityTruth],
    call_rows: list[dict[str, str]],
) -> tuple[
    list[dict[str, object]],
    list[dict[str, object]],
    dict[str, object],
]:
    """Summarise expected and unexpected species for one benchmark sample.

    Args:
        task: Task manifest row.
        truth_records: Validated truth records.
        call_rows: KmerSutra species-call rows.

    Returns:
        Expected-species records, off-target records and sample summary.
    """
    metadata = _task_metadata(task)
    grouped_calls = _call_rows_by_species(call_rows)
    truth_names = {
        name for truth in truth_records for name in truth.accepted_species_names
    }
    expected_results: list[dict[str, object]] = []
    for truth in truth_records:
        if not truth.expected:
            continue
        matches = _matching_rows(truth=truth, grouped_calls=grouped_calls)
        raw_detected = any(row_has_evidence(record=row) for row in matches)
        reportable = any(row_is_reportable(record=row) for row in matches)
        expected_results.append(
            {
                **metadata,
                "benchmark_id": truth.benchmark_id,
                "truth_organism_id": truth.organism_id,
                "species_name": truth.species_name,
                "expected_abundance_fraction": (
                    f"{truth.expected_abundance_fraction:.8f}"
                ),
                "abundance_tier": truth.abundance_tier,
                "raw_evidence_detected": str(raw_detected),
                "species_reportable": str(reportable),
                "n_matching_call_rows": len(matches),
                "max_unique_kmers": max(
                    (safe_int(row.get("n_unique_kmers")) for row in matches),
                    default=0,
                ),
                "max_positive_sequences": max(
                    (
                        safe_int(row.get("n_positive_sequences"))
                        for row in matches
                    ),
                    default=0,
                ),
                "best_k": max(
                    (safe_int(row.get("best_k")) for row in matches),
                    default=0,
                ),
            }
        )

    off_targets: list[dict[str, object]] = []
    for species_key, rows in sorted(grouped_calls.items()):
        if species_key in truth_names:
            continue
        evidence_rows = [row for row in rows if row_has_evidence(record=row)]
        if not evidence_rows:
            continue
        representative = max(
            evidence_rows,
            key=lambda row: (
                safe_int(row.get("n_unique_kmers")),
                safe_int(row.get("n_positive_sequences")),
                safe_int(row.get("best_k")),
            ),
        )
        off_targets.append(
            {
                **metadata,
                "species_name": infer_species_name(record=representative),
                "raw_evidence_detected": "True",
                "species_reportable": str(
                    any(row_is_reportable(record=row) for row in rows)
                ),
                "n_unique_kmers": safe_int(
                    representative.get("n_unique_kmers")
                ),
                "n_positive_sequences": safe_int(
                    representative.get("n_positive_sequences")
                ),
                "best_k": safe_int(representative.get("best_k")),
                "call": str(
                    representative.get(
                        "call",
                        representative.get("call_status", ""),
                    )
                ),
            }
        )

    n_expected = len(expected_results)
    n_expected_raw = sum(
        row["raw_evidence_detected"] == "True" for row in expected_results
    )
    n_expected_reportable = sum(
        row["species_reportable"] == "True" for row in expected_results
    )
    n_off_target_raw = len(off_targets)
    n_off_target_reportable = sum(
        row["species_reportable"] == "True" for row in off_targets
    )
    reportable_total = n_expected_reportable + n_off_target_reportable
    sample_summary = {
        **metadata,
        "n_expected_species": n_expected,
        "n_expected_raw_detected": n_expected_raw,
        "n_expected_reportable": n_expected_reportable,
        "expected_raw_recall": (
            f"{n_expected_raw / n_expected:.8f}" if n_expected else ""
        ),
        "expected_reportable_recall": (
            f"{n_expected_reportable / n_expected:.8f}" if n_expected else ""
        ),
        "n_unexpected_raw_evidence_species": n_off_target_raw,
        "n_reportable_off_target_species": n_off_target_reportable,
        "reportable_precision": (
            f"{n_expected_reportable / reportable_total:.8f}"
            if reportable_total
            else ""
        ),
        "strictly_clean": str(n_off_target_reportable == 0),
    }
    return expected_results, off_targets, sample_summary


def summarise_tiers(
    expected_results: list[dict[str, object]],
) -> list[dict[str, object]]:
    """Aggregate expected-species outcomes by abundance tier.

    Args:
        expected_results: Expected-species result rows.

    Returns:
        Tier-level recall records.
    """
    groups: dict[tuple[str, str, str, str], list[dict[str, object]]] = defaultdict(
        list
    )
    for row in expected_results:
        key = (
            str(row.get("analysis_type", "")),
            str(row.get("fraction", "")),
            str(row.get("k_value", "")),
            str(row.get("abundance_tier", "")),
        )
        groups[key].append(row)
    output: list[dict[str, object]] = []
    for key, rows in sorted(groups.items()):
        n_rows = len(rows)
        n_raw = sum(row["raw_evidence_detected"] == "True" for row in rows)
        n_reportable = sum(
            row["species_reportable"] == "True" for row in rows
        )
        output.append(
            {
                "analysis_type": key[0],
                "fraction": key[1],
                "k_value": key[2],
                "abundance_tier": key[3],
                "n_sample_species_opportunities": n_rows,
                "n_raw_detected": n_raw,
                "n_reportable": n_reportable,
                "raw_recall": f"{n_raw / n_rows:.8f}",
                "reportable_recall": f"{n_reportable / n_rows:.8f}",
            }
        )
    return output


def summarise_mock_benchmark(
    *,
    task_manifest: str | Path,
    truth_manifest: str | Path,
    output_dir: str | Path,
    logger: logging.Logger | None = None,
) -> dict[str, Path]:
    """Summarise every completed task in a mock-community benchmark.

    Args:
        task_manifest: Table listing sample metadata and calls-table paths.
        truth_manifest: Mock-community truth manifest.
        output_dir: Summary output directory.
        logger: Optional logger.

    Returns:
        Named output paths.

    Raises:
        ValueError: If no task rows can be summarised.
        FileNotFoundError: If a declared calls table is missing.
    """
    tasks = read_records_table(
        input_path=task_manifest,
        required_columns=["sample_id", "analysis_type", "calls_table"],
        logger=logger,
    )
    if not tasks:
        raise ValueError(f"Task manifest contains no records: {task_manifest}")
    truth_records = load_truth_manifest(
        manifest_path=truth_manifest,
        logger=logger,
    )
    expected_results: list[dict[str, object]] = []
    off_targets: list[dict[str, object]] = []
    sample_summaries: list[dict[str, object]] = []
    for task in tasks:
        calls_path = Path(str(task.get("calls_table", ""))).expanduser()
        if not calls_path.is_file() or calls_path.stat().st_size <= 0:
            raise FileNotFoundError(
                f"Calls table is missing or empty for {task.get('sample_id')}: "
                f"{calls_path}"
            )
        call_rows = read_records_table(input_path=calls_path, logger=logger)
        expected, unexpected, sample = summarise_one_sample(
            task=task,
            truth_records=truth_records,
            call_rows=call_rows,
        )
        expected_results.extend(expected)
        off_targets.extend(unexpected)
        sample_summaries.append(sample)

    destination = Path(output_dir).expanduser().resolve()
    destination.mkdir(parents=True, exist_ok=True)
    paths = {
        "expected_species_results": destination
        / "expected_species_results.tsv.gz",
        "off_target_results": destination / "off_target_results.tsv.gz",
        "sample_summary": destination / "sample_summary.tsv",
        "tier_summary": destination / "abundance_tier_summary.tsv",
    }
    write_records_table(
        records=expected_results,
        output_path=paths["expected_species_results"],
        fieldnames=EXPECTED_RESULT_FIELDS,
        logger=logger,
    )
    write_records_table(
        records=off_targets,
        output_path=paths["off_target_results"],
        fieldnames=OFF_TARGET_FIELDS,
        logger=logger,
    )
    write_records_table(
        records=sample_summaries,
        output_path=paths["sample_summary"],
        fieldnames=SAMPLE_SUMMARY_FIELDS,
        logger=logger,
    )
    write_records_table(
        records=summarise_tiers(expected_results),
        output_path=paths["tier_summary"],
        fieldnames=TIER_SUMMARY_FIELDS,
        logger=logger,
    )
    return paths
