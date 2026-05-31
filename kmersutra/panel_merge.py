"""Merge and validate KmerSutra k-mer panels.

This module supports the scalable KmerSutra workflow where panels are built
separately for different taxonomic modules and then merged into a globally
validated master panel. During merging, shared k-mers are either retained at
an appropriate taxonomic evidence level, downgraded to a broader level, or
removed when they are not diagnostically useful.
"""

from __future__ import annotations

import logging
from collections import Counter, defaultdict
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path

from kmersutra.io import read_tsv, write_json, write_tsv
from kmersutra.taxonomy import CORE_RANK_ORDER, TaxonomyDatabase

PANEL_REQUIRED_COLUMNS = [
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

PANEL_OUTPUT_COLUMNS = PANEL_REQUIRED_COLUMNS + ["source_panels", "n_panel_records"]

VALIDATION_ISSUE_COLUMNS = [
    "issue_type",
    "severity",
    "k",
    "kmer",
    "n_records",
    "details",
]


@dataclass(frozen=True)
class PanelMergeResult:
    """Container for merged panel outputs.

    Attributes
    ----------
    master_records : list[dict[str, object]]
        Records retained in the validated master panel.
    removed_records : list[dict[str, object]]
        Records removed because no useful global evidence level was found.
    downgraded_records : list[dict[str, object]]
        Records whose evidence rank became broader during merging.
    validation_rows : list[dict[str, object]]
        Summary rows describing the merge process.
    taxonomic_level_rows : list[dict[str, object]]
        Counts by k, panel type, evidence rank and evidence name.
    """

    master_records: list[dict[str, object]]
    removed_records: list[dict[str, object]]
    downgraded_records: list[dict[str, object]]
    validation_rows: list[dict[str, object]]
    taxonomic_level_rows: list[dict[str, object]]


def _split_semicolon(value: object) -> set[str]:
    """Split a semicolon-separated value into non-empty tokens.

    Parameters
    ----------
    value : object
        Value to split.

    Returns
    -------
    set[str]
        Non-empty stripped tokens.
    """
    if value is None:
        return set()
    text = str(value).strip()
    if not text:
        return set()
    return {item.strip() for item in text.split(";") if item.strip()}


def _join_sorted(values: Iterable[object]) -> str:
    """Join non-empty values as a stable semicolon-separated string.

    Parameters
    ----------
    values : iterable of object
        Values to join.

    Returns
    -------
    str
        Sorted semicolon-separated values.
    """
    clean_values = sorted({str(value).strip() for value in values if str(value).strip()})
    return ";".join(clean_values)


def _normalise_panel_record(
    *,
    record: dict[str, str],
    source_panel: str,
) -> dict[str, object]:
    """Normalise one panel record to the current master-panel schema.

    Parameters
    ----------
    record : dict[str, str]
        Raw input panel record.
    source_panel : str
        Path or label for the panel from which the record was read.

    Returns
    -------
    dict[str, object]
        Normalised panel record.
    """
    output: dict[str, object] = {column: record.get(column, "") for column in PANEL_REQUIRED_COLUMNS}
    output["k"] = int(str(output["k"])) if str(output["k"]).strip() else 0
    output["example_position"] = (
        int(str(output["example_position"]))
        if str(output["example_position"]).strip()
        else 0
    )
    output["source_panels"] = source_panel
    output["n_panel_records"] = 1
    return output


def load_panel_records(
    *,
    panel_paths: Sequence[str | Path],
    logger: logging.Logger | None = None,
) -> list[dict[str, object]]:
    """Load one or more KmerSutra panel files.

    Parameters
    ----------
    panel_paths : sequence of str or pathlib.Path
        Panel paths to load. Paths may point to TSV or TSV.GZ files.
    logger : logging.Logger | None, optional
        Logger for progress messages.

    Returns
    -------
    list[dict[str, object]]
        Normalised records with a ``source_panels`` column.
    """
    records: list[dict[str, object]] = []
    for panel_path in panel_paths:
        path = Path(panel_path)
        if logger:
            logger.info("Loading panel: %s", path)
        raw_records = read_tsv(input_path=path)
        if logger:
            logger.info("Loaded %d rows from %s", len(raw_records), path)
        for record in raw_records:
            records.append(
                _normalise_panel_record(
                    record=record,
                    source_panel=str(path),
                )
            )
    return records


def _record_source_taxids(record: dict[str, object]) -> set[str]:
    """Return all source taxids represented by a panel record.

    Parameters
    ----------
    record : dict[str, object]
        Panel record.

    Returns
    -------
    set[str]
        Taxids from ``source_taxids`` plus fallback evidence taxid.
    """
    taxids = _split_semicolon(record.get("source_taxids", ""))
    evidence_taxid = str(record.get("evidence_taxid", "")).strip()
    if not taxids and evidence_taxid:
        taxids.add(evidence_taxid)
    return taxids


def _choose_without_taxonomy(
    *,
    group: list[dict[str, object]],
) -> tuple[dict[str, object] | None, str]:
    """Assign a merged record without a taxonomy database.

    This fallback is intentionally conservative. It keeps k-mers that remain
    unique to one species or one clade and removes k-mers that span multiple
    clades without taxonomic context.

    Parameters
    ----------
    group : list[dict[str, object]]
        Records sharing the same k and k-mer.

    Returns
    -------
    tuple[dict[str, object] | None, str]
        Merged record and status string. The record is ``None`` when removed.
    """
    species = {str(row.get("species_name", "")).strip() for row in group}
    species.discard("")
    clades = {str(row.get("clade", "")).strip() for row in group}
    clades.discard("")

    first = group[0]
    if len(species) == 1:
        species_name = next(iter(species))
        return (
            _build_merged_record(
                group=group,
                panel_type="species_unique",
                species_name=species_name,
                clade=next(iter(clades)) if clades else "",
                evidence_taxid=str(first.get("evidence_taxid", "")),
                evidence_name=species_name,
                evidence_rank="species",
                lineage_taxids=str(first.get("lineage_taxids", "")),
            ),
            "retained_species_unique",
        )
    if len(clades) == 1:
        clade = next(iter(clades))
        return (
            _build_merged_record(
                group=group,
                panel_type="clade_core",
                species_name="",
                clade=clade,
                evidence_taxid="",
                evidence_name=clade,
                evidence_rank="clade",
                lineage_taxids="",
            ),
            "downgraded_clade_core",
        )
    return None, "removed_cross_clade_without_taxonomy"


def _build_merged_record(
    *,
    group: list[dict[str, object]],
    panel_type: str,
    species_name: str,
    clade: str,
    evidence_taxid: str,
    evidence_name: str,
    evidence_rank: str,
    lineage_taxids: str,
) -> dict[str, object]:
    """Build a single merged panel record from grouped records.

    Parameters
    ----------
    group : list[dict[str, object]]
        Records sharing the same k and k-mer.
    panel_type : str
        Assigned panel type.
    species_name : str
        Species name if species-level evidence is retained.
    clade : str
        Clade or broad group label.
    evidence_taxid : str
        Assigned evidence taxid.
    evidence_name : str
        Assigned taxonomic evidence name.
    evidence_rank : str
        Assigned evidence rank.
    lineage_taxids : str
        Assigned lineage taxids.

    Returns
    -------
    dict[str, object]
        Merged panel record.
    """
    first = group[0]
    source_taxids: set[str] = set()
    for record in group:
        source_taxids.update(_record_source_taxids(record))

    return {
        "kmer": first["kmer"],
        "k": first["k"],
        "panel_type": panel_type,
        "species_name": species_name,
        "clade": clade,
        "source_genomes": _join_sorted(row.get("source_genomes", "") for row in group),
        "source_contigs": _join_sorted(row.get("source_contigs", "") for row in group),
        "example_position": first.get("example_position", 0),
        "evidence_taxid": evidence_taxid,
        "evidence_name": evidence_name,
        "evidence_rank": evidence_rank,
        "lineage_taxids": lineage_taxids,
        "source_taxids": _join_sorted(source_taxids),
        "source_panels": _join_sorted(row.get("source_panels", "") for row in group),
        "n_panel_records": len(group),
    }


def _choose_with_taxonomy(
    *,
    group: list[dict[str, object]],
    taxonomy_db: TaxonomyDatabase,
    preferred_ranks: list[str],
) -> tuple[dict[str, object] | None, str]:
    """Assign a merged record using NCBI taxonomy.

    Parameters
    ----------
    group : list[dict[str, object]]
        Records sharing the same k and k-mer.
    taxonomy_db : TaxonomyDatabase
        Parsed taxonomy database.
    preferred_ranks : list[str]
        Ranks allowed as diagnostic evidence levels.

    Returns
    -------
    tuple[dict[str, object] | None, str]
        Merged record and status string. The record is ``None`` when removed.
    """
    source_taxids: set[str] = set()
    for record in group:
        for taxid in _record_source_taxids(record):
            normalised = taxonomy_db.normalise_taxid(taxid)
            if normalised:
                source_taxids.add(normalised)

    if not source_taxids:
        return None, "removed_no_source_taxids"

    evidence_node = taxonomy_db.best_named_ancestor(
        taxids=source_taxids,
        preferred_ranks=preferred_ranks,
    )
    if evidence_node is None:
        return None, "removed_no_supported_taxonomic_rank"

    if evidence_node.rank not in preferred_ranks:
        return None, "removed_unsupported_taxonomic_rank"

    panel_type = "species_unique" if evidence_node.rank == "species" else f"{evidence_node.rank}_core"
    species_name = evidence_node.name if evidence_node.rank == "species" else ""
    lineage_taxids = ";".join(taxonomy_db.get_lineage(evidence_node.taxid))
    record = _build_merged_record(
        group=group,
        panel_type=panel_type,
        species_name=species_name,
        clade=evidence_node.name,
        evidence_taxid=evidence_node.taxid,
        evidence_name=evidence_node.name,
        evidence_rank=evidence_node.rank,
        lineage_taxids=lineage_taxids,
    )

    original_ranks = {str(row.get("evidence_rank", "")).strip() for row in group}
    original_ranks.discard("")
    status = "retained_taxonomic_evidence"
    if original_ranks and evidence_node.rank not in original_ranks:
        status = "downgraded_taxonomic_evidence"
    if len(group) > 1 and evidence_node.rank != "species":
        status = "downgraded_taxonomic_evidence"
    return record, status


def _summarise_taxonomic_levels(records: Iterable[dict[str, object]]) -> list[dict[str, object]]:
    """Summarise merged records by taxonomic evidence level.

    Parameters
    ----------
    records : iterable of dict[str, object]
        Panel records.

    Returns
    -------
    list[dict[str, object]]
        Count summaries by k, panel type, evidence rank and evidence name.
    """
    counter: Counter[tuple[str, str, str, str]] = Counter()
    for record in records:
        key = (
            str(record.get("k", "")),
            str(record.get("panel_type", "")),
            str(record.get("evidence_rank", "")),
            str(record.get("evidence_name", "")),
        )
        counter[key] += 1

    rows = [
        {
            "k": key[0],
            "panel_type": key[1],
            "evidence_rank": key[2],
            "evidence_name": key[3],
            "n_kmers": value,
        }
        for key, value in sorted(counter.items())
    ]
    return rows


def merge_panel_records(
    *,
    records: Iterable[dict[str, object]],
    taxonomy_db: TaxonomyDatabase | None = None,
    preferred_ranks: list[str] | None = None,
    logger: logging.Logger | None = None,
) -> PanelMergeResult:
    """Merge panel records into a globally validated master panel.

    Parameters
    ----------
    records : iterable of dict[str, object]
        Input panel records.
    taxonomy_db : TaxonomyDatabase | None, optional
        Taxonomy database used to assign globally valid evidence levels.
    preferred_ranks : list[str] | None, optional
        Diagnostic evidence ranks to retain.
    logger : logging.Logger | None, optional
        Logger for progress messages.

    Returns
    -------
    PanelMergeResult
        Merged records and validation summaries.
    """
    input_records = list(records)
    ranks = preferred_ranks or CORE_RANK_ORDER
    grouped: dict[tuple[int, str], list[dict[str, object]]] = defaultdict(list)
    for record in input_records:
        grouped[(int(record.get("k", 0)), str(record.get("kmer", "")))].append(record)

    if logger:
        logger.info(
            "Merging %d panel rows across %d distinct k-mer keys",
            len(input_records),
            len(grouped),
        )
        if taxonomy_db is not None:
            logger.info("Taxonomy-aware merging enabled with ranks: %s", ", ".join(ranks))
        else:
            logger.info("Taxonomy-aware merging disabled; using conservative label-based fallback")

    master_records: list[dict[str, object]] = []
    removed_records: list[dict[str, object]] = []
    downgraded_records: list[dict[str, object]] = []
    status_counter: Counter[str] = Counter()

    for (_k, _kmer), group in grouped.items():
        if taxonomy_db is not None:
            merged_record, status = _choose_with_taxonomy(
                group=group,
                taxonomy_db=taxonomy_db,
                preferred_ranks=ranks,
            )
        else:
            merged_record, status = _choose_without_taxonomy(group=group)

        status_counter[status] += 1
        if merged_record is None:
            for record in group:
                removed = dict(record)
                removed["removal_reason"] = status
                removed_records.append(removed)
            continue

        master_records.append(merged_record)
        if status.startswith("downgraded"):
            downgraded = dict(merged_record)
            downgraded["downgrade_reason"] = status
            downgraded["original_panel_types"] = _join_sorted(row.get("panel_type", "") for row in group)
            downgraded["original_evidence_ranks"] = _join_sorted(row.get("evidence_rank", "") for row in group)
            downgraded_records.append(downgraded)

    validation_rows = [
        {"metric": "input_records", "value": len(input_records)},
        {"metric": "distinct_input_kmer_keys", "value": len(grouped)},
        {"metric": "retained_master_records", "value": len(master_records)},
        {"metric": "removed_records", "value": len(removed_records)},
        {"metric": "downgraded_records", "value": len(downgraded_records)},
    ]
    validation_rows.extend(
        {"metric": f"status_{status}", "value": count}
        for status, count in sorted(status_counter.items())
    )

    taxonomic_level_rows = _summarise_taxonomic_levels(master_records)

    if logger:
        logger.info(
            "Master panel retained %d records; removed %d records; downgraded %d records",
            len(master_records),
            len(removed_records),
            len(downgraded_records),
        )

    return PanelMergeResult(
        master_records=master_records,
        removed_records=removed_records,
        downgraded_records=downgraded_records,
        validation_rows=validation_rows,
        taxonomic_level_rows=taxonomic_level_rows,
    )


def merge_panel_files(
    *,
    panel_paths: Sequence[str | Path],
    out_dir: str | Path,
    taxonomy_db: TaxonomyDatabase | None = None,
    preferred_ranks: list[str] | None = None,
    logger: logging.Logger | None = None,
) -> PanelMergeResult:
    """Merge panel files and write a master-panel directory.

    Parameters
    ----------
    panel_paths : sequence of str or pathlib.Path
        Input panel files.
    out_dir : str or pathlib.Path
        Output directory.
    taxonomy_db : TaxonomyDatabase | None, optional
        Optional taxonomy database for taxonomic evidence reassignment.
    preferred_ranks : list[str] | None, optional
        Evidence ranks to retain.
    logger : logging.Logger | None, optional
        Logger for progress messages.

    Returns
    -------
    PanelMergeResult
        Result object containing all written table data.
    """
    output_dir = Path(out_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    records = load_panel_records(panel_paths=panel_paths, logger=logger)
    result = merge_panel_records(
        records=records,
        taxonomy_db=taxonomy_db,
        preferred_ranks=preferred_ranks,
        logger=logger,
    )

    write_tsv(
        records=result.master_records,
        output_path=output_dir / "master_kmer_panel.tsv.gz",
        fieldnames=PANEL_OUTPUT_COLUMNS,
    )
    write_tsv(
        records=result.removed_records,
        output_path=output_dir / "removed_conflicting_kmers.tsv.gz",
        fieldnames=PANEL_OUTPUT_COLUMNS + ["removal_reason"],
    )
    write_tsv(
        records=result.downgraded_records,
        output_path=output_dir / "downgraded_kmers.tsv.gz",
        fieldnames=PANEL_OUTPUT_COLUMNS
        + ["downgrade_reason", "original_panel_types", "original_evidence_ranks"],
    )
    write_tsv(
        records=result.validation_rows,
        output_path=output_dir / "master_validation_summary.tsv",
        fieldnames=["metric", "value"],
    )
    write_tsv(
        records=result.taxonomic_level_rows,
        output_path=output_dir / "taxonomic_level_summary.tsv",
        fieldnames=["k", "panel_type", "evidence_rank", "evidence_name", "n_kmers"],
    )
    write_json(
        data={
            "panel_paths": [str(Path(path)) for path in panel_paths],
            "n_input_panels": len(panel_paths),
            "n_master_records": len(result.master_records),
            "n_removed_records": len(result.removed_records),
            "n_downgraded_records": len(result.downgraded_records),
        },
        output_path=output_dir / "master_panel_metadata.json",
    )
    return result


def validate_panel_records(
    *,
    records: Iterable[dict[str, object]],
    logger: logging.Logger | None = None,
) -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]:
    """Validate one panel and return summary, issue and level tables.

    Parameters
    ----------
    records : iterable of dict[str, object]
        Panel records to validate.
    logger : logging.Logger | None, optional
        Logger for progress messages.

    Returns
    -------
    tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]
        Validation summary rows, issue rows and taxonomic-level summaries.
    """
    record_list = list(records)
    issues: list[dict[str, object]] = []
    grouped: dict[tuple[int, str], list[dict[str, object]]] = defaultdict(list)
    for record in record_list:
        kmer = str(record.get("kmer", ""))
        k = int(record.get("k", 0) or 0)
        grouped[(k, kmer)].append(record)
        for column in PANEL_REQUIRED_COLUMNS:
            if column not in record:
                issues.append(
                    {
                        "issue_type": "missing_column",
                        "severity": "error",
                        "k": k,
                        "kmer": kmer,
                        "n_records": 1,
                        "details": column,
                    }
                )
        if k and len(kmer) != k:
            issues.append(
                {
                    "issue_type": "kmer_length_mismatch",
                    "severity": "error",
                    "k": k,
                    "kmer": kmer,
                    "n_records": 1,
                    "details": f"observed_length={len(kmer)}",
                }
            )
        if not str(record.get("evidence_rank", "")).strip():
            issues.append(
                {
                    "issue_type": "missing_evidence_rank",
                    "severity": "warning",
                    "k": k,
                    "kmer": kmer,
                    "n_records": 1,
                    "details": "evidence_rank is empty",
                }
            )

    duplicate_keys = 0
    conflicting_duplicate_keys = 0
    for (k, kmer), group in grouped.items():
        if len(group) <= 1:
            continue
        duplicate_keys += 1
        evidence_keys = {
            (
                str(row.get("evidence_taxid", "")),
                str(row.get("evidence_rank", "")),
                str(row.get("evidence_name", "")),
            )
            for row in group
        }
        severity = "warning" if len(evidence_keys) == 1 else "error"
        issue_type = "duplicate_kmer_key" if len(evidence_keys) == 1 else "conflicting_duplicate_kmer_key"
        if len(evidence_keys) > 1:
            conflicting_duplicate_keys += 1
        issues.append(
            {
                "issue_type": issue_type,
                "severity": severity,
                "k": k,
                "kmer": kmer,
                "n_records": len(group),
                "details": ";".join("|".join(key) for key in sorted(evidence_keys)),
            }
        )

    summary_rows = [
        {"metric": "n_records", "value": len(record_list)},
        {"metric": "n_distinct_kmer_keys", "value": len(grouped)},
        {"metric": "n_duplicate_kmer_keys", "value": duplicate_keys},
        {"metric": "n_conflicting_duplicate_kmer_keys", "value": conflicting_duplicate_keys},
        {"metric": "n_issues", "value": len(issues)},
        {"metric": "n_error_issues", "value": sum(row["severity"] == "error" for row in issues)},
        {"metric": "n_warning_issues", "value": sum(row["severity"] == "warning" for row in issues)},
    ]
    level_rows = _summarise_taxonomic_levels(record_list)

    if logger:
        logger.info(
            "Validated %d panel records: %d issue(s), %d duplicate key(s)",
            len(record_list),
            len(issues),
            duplicate_keys,
        )
    return summary_rows, issues, level_rows


def validate_panel_file(
    *,
    panel_path: str | Path,
    out_dir: str | Path,
    logger: logging.Logger | None = None,
) -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]:
    """Validate a panel file and write validation tables.

    Parameters
    ----------
    panel_path : str or pathlib.Path
        Panel file to validate.
    out_dir : str or pathlib.Path
        Output directory.
    logger : logging.Logger | None, optional
        Logger for progress messages.

    Returns
    -------
    tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]
        Validation summary rows, issue rows and taxonomic-level summaries.
    """
    output_dir = Path(out_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    records = load_panel_records(panel_paths=[panel_path], logger=logger)
    summary_rows, issue_rows, level_rows = validate_panel_records(
        records=records,
        logger=logger,
    )
    write_tsv(
        records=summary_rows,
        output_path=output_dir / "panel_validation_summary.tsv",
        fieldnames=["metric", "value"],
    )
    write_tsv(
        records=issue_rows,
        output_path=output_dir / "panel_validation_issues.tsv",
        fieldnames=VALIDATION_ISSUE_COLUMNS,
    )
    write_tsv(
        records=level_rows,
        output_path=output_dir / "taxonomic_level_summary.tsv",
        fieldnames=["k", "panel_type", "evidence_rank", "evidence_name", "n_kmers"],
    )
    return summary_rows, issue_rows, level_rows
