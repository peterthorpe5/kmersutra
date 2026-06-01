"""Lowest-common-ancestor reporting helpers for KmerSutra outputs.

The LCA layer is intentionally conservative. It does not replace the existing
species, neighbour-lineage or background-candidate calls. Instead, it places
supported KmerSutra evidence onto the NCBI taxonomy so that unknown or partially
represented samples can be reported at the most specific defensible taxonomic
level.
"""

from __future__ import annotations

import logging
from collections import defaultdict
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path

from kmersutra.table_io import read_records_table, write_records_table
from kmersutra.taxonomy import TaxonomyDatabase, TaxonomyNode

NAME_COLUMN_CANDIDATES = [
    "species_name",
    "report_label",
    "normalised_report_label",
    "taxon_name",
    "evidence_name",
    "name",
]
TAXID_COLUMN_CANDIDATES = [
    "taxid",
    "species_taxid",
    "ncbi_taxid",
    "evidence_taxid",
    "report_taxid",
]
LCA_REPORT_FIELDNAMES = [
    "sample_id",
    "benchmark_family",
    "panel",
    "replicate",
    "spike_n",
    "expected_targets",
    "lca_scope",
    "lca_taxid",
    "lca_rank",
    "lca_name",
    "lca_lineage",
    "lca_interpretation",
    "n_input_rows",
    "n_taxa",
    "support_taxids",
    "support_names",
    "total_unique_kmers",
    "total_positive_sequences",
    "max_best_k",
    "max_k_values_positive",
    "top_taxid",
    "top_name",
    "top_unique_kmers",
    "top_positive_sequences",
    "top_best_k",
    "top_score",
]
SCOPE_COLUMN_RULES = {
    "all_supported_evidence": [],
    "dominant_lineage": [],
    "background_candidate": [
        "is_background_candidate_signal",
        "is_positive_background_candidate",
        "is_background_candidate_call",
    ],
    "neighbour_lineage": [
        "is_positive_neighbour_lineage",
        "is_expected_genus_neighbour",
    ],
    "reportable_positive": [
        "is_positive_call",
        "is_positive_expected",
    ],
}
METADATA_COLUMNS = [
    "benchmark_family",
    "panel",
    "replicate",
    "spike_n",
    "expected_targets",
]


def _as_int(*, value: object, default: int = 0) -> int:
    """Return an integer parsed from a table value.

    Parameters
    ----------
    value : object
        Input value from a KmerSutra table row.
    default : int, optional
        Value returned when parsing fails or the input is blank.

    Returns
    -------
    int
        Parsed integer.
    """
    if value in (None, ""):
        return default
    try:
        return int(float(str(value)))
    except (TypeError, ValueError):
        return default


def _as_bool(*, value: object) -> bool:
    """Return a boolean parsed from common table values.

    Parameters
    ----------
    value : object
        Input table value.

    Returns
    -------
    bool
        Parsed boolean.
    """
    return str(value).strip().lower() in {"1", "true", "t", "yes", "y"}


def _normalise_name(*, value: object) -> str:
    """Return a normalised taxon-name key.

    Parameters
    ----------
    value : object
        Raw taxon name.

    Returns
    -------
    str
        Lower-case, single-spaced name.
    """
    text = str(value or "").strip().lower()
    return " ".join(text.split())


def _first_existing_column(
    *,
    records: Sequence[Mapping[str, object]],
    candidates: Sequence[str],
) -> str | None:
    """Return the first candidate column present in records.

    Parameters
    ----------
    records : sequence of mappings
        Records to inspect.
    candidates : sequence of str
        Candidate column names.

    Returns
    -------
    str or None
        First matching column, or ``None``.
    """
    if not records:
        return None
    columns = set(records[0])
    for column in candidates:
        if column in columns:
            return column
    return None


def build_taxon_name_map(
    *,
    taxon_map_records: Iterable[Mapping[str, object]],
    name_column: str | None = None,
    taxid_column: str | None = None,
) -> dict[str, str]:
    """Build a normalised taxon-name to taxid mapping.

    Parameters
    ----------
    taxon_map_records : iterable of mappings
        Records containing taxon names and taxids.
    name_column : str or None, optional
        Name column. If omitted, a common name column is inferred.
    taxid_column : str or None, optional
        Taxid column. If omitted, a common taxid column is inferred.

    Returns
    -------
    dict[str, str]
        Mapping from normalised names to taxids.

    Raises
    ------
    ValueError
        If required columns cannot be inferred or if a name maps to multiple
        different taxids.
    """
    records = [dict(record) for record in taxon_map_records]
    if not records:
        return {}
    resolved_name_column = name_column or _first_existing_column(
        records=records,
        candidates=NAME_COLUMN_CANDIDATES,
    )
    resolved_taxid_column = taxid_column or _first_existing_column(
        records=records,
        candidates=TAXID_COLUMN_CANDIDATES,
    )
    if resolved_name_column is None:
        raise ValueError(
            "Could not infer a taxon-name column in the taxon map. "
            f"Tried: {', '.join(NAME_COLUMN_CANDIDATES)}"
        )
    if resolved_taxid_column is None:
        raise ValueError(
            "Could not infer a taxid column in the taxon map. "
            f"Tried: {', '.join(TAXID_COLUMN_CANDIDATES)}"
        )

    mapping: dict[str, str] = {}
    for row in records:
        name = _normalise_name(value=row.get(resolved_name_column, ""))
        taxid = str(row.get(resolved_taxid_column, "")).strip()
        if not name or not taxid:
            continue
        previous = mapping.get(name)
        if previous is not None and previous != taxid:
            raise ValueError(
                "Taxon map contains conflicting taxids for "
                f"{row.get(resolved_name_column)!r}: {previous} and {taxid}"
            )
        mapping[name] = taxid
    return mapping


def _resolve_row_taxid(
    *,
    row: Mapping[str, object],
    name_to_taxid: Mapping[str, str],
    taxid_column: str | None = None,
    taxon_name_column: str | None = None,
) -> str:
    """Resolve a taxid for an evidence row.

    Parameters
    ----------
    row : mapping
        KmerSutra evidence/call row.
    name_to_taxid : mapping
        Optional normalised name-to-taxid lookup.
    taxid_column : str or None, optional
        Explicit taxid column to inspect first.
    taxon_name_column : str or None, optional
        Explicit name column used for name lookup.

    Returns
    -------
    str
        Resolved taxid or an empty string.
    """
    taxid_columns = [taxid_column] if taxid_column else []
    taxid_columns.extend(column for column in TAXID_COLUMN_CANDIDATES if column not in taxid_columns)
    for column in taxid_columns:
        if column and str(row.get(column, "")).strip():
            return str(row.get(column, "")).strip()

    name_columns = [taxon_name_column] if taxon_name_column else []
    name_columns.extend(column for column in NAME_COLUMN_CANDIDATES if column not in name_columns)
    for column in name_columns:
        name = _normalise_name(value=row.get(column, ""))
        if name in name_to_taxid:
            return str(name_to_taxid[name])
    return ""


def _evidence_score(*, row: Mapping[str, object]) -> float:
    """Return a simple ranking score for one evidence row.

    Parameters
    ----------
    row : mapping
        KmerSutra evidence/call row.

    Returns
    -------
    float
        Score used only for ordering and dominant-lineage selection.
    """
    unique = _as_int(value=row.get("n_unique_kmers", 0))
    sequences = _as_int(value=row.get("n_positive_sequences", 0))
    best_k = _as_int(value=row.get("best_k", 0))
    k_values = _as_int(value=row.get("n_k_values_positive", 0))
    return float(unique) + 0.20 * sequences + 0.05 * best_k + 5.0 * k_values


def _row_passes_minima(
    *,
    row: Mapping[str, object],
    min_unique_kmers: int,
    min_positive_sequences: int,
    min_best_k: int,
) -> bool:
    """Return whether a row has enough support for LCA reporting.

    Parameters
    ----------
    row : mapping
        KmerSutra evidence/call row.
    min_unique_kmers : int
        Minimum unique k-mers.
    min_positive_sequences : int
        Minimum positive sequences.
    min_best_k : int
        Minimum longest supported k-mer.

    Returns
    -------
    bool
        True if row passes all minima.
    """
    return (
        _as_int(value=row.get("n_unique_kmers", 0)) >= min_unique_kmers
        and _as_int(value=row.get("n_positive_sequences", 0)) >= min_positive_sequences
        and _as_int(value=row.get("best_k", 0)) >= min_best_k
    )


def _scope_rows(
    *,
    rows: Sequence[dict[str, object]],
    scope: str,
) -> list[dict[str, object]]:
    """Return rows belonging to an LCA scope.

    Parameters
    ----------
    rows : sequence of dict
        Supported rows for one sample.
    scope : str
        Scope name.

    Returns
    -------
    list[dict[str, object]]
        Rows selected for the scope.
    """
    if scope in {"all_supported_evidence", "dominant_lineage"}:
        return list(rows)
    columns = SCOPE_COLUMN_RULES.get(scope, [])
    selected: list[dict[str, object]] = []
    for row in rows:
        if any(_as_bool(value=row.get(column, False)) for column in columns):
            selected.append(row)
        elif scope == "neighbour_lineage" and str(row.get("benchmark_report_layer", "")) == scope:
            selected.append(row)
        elif scope == "background_candidate" and str(row.get("benchmark_report_layer", "")) == "background_candidate_signal":
            selected.append(row)
    return selected


def _common_named_ancestor_rank(
    *,
    taxonomy: TaxonomyDatabase,
    taxids: Sequence[str],
    preferred_ranks: Sequence[str],
) -> str:
    """Return the rank of the best named LCA ancestor.

    Parameters
    ----------
    taxonomy : TaxonomyDatabase
        Taxonomy lookup.
    taxids : sequence of str
        Taxids to place.
    preferred_ranks : sequence of str
        Preferred named ranks.

    Returns
    -------
    str
        Rank name or an empty string.
    """
    node = taxonomy.best_named_ancestor(taxids=list(taxids), preferred_ranks=list(preferred_ranks))
    return node.rank if node is not None else ""


def _select_dominant_lineage_rows(
    *,
    rows: Sequence[dict[str, object]],
    taxonomy: TaxonomyDatabase,
    min_score_fraction: float,
    compatible_ranks: Sequence[str],
) -> list[dict[str, object]]:
    """Select a dominant lineage without allowing weak distant taxa to collapse LCA.

    Parameters
    ----------
    rows : sequence of dict
        Supported rows for one sample.
    taxonomy : TaxonomyDatabase
        Taxonomy lookup.
    min_score_fraction : float
        Minimum fraction of the top score required for non-top rows.
    compatible_ranks : sequence of str
        Named ranks considered compatible with the top taxon.

    Returns
    -------
    list[dict[str, object]]
        Dominant-lineage rows.
    """
    if not rows:
        return []
    ranked = sorted(rows, key=lambda item: (-float(item["_lca_score"]), item["_taxon_name"]))
    top = ranked[0]
    top_score = max(float(top["_lca_score"]), 1.0)
    selected = [top]
    for row in ranked[1:]:
        score = float(row["_lca_score"])
        if score < top_score * min_score_fraction:
            continue
        rank = _common_named_ancestor_rank(
            taxonomy=taxonomy,
            taxids=[str(top["_taxid"]), str(row["_taxid"])],
            preferred_ranks=compatible_ranks,
        )
        if rank in compatible_ranks:
            selected.append(row)
    return selected


def _lineage_summary(
    *,
    taxonomy: TaxonomyDatabase,
    taxid: str,
) -> str:
    """Return a compact named lineage string for a taxid.

    Parameters
    ----------
    taxonomy : TaxonomyDatabase
        Taxonomy lookup.
    taxid : str
        Taxid of interest.

    Returns
    -------
    str
        Semicolon-delimited ranked lineage.
    """
    parts: list[str] = []
    for node in taxonomy.get_ranked_lineage(taxid):
        if node.rank and node.rank != "no rank" and node.name:
            parts.append(f"{node.rank}:{node.name}")
    return ";".join(parts)


def _interpret_lca_node(*, node: TaxonomyNode | None) -> str:
    """Return a high-level interpretation for an LCA node.

    Parameters
    ----------
    node : TaxonomyNode or None
        LCA node.

    Returns
    -------
    str
        Interpretation label.
    """
    if node is None:
        return "no_supported_taxonomic_placement"
    rank = node.rank.lower()
    if rank == "species":
        return "species_resolved"
    if rank == "genus":
        return "genus_resolved"
    if rank in {"family", "order", "class", "phylum"}:
        return "broad_lineage_resolved"
    if rank in {"superkingdom", "no rank"}:
        return "very_broad_or_root_only"
    return "rank_resolved"


def _make_lca_record(
    *,
    sample_id: str,
    sample_rows: Sequence[dict[str, object]],
    scoped_rows: Sequence[dict[str, object]],
    scope: str,
    taxonomy: TaxonomyDatabase,
) -> dict[str, object]:
    """Create one LCA report record for one sample/scope.

    Parameters
    ----------
    sample_id : str
        Sample identifier.
    sample_rows : sequence of dict
        All rows for the sample.
    scoped_rows : sequence of dict
        Rows selected for the requested scope.
    scope : str
        LCA scope name.
    taxonomy : TaxonomyDatabase
        Taxonomy lookup.

    Returns
    -------
    dict[str, object]
        LCA report row.
    """
    metadata_source = sample_rows[0] if sample_rows else {}
    record: dict[str, object] = {"sample_id": sample_id}
    for column in METADATA_COLUMNS:
        record[column] = metadata_source.get(column, "")
    record["lca_scope"] = scope

    if not scoped_rows:
        record.update(
            {
                "lca_taxid": "",
                "lca_rank": "",
                "lca_name": "",
                "lca_lineage": "",
                "lca_interpretation": "no_supported_taxonomic_placement",
                "n_input_rows": 0,
                "n_taxa": 0,
                "support_taxids": "",
                "support_names": "",
                "total_unique_kmers": 0,
                "total_positive_sequences": 0,
                "max_best_k": 0,
                "max_k_values_positive": 0,
                "top_taxid": "",
                "top_name": "",
                "top_unique_kmers": 0,
                "top_positive_sequences": 0,
                "top_best_k": 0,
                "top_score": 0.0,
            }
        )
        return record

    unique_taxids = sorted({str(row["_taxid"]) for row in scoped_rows if row.get("_taxid")})
    lca_taxid = taxonomy.lowest_common_ancestor(unique_taxids)
    node = taxonomy.get_node(lca_taxid)
    top = sorted(
        scoped_rows,
        key=lambda item: (-float(item["_lca_score"]), item["_taxon_name"]),
    )[0]
    support_names = sorted({str(row["_taxon_name"]) for row in scoped_rows})

    record.update(
        {
            "lca_taxid": lca_taxid,
            "lca_rank": node.rank if node is not None else "",
            "lca_name": node.name if node is not None else "",
            "lca_lineage": _lineage_summary(taxonomy=taxonomy, taxid=lca_taxid),
            "lca_interpretation": _interpret_lca_node(node=node),
            "n_input_rows": len(scoped_rows),
            "n_taxa": len(unique_taxids),
            "support_taxids": ";".join(unique_taxids),
            "support_names": ";".join(support_names),
            "total_unique_kmers": sum(
                _as_int(value=row.get("n_unique_kmers", 0)) for row in scoped_rows
            ),
            "total_positive_sequences": sum(
                _as_int(value=row.get("n_positive_sequences", 0)) for row in scoped_rows
            ),
            "max_best_k": max(_as_int(value=row.get("best_k", 0)) for row in scoped_rows),
            "max_k_values_positive": max(
                _as_int(value=row.get("n_k_values_positive", 0)) for row in scoped_rows
            ),
            "top_taxid": top.get("_taxid", ""),
            "top_name": top.get("_taxon_name", ""),
            "top_unique_kmers": _as_int(value=top.get("n_unique_kmers", 0)),
            "top_positive_sequences": _as_int(value=top.get("n_positive_sequences", 0)),
            "top_best_k": _as_int(value=top.get("best_k", 0)),
            "top_score": round(float(top.get("_lca_score", 0.0)), 4),
        }
    )
    return record


def summarise_lca_by_sample(
    *,
    evidence_records: Iterable[Mapping[str, object]],
    taxonomy: TaxonomyDatabase,
    name_to_taxid: Mapping[str, str] | None = None,
    sample_column: str = "sample_id",
    taxid_column: str | None = None,
    taxon_name_column: str | None = None,
    min_unique_kmers: int = 1,
    min_positive_sequences: int = 1,
    min_best_k: int = 0,
    dominant_min_score_fraction: float = 0.25,
    dominant_compatible_ranks: Sequence[str] = ("species", "genus", "family"),
    scopes: Sequence[str] = (
        "dominant_lineage",
        "all_supported_evidence",
        "background_candidate",
        "neighbour_lineage",
        "reportable_positive",
    ),
    logger: logging.Logger | None = None,
) -> list[dict[str, object]]:
    """Summarise KmerSutra evidence using sample-level LCA placement.

    Parameters
    ----------
    evidence_records : iterable of mappings
        KmerSutra detection or evidence rows.
    taxonomy : TaxonomyDatabase
        NCBI taxonomy lookup used for LCA placement.
    name_to_taxid : mapping or None, optional
        Optional taxon-name to taxid lookup used when evidence rows do not
        already contain taxid columns.
    sample_column : str, optional
        Sample identifier column.
    taxid_column : str or None, optional
        Explicit taxid column in the evidence table.
    taxon_name_column : str or None, optional
        Explicit taxon-name column in the evidence table.
    min_unique_kmers : int, optional
        Minimum unique k-mers required for a row to contribute to LCA.
    min_positive_sequences : int, optional
        Minimum positive sequences required for a row to contribute to LCA.
    min_best_k : int, optional
        Minimum longest supported k required for a row to contribute to LCA.
    dominant_min_score_fraction : float, optional
        Fraction of the top row score required for other rows to enter the
        dominant-lineage LCA.
    dominant_compatible_ranks : sequence of str, optional
        Named ranks considered compatible with the top taxon for dominant LCA.
    scopes : sequence of str, optional
        LCA scopes to report.
    logger : logging.Logger or None, optional
        Logger for progress messages.

    Returns
    -------
    list[dict[str, object]]
        LCA summary records.
    """
    lookup = name_to_taxid or {}
    rows_by_sample: dict[str, list[dict[str, object]]] = defaultdict(list)
    n_rows = 0
    n_supported = 0
    n_without_taxid = 0

    for input_row in evidence_records:
        n_rows += 1
        sample_id = str(input_row.get(sample_column, "")).strip()
        if not sample_id:
            continue
        row = dict(input_row)
        if not _row_passes_minima(
            row=row,
            min_unique_kmers=min_unique_kmers,
            min_positive_sequences=min_positive_sequences,
            min_best_k=min_best_k,
        ):
            rows_by_sample.setdefault(sample_id, [])
            continue
        taxid = _resolve_row_taxid(
            row=row,
            name_to_taxid=lookup,
            taxid_column=taxid_column,
            taxon_name_column=taxon_name_column,
        )
        taxid = taxonomy.normalise_taxid(taxid)
        if not taxid or taxonomy.get_node(taxid) is None:
            n_without_taxid += 1
            rows_by_sample.setdefault(sample_id, [])
            continue
        name = ""
        if taxon_name_column:
            name = str(row.get(taxon_name_column, "")).strip()
        if not name:
            for column in NAME_COLUMN_CANDIDATES:
                if str(row.get(column, "")).strip():
                    name = str(row.get(column, "")).strip()
                    break
        row["_taxid"] = taxid
        row["_taxon_name"] = name or taxonomy.get_name(taxid) or taxid
        row["_lca_score"] = _evidence_score(row=row)
        rows_by_sample[sample_id].append(row)
        n_supported += 1

    if logger:
        logger.info("Read %d input rows for LCA reporting", n_rows)
        logger.info("Rows passing support and taxid filters: %d", n_supported)
        logger.info("Supported rows without resolvable taxids: %d", n_without_taxid)
        logger.info("Samples represented in LCA input: %d", len(rows_by_sample))

    records: list[dict[str, object]] = []
    for sample_id in sorted(rows_by_sample):
        sample_rows = rows_by_sample[sample_id]
        for scope in scopes:
            scoped = _scope_rows(rows=sample_rows, scope=scope)
            if scope == "dominant_lineage":
                scoped = _select_dominant_lineage_rows(
                    rows=scoped,
                    taxonomy=taxonomy,
                    min_score_fraction=dominant_min_score_fraction,
                    compatible_ranks=dominant_compatible_ranks,
                )
            records.append(
                _make_lca_record(
                    sample_id=sample_id,
                    sample_rows=sample_rows,
                    scoped_rows=scoped,
                    scope=scope,
                    taxonomy=taxonomy,
                )
            )
    return records


def summarise_lca_table(
    *,
    evidence_table: str | Path,
    output_table: str | Path,
    taxonomy_dir: str | Path,
    taxon_map_table: str | Path | None = None,
    download_taxonomy_if_missing: bool = False,
    sample_column: str = "sample_id",
    taxid_column: str | None = None,
    taxon_name_column: str | None = None,
    taxon_map_name_column: str | None = None,
    taxon_map_taxid_column: str | None = None,
    min_unique_kmers: int = 1,
    min_positive_sequences: int = 1,
    min_best_k: int = 0,
    dominant_min_score_fraction: float = 0.25,
    scopes: Sequence[str] = (
        "dominant_lineage",
        "all_supported_evidence",
        "background_candidate",
        "neighbour_lineage",
        "reportable_positive",
    ),
    logger: logging.Logger | None = None,
) -> list[dict[str, object]]:
    """Read KmerSutra evidence and write an LCA report table.

    Parameters
    ----------
    evidence_table : str or pathlib.Path
        Input KmerSutra detection/evidence table.
    output_table : str or pathlib.Path
        Output LCA report table.
    taxonomy_dir : str or pathlib.Path
        Directory containing NCBI taxdump files.
    taxon_map_table : str or pathlib.Path or None, optional
        Optional table mapping taxon names to taxids.
    download_taxonomy_if_missing : bool, optional
        Download NCBI taxdump if missing.
    sample_column : str, optional
        Sample identifier column.
    taxid_column : str or None, optional
        Taxid column in the evidence table.
    taxon_name_column : str or None, optional
        Taxon-name column in the evidence table.
    taxon_map_name_column : str or None, optional
        Name column in the taxon map.
    taxon_map_taxid_column : str or None, optional
        Taxid column in the taxon map.
    min_unique_kmers : int, optional
        Minimum unique k-mers for row inclusion.
    min_positive_sequences : int, optional
        Minimum positive sequences for row inclusion.
    min_best_k : int, optional
        Minimum best k value for row inclusion.
    dominant_min_score_fraction : float, optional
        Minimum score fraction for dominant-lineage non-top rows.
    scopes : sequence of str, optional
        LCA scopes to report.
    logger : logging.Logger or None, optional
        Logger.

    Returns
    -------
    list[dict[str, object]]
        Written LCA records.
    """
    taxonomy = TaxonomyDatabase.from_taxdump(
        taxonomy_dir=taxonomy_dir,
        download_if_missing=download_taxonomy_if_missing,
        logger=logger,
    )
    evidence_records = read_records_table(
        input_path=evidence_table,
        required_columns=[sample_column],
        logger=logger,
    )
    name_to_taxid: dict[str, str] = {}
    if taxon_map_table is not None:
        taxon_map_records = read_records_table(input_path=taxon_map_table, logger=logger)
        name_to_taxid = build_taxon_name_map(
            taxon_map_records=taxon_map_records,
            name_column=taxon_map_name_column,
            taxid_column=taxon_map_taxid_column,
        )
        if logger:
            logger.info("Loaded %d name-to-taxid mappings", len(name_to_taxid))

    records = summarise_lca_by_sample(
        evidence_records=evidence_records,
        taxonomy=taxonomy,
        name_to_taxid=name_to_taxid,
        sample_column=sample_column,
        taxid_column=taxid_column,
        taxon_name_column=taxon_name_column,
        min_unique_kmers=min_unique_kmers,
        min_positive_sequences=min_positive_sequences,
        min_best_k=min_best_k,
        dominant_min_score_fraction=dominant_min_score_fraction,
        scopes=scopes,
        logger=logger,
    )
    write_records_table(
        records=records,
        output_path=output_table,
        fieldnames=LCA_REPORT_FIELDNAMES,
        logger=logger,
    )
    return records
