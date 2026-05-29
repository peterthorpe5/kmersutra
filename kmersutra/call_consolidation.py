"""Consolidate KmerSutra species calls into reportable interpretation layers.

This module sits after raw species thresholding. It does not discard raw
species evidence. Instead, it labels dominated same-genus neighbours as
neighbour-lineage evidence and separates user-specified empirical-background
candidate taxa from ordinary off-target species calls.
"""

from __future__ import annotations

import logging
import math
import re
from collections import defaultdict
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path

REPORTABLE_SPECIES_CALLS = {
    "present_high_confidence",
    "present_in_mixed_sample",
    "mixed_species_present",
    "present",
}
BACKGROUND_CANDIDATE_CALL = "background_candidate_signal"
NEIGHBOUR_LINEAGE_CALL = "neighbour_lineage_evidence"
CONSOLIDATED_CALL_FIELDNAMES = [
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
    "pre_consolidation_call",
    "call",
    "report_layer",
    "consolidation_reason",
    "taxon_genus",
    "primary_species",
    "primary_unique_kmers",
    "primary_to_candidate_unique_margin",
    "primary_to_candidate_unique_ratio",
    "is_background_candidate",
]


def normalise_taxon_name(*, value: object) -> str:
    """Normalise a taxon label for robust matching.

    Parameters
    ----------
    value : object
        Raw taxon label.

    Returns
    -------
    str
        Lower-case, whitespace-normalised taxon label.
    """
    if value is None:
        return ""
    text = str(value).replace("_", " ").strip().lower()
    text = re.sub(r"\s+", " ", text)
    return text


def extract_genus(*, species_name: object) -> str:
    """Extract the first binomial token as a genus proxy.

    Parameters
    ----------
    species_name : object
        Raw species name.

    Returns
    -------
    str
        Normalised genus token, or an empty string when unavailable.
    """
    label = normalise_taxon_name(value=species_name)
    if not label:
        return ""
    return label.split(" ", maxsplit=1)[0]


def numeric_value(*, row: Mapping[str, object], key: str, default: float = 0.0) -> float:
    """Return a numeric field from a row with defensive coercion.

    Parameters
    ----------
    row : Mapping[str, object]
        Input record.
    key : str
        Field to parse.
    default : float, optional
        Value returned for missing, blank or unparsable fields.

    Returns
    -------
    float
        Parsed value.
    """
    value = row.get(key, default)
    if value in (None, ""):
        return default
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return default
    if math.isnan(parsed):
        return default
    return parsed


def read_taxon_list(*, path: str | Path | None) -> list[str]:
    """Read one taxon label per line from a text file.

    Parameters
    ----------
    path : str, pathlib.Path or None
        Optional file path.

    Returns
    -------
    list[str]
        Normalised taxon labels, preserving first-seen order.

    Raises
    ------
    FileNotFoundError
        If a non-empty path is supplied but does not exist.
    """
    if path is None or str(path).strip() == "":
        return []
    input_path = Path(path)
    if not input_path.is_file():
        raise FileNotFoundError(f"Taxon list file not found: {input_path}")
    labels: list[str] = []
    seen: set[str] = set()
    for line in input_path.read_text(encoding="utf-8").splitlines():
        if line.strip().startswith("#"):
            continue
        label = normalise_taxon_name(value=line)
        if not label or label in seen:
            continue
        seen.add(label)
        labels.append(label)
    return labels


def merge_background_taxa(
    *,
    background_candidate_taxa: Sequence[str] | None = None,
    background_candidate_file: str | Path | None = None,
) -> set[str]:
    """Merge background-candidate labels from CLI values and a file.

    Parameters
    ----------
    background_candidate_taxa : sequence of str or None, optional
        Taxon labels supplied directly.
    background_candidate_file : str, pathlib.Path or None, optional
        Optional text file with one taxon label per line.

    Returns
    -------
    set[str]
        Normalised taxon labels.
    """
    labels = {
        normalise_taxon_name(value=label)
        for label in (background_candidate_taxa or [])
        if normalise_taxon_name(value=label)
    }
    labels.update(read_taxon_list(path=background_candidate_file))
    return labels


def is_reportable_species_call(*, call: object) -> bool:
    """Return whether a call label represents reportable species presence.

    Parameters
    ----------
    call : object
        Call label.

    Returns
    -------
    bool
        True for reportable species calls.
    """
    return str(call) in REPORTABLE_SPECIES_CALLS


def _rank_species_rows(*, rows: Sequence[dict[str, object]]) -> list[dict[str, object]]:
    """Sort species rows by descending support.

    Parameters
    ----------
    rows : sequence of dict
        Species-call rows.

    Returns
    -------
    list[dict[str, object]]
        Ranked rows.
    """
    return sorted(
        rows,
        key=lambda row: (
            -numeric_value(row=row, key="n_unique_kmers"),
            -numeric_value(row=row, key="n_positive_sequences"),
            -numeric_value(row=row, key="best_k"),
            str(row.get("species_name", "")),
        ),
    )


def _dominates_candidate(
    *,
    primary: Mapping[str, object],
    candidate: Mapping[str, object],
    min_margin: int,
    min_ratio: float,
) -> tuple[bool, float, float]:
    """Return whether a primary species dominates a candidate species.

    Parameters
    ----------
    primary : Mapping[str, object]
        Highest-supported same-genus species row.
    candidate : Mapping[str, object]
        Candidate row being evaluated.
    min_margin : int
        Minimum unique-k-mer margin for demotion.
    min_ratio : float
        Minimum primary/candidate unique-k-mer ratio for demotion.

    Returns
    -------
    tuple[bool, float, float]
        Whether the candidate is dominated, the margin and the ratio.
    """
    primary_unique = numeric_value(row=primary, key="n_unique_kmers")
    candidate_unique = numeric_value(row=candidate, key="n_unique_kmers")
    margin = primary_unique - candidate_unique
    if candidate_unique <= 0:
        ratio = math.inf if primary_unique > 0 else 0.0
    else:
        ratio = primary_unique / candidate_unique
    return margin >= min_margin and ratio >= min_ratio, margin, ratio


def consolidate_species_calls(
    *,
    species_calls: Iterable[Mapping[str, object]],
    background_candidate_taxa: Sequence[str] | set[str] | None = None,
    demote_same_genus_neighbours: bool = True,
    dominant_species_min_margin: int = 25,
    dominant_species_min_ratio: float = 2.0,
    logger: logging.Logger | None = None,
) -> list[dict[str, object]]:
    """Consolidate raw species calls into reportable interpretation layers.

    Parameters
    ----------
    species_calls : iterable of mappings
        Raw species-call rows from thresholding.
    background_candidate_taxa : sequence of str, set of str or None, optional
        Taxa treated as plausible empirical-background candidates rather than
        ordinary off-target species calls.
    demote_same_genus_neighbours : bool, optional
        Whether to demote dominated same-genus reportable species calls to
        neighbouring-lineage evidence.
    dominant_species_min_margin : int, optional
        Minimum unique-k-mer margin required for same-genus demotion.
    dominant_species_min_ratio : float, optional
        Minimum primary/candidate unique-k-mer ratio required for demotion.
    logger : logging.Logger or None, optional
        Optional logger for summary messages.

    Returns
    -------
    list[dict[str, object]]
        Consolidated call rows. The original call is retained in
        ``pre_consolidation_call`` and the reportable call is written to
        ``call``.
    """
    background_taxa = {
        normalise_taxon_name(value=taxon) for taxon in (background_candidate_taxa or [])
    }
    rows = [dict(row) for row in species_calls]
    by_sample_genus: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        sample_id = str(row.get("sample_id", ""))
        genus = extract_genus(species_name=row.get("species_name", ""))
        if sample_id and genus:
            by_sample_genus[(sample_id, genus)].append(row)

    primary_by_sample_genus: dict[tuple[str, str], dict[str, object]] = {}
    for key, group in by_sample_genus.items():
        reportable = [row for row in group if is_reportable_species_call(call=row.get("call"))]
        if reportable:
            primary_by_sample_genus[key] = _rank_species_rows(rows=reportable)[0]

    output: list[dict[str, object]] = []
    demoted = 0
    background = 0
    for row in rows:
        species_name = str(row.get("species_name", ""))
        taxon_norm = normalise_taxon_name(value=species_name)
        sample_id = str(row.get("sample_id", ""))
        genus = extract_genus(species_name=species_name)
        original_call = str(row.get("call", ""))
        call = original_call
        report_layer = "raw_species_evidence"
        reason = "unchanged"
        primary_species = ""
        primary_unique = ""
        margin_value: float | str = ""
        ratio_value: float | str = ""
        is_background = taxon_norm in background_taxa

        if is_background and is_reportable_species_call(call=original_call):
            call = BACKGROUND_CANDIDATE_CALL
            report_layer = "background_candidate"
            reason = "taxon_supplied_as_empirical_background_candidate"
            background += 1
        elif demote_same_genus_neighbours and is_reportable_species_call(call=original_call):
            primary = primary_by_sample_genus.get((sample_id, genus))
            if primary is not None and str(primary.get("species_name", "")) != species_name:
                dominated, margin, ratio = _dominates_candidate(
                    primary=primary,
                    candidate=row,
                    min_margin=dominant_species_min_margin,
                    min_ratio=dominant_species_min_ratio,
                )
                primary_species = str(primary.get("species_name", ""))
                primary_unique = int(numeric_value(row=primary, key="n_unique_kmers"))
                margin_value = round(margin, 4)
                ratio_value = "inf" if math.isinf(ratio) else round(ratio, 4)
                if dominated:
                    call = NEIGHBOUR_LINEAGE_CALL
                    report_layer = "neighbour_lineage_evidence"
                    reason = "dominated_same_genus_neighbour"
                    demoted += 1
                else:
                    report_layer = "reportable_species"
                    reason = "co_dominant_or_insufficiently_dominated_species"
            else:
                report_layer = "reportable_species"
                reason = "primary_species_for_sample_genus"
        elif call in {"neighbour_lineage_evidence", "observed_below_threshold", "present_low_confidence"}:
            report_layer = "lineage_context"
            reason = "non_reportable_species_context"
        elif call == "not_detected":
            report_layer = "not_detected"
            reason = "no_reportable_signal"

        row["pre_consolidation_call"] = original_call
        row["call"] = call
        row["report_layer"] = report_layer
        row["consolidation_reason"] = reason
        row["taxon_genus"] = genus
        row["primary_species"] = primary_species
        row["primary_unique_kmers"] = primary_unique
        row["primary_to_candidate_unique_margin"] = margin_value
        row["primary_to_candidate_unique_ratio"] = ratio_value
        row["is_background_candidate"] = int(is_background)
        output.append(row)

    if logger is not None:
        logger.info(
            "Consolidated %d species calls; demoted_same_genus=%d; background_candidates=%d",
            len(output),
            demoted,
            background,
        )
    return output
