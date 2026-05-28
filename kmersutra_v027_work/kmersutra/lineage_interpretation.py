"""Interpret KmerSutra evidence at species and taxonomic-lineage levels."""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterable


REPORTABLE_SPECIES_CALLS = {
    "present_high_confidence",
    "present_in_mixed_sample",
}
WEAK_OR_CONFLICTING_SPECIES_CALLS = {
    "present_low_confidence",
    "observed_below_threshold",
    "ambiguous_conflicting_signal",
    "neighbour_lineage_evidence",
}
TAXONOMIC_RANK_ORDER = {
    "species": 0,
    "genus": 1,
    "family": 2,
    "order": 3,
    "class": 4,
    "phylum": 5,
    "superkingdom": 6,
    "clade": 7,
}
UNRESOLVED_REPORT_RANKS = {
    "genus",
    "family",
    "order",
    "class",
    "phylum",
    "superkingdom",
    "clade",
}
LINEAGE_INTERPRETATION_FIELDNAMES = [
    "sample_id",
    "lineage_call",
    "possible_novel_lineage",
    "novelty_reason",
    "report_rank",
    "report_name",
    "report_taxid",
    "n_reportable_species",
    "n_observed_species",
    "n_low_evidence_species",
    "n_conflicting_species",
    "n_neighbour_lineage_species",
    "n_neighbour_species",
    "best_species",
    "best_species_call",
    "best_species_unique_kmers",
    "second_species",
    "second_species_unique_kmers",
    "best_species_margin",
    "best_species_ratio",
    "best_taxonomic_rank",
    "best_taxonomic_name",
    "best_taxonomic_taxid",
    "best_taxonomic_unique_kmers",
    "best_taxonomic_positive_sequences",
    "best_taxonomic_k_values",
    "best_taxonomic_best_k",
    "best_taxonomic_confidence_score",
]


def _integer_from_row(*, row: dict[str, object], key: str) -> int:
    """Return an integer value from a dictionary row.

    Parameters
    ----------
    row : dict[str, object]
        Input row.
    key : str
        Field name to parse.

    Returns
    -------
    int
        Parsed integer, or zero when the field is blank or missing.
    """
    value = row.get(key, 0)
    if value in (None, ""):
        return 0
    return int(float(value))


def _float_from_row(*, row: dict[str, object], key: str) -> float:
    """Return a floating-point value from a dictionary row.

    Parameters
    ----------
    row : dict[str, object]
        Input row.
    key : str
        Field name to parse.

    Returns
    -------
    float
        Parsed float, or zero when the field is blank or missing.
    """
    value = row.get(key, 0.0)
    if value in (None, ""):
        return 0.0
    return float(value)


def _rank_priority(*, rank: str) -> int:
    """Return ordering priority for a taxonomic rank.

    Parameters
    ----------
    rank : str
        Taxonomic rank.

    Returns
    -------
    int
        Lower values are more specific.
    """
    return TAXONOMIC_RANK_ORDER.get(rank.lower(), 99)


def _passes_taxonomic_thresholds(
    *,
    row: dict[str, object],
    min_taxonomic_unique_kmers: int,
    min_taxonomic_positive_sequences: int,
    min_taxonomic_k_values: int,
    min_taxonomic_best_k: int,
) -> bool:
    """Return whether a taxonomic evidence row passes unresolved-lineage minima.

    Parameters
    ----------
    row : dict[str, object]
        Taxonomic evidence row.
    min_taxonomic_unique_kmers : int
        Minimum unique k-mers.
    min_taxonomic_positive_sequences : int
        Minimum independent positive sequences.
    min_taxonomic_k_values : int
        Minimum number of positive k values.
    min_taxonomic_best_k : int
        Minimum longest supported k value.

    Returns
    -------
    bool
        True when all taxonomic-evidence requirements are met.
    """
    return (
        _integer_from_row(row=row, key="n_unique_kmers") >= min_taxonomic_unique_kmers
        and _integer_from_row(row=row, key="n_positive_sequences")
        >= min_taxonomic_positive_sequences
        and _integer_from_row(row=row, key="n_k_values_positive") >= min_taxonomic_k_values
        and _integer_from_row(row=row, key="best_k") >= min_taxonomic_best_k
    )


def calculate_taxonomic_confidence_score(
    *,
    n_unique_kmers: int,
    n_positive_sequences: int,
    n_k_values_positive: int,
    best_k: int,
) -> float:
    """Calculate a bounded heuristic score for unresolved lineage evidence.

    This score is intentionally conservative and should not be interpreted as a
    calibrated probability. It is used to help rank genus-level or broader
    evidence when a species-level call is not supported.

    Parameters
    ----------
    n_unique_kmers : int
        Number of unique diagnostic k-mers.
    n_positive_sequences : int
        Number of independent positive reads or contigs.
    n_k_values_positive : int
        Number of positive k-mer lengths.
    best_k : int
        Longest supported k-mer length.

    Returns
    -------
    float
        Score between zero and one.
    """
    kmer_score = min(max(n_unique_kmers, 0) / 20.0, 1.0)
    sequence_score = min(max(n_positive_sequences, 0) / 5.0, 1.0)
    multi_k_score = min(max(n_k_values_positive, 0) / 2.0, 1.0)
    long_k_score = min(max(best_k, 0) / 101.0, 1.0)
    score = 0.35 * kmer_score + 0.25 * sequence_score + 0.20 * multi_k_score + 0.20 * long_k_score
    return round(max(0.0, min(1.0, score)), 4)


def _best_taxonomic_row(
    *,
    rows: Iterable[dict[str, object]],
    unresolved_report_ranks: set[str],
) -> dict[str, object] | None:
    """Return the strongest most-specific unresolved taxonomic evidence row.

    Parameters
    ----------
    rows : iterable of dict[str, object]
        Taxonomic evidence rows.
    unresolved_report_ranks : set[str]
        Ranks that may be used as unresolved-lineage reports.

    Returns
    -------
    dict[str, object] or None
        Best row, or ``None`` when no eligible row exists.
    """
    eligible = [
        row
        for row in rows
        if str(row.get("evidence_rank", "")).lower() in unresolved_report_ranks
        and _integer_from_row(row=row, key="n_unique_kmers") > 0
    ]
    if not eligible:
        return None
    return sorted(
        eligible,
        key=lambda row: (
            _rank_priority(rank=str(row.get("evidence_rank", ""))),
            -_integer_from_row(row=row, key="n_unique_kmers"),
            -_integer_from_row(row=row, key="n_positive_sequences"),
            -_integer_from_row(row=row, key="best_k"),
            str(row.get("evidence_name", "")),
        ),
    )[0]


def _species_rank_rows(*, rows: Iterable[dict[str, object]]) -> list[dict[str, object]]:
    """Return species call rows ordered by unique k-mer support.

    Parameters
    ----------
    rows : iterable of dict[str, object]
        Species detection-call rows.

    Returns
    -------
    list[dict[str, object]]
        Rows sorted by descending species-level support.
    """
    return sorted(
        [row for row in rows if str(row.get("species_name", ""))],
        key=lambda row: (
            -_integer_from_row(row=row, key="n_unique_kmers"),
            -_integer_from_row(row=row, key="n_positive_sequences"),
            -_integer_from_row(row=row, key="best_k"),
            str(row.get("species_name", "")),
        ),
    )


def _summarise_species_context(
    *,
    species_rows: Iterable[dict[str, object]],
) -> dict[str, object]:
    """Summarise species-level context for a sample.

    Parameters
    ----------
    species_rows : iterable of dict[str, object]
        Species detection-call rows for one sample.

    Returns
    -------
    dict[str, object]
        Species context values used by lineage interpretation.
    """
    ranked = _species_rank_rows(rows=species_rows)
    best = ranked[0] if ranked else {}
    second = ranked[1] if len(ranked) > 1 else {}
    best_unique = _integer_from_row(row=best, key="n_unique_kmers") if best else 0
    second_unique = _integer_from_row(row=second, key="n_unique_kmers") if second else 0
    ratio = ""
    if second_unique > 0:
        ratio = round(best_unique / second_unique, 4)
    elif best_unique > 0:
        ratio = "inf"

    return {
        "n_observed_species": sum(
            _integer_from_row(row=row, key="n_unique_kmers") > 0 for row in species_rows
        ),
        "n_low_evidence_species": sum(
            str(row.get("call", "")) in {"present_low_confidence", "observed_below_threshold"}
            for row in species_rows
        ),
        "n_conflicting_species": sum(
            str(row.get("call", "")) == "ambiguous_conflicting_signal" for row in species_rows
        ),
        "n_neighbour_lineage_species": sum(
            str(row.get("call", "")) == "neighbour_lineage_evidence"
            for row in species_rows
        ),
        "best_species": best.get("species_name", ""),
        "best_species_call": best.get("call", ""),
        "best_species_unique_kmers": best_unique,
        "second_species": second.get("species_name", ""),
        "second_species_unique_kmers": second_unique,
        "best_species_margin": best_unique - second_unique,
        "best_species_ratio": ratio,
    }


def interpret_lineage_evidence(
    *,
    species_calls: Iterable[dict[str, object]],
    taxonomic_evidence: Iterable[dict[str, object]],
    min_taxonomic_unique_kmers: int = 20,
    min_taxonomic_positive_sequences: int = 5,
    min_taxonomic_k_values: int = 1,
    min_taxonomic_best_k: int = 77,
    min_taxonomic_confidence_score: float = 0.40,
    min_neighbour_species_for_novelty: int = 2,
    unresolved_report_ranks: Iterable[str] | None = None,
) -> list[dict[str, object]]:
    """Interpret species and lineage evidence for unresolved or novel signals.

    The function does not replace species-level detection calls. Instead, it
    adds a conservative sample-level interpretation. Reportable species calls
    remain species calls; weak neighbouring-species evidence is retained as
    context and may support an unresolved or possible-novel lineage call when
    genus-level or broader evidence is strong enough.

    Parameters
    ----------
    species_calls : iterable of dict[str, object]
        Rows from ``species_detection_calls.tsv``.
    taxonomic_evidence : iterable of dict[str, object]
        Rows from ``sample_taxonomic_kmer_evidence.tsv``.
    min_taxonomic_unique_kmers : int, optional
        Minimum unique k-mers for unresolved taxonomic support.
    min_taxonomic_positive_sequences : int, optional
        Minimum independent positive sequences for unresolved support.
    min_taxonomic_k_values : int, optional
        Minimum positive k values for unresolved support.
    min_taxonomic_best_k : int, optional
        Minimum longest k value for unresolved support.
    min_taxonomic_confidence_score : float, optional
        Minimum heuristic taxonomic confidence score for possible novelty.
    min_neighbour_species_for_novelty : int, optional
        Minimum weak neighbouring species count that strengthens a possible
        novel/unsampled lineage interpretation.
    unresolved_report_ranks : iterable of str or None, optional
        Taxonomic ranks eligible for unresolved-lineage reporting. Species rank
        is excluded by default so weak species evidence is not overcalled.

    Returns
    -------
    list[dict[str, object]]
        One lineage-interpretation record per sample.
    """
    allowed_ranks = {
        rank.lower()
        for rank in (unresolved_report_ranks or UNRESOLVED_REPORT_RANKS)
    }
    calls_by_sample: dict[str, list[dict[str, object]]] = defaultdict(list)
    taxonomic_by_sample: dict[str, list[dict[str, object]]] = defaultdict(list)
    sample_ids: set[str] = set()

    for row in species_calls:
        sample_id = str(row.get("sample_id", ""))
        if not sample_id:
            continue
        calls_by_sample[sample_id].append(row)
        sample_ids.add(sample_id)

    for row in taxonomic_evidence:
        sample_id = str(row.get("sample_id", ""))
        if not sample_id:
            continue
        taxonomic_by_sample[sample_id].append(row)
        sample_ids.add(sample_id)

    records: list[dict[str, object]] = []
    for sample_id in sorted(sample_ids):
        sample_species_rows = calls_by_sample.get(sample_id, [])
        sample_taxonomic_rows = taxonomic_by_sample.get(sample_id, [])
        species_context = _summarise_species_context(species_rows=sample_species_rows)
        reportable_species = [
            row
            for row in sample_species_rows
            if str(row.get("call", "")) in REPORTABLE_SPECIES_CALLS
        ]
        weak_neighbour_rows = [
            row
            for row in sample_species_rows
            if str(row.get("call", "")) in WEAK_OR_CONFLICTING_SPECIES_CALLS
            and _integer_from_row(row=row, key="n_unique_kmers") > 0
        ]
        best_taxonomic = _best_taxonomic_row(
            rows=sample_taxonomic_rows,
            unresolved_report_ranks=allowed_ranks,
        )
        tax_score = 0.0
        tax_passes = False
        if best_taxonomic is not None:
            tax_score = calculate_taxonomic_confidence_score(
                n_unique_kmers=_integer_from_row(row=best_taxonomic, key="n_unique_kmers"),
                n_positive_sequences=_integer_from_row(
                    row=best_taxonomic,
                    key="n_positive_sequences",
                ),
                n_k_values_positive=_integer_from_row(
                    row=best_taxonomic,
                    key="n_k_values_positive",
                ),
                best_k=_integer_from_row(row=best_taxonomic, key="best_k"),
            )
            tax_passes = (
                _passes_taxonomic_thresholds(
                    row=best_taxonomic,
                    min_taxonomic_unique_kmers=min_taxonomic_unique_kmers,
                    min_taxonomic_positive_sequences=min_taxonomic_positive_sequences,
                    min_taxonomic_k_values=min_taxonomic_k_values,
                    min_taxonomic_best_k=min_taxonomic_best_k,
                )
                and tax_score >= min_taxonomic_confidence_score
            )

        if len(reportable_species) > 1:
            lineage_call = "mixed_species_detected"
            possible_novel = False
            novelty_reason = "reportable_species_present"
            report_rank = "species"
            report_name = ";".join(str(row.get("species_name", "")) for row in reportable_species)
            report_taxid = ""
        elif len(reportable_species) == 1:
            lineage_call = "species_detected"
            possible_novel = False
            novelty_reason = "reportable_species_present"
            report_rank = "species"
            report_name = str(reportable_species[0].get("species_name", ""))
            report_taxid = ""
        elif tax_passes and len(weak_neighbour_rows) >= min_neighbour_species_for_novelty:
            lineage_call = "possible_novel_or_unsampled_lineage"
            possible_novel = True
            novelty_reason = "strong_taxonomic_signal_with_multiple_weak_neighbours"
            report_rank = str(best_taxonomic.get("evidence_rank", ""))
            report_name = str(best_taxonomic.get("evidence_name", ""))
            report_taxid = str(best_taxonomic.get("evidence_taxid", ""))
        elif tax_passes:
            lineage_call = "unresolved_taxonomic_signal"
            possible_novel = False
            novelty_reason = "strong_taxonomic_signal_without_species_call"
            report_rank = str(best_taxonomic.get("evidence_rank", ""))
            report_name = str(best_taxonomic.get("evidence_name", ""))
            report_taxid = str(best_taxonomic.get("evidence_taxid", ""))
        elif weak_neighbour_rows:
            lineage_call = "weak_unresolved_neighbour_signal"
            possible_novel = False
            novelty_reason = "weak_species_neighbour_evidence_below_taxonomic_threshold"
            report_rank = "unresolved"
            report_name = "weak_neighbour_signal"
            report_taxid = ""
        else:
            lineage_call = "no_supported_signal"
            possible_novel = False
            novelty_reason = "no_species_or_taxonomic_evidence"
            report_rank = "none"
            report_name = "not_detected"
            report_taxid = ""

        best_taxonomic_rank = ""
        best_taxonomic_name = ""
        best_taxonomic_taxid = ""
        best_taxonomic_unique = 0
        best_taxonomic_sequences = 0
        best_taxonomic_k_values = 0
        best_taxonomic_best_k = 0
        if best_taxonomic is not None:
            best_taxonomic_rank = str(best_taxonomic.get("evidence_rank", ""))
            best_taxonomic_name = str(best_taxonomic.get("evidence_name", ""))
            best_taxonomic_taxid = str(best_taxonomic.get("evidence_taxid", ""))
            best_taxonomic_unique = _integer_from_row(
                row=best_taxonomic,
                key="n_unique_kmers",
            )
            best_taxonomic_sequences = _integer_from_row(
                row=best_taxonomic,
                key="n_positive_sequences",
            )
            best_taxonomic_k_values = _integer_from_row(
                row=best_taxonomic,
                key="n_k_values_positive",
            )
            best_taxonomic_best_k = _integer_from_row(row=best_taxonomic, key="best_k")

        records.append(
            {
                "sample_id": sample_id,
                "lineage_call": lineage_call,
                "possible_novel_lineage": possible_novel,
                "novelty_reason": novelty_reason,
                "report_rank": report_rank,
                "report_name": report_name,
                "report_taxid": report_taxid,
                "n_reportable_species": len(reportable_species),
                "n_observed_species": species_context["n_observed_species"],
                "n_low_evidence_species": species_context["n_low_evidence_species"],
                "n_conflicting_species": species_context["n_conflicting_species"],
                "n_neighbour_lineage_species": species_context[
                    "n_neighbour_lineage_species"
                ],
                "n_neighbour_species": len(weak_neighbour_rows),
                "best_species": species_context["best_species"],
                "best_species_call": species_context["best_species_call"],
                "best_species_unique_kmers": species_context["best_species_unique_kmers"],
                "second_species": species_context["second_species"],
                "second_species_unique_kmers": species_context["second_species_unique_kmers"],
                "best_species_margin": species_context["best_species_margin"],
                "best_species_ratio": species_context["best_species_ratio"],
                "best_taxonomic_rank": best_taxonomic_rank,
                "best_taxonomic_name": best_taxonomic_name,
                "best_taxonomic_taxid": best_taxonomic_taxid,
                "best_taxonomic_unique_kmers": best_taxonomic_unique,
                "best_taxonomic_positive_sequences": best_taxonomic_sequences,
                "best_taxonomic_k_values": best_taxonomic_k_values,
                "best_taxonomic_best_k": best_taxonomic_best_k,
                "best_taxonomic_confidence_score": round(tax_score, 4),
            }
        )
    return records
