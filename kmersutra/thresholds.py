"""Detection-call and confidence scoring logic for KmerSutra."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from typing import Literal


WEAK_SIGNAL_CALLS = {"present_low_confidence", "observed_below_threshold"}
NON_REPORTABLE_LINEAGE_CALLS = {"neighbour_lineage_evidence"}


@dataclass(frozen=True)
class SpeciesCallPreset:
    """Threshold settings for species-level detection calls.

    Attributes
    ----------
    min_unique_kmers : int
        Minimum number of unique diagnostic k-mers required.
    min_positive_sequences : int
        Minimum number of independent positive sequences required.
    min_k_values_positive : int
        Minimum number of k-mer lengths with positive evidence.
    max_conflict_ratio : float
        Maximum conflicting-evidence ratio for a high-confidence call.
    min_best_k : int
        Minimum longest supported k-mer length required.
    min_exact_hits : int
        Minimum number of exact k-mer hits required.
    min_confidence_score : float
        Minimum heuristic confidence score required.
    low_evidence_call : str
        Call label used when evidence exists but does not pass thresholds.
    min_mixed_species_fraction : float
        Minimum fraction of the strongest species-level support required for a
        passing species to be promoted to a reportable mixed-species call.
        Species below this fraction are retained as neighbouring-lineage
        evidence rather than over-reported as present species.
    """

    min_unique_kmers: int
    min_positive_sequences: int
    min_k_values_positive: int
    max_conflict_ratio: float
    min_best_k: int
    min_exact_hits: int
    min_confidence_score: float
    low_evidence_call: str
    min_mixed_species_fraction: float


CALL_PRESETS: dict[str, SpeciesCallPreset] = {
    "legacy": SpeciesCallPreset(
        min_unique_kmers=3,
        min_positive_sequences=2,
        min_k_values_positive=1,
        max_conflict_ratio=0.10,
        min_best_k=0,
        min_exact_hits=0,
        min_confidence_score=0.0,
        low_evidence_call="present_low_confidence",
        min_mixed_species_fraction=0.0,
    ),
    "conservative": SpeciesCallPreset(
        min_unique_kmers=20,
        min_positive_sequences=5,
        min_k_values_positive=2,
        max_conflict_ratio=0.10,
        min_best_k=101,
        min_exact_hits=20,
        min_confidence_score=0.50,
        low_evidence_call="observed_below_threshold",
        min_mixed_species_fraction=0.0,
    ),
    "lineage_aware": SpeciesCallPreset(
        min_unique_kmers=20,
        min_positive_sequences=5,
        min_k_values_positive=2,
        max_conflict_ratio=0.10,
        min_best_k=101,
        min_exact_hits=20,
        min_confidence_score=0.50,
        low_evidence_call="observed_below_threshold",
        min_mixed_species_fraction=0.25,
    ),
    "strict": SpeciesCallPreset(
        min_unique_kmers=50,
        min_positive_sequences=10,
        min_k_values_positive=2,
        max_conflict_ratio=0.05,
        min_best_k=101,
        min_exact_hits=50,
        min_confidence_score=0.70,
        low_evidence_call="observed_below_threshold",
        min_mixed_species_fraction=0.50,
    ),
}


def get_species_call_preset(*, preset_name: str) -> SpeciesCallPreset:
    """Return a named species-call preset.

    Parameters
    ----------
    preset_name : str
        Preset name. Supported values are ``legacy``, ``conservative``,
        ``lineage_aware`` and ``strict``.

    Returns
    -------
    SpeciesCallPreset
        Preset thresholds.

    Raises
    ------
    ValueError
        If the preset name is unsupported.
    """
    try:
        return CALL_PRESETS[preset_name]
    except KeyError as exc:
        raise ValueError(f"Unsupported species-call preset: {preset_name}") from exc


def calculate_conflict_ratio(
    *,
    target_unique_kmers: int,
    conflicting_unique_kmers: int,
) -> float:
    """Calculate the conflicting-evidence ratio.

    Parameters
    ----------
    target_unique_kmers : int
        Unique diagnostic k-mers for the focal species.
    conflicting_unique_kmers : int
        Unique diagnostic k-mers supporting other species.

    Returns
    -------
    float
        Conflict ratio in the range 0 to 1 where possible.
    """
    denominator = target_unique_kmers + conflicting_unique_kmers
    if denominator == 0:
        return 0.0
    return conflicting_unique_kmers / denominator


def calculate_confidence_score(
    *,
    n_unique_kmers: int,
    n_positive_sequences: int,
    n_k_values_positive: int,
    best_k: int,
    conflict_ratio: float,
) -> float:
    """Calculate a bounded heuristic confidence score.

    Parameters
    ----------
    n_unique_kmers : int
        Number of unique diagnostic k-mers.
    n_positive_sequences : int
        Number of independent positive reads or contigs.
    n_k_values_positive : int
        Number of positive k values.
    best_k : int
        Longest positive k value.
    conflict_ratio : float
        Conflicting-evidence ratio used as a penalty.

    Returns
    -------
    float
        Score between 0 and 1. This is not yet a calibrated probability.
    """
    kmer_score = min(max(n_unique_kmers, 0) / 10.0, 1.0)
    sequence_score = min(max(n_positive_sequences, 0) / 5.0, 1.0)
    multi_k_score = min(max(n_k_values_positive, 0) / 3.0, 1.0)
    long_k_score = min(max(best_k, 0) / 101.0, 1.0)
    raw_score = (
        0.35 * kmer_score
        + 0.25 * sequence_score
        + 0.20 * multi_k_score
        + 0.20 * long_k_score
    )
    penalty = min(max(conflict_ratio, 0.0), 1.0)
    penalised = raw_score * (1.0 - penalty)
    return round(max(0.0, min(1.0, penalised)), 4)


def _integer_from_row(*, row: dict[str, object], key: str) -> int:
    """Read an integer value from an evidence row.

    Parameters
    ----------
    row : dict[str, object]
        Evidence row.
    key : str
        Column name.

    Returns
    -------
    int
        Parsed integer, or zero when missing/blank.
    """
    value = row.get(key, 0)
    if value in (None, ""):
        return 0
    return int(value)


def _passes_species_evidence(
    *,
    row: dict[str, object],
    min_unique_kmers: int,
    min_positive_sequences: int,
    min_k_values_positive: int,
    min_best_k: int,
    min_exact_hits: int,
) -> bool:
    """Return whether a species evidence row passes detection minima.

    Parameters
    ----------
    row : dict[str, object]
        Species evidence row.
    min_unique_kmers : int
        Minimum unique k-mers.
    min_positive_sequences : int
        Minimum independent positive sequences.
    min_k_values_positive : int
        Minimum number of positive k values.
    min_best_k : int
        Minimum longest k value supported.
    min_exact_hits : int
        Minimum number of exact k-mer hits.

    Returns
    -------
    bool
        True if all evidence thresholds are met.
    """
    return (
        _integer_from_row(row=row, key="n_unique_kmers") >= min_unique_kmers
        and _integer_from_row(row=row, key="n_positive_sequences") >= min_positive_sequences
        and _integer_from_row(row=row, key="n_k_values_positive") >= min_k_values_positive
        and _integer_from_row(row=row, key="best_k") >= min_best_k
        and _integer_from_row(row=row, key="n_exact_hits") >= min_exact_hits
    )


def _passes_relative_support(
    *,
    target_unique_kmers: int,
    second_best_unique_kmers: int,
    min_unique_kmer_margin: int,
    min_unique_kmer_ratio: float,
) -> bool:
    """Return whether focal evidence is separated from the next best species.

    Parameters
    ----------
    target_unique_kmers : int
        Unique k-mers supporting the focal species.
    second_best_unique_kmers : int
        Unique k-mers supporting the strongest alternative species.
    min_unique_kmer_margin : int
        Required absolute margin over the second-best species.
    min_unique_kmer_ratio : float
        Required focal-to-second-best ratio. Values <= 0 disable this filter.

    Returns
    -------
    bool
        True if relative support requirements are met.
    """
    if min_unique_kmer_margin > 0:
        if target_unique_kmers - second_best_unique_kmers < min_unique_kmer_margin:
            return False
    if min_unique_kmer_ratio > 0 and second_best_unique_kmers > 0:
        if target_unique_kmers / second_best_unique_kmers < min_unique_kmer_ratio:
            return False
    return True


def call_species_presence(
    *,
    evidence_records: Iterable[dict[str, object]],
    min_unique_kmers: int = 3,
    min_positive_sequences: int = 2,
    min_k_values_positive: int = 1,
    max_conflict_ratio: float = 0.10,
    allow_mixed_species: bool = True,
    min_best_k: int = 0,
    min_exact_hits: int = 0,
    min_confidence_score: float = 0.0,
    min_unique_kmer_margin: int = 0,
    min_unique_kmer_ratio: float = 0.0,
    low_evidence_call: str = "present_low_confidence",
    min_mixed_species_fraction: float = 0.0,
) -> list[dict[str, object]]:
    """Call species presence from summarised evidence.

    Parameters
    ----------
    evidence_records : iterable of dict[str, object]
        Species-level evidence records.
    min_unique_kmers : int, optional
        Minimum unique diagnostic k-mers for a reportable call.
    min_positive_sequences : int, optional
        Minimum independent reads or contigs.
    min_k_values_positive : int, optional
        Minimum positive k values.
    max_conflict_ratio : float, optional
        Maximum tolerated conflict ratio for a single-species high-confidence
        call.
    allow_mixed_species : bool, optional
        If true, multiple species that independently pass evidence thresholds
        are called as ``present_in_mixed_sample``. If false, multiple passing
        species are labelled as conflicting.
    min_best_k : int, optional
        Minimum longest k-mer length required for a reportable species call.
    min_exact_hits : int, optional
        Minimum exact-hit count required for a reportable species call.
    min_confidence_score : float, optional
        Minimum heuristic confidence score for a reportable species call.
    min_unique_kmer_margin : int, optional
        Minimum absolute margin over the next-best species.
    min_unique_kmer_ratio : float, optional
        Minimum focal-to-next-best unique-k-mer ratio. Values <= 0 disable this.
    low_evidence_call : str, optional
        Call label for observed evidence below reportable thresholds.
    min_mixed_species_fraction : float, optional
        Minimum fraction of the strongest passing species support required for
        promotion to a reportable mixed-species call. Values <= 0 preserve the
        legacy behaviour where all passing species are reported in mixed
        samples.

    Returns
    -------
    list[dict[str, object]]
        Detection-call records.

    Raises
    ------
    ValueError
        If ``low_evidence_call`` is unsupported.
    """
    if low_evidence_call not in WEAK_SIGNAL_CALLS:
        raise ValueError(
            "low_evidence_call must be one of: " + ", ".join(sorted(WEAK_SIGNAL_CALLS))
        )

    records = list(evidence_records)
    calls: list[dict[str, object]] = []
    by_sample: dict[str, list[dict[str, object]]] = {}
    for record in records:
        by_sample.setdefault(str(record["sample_id"]), []).append(record)

    for sample_id, sample_records in sorted(by_sample.items()):
        unique_by_species = {
            str(row["species_name"]): _integer_from_row(row=row, key="n_unique_kmers")
            for row in sample_records
        }
        passing_species: set[str] = set()
        intermediate_rows: list[dict[str, object]] = []

        for row in sample_records:
            species_name = str(row["species_name"])
            target_unique = _integer_from_row(row=row, key="n_unique_kmers")
            alternative_unique_values = [
                value for species, value in unique_by_species.items() if species != species_name
            ]
            second_best_unique = max(alternative_unique_values) if alternative_unique_values else 0
            conflict_unique = sum(alternative_unique_values)
            conflict_ratio = calculate_conflict_ratio(
                target_unique_kmers=target_unique,
                conflicting_unique_kmers=conflict_unique,
            )
            passes_basic_evidence = _passes_species_evidence(
                row=row,
                min_unique_kmers=min_unique_kmers,
                min_positive_sequences=min_positive_sequences,
                min_k_values_positive=min_k_values_positive,
                min_best_k=min_best_k,
                min_exact_hits=min_exact_hits,
            )
            passes_relative_support = _passes_relative_support(
                target_unique_kmers=target_unique,
                second_best_unique_kmers=second_best_unique,
                min_unique_kmer_margin=min_unique_kmer_margin,
                min_unique_kmer_ratio=min_unique_kmer_ratio,
            )

            signal_confidence_score = calculate_confidence_score(
                n_unique_kmers=target_unique,
                n_positive_sequences=_integer_from_row(row=row, key="n_positive_sequences"),
                n_k_values_positive=_integer_from_row(row=row, key="n_k_values_positive"),
                best_k=_integer_from_row(row=row, key="best_k"),
                conflict_ratio=0.0,
            )
            confidence_score = calculate_confidence_score(
                n_unique_kmers=target_unique,
                n_positive_sequences=_integer_from_row(row=row, key="n_positive_sequences"),
                n_k_values_positive=_integer_from_row(row=row, key="n_k_values_positive"),
                best_k=_integer_from_row(row=row, key="best_k"),
                conflict_ratio=conflict_ratio,
            )
            passes_evidence = (
                passes_basic_evidence
                and passes_relative_support
                and signal_confidence_score >= min_confidence_score
            )
            if passes_evidence:
                passing_species.add(species_name)
            intermediate_rows.append(
                {
                    "row": row,
                    "species_name": species_name,
                    "target_unique": target_unique,
                    "conflict_unique": conflict_unique,
                    "conflict_ratio": conflict_ratio,
                    "confidence_score": confidence_score,
                    "signal_confidence_score": signal_confidence_score,
                    "passes_evidence": passes_evidence,
                }
            )

        passing_items = [
            item for item in intermediate_rows if bool(item["passes_evidence"])
        ]
        strongest_unique = max(
            (int(item["target_unique"]) for item in passing_items),
            default=0,
        )
        reportable_species: set[str] = set(passing_species)
        lineage_neighbour_species: set[str] = set()
        if (
            allow_mixed_species
            and len(passing_species) > 1
            and min_mixed_species_fraction > 0
            and strongest_unique > 0
        ):
            reportable_species = {
                str(item["species_name"])
                for item in passing_items
                if int(item["target_unique"])
                >= strongest_unique * min_mixed_species_fraction
            }
            lineage_neighbour_species = passing_species - reportable_species

        sample_is_mixed = allow_mixed_species and len(reportable_species) > 1
        for item in intermediate_rows:
            row = item["row"]
            species_name = str(item["species_name"])
            target_unique = int(item["target_unique"])
            conflict_ratio = float(item["conflict_ratio"])
            passes_evidence = bool(item["passes_evidence"])
            confidence_score = float(item["confidence_score"])
            signal_confidence_score = float(item["signal_confidence_score"])
            reportable_fraction = 0.0
            if strongest_unique > 0 and passes_evidence:
                reportable_fraction = target_unique / strongest_unique

            is_reportable_species = species_name in reportable_species
            is_neighbour_lineage = species_name in lineage_neighbour_species
            reportable_conflict_unique = sum(
                int(other["target_unique"])
                for other in passing_items
                if str(other["species_name"]) != species_name
                and str(other["species_name"]) in reportable_species
            )
            reportable_conflict_ratio = calculate_conflict_ratio(
                target_unique_kmers=target_unique,
                conflicting_unique_kmers=reportable_conflict_unique,
            )

            if passes_evidence and is_neighbour_lineage:
                call = "neighbour_lineage_evidence"
            elif passes_evidence and sample_is_mixed and is_reportable_species:
                confidence_score = calculate_confidence_score(
                    n_unique_kmers=target_unique,
                    n_positive_sequences=_integer_from_row(
                        row=row,
                        key="n_positive_sequences",
                    ),
                    n_k_values_positive=_integer_from_row(
                        row=row,
                        key="n_k_values_positive",
                    ),
                    best_k=_integer_from_row(row=row, key="best_k"),
                    conflict_ratio=0.0,
                )
                call = "present_in_mixed_sample"
            elif passes_evidence and is_reportable_species:
                if reportable_conflict_ratio <= max_conflict_ratio:
                    call = "present_high_confidence"
                else:
                    call = "ambiguous_conflicting_signal"
            elif passes_evidence:
                call = "neighbour_lineage_evidence"
            elif target_unique > 0:
                call = low_evidence_call
            else:
                call = "not_detected"

            calls.append(
                {
                    **row,
                    "conflicting_unique_kmers": int(item["conflict_unique"]),
                    "conflict_ratio": round(conflict_ratio, 4),
                    "reportable_conflicting_unique_kmers": int(reportable_conflict_unique),
                    "reportable_conflict_ratio": round(reportable_conflict_ratio, 4),
                    "mixed_species_support_fraction": round(reportable_fraction, 4),
                    "confidence_score": round(confidence_score, 4),
                    "signal_confidence_score": round(signal_confidence_score, 4),
                    "call": call,
                }
            )
    return calls


def apply_species_call_preset(
    *,
    preset_name: str,
    min_unique_kmers: int | None = None,
    min_positive_sequences: int | None = None,
    min_k_values_positive: int | None = None,
    max_conflict_ratio: float | None = None,
    min_best_k: int | None = None,
    min_exact_hits: int | None = None,
    min_confidence_score: float | None = None,
    low_evidence_call: Literal[
        "present_low_confidence", "observed_below_threshold"
    ] | None = None,
    min_mixed_species_fraction: float | None = None,
) -> dict[str, object]:
    """Resolve species-call thresholds from a preset and optional overrides.

    Parameters
    ----------
    preset_name : str
        Preset name.
    min_unique_kmers : int or None, optional
        Optional unique-k-mer override.
    min_positive_sequences : int or None, optional
        Optional positive-sequence override.
    min_k_values_positive : int or None, optional
        Optional multi-k override.
    max_conflict_ratio : float or None, optional
        Optional conflict-ratio override.
    min_best_k : int or None, optional
        Optional longest-k override.
    min_exact_hits : int or None, optional
        Optional exact-hit override.
    min_confidence_score : float or None, optional
        Optional confidence-score override.
    low_evidence_call : str or None, optional
        Optional low-evidence call label override.
    min_mixed_species_fraction : float or None, optional
        Optional mixed-species co-dominance fraction override.

    Returns
    -------
    dict[str, object]
        Resolved call-threshold settings.
    """
    preset = get_species_call_preset(preset_name=preset_name)
    return {
        "min_unique_kmers": preset.min_unique_kmers if min_unique_kmers is None else min_unique_kmers,
        "min_positive_sequences": (
            preset.min_positive_sequences
            if min_positive_sequences is None
            else min_positive_sequences
        ),
        "min_k_values_positive": (
            preset.min_k_values_positive
            if min_k_values_positive is None
            else min_k_values_positive
        ),
        "max_conflict_ratio": (
            preset.max_conflict_ratio if max_conflict_ratio is None else max_conflict_ratio
        ),
        "min_best_k": preset.min_best_k if min_best_k is None else min_best_k,
        "min_exact_hits": preset.min_exact_hits if min_exact_hits is None else min_exact_hits,
        "min_confidence_score": (
            preset.min_confidence_score
            if min_confidence_score is None
            else min_confidence_score
        ),
        "low_evidence_call": (
            preset.low_evidence_call if low_evidence_call is None else low_evidence_call
        ),
        "min_mixed_species_fraction": (
            preset.min_mixed_species_fraction
            if min_mixed_species_fraction is None
            else min_mixed_species_fraction
        ),
    }
