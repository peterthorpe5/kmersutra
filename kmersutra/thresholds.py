"""Detection-call and confidence scoring logic for KmerSutra."""

from __future__ import annotations

from collections.abc import Iterable


def calculate_conflict_ratio(
    *,
    target_unique_kmers: int,
    conflicting_unique_kmers: int,
) -> float:
    """Calculate the conflicting-evidence ratio.

    Parameters
    ----------
    target_unique_kmers : int
        Unique diagnostic k-mers for the strongest species.
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
        Conflicting-evidence ratio.

    Returns
    -------
    float
        Score between 0 and 1. This is not yet a calibrated probability.
    """
    kmer_score = min(n_unique_kmers / 10.0, 1.0)
    sequence_score = min(n_positive_sequences / 5.0, 1.0)
    multi_k_score = min(n_k_values_positive / 3.0, 1.0)
    long_k_score = min(best_k / 101.0, 1.0)
    raw_score = (
        0.35 * kmer_score
        + 0.25 * sequence_score
        + 0.20 * multi_k_score
        + 0.20 * long_k_score
    )
    penalised = raw_score * (1.0 - min(max(conflict_ratio, 0.0), 1.0))
    return round(max(0.0, min(1.0, penalised)), 4)


def call_species_presence(
    *,
    evidence_records: Iterable[dict[str, object]],
    min_unique_kmers: int = 3,
    min_positive_sequences: int = 2,
    min_k_values_positive: int = 1,
    max_conflict_ratio: float = 0.10,
) -> list[dict[str, object]]:
    """Call species presence from summarised evidence.

    Parameters
    ----------
    evidence_records : iterable of dict[str, object]
        Species-level evidence records.
    min_unique_kmers : int, optional
        Minimum unique diagnostic k-mers for a high-confidence call.
    min_positive_sequences : int, optional
        Minimum independent reads or contigs.
    min_k_values_positive : int, optional
        Minimum positive k values.
    max_conflict_ratio : float, optional
        Maximum tolerated conflict ratio.

    Returns
    -------
    list[dict[str, object]]
        Detection-call records.
    """
    records = list(evidence_records)
    calls: list[dict[str, object]] = []
    by_sample: dict[str, list[dict[str, object]]] = {}
    for record in records:
        by_sample.setdefault(str(record["sample_id"]), []).append(record)

    for sample_id, sample_records in sorted(by_sample.items()):
        unique_by_species = {
            str(row["species_name"]): int(row["n_unique_kmers"])
            for row in sample_records
        }
        for row in sample_records:
            species_name = str(row["species_name"])
            target_unique = int(row["n_unique_kmers"])
            conflict_unique = sum(
                value for species, value in unique_by_species.items() if species != species_name
            )
            conflict_ratio = calculate_conflict_ratio(
                target_unique_kmers=target_unique,
                conflicting_unique_kmers=conflict_unique,
            )
            confidence_score = calculate_confidence_score(
                n_unique_kmers=target_unique,
                n_positive_sequences=int(row["n_positive_sequences"]),
                n_k_values_positive=int(row["n_k_values_positive"]),
                best_k=int(row["best_k"]),
                conflict_ratio=conflict_ratio,
            )
            passes_evidence = (
                target_unique >= min_unique_kmers
                and int(row["n_positive_sequences"]) >= min_positive_sequences
                and int(row["n_k_values_positive"]) >= min_k_values_positive
            )
            if passes_evidence and conflict_ratio <= max_conflict_ratio:
                call = "present_high_confidence"
            elif passes_evidence and conflict_ratio > max_conflict_ratio:
                call = "ambiguous_mixed_signal"
            elif target_unique > 0:
                call = "present_low_confidence"
            else:
                call = "not_detected"

            calls.append(
                {
                    **row,
                    "conflicting_unique_kmers": conflict_unique,
                    "conflict_ratio": round(conflict_ratio, 4),
                    "confidence_score": confidence_score,
                    "call": call,
                }
            )
    return calls
