"""Summarise KmerSutra diagnostic k-mer hits."""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterable

from kmersutra.screen_reads import KmerHit


def summarise_species_hits(*, hits: Iterable[KmerHit]) -> list[dict[str, object]]:
    """Summarise hits by sample, species or clade, panel type and k.

    Parameters
    ----------
    hits : iterable of KmerHit
        Diagnostic k-mer hits.

    Returns
    -------
    list[dict[str, object]]
        Summary records.
    """
    grouped: dict[tuple[str, str, str, str, int], list[KmerHit]] = defaultdict(list)
    for hit in hits:
        label = hit.species_name or hit.clade
        key = (hit.sample_id, hit.panel_type, label, hit.clade, hit.k)
        grouped[key].append(hit)

    records: list[dict[str, object]] = []
    for (sample_id, panel_type, label, clade, k), group in sorted(grouped.items()):
        records.append(
            {
                "sample_id": sample_id,
                "panel_type": panel_type,
                "label": label,
                "clade": clade,
                "k": k,
                "n_hits": len(group),
                "n_unique_kmers": len({hit.matched_kmer for hit in group}),
                "n_positive_sequences": len({hit.sequence_id for hit in group}),
                "n_exact_hits": sum(hit.mismatches == 0 for hit in group),
                "n_fuzzy_hits": sum(hit.mismatches > 0 for hit in group),
                "min_mismatches": min(hit.mismatches for hit in group),
                "max_mismatches": max(hit.mismatches for hit in group),
            }
        )
    return records


def summarise_sample_species_evidence(
    *,
    species_summary: Iterable[dict[str, object]],
) -> list[dict[str, object]]:
    """Collapse k-specific summary rows into species-level evidence.

    Parameters
    ----------
    species_summary : iterable of dict[str, object]
        Output from :func:`summarise_species_hits`.

    Returns
    -------
    list[dict[str, object]]
        Collapsed sample/species evidence records.
    """
    grouped: dict[tuple[str, str, str], list[dict[str, object]]] = defaultdict(list)
    for row in species_summary:
        if row["panel_type"] != "species_unique":
            continue
        key = (str(row["sample_id"]), str(row["label"]), str(row["clade"]))
        grouped[key].append(row)

    records: list[dict[str, object]] = []
    for (sample_id, species_name, clade), rows in sorted(grouped.items()):
        records.append(
            {
                "sample_id": sample_id,
                "species_name": species_name,
                "clade": clade,
                "n_hits": sum(int(row["n_hits"]) for row in rows),
                "n_unique_kmers": sum(int(row["n_unique_kmers"]) for row in rows),
                "n_positive_sequences": max(int(row["n_positive_sequences"]) for row in rows),
                "n_k_values_positive": len({int(row["k"]) for row in rows}),
                "best_k": max(int(row["k"]) for row in rows),
                "n_exact_hits": sum(int(row["n_exact_hits"]) for row in rows),
                "n_fuzzy_hits": sum(int(row["n_fuzzy_hits"]) for row in rows),
            }
        )
    return records
