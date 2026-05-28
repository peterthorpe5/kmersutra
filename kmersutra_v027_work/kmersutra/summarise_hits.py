"""Summarise KmerSutra diagnostic k-mer hits."""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterable
from pathlib import Path

from kmersutra.io import read_tsv
from kmersutra.screen_reads import KmerHit


SPECIES_EVIDENCE_FIELDNAMES = [
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
]


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


def load_panel_species_metadata(*, panel_path: str | Path) -> list[dict[str, str]]:
    """Load species labels represented in a KmerSutra panel.

    Parameters
    ----------
    panel_path : str or pathlib.Path
        KmerSutra diagnostic k-mer panel TSV or TSV.GZ.

    Returns
    -------
    list[dict[str, str]]
        One record per species with ``species_name`` and ``clade`` fields.
    """
    records = read_tsv(input_path=panel_path)
    species_to_clade: dict[str, str] = {}
    for record in records:
        if record.get("panel_type") != "species_unique":
            continue
        species_name = record.get("species_name", "")
        if not species_name:
            continue
        species_to_clade.setdefault(species_name, record.get("clade", ""))
    return [
        {"species_name": species_name, "clade": species_to_clade[species_name]}
        for species_name in sorted(species_to_clade)
    ]


def complete_sample_species_evidence(
    *,
    evidence_records: Iterable[dict[str, object]],
    expected_species: Iterable[dict[str, str]],
    sample_id: str,
) -> list[dict[str, object]]:
    """Add explicit zero-evidence rows for expected species not observed.

    Parameters
    ----------
    evidence_records : iterable of dict[str, object]
        Observed sample/species evidence records.
    expected_species : iterable of dict[str, str]
        Expected species metadata, usually from the diagnostic panel.
    sample_id : str
        Sample identifier to use for zero-evidence records.

    Returns
    -------
    list[dict[str, object]]
        Completed evidence table with one row per expected species.
    """
    completed = [dict(record) for record in evidence_records]
    observed = {str(record.get("species_name", "")) for record in completed}

    for species_record in expected_species:
        species_name = str(species_record.get("species_name", ""))
        if not species_name or species_name in observed:
            continue
        completed.append(
            {
                "sample_id": sample_id,
                "species_name": species_name,
                "clade": species_record.get("clade", ""),
                "n_hits": 0,
                "n_unique_kmers": 0,
                "n_positive_sequences": 0,
                "n_k_values_positive": 0,
                "best_k": 0,
                "n_exact_hits": 0,
                "n_fuzzy_hits": 0,
            }
        )

    completed.sort(key=lambda row: (str(row.get("sample_id", "")), str(row.get("species_name", ""))))
    return completed

TAXONOMIC_EVIDENCE_FIELDNAMES = [
    "sample_id",
    "evidence_rank",
    "evidence_name",
    "evidence_taxid",
    "panel_type",
    "clade",
    "n_hits",
    "n_unique_kmers",
    "n_positive_sequences",
    "n_k_values_positive",
    "best_k",
    "n_exact_hits",
    "n_fuzzy_hits",
]


def _taxonomic_label_from_hit(*, hit: KmerHit) -> tuple[str, str, str]:
    """Return evidence-rank, evidence-name and taxid labels for one hit.

    Parameters
    ----------
    hit : KmerHit
        Diagnostic k-mer hit.

    Returns
    -------
    tuple[str, str, str]
        Evidence rank, evidence name and evidence taxid.
    """
    rank = hit.evidence_rank or ("species" if hit.species_name else "clade")
    name = hit.evidence_name or hit.species_name or hit.clade
    taxid = hit.evidence_taxid
    return rank, name, taxid


def summarise_taxonomic_hits(*, hits: Iterable[KmerHit]) -> list[dict[str, object]]:
    """Summarise hits by sample, evidence rank/name and k value.

    This retains genus, family or broader evidence that is intentionally not
    forced into a species-level call.

    Parameters
    ----------
    hits : iterable of KmerHit
        Diagnostic k-mer hits.

    Returns
    -------
    list[dict[str, object]]
        K-specific taxonomic evidence records.
    """
    grouped: dict[tuple[str, str, str, str, str, str, int], list[KmerHit]] = defaultdict(list)
    for hit in hits:
        evidence_rank, evidence_name, evidence_taxid = _taxonomic_label_from_hit(hit=hit)
        if not evidence_name:
            continue
        key = (
            hit.sample_id,
            evidence_rank,
            evidence_name,
            evidence_taxid,
            hit.panel_type,
            hit.clade,
            hit.k,
        )
        grouped[key].append(hit)

    records: list[dict[str, object]] = []
    for (
        sample_id,
        evidence_rank,
        evidence_name,
        evidence_taxid,
        panel_type,
        clade,
        k,
    ), group in sorted(grouped.items()):
        records.append(
            {
                "sample_id": sample_id,
                "evidence_rank": evidence_rank,
                "evidence_name": evidence_name,
                "evidence_taxid": evidence_taxid,
                "panel_type": panel_type,
                "clade": clade,
                "k": k,
                "n_hits": len(group),
                "n_unique_kmers": len({hit.matched_kmer for hit in group}),
                "n_positive_sequences": len({hit.sequence_id for hit in group}),
                "n_exact_hits": sum(hit.mismatches == 0 for hit in group),
                "n_fuzzy_hits": sum(hit.mismatches > 0 for hit in group),
            }
        )
    return records


def summarise_sample_taxonomic_evidence(
    *,
    taxonomic_summary: Iterable[dict[str, object]],
) -> list[dict[str, object]]:
    """Collapse k-specific taxonomic records into evidence-level records.

    Parameters
    ----------
    taxonomic_summary : iterable of dict[str, object]
        K-specific records from :func:`summarise_taxonomic_hits`.

    Returns
    -------
    list[dict[str, object]]
        Collapsed taxonomic evidence table.
    """
    grouped: dict[tuple[str, str, str, str, str, str], list[dict[str, object]]] = defaultdict(list)
    for row in taxonomic_summary:
        key = (
            str(row.get("sample_id", "")),
            str(row.get("evidence_rank", "")),
            str(row.get("evidence_name", "")),
            str(row.get("evidence_taxid", "")),
            str(row.get("panel_type", "")),
            str(row.get("clade", "")),
        )
        grouped[key].append(row)

    records: list[dict[str, object]] = []
    for (
        sample_id,
        evidence_rank,
        evidence_name,
        evidence_taxid,
        panel_type,
        clade,
    ), rows in sorted(grouped.items()):
        records.append(
            {
                "sample_id": sample_id,
                "evidence_rank": evidence_rank,
                "evidence_name": evidence_name,
                "evidence_taxid": evidence_taxid,
                "panel_type": panel_type,
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
