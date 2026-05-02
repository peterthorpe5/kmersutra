"""Feature extraction for KmerSutra-ML workflows.

This module converts diagnostic k-mer hits into interpretable feature rows
that can be used by rule-based thresholds or by the lightweight open-set
classifier in :mod:`kmersutra.ml`.
"""

from __future__ import annotations

import logging
from collections import Counter, defaultdict
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

from kmersutra.io import read_tsv, write_tsv
from kmersutra.screen_reads import KmerHit


METADATA_COLUMNS = {
    "sample_id",
    "sequence_id",
    "sequence_type",
    "true_label",
    "true_species",
    "true_clade",
    "source_genome",
    "split",
}


@dataclass(frozen=True)
class SequenceFeatureRecord:
    """KmerSutra feature summary for one read or contig.

    Attributes
    ----------
    sample_id : str
        Sample identifier.
    sequence_id : str
        Read or contig identifier.
    sequence_type : str
        Sequence type, usually ``read`` or ``contig``.
    n_total_hits : int
        Total diagnostic k-mer hits.
    n_exact_hits : int
        Number of exact diagnostic k-mer hits.
    n_fuzzy_hits : int
        Number of fuzzy diagnostic k-mer hits.
    n_species_unique_hits : int
        Number of species-unique hit observations.
    n_clade_core_hits : int
        Number of clade-core hit observations.
    n_unique_species_kmers : int
        Number of unique species-level diagnostic k-mers.
    n_unique_clade_kmers : int
        Number of unique clade-level diagnostic k-mers.
    n_positive_k_values : int
        Number of distinct k values with at least one hit.
    longest_k : int
        Longest positive k value.
    best_species : str
        Species receiving the highest hit count.
    best_species_hits : int
        Hit count for the best species.
    second_species : str
        Species receiving the second-highest hit count.
    second_species_hits : int
        Hit count for the second species.
    species_hit_margin : int
        Difference between best and second species hit counts.
    species_hit_ratio : float
        Ratio of best species hits to best plus second species hits.
    conflict_ratio : float
        Fraction of species-level hits not supporting the best species.
    best_clade : str
        Clade receiving the highest hit count.
    best_clade_hits : int
        Hit count for the best clade.
    unresolved_clade_signal : int
        One when clade evidence exists but species evidence is weak or absent.
    """

    sample_id: str
    sequence_id: str
    sequence_type: str
    n_total_hits: int
    n_exact_hits: int
    n_fuzzy_hits: int
    n_species_unique_hits: int
    n_clade_core_hits: int
    n_unique_species_kmers: int
    n_unique_clade_kmers: int
    n_positive_k_values: int
    longest_k: int
    best_species: str
    best_species_hits: int
    second_species: str
    second_species_hits: int
    species_hit_margin: int
    species_hit_ratio: float
    conflict_ratio: float
    best_clade: str
    best_clade_hits: int
    unresolved_clade_signal: int

    def to_record(self) -> dict[str, object]:
        """Return a serialisable dictionary.

        Returns
        -------
        dict[str, object]
            Feature record.
        """
        return {
            "sample_id": self.sample_id,
            "sequence_id": self.sequence_id,
            "sequence_type": self.sequence_type,
            "n_total_hits": self.n_total_hits,
            "n_exact_hits": self.n_exact_hits,
            "n_fuzzy_hits": self.n_fuzzy_hits,
            "n_species_unique_hits": self.n_species_unique_hits,
            "n_clade_core_hits": self.n_clade_core_hits,
            "n_unique_species_kmers": self.n_unique_species_kmers,
            "n_unique_clade_kmers": self.n_unique_clade_kmers,
            "n_positive_k_values": self.n_positive_k_values,
            "longest_k": self.longest_k,
            "best_species": self.best_species,
            "best_species_hits": self.best_species_hits,
            "second_species": self.second_species,
            "second_species_hits": self.second_species_hits,
            "species_hit_margin": self.species_hit_margin,
            "species_hit_ratio": round(self.species_hit_ratio, 6),
            "conflict_ratio": round(self.conflict_ratio, 6),
            "best_clade": self.best_clade,
            "best_clade_hits": self.best_clade_hits,
            "unresolved_clade_signal": self.unresolved_clade_signal,
        }


FEATURE_FIELDNAMES = list(SequenceFeatureRecord(
    sample_id="",
    sequence_id="",
    sequence_type="",
    n_total_hits=0,
    n_exact_hits=0,
    n_fuzzy_hits=0,
    n_species_unique_hits=0,
    n_clade_core_hits=0,
    n_unique_species_kmers=0,
    n_unique_clade_kmers=0,
    n_positive_k_values=0,
    longest_k=0,
    best_species="",
    best_species_hits=0,
    second_species="",
    second_species_hits=0,
    species_hit_margin=0,
    species_hit_ratio=0.0,
    conflict_ratio=0.0,
    best_clade="",
    best_clade_hits=0,
    unresolved_clade_signal=0,
).to_record().keys())


def hit_from_record(record: dict[str, str]) -> KmerHit:
    """Convert a TSV hit row into a :class:`KmerHit`.

    Parameters
    ----------
    record : dict[str, str]
        Parsed row from ``read_level_species_kmer_hits.tsv.gz``.

    Returns
    -------
    KmerHit
        Hit object.
    """
    return KmerHit(
        sample_id=record.get("sample_id", ""),
        sequence_id=record.get("sequence_id", ""),
        sequence_type=record.get("sequence_type", "read"),
        k=int(record.get("k", "0") or 0),
        query_position=int(record.get("query_position", "0") or 0),
        matched_kmer=record.get("matched_kmer", ""),
        query_kmer=record.get("query_kmer", ""),
        mismatches=int(record.get("mismatches", "0") or 0),
        panel_type=record.get("panel_type", ""),
        species_name=record.get("species_name", ""),
        clade=record.get("clade", ""),
    )


def _best_and_second(counter: Counter[str]) -> tuple[str, int, str, int]:
    """Return the top two labels and counts from a counter.

    Parameters
    ----------
    counter : collections.Counter[str]
        Label counts.

    Returns
    -------
    tuple[str, int, str, int]
        Best label, best count, second label, second count.
    """
    if not counter:
        return "", 0, "", 0
    ranked = sorted(counter.items(), key=lambda item: (-item[1], item[0]))
    best_label, best_count = ranked[0]
    if len(ranked) == 1:
        return best_label, best_count, "", 0
    second_label, second_count = ranked[1]
    return best_label, best_count, second_label, second_count


def summarise_sequence_features(
    *,
    hits: Iterable[KmerHit],
    logger: logging.Logger | None = None,
) -> list[dict[str, object]]:
    """Summarise diagnostic hit evidence per read or contig.

    Parameters
    ----------
    hits : iterable of KmerHit
        Diagnostic k-mer hits.
    logger : logging.Logger | None, optional
        Logger for progress messages.

    Returns
    -------
    list[dict[str, object]]
        Per-sequence feature records.
    """
    grouped: dict[tuple[str, str, str], list[KmerHit]] = defaultdict(list)
    for hit in hits:
        grouped[(hit.sample_id, hit.sequence_id, hit.sequence_type)].append(hit)

    if logger:
        logger.info("Summarising KmerSutra-ML features for %d positive sequences", len(grouped))

    features: list[dict[str, object]] = []
    for (sample_id, sequence_id, sequence_type), group in sorted(grouped.items()):
        species_hits = Counter(
            hit.species_name for hit in group
            if hit.panel_type == "species_unique" and hit.species_name
        )
        clade_hits = Counter(hit.clade for hit in group if hit.clade)
        best_species, best_species_hits, second_species, second_species_hits = _best_and_second(species_hits)
        best_clade, best_clade_hits, _, _ = _best_and_second(clade_hits)
        n_species_hits = sum(species_hits.values())
        non_best_species_hits = n_species_hits - best_species_hits
        species_denominator = n_species_hits or 0
        conflict_ratio = (
            non_best_species_hits / species_denominator
            if species_denominator > 0
            else 0.0
        )
        ratio_denominator = best_species_hits + second_species_hits
        species_hit_ratio = (
            best_species_hits / ratio_denominator
            if ratio_denominator > 0
            else 0.0
        )
        n_clade_core_hits = sum(hit.panel_type == "clade_core" for hit in group)
        n_species_unique_hits = sum(hit.panel_type == "species_unique" for hit in group)
        unresolved_clade_signal = int(n_clade_core_hits > 0 and best_species_hits == 0)

        record = SequenceFeatureRecord(
            sample_id=sample_id,
            sequence_id=sequence_id,
            sequence_type=sequence_type,
            n_total_hits=len(group),
            n_exact_hits=sum(hit.mismatches == 0 for hit in group),
            n_fuzzy_hits=sum(hit.mismatches > 0 for hit in group),
            n_species_unique_hits=n_species_unique_hits,
            n_clade_core_hits=n_clade_core_hits,
            n_unique_species_kmers=len(
                {
                    hit.matched_kmer
                    for hit in group
                    if hit.panel_type == "species_unique"
                }
            ),
            n_unique_clade_kmers=len(
                {
                    hit.matched_kmer
                    for hit in group
                    if hit.panel_type == "clade_core"
                }
            ),
            n_positive_k_values=len({hit.k for hit in group}),
            longest_k=max(hit.k for hit in group),
            best_species=best_species,
            best_species_hits=best_species_hits,
            second_species=second_species,
            second_species_hits=second_species_hits,
            species_hit_margin=best_species_hits - second_species_hits,
            species_hit_ratio=species_hit_ratio,
            conflict_ratio=conflict_ratio,
            best_clade=best_clade,
            best_clade_hits=best_clade_hits,
            unresolved_clade_signal=unresolved_clade_signal,
        )
        features.append(record.to_record())

    return features


def load_hits_as_features(
    *,
    hits_tsv: str | Path,
    logger: logging.Logger | None = None,
) -> list[dict[str, object]]:
    """Load read-level hit TSV and return feature records.

    Parameters
    ----------
    hits_tsv : str or pathlib.Path
        KmerSutra read-level hit table.
    logger : logging.Logger | None, optional
        Logger for progress messages.

    Returns
    -------
    list[dict[str, object]]
        Per-sequence feature records.
    """
    records = read_tsv(input_path=hits_tsv)
    hits = [hit_from_record(record) for record in records]
    if logger:
        logger.info("Loaded %d read-level hit rows from %s", len(hits), hits_tsv)
    return summarise_sequence_features(hits=hits, logger=logger)


def write_sequence_features(
    *,
    features: Iterable[dict[str, object]],
    output_path: str | Path,
) -> None:
    """Write sequence-level features to TSV.

    Parameters
    ----------
    features : iterable of dict[str, object]
        Feature records.
    output_path : str or pathlib.Path
        Output TSV path.
    """
    write_tsv(records=features, output_path=output_path, fieldnames=FEATURE_FIELDNAMES)


def infer_numeric_feature_columns(
    *,
    records: Iterable[dict[str, object]],
    label_column: str | None = None,
) -> list[str]:
    """Infer numeric feature columns from feature records.

    Parameters
    ----------
    records : iterable of dict[str, object]
        Feature table rows.
    label_column : str | None, optional
        Label column to exclude.

    Returns
    -------
    list[str]
        Numeric feature column names.
    """
    rows = list(records)
    if not rows:
        return []
    excluded = set(METADATA_COLUMNS)
    if label_column:
        excluded.add(label_column)
    columns = sorted(set().union(*(row.keys() for row in rows)) - excluded)
    numeric_columns: list[str] = []
    for column in columns:
        all_numeric = True
        for row in rows:
            value = row.get(column, "")
            if value in {"", None}:
                continue
            try:
                float(value)
            except (TypeError, ValueError):
                all_numeric = False
                break
        if all_numeric:
            numeric_columns.append(column)
    return numeric_columns
