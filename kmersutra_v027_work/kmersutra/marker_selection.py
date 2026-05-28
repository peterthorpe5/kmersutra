"""Genome-aware diagnostic marker selection for KmerSutra panels.

This module contains deterministic selection helpers used when a diagnostic
panel has more valid k-mers than should be retained for a taxon/k bucket. The
main aim is to avoid positional bias caused by retaining dense, overlapping
sliding-window k-mers from one small genomic region.
"""

from __future__ import annotations

import hashlib
import heapq
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Iterable, Iterator

from kmersutra.build_panel import DiagnosticKmer


@dataclass(frozen=True)
class MarkerSelectionConfig:
    """Configuration for marker thinning/selection.

    Attributes
    ----------
    strategy : str
        Selection strategy. Supported values are ``first_seen`` and
        ``genome_spread``.
    max_per_bucket : int | None
        Maximum retained markers per evidence bucket/k value. ``None`` means
        no overall cap.
    genome_bin_size : int
        Number of reference bases per positional bin for genome-spread mode.
    max_per_genome_bin : int
        Maximum retained markers from one genome/contig/bin within an evidence
        bucket.
    """

    strategy: str = "first_seen"
    max_per_bucket: int | None = None
    genome_bin_size: int = 10000
    max_per_genome_bin: int = 10

    def validate(self) -> None:
        """Validate marker-selection settings.

        Raises
        ------
        ValueError
            If a supplied setting is unsupported or non-positive where a
            positive value is required.
        """
        if self.strategy not in {"first_seen", "genome_spread"}:
            raise ValueError("marker_selection must be 'first_seen' or 'genome_spread'")
        if self.max_per_bucket is not None and self.max_per_bucket <= 0:
            raise ValueError("max_per_bucket must be positive when supplied")
        if self.genome_bin_size <= 0:
            raise ValueError("genome_bin_size must be positive")
        if self.max_per_genome_bin <= 0:
            raise ValueError("max_per_genome_bin must be positive")


@dataclass
class _SelectedMarker:
    """Internal selected marker record."""

    item: DiagnosticKmer
    score: int
    bin_key: tuple[str, str, int]


@dataclass
class _BucketState:
    """Internal marker-selection state for one evidence bucket."""

    selected: list[_SelectedMarker] = field(default_factory=list)
    bin_counts: dict[tuple[str, str, int], int] = field(
        default_factory=lambda: defaultdict(int)
    )


def diagnostic_retention_key(item: DiagnosticKmer) -> tuple[str, str, str, str, int]:
    """Return the evidence bucket key for one diagnostic k-mer.

    Parameters
    ----------
    item : DiagnosticKmer
        Diagnostic k-mer record.

    Returns
    -------
    tuple[str, str, str, str, int]
        Key based on evidence class, taxon and k value.
    """
    return (
        item.panel_type,
        item.evidence_taxid,
        item.species_name,
        item.clade,
        int(item.k),
    )


def first_semicolon_value(value: str) -> str:
    """Return the first non-empty item in a semicolon-separated field.

    Parameters
    ----------
    value : str
        Semicolon-separated field value.

    Returns
    -------
    str
        First non-empty value, or an empty string.
    """
    for item in str(value or "").split(";"):
        if item:
            return item
    return ""


def genome_bin_key(
    *,
    item: DiagnosticKmer,
    genome_bin_size: int,
) -> tuple[str, str, int]:
    """Return a positional genome-bin key for one diagnostic k-mer.

    Parameters
    ----------
    item : DiagnosticKmer
        Diagnostic k-mer record.
    genome_bin_size : int
        Number of reference bases per bin.

    Returns
    -------
    tuple[str, str, int]
        Source-genome, source-contig and positional-bin key.
    """
    if genome_bin_size <= 0:
        raise ValueError("genome_bin_size must be positive")
    genome_id = first_semicolon_value(item.source_genomes)
    contig_id = first_semicolon_value(item.source_contigs)
    position = max(0, int(item.example_position))
    return (genome_id, contig_id, position // genome_bin_size)


def marker_score(item: DiagnosticKmer) -> int:
    """Return a deterministic pseudo-random score for one marker.

    Lower scores are preferred. The score is based on stable marker metadata and
    is independent of input order, making genome-spread selection reproducible.

    Parameters
    ----------
    item : DiagnosticKmer
        Diagnostic k-mer record.

    Returns
    -------
    int
        Deterministic integer score.
    """
    payload = "\t".join(
        [
            str(item.k),
            item.kmer,
            item.panel_type,
            item.evidence_taxid,
            item.evidence_rank,
            item.species_name,
            item.clade,
            item.source_genomes,
            item.source_contigs,
            str(item.example_position),
        ]
    )
    digest = hashlib.sha256(payload.encode("utf-8")).hexdigest()
    return int(digest[:16], 16)


def _find_worst_index(selected: list[_SelectedMarker]) -> int:
    """Return the index of the least desirable selected marker.

    Parameters
    ----------
    selected : list[_SelectedMarker]
        Selected marker records.

    Returns
    -------
    int
        Index of the marker with the highest score.
    """
    if not selected:
        raise ValueError("Cannot find worst marker from an empty selection")
    return max(range(len(selected)), key=lambda index: selected[index].score)


def _find_worst_index_in_bin(
    selected: list[_SelectedMarker],
    bin_key: tuple[str, str, int],
) -> int | None:
    """Return the least desirable selected marker index for one bin.

    Parameters
    ----------
    selected : list[_SelectedMarker]
        Selected marker records.
    bin_key : tuple[str, str, int]
        Genome-bin key.

    Returns
    -------
    int or None
        Matching marker index, or ``None`` if the bin is absent.
    """
    candidate_indices = [
        index for index, selected_item in enumerate(selected)
        if selected_item.bin_key == bin_key
    ]
    if not candidate_indices:
        return None
    return max(candidate_indices, key=lambda index: selected[index].score)


def _replace_marker(
    *,
    state: _BucketState,
    index: int,
    candidate: _SelectedMarker,
) -> None:
    """Replace one selected marker and update bin counts.

    Parameters
    ----------
    state : _BucketState
        Bucket state to update.
    index : int
        Selected marker index to replace.
    candidate : _SelectedMarker
        Replacement marker.
    """
    previous = state.selected[index]
    state.bin_counts[previous.bin_key] -= 1
    if state.bin_counts[previous.bin_key] <= 0:
        del state.bin_counts[previous.bin_key]
    state.selected[index] = candidate
    state.bin_counts[candidate.bin_key] += 1


def _add_or_replace_genome_spread_marker(
    *,
    state: _BucketState,
    candidate: _SelectedMarker,
    config: MarkerSelectionConfig,
) -> None:
    """Add or replace a marker under genome-spread constraints.

    Parameters
    ----------
    state : _BucketState
        Selection state for one evidence bucket.
    candidate : _SelectedMarker
        Candidate marker.
    config : MarkerSelectionConfig
        Validated selection configuration.
    """
    max_per_bucket = config.max_per_bucket
    bin_count = state.bin_counts.get(candidate.bin_key, 0)

    if bin_count >= config.max_per_genome_bin:
        worst_same_bin = _find_worst_index_in_bin(state.selected, candidate.bin_key)
        if worst_same_bin is not None and candidate.score < state.selected[worst_same_bin].score:
            _replace_marker(state=state, index=worst_same_bin, candidate=candidate)
        return

    if max_per_bucket is None or len(state.selected) < max_per_bucket:
        state.selected.append(candidate)
        state.bin_counts[candidate.bin_key] += 1
        return

    # At the global cap, favour adding evidence from a previously unrepresented
    # genome bin; otherwise replace only if the deterministic score improves.
    if candidate.bin_key not in state.bin_counts:
        replace_index = _find_worst_index(state.selected)
        _replace_marker(state=state, index=replace_index, candidate=candidate)
        return

    replace_index = _find_worst_index(state.selected)
    if candidate.score < state.selected[replace_index].score:
        _replace_marker(state=state, index=replace_index, candidate=candidate)


def _marker_sort_key(item: DiagnosticKmer) -> tuple[object, ...]:
    """Return a deterministic output sort key for a marker.

    Parameters
    ----------
    item : DiagnosticKmer
        Marker to sort.

    Returns
    -------
    tuple[object, ...]
        Stable sort key.
    """
    bin_key = genome_bin_key(item=item, genome_bin_size=1)
    return (
        item.k,
        item.evidence_rank,
        item.evidence_name,
        bin_key[0],
        item.source_contigs,
        int(item.example_position),
        item.kmer,
    )


def _choose_bucket_markers(
    *,
    bin_markers: dict[tuple[str, str, int], list[DiagnosticKmer]],
    config: MarkerSelectionConfig,
) -> list[DiagnosticKmer]:
    """Choose final markers for one evidence bucket.

    Parameters
    ----------
    bin_markers : dict[tuple[str, str, int], list[DiagnosticKmer]]
        Candidate markers retained within each source genome bin.
    config : MarkerSelectionConfig
        Marker-selection settings.

    Returns
    -------
    list[DiagnosticKmer]
        Selected markers for the evidence bucket.
    """
    candidates: list[DiagnosticKmer] = []
    for markers in bin_markers.values():
        candidates.extend(markers)
    if config.max_per_bucket is None or len(candidates) <= config.max_per_bucket:
        return sorted(candidates, key=_marker_sort_key)

    # Prefer broad bin coverage by taking the best marker from as many bins as
    # possible before filling remaining capacity with the next-best markers.
    best_per_bin: list[tuple[int, tuple[str, str, int], DiagnosticKmer]] = []
    remaining: list[tuple[int, tuple[str, str, int], DiagnosticKmer]] = []
    for bin_key, markers in bin_markers.items():
        ordered = sorted(markers, key=lambda item: (marker_score(item), _marker_sort_key(item)))
        if ordered:
            best_per_bin.append((marker_score(ordered[0]), bin_key, ordered[0]))
            for item in ordered[1:]:
                remaining.append((marker_score(item), bin_key, item))

    selected_entries = sorted(best_per_bin, key=lambda value: (value[0], value[1]))[
        : config.max_per_bucket
    ]
    selected = [entry[2] for entry in selected_entries]
    if len(selected) < config.max_per_bucket:
        selected_ids = {id(item) for item in selected}
        for _, _, item in sorted(remaining, key=lambda value: (value[0], value[1])):
            if id(item) in selected_ids:
                continue
            selected.append(item)
            selected_ids.add(id(item))
            if len(selected) >= config.max_per_bucket:
                break
    return sorted(selected, key=_marker_sort_key)


def select_genome_spread_markers(
    *,
    diagnostic_kmers: Iterable[DiagnosticKmer],
    config: MarkerSelectionConfig,
) -> Iterator[DiagnosticKmer]:
    """Select a deterministic genome-spread subset of diagnostic k-mers.

    Parameters
    ----------
    diagnostic_kmers : iterable of DiagnosticKmer
        Candidate diagnostic k-mers.
    config : MarkerSelectionConfig
        Marker-selection settings. ``strategy`` must be ``genome_spread``.

    Yields
    ------
    DiagnosticKmer
        Selected diagnostic k-mer records.

    Notes
    -----
    The implementation keeps a small bounded heap for each evidence-bucket and
    genome-bin pair. This avoids the earlier O(N * max_per_bucket) replacement
    scan when a bucket contained millions of candidate markers.
    """
    config.validate()
    if config.strategy != "genome_spread":
        raise ValueError("select_genome_spread_markers requires genome_spread strategy")

    # bucket -> bin -> heap of (-score, stable_order, item). The negative score
    # makes heap[0] the current worst marker in the bin, so better candidates can
    # replace it cheaply.
    states: dict[
        tuple[str, str, str, str, int],
        dict[tuple[str, str, int], list[tuple[int, int, DiagnosticKmer]]],
    ] = defaultdict(lambda: defaultdict(list))

    stable_order = 0
    for item in diagnostic_kmers:
        key = diagnostic_retention_key(item)
        bin_key = genome_bin_key(item=item, genome_bin_size=config.genome_bin_size)
        score = marker_score(item)
        heap = states[key][bin_key]
        entry = (-score, stable_order, item)
        stable_order += 1
        if len(heap) < config.max_per_genome_bin:
            heapq.heappush(heap, entry)
            continue
        if score < -heap[0][0]:
            heapq.heapreplace(heap, entry)

    for key in sorted(states):
        bin_markers: dict[tuple[str, str, int], list[DiagnosticKmer]] = {}
        for bin_key, heap in states[key].items():
            bin_markers[bin_key] = [entry[2] for entry in heap]
        for item in _choose_bucket_markers(bin_markers=bin_markers, config=config):
            yield item
