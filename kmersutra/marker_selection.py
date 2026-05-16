"""Genome-aware diagnostic marker selection for KmerSutra panels.

This module contains deterministic selection helpers used when a diagnostic
panel has more valid k-mers than should be retained for a taxon/k bucket. The
main aim is to avoid positional bias caused by retaining dense, overlapping
sliding-window k-mers from one small genomic region.
"""

from __future__ import annotations

import hashlib
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
    """
    config.validate()
    if config.strategy != "genome_spread":
        raise ValueError("select_genome_spread_markers requires genome_spread strategy")

    states: dict[tuple[str, str, str, str, int], _BucketState] = defaultdict(_BucketState)
    for item in diagnostic_kmers:
        key = diagnostic_retention_key(item)
        candidate = _SelectedMarker(
            item=item,
            score=marker_score(item),
            bin_key=genome_bin_key(item=item, genome_bin_size=config.genome_bin_size),
        )
        _add_or_replace_genome_spread_marker(
            state=states[key],
            candidate=candidate,
            config=config,
        )

    selected_markers: list[DiagnosticKmer] = []
    for key in sorted(states):
        state = states[key]
        selected = sorted(
            state.selected,
            key=lambda marker: (
                marker.item.k,
                marker.item.evidence_rank,
                marker.item.evidence_name,
                marker.bin_key[0],
                marker.bin_key[1],
                marker.bin_key[2],
                int(marker.item.example_position),
                marker.item.kmer,
            ),
        )
        selected_markers.extend(marker.item for marker in selected)

    for item in selected_markers:
        yield item
