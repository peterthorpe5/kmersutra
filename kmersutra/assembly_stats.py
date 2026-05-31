"""Assembly-aware helpers for KmerSutra marker sampling.

Candidate-universe sampling should not let thousands of short contigs create
thousands of independent local bin-zero opportunities, and small viral-scale
assemblies should not be sampled as though they were large eukaryotic genomes.
These helpers calculate simple assembly statistics and choose deterministic
sampling bins based on effective assembly length.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

from kmersutra.fasta import read_fasta_records


@dataclass(frozen=True)
class AssemblyStats:
    """Summary statistics used for assembly-aware marker sampling.

    Attributes
    ----------
    assembly_id : str
        Assembly or genome identifier.
    total_length : int
        Total number of bases across all contigs.
    effective_length : int
        Total bases in contigs long enough to contain at least one requested
        k-mer at the smallest k value.
    n_contigs : int
        Number of contigs in the assembly.
    n_effective_contigs : int
        Number of contigs contributing to ``effective_length``.
    n_contigs_ge_max_k : int
        Number of contigs long enough to contain the longest requested k-mer.
    max_contig_length : int
        Length of the longest contig.
    mean_contig_length : float
        Arithmetic mean contig length.
    median_contig_length : float
        Median contig length.
    n50 : int
        Assembly N50 calculated from all contigs.
    min_k : int
        Smallest requested k-mer length.
    max_k : int
        Largest requested k-mer length.
    """

    assembly_id: str
    total_length: int
    effective_length: int
    n_contigs: int
    n_effective_contigs: int
    n_contigs_ge_max_k: int
    max_contig_length: int
    mean_contig_length: float
    median_contig_length: float
    n50: int
    min_k: int
    max_k: int


@dataclass(frozen=True)
class AssemblyBinPlan:
    """Assembly-aware binning plan for candidate marker sampling.

    Attributes
    ----------
    requested_bin_size : int
        User-requested base bin size.
    effective_bin_size : int
        Bin size actually used for sampling.
    estimated_global_bins : int
        Estimated number of bins from effective assembly length.
    use_global_assembly_bins : bool
        Whether contigs should be laid onto a cumulative assembly coordinate
        before assigning bins.
    is_small_assembly : bool
        Whether the assembly is below the configured small-genome threshold.
    is_fragmented_assembly : bool
        Whether the assembly meets the configured fragmentation heuristic.
    reason : str
        Human-readable rationale for audit tables and logs.
    """

    requested_bin_size: int
    effective_bin_size: int
    estimated_global_bins: int
    use_global_assembly_bins: bool
    is_small_assembly: bool
    is_fragmented_assembly: bool
    reason: str


def calculate_n50(*, lengths: list[int]) -> int:
    """Calculate N50 from contig lengths.

    Parameters
    ----------
    lengths : list of int
        Contig lengths.

    Returns
    -------
    int
        N50 length, or zero when no positive lengths are supplied.
    """
    positive = sorted((int(length) for length in lengths if int(length) > 0), reverse=True)
    if not positive:
        return 0
    half_total = sum(positive) / 2.0
    cumulative = 0
    for length in positive:
        cumulative += length
        if cumulative >= half_total:
            return length
    return positive[-1]


def calculate_median_length(*, lengths: list[int]) -> float:
    """Return the median contig length.

    Parameters
    ----------
    lengths : list of int
        Contig lengths.

    Returns
    -------
    float
        Median length, or zero for an empty list.
    """
    positive = sorted(int(length) for length in lengths if int(length) > 0)
    if not positive:
        return 0.0
    midpoint = len(positive) // 2
    if len(positive) % 2:
        return float(positive[midpoint])
    return (positive[midpoint - 1] + positive[midpoint]) / 2.0


def calculate_assembly_stats(
    *,
    contig_lengths: dict[str, int],
    assembly_id: str,
    k_values: list[int],
) -> AssemblyStats:
    """Calculate assembly statistics from contig lengths.

    Parameters
    ----------
    contig_lengths : dict[str, int]
        Mapping from contig identifier to contig length.
    assembly_id : str
        Assembly or genome identifier.
    k_values : list of int
        Requested k-mer lengths.

    Returns
    -------
    AssemblyStats
        Assembly summary statistics.

    Raises
    ------
    ValueError
        If no positive k value is supplied.
    """
    positive_k = sorted({int(k_value) for k_value in k_values if int(k_value) > 0})
    if not positive_k:
        raise ValueError("At least one positive k value is required")
    lengths = [max(0, int(length)) for length in contig_lengths.values()]
    min_k = positive_k[0]
    max_k = positive_k[-1]
    total_length = sum(lengths)
    effective_lengths = [length for length in lengths if length >= min_k]
    return AssemblyStats(
        assembly_id=assembly_id,
        total_length=total_length,
        effective_length=sum(effective_lengths),
        n_contigs=len(lengths),
        n_effective_contigs=len(effective_lengths),
        n_contigs_ge_max_k=sum(1 for length in lengths if length >= max_k),
        max_contig_length=max(lengths) if lengths else 0,
        mean_contig_length=(total_length / len(lengths)) if lengths else 0.0,
        median_contig_length=calculate_median_length(lengths=lengths),
        n50=calculate_n50(lengths=lengths),
        min_k=min_k,
        max_k=max_k,
    )


def collect_fasta_contig_lengths(*, fasta_path: str | Path) -> dict[str, int]:
    """Collect contig lengths from a FASTA file without storing sequences.

    Parameters
    ----------
    fasta_path : str or pathlib.Path
        FASTA file path.

    Returns
    -------
    dict[str, int]
        Mapping from contig identifier to sequence length.
    """
    lengths: dict[str, int] = {}
    for record in read_fasta_records(fasta_path=fasta_path):
        lengths[record.identifier] = len(record.sequence)
    return lengths


def build_contig_offsets(
    *,
    contig_lengths: dict[str, int],
    min_contig_length: int,
) -> dict[str, int]:
    """Assign deterministic cumulative assembly offsets to contigs.

    Parameters
    ----------
    contig_lengths : dict[str, int]
        Mapping from contig identifier to contig length.
    min_contig_length : int
        Minimum contig length retained in the effective coordinate system.

    Returns
    -------
    dict[str, int]
        Mapping from contig identifier to cumulative zero-based offset.
    """
    offsets: dict[str, int] = {}
    cumulative = 0
    for contig_id, length in contig_lengths.items():
        clean_length = max(0, int(length))
        if clean_length < min_contig_length:
            continue
        offsets[contig_id] = cumulative
        cumulative += clean_length
    return offsets


def choose_assembly_aware_bin_plan(
    *,
    stats: AssemblyStats,
    requested_bin_size: int,
    small_assembly_length: int = 250000,
    small_assembly_min_bin_size: int = 10000,
    small_assembly_target_bins: int = 25,
    fragmented_contig_count: int = 500,
    fragmented_n50_multiplier: float = 2.0,
    fragmented_max_global_bins: int = 1000,
) -> AssemblyBinPlan:
    """Choose an assembly-aware candidate-sampling bin plan.

    Parameters
    ----------
    stats : AssemblyStats
        Assembly statistics.
    requested_bin_size : int
        User-requested bin size.
    small_assembly_length : int, optional
        Assemblies at or below this effective length use small-genome logic.
    small_assembly_min_bin_size : int, optional
        Minimum bin size for small assemblies.
    small_assembly_target_bins : int, optional
        Approximate maximum number of bins targeted for small assemblies.
    fragmented_contig_count : int, optional
        Contig-count threshold for fragmented-assembly logic.
    fragmented_n50_multiplier : float, optional
        Assembly is considered fragmented when ``n50`` is less than this
        multiple of the requested bin size and contig count is high.
    fragmented_max_global_bins : int, optional
        Approximate maximum number of cumulative global bins for fragmented
        assemblies.

    Returns
    -------
    AssemblyBinPlan
        Deterministic binning plan.

    Raises
    ------
    ValueError
        If numeric parameters are invalid.
    """
    if requested_bin_size <= 0:
        raise ValueError("requested_bin_size must be positive")
    if small_assembly_target_bins <= 0:
        raise ValueError("small_assembly_target_bins must be positive")
    if fragmented_max_global_bins <= 0:
        raise ValueError("fragmented_max_global_bins must be positive")
    effective_length = max(0, int(stats.effective_length))
    is_small = 0 < effective_length <= int(small_assembly_length)
    is_fragmented = (
        stats.n_effective_contigs >= int(fragmented_contig_count)
        and stats.n50 < (float(fragmented_n50_multiplier) * requested_bin_size)
    )
    effective_bin_size = int(requested_bin_size)
    reasons = ["requested_bin_size"]
    if is_small:
        small_target_size = math.ceil(effective_length / small_assembly_target_bins)
        effective_bin_size = max(
            effective_bin_size,
            int(small_assembly_min_bin_size),
            int(small_target_size),
        )
        reasons.append("small_assembly")
    if is_fragmented:
        fragmented_target_size = math.ceil(effective_length / fragmented_max_global_bins)
        effective_bin_size = max(effective_bin_size, int(fragmented_target_size))
        reasons.append("fragmented_assembly")
    estimated_bins = (
        math.ceil(effective_length / effective_bin_size)
        if effective_length and effective_bin_size
        else 0
    )
    return AssemblyBinPlan(
        requested_bin_size=int(requested_bin_size),
        effective_bin_size=effective_bin_size,
        estimated_global_bins=int(estimated_bins),
        use_global_assembly_bins=True,
        is_small_assembly=is_small,
        is_fragmented_assembly=is_fragmented,
        reason=";".join(reasons),
    )
