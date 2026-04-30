"""K-mer operations for KmerSutra."""

from __future__ import annotations

from collections import Counter
from collections.abc import Iterator

VALID_BASES = frozenset("ACGT")
_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def reverse_complement(sequence: str) -> str:
    """Return the reverse complement of a normalised DNA sequence.

    Parameters
    ----------
    sequence : str
        DNA sequence containing A, C, G and T.

    Returns
    -------
    str
        Reverse-complement sequence.
    """
    return sequence.translate(_COMPLEMENT)[::-1]


def canonical_kmer(kmer: str) -> str:
    """Return the canonical representation of a k-mer.

    Parameters
    ----------
    kmer : str
        Input k-mer.

    Returns
    -------
    str
        Lexicographically smaller of the forward k-mer and its reverse
        complement.
    """
    revcomp = reverse_complement(kmer)
    return min(kmer, revcomp)


def is_valid_kmer(kmer: str) -> bool:
    """Return whether a k-mer contains only A, C, G and T.

    Parameters
    ----------
    kmer : str
        Input k-mer.

    Returns
    -------
    bool
        True if all bases are unambiguous DNA bases.
    """
    return bool(kmer) and set(kmer).issubset(VALID_BASES)


def iter_kmers(
    *,
    sequence: str,
    k: int,
    canonical: bool = True,
    skip_ambiguous: bool = True,
) -> Iterator[tuple[int, str]]:
    """Yield k-mers from a sequence.

    Parameters
    ----------
    sequence : str
        Normalised nucleotide sequence.
    k : int
        K-mer length.
    canonical : bool, optional
        If true, return canonical k-mers.
    skip_ambiguous : bool, optional
        If true, skip k-mers containing non-ACGT bases.

    Yields
    ------
    tuple[int, str]
        Zero-based start position and k-mer.
    """
    if k <= 0:
        raise ValueError("k must be a positive integer")
    if len(sequence) < k:
        return
    for start in range(0, len(sequence) - k + 1):
        kmer = sequence[start : start + k]
        if skip_ambiguous and not is_valid_kmer(kmer):
            continue
        yield start, canonical_kmer(kmer) if canonical else kmer


def count_kmers(*, sequence: str, k: int) -> Counter[str]:
    """Count canonical k-mers in a sequence.

    Parameters
    ----------
    sequence : str
        Normalised sequence.
    k : int
        K-mer length.

    Returns
    -------
    collections.Counter[str]
        K-mer counts.
    """
    return Counter(kmer for _, kmer in iter_kmers(sequence=sequence, k=k))


def hamming_distance(*, left: str, right: str) -> int:
    """Calculate Hamming distance between equal-length strings.

    Parameters
    ----------
    left : str
        First sequence.
    right : str
        Second sequence.

    Returns
    -------
    int
        Number of differing positions.
    """
    if len(left) != len(right):
        raise ValueError("Hamming distance requires equal-length strings")
    return sum(a != b for a, b in zip(left, right))
