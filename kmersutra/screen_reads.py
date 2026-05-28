"""Screen reads or assemblies against a KmerSutra panel."""

from __future__ import annotations

import logging
from collections.abc import Iterable, Iterator
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
import time
from pathlib import Path

from kmersutra.build_panel import DiagnosticKmer
from kmersutra.panel_cache import load_panel_with_cache
from kmersutra.fasta import SequenceRecord, read_fasta_records, read_fastq_records
from kmersutra.kmers import VALID_BASES, hamming_distance, iter_kmers, reverse_complement

_GLOBAL_PANEL_INDEX: dict[int, dict[str, list[DiagnosticKmer]]] | None = None
_GLOBAL_SAMPLE_ID = ""
_GLOBAL_SEQUENCE_TYPE = ""
_GLOBAL_MAX_MISMATCHES = 0
_GLOBAL_FUZZY_MIN_K = 71
_GLOBAL_EXACT_ORIENTATION_INDEX = False


@dataclass(frozen=True)
class KmerHit:
    """A diagnostic k-mer hit in a query sequence.

    Attributes
    ----------
    sample_id : str
        Sample identifier.
    sequence_id : str
        Read or contig identifier.
    sequence_type : str
        Query type, usually read or contig.
    k : int
        K-mer length.
    query_position : int
        Zero-based query position.
    matched_kmer : str
        Diagnostic panel k-mer matched.
    query_kmer : str
        Query k-mer.
    mismatches : int
        Number of mismatches between query and panel k-mer.
    panel_type : str
        Panel class.
    species_name : str
        Species label.
    clade : str
        Clade label.
    evidence_taxid : str
        Taxid represented by this evidence level, when available.
    evidence_name : str
        Taxonomic name represented by this evidence level, when available.
    evidence_rank : str
        Taxonomic rank represented by this evidence level, when available.
    """

    sample_id: str
    sequence_id: str
    sequence_type: str
    k: int
    query_position: int
    matched_kmer: str
    query_kmer: str
    mismatches: int
    panel_type: str
    species_name: str
    clade: str
    evidence_taxid: str = ""
    evidence_name: str = ""
    evidence_rank: str = ""

    def to_record(self) -> dict[str, object]:
        """Convert the hit to a serialisable record.

        Returns
        -------
        dict[str, object]
            Dictionary representation.
        """
        return {
            "sample_id": self.sample_id,
            "sequence_id": self.sequence_id,
            "sequence_type": self.sequence_type,
            "k": self.k,
            "query_position": self.query_position,
            "matched_kmer": self.matched_kmer,
            "query_kmer": self.query_kmer,
            "mismatches": self.mismatches,
            "panel_type": self.panel_type,
            "species_name": self.species_name,
            "clade": self.clade,
            "evidence_taxid": self.evidence_taxid,
            "evidence_name": self.evidence_name,
            "evidence_rank": self.evidence_rank,
        }



def build_orientation_aware_exact_index(
    *,
    panel_index: dict[int, dict[str, list[DiagnosticKmer]]],
) -> dict[int, dict[str, list[DiagnosticKmer]]]:
    """Build a forward-or-reverse-complement exact-match index.

    Parameters
    ----------
    panel_index : dict[int, dict[str, list[DiagnosticKmer]]]
        Canonical panel index loaded from a KmerSutra panel.

    Returns
    -------
    dict[int, dict[str, list[DiagnosticKmer]]]
        Orientation-aware exact-match index. For each diagnostic k-mer, both
        the stored k-mer and its reverse complement are indexed. This allows
        exact screening to avoid canonicalising every query k-mer, which is
        expensive for long k values and large FASTQ files.
    """
    oriented: dict[int, dict[str, list[DiagnosticKmer]]] = {}
    for k, kmer_map in panel_index.items():
        oriented_map: dict[str, list[DiagnosticKmer]] = {}
        for panel_kmer, diagnostics in kmer_map.items():
            if not panel_kmer:
                continue
            panel_kmer_upper = panel_kmer.upper()
            oriented_map.setdefault(panel_kmer_upper, []).extend(diagnostics)
            reverse = reverse_complement(panel_kmer_upper)
            if reverse != panel_kmer_upper:
                oriented_map.setdefault(reverse, []).extend(diagnostics)
        oriented[int(k)] = oriented_map
    return oriented


def iter_unambiguous_windows(
    *,
    sequence: str,
    k: int,
) -> Iterator[tuple[int, str]]:
    """Yield raw k-length windows from unambiguous sequence segments.

    Parameters
    ----------
    sequence : str
        Normalised nucleotide sequence.
    k : int
        K-mer length.

    Yields
    ------
    tuple[int, str]
        Zero-based start position and raw forward-orientation k-mer.

    Notes
    -----
    This iterator avoids the per-window ``set(kmer)`` validity check and
    reverse-complement canonicalisation used by the general-purpose iterator.
    It is intended for exact matching against an orientation-aware panel index.
    """
    if k <= 0:
        raise ValueError("k must be a positive integer")
    sequence_length = len(sequence)
    if sequence_length < k:
        return

    segment_start = 0
    while segment_start < sequence_length:
        while segment_start < sequence_length and sequence[segment_start] not in VALID_BASES:
            segment_start += 1
        segment_end = segment_start
        while segment_end < sequence_length and sequence[segment_end] in VALID_BASES:
            segment_end += 1

        if segment_end - segment_start >= k:
            last_start = segment_end - k
            for start in range(segment_start, last_start + 1):
                yield start, sequence[start : start + k]

        segment_start = segment_end + 1

def _exact_hits(
    *,
    query_kmer: str,
    diagnostics: dict[str, list[DiagnosticKmer]],
) -> list[tuple[int, DiagnosticKmer]]:
    """Find exact diagnostic k-mer hits."""
    return [(0, item) for item in diagnostics.get(query_kmer, [])]


def iter_mismatch_neighbourhood(*, kmer: str, max_mismatches: int) -> Iterator[str]:
    """Yield sequence neighbours within one or two substitutions.

    Parameters
    ----------
    kmer : str
        Query k-mer.
    max_mismatches : int
        Maximum number of substitutions to generate.

    Yields
    ------
    str
        Candidate neighbouring k-mer.
    """
    if max_mismatches <= 0:
        return
    if max_mismatches > 2:
        raise ValueError("Fuzzy matching currently supports at most two mismatches")

    bases = tuple(sorted(VALID_BASES))
    kmer_list = list(kmer)
    seen: set[str] = set()

    for first_index, original_first in enumerate(kmer_list):
        for first_base in bases:
            if first_base == original_first:
                continue
            mutated = kmer_list.copy()
            mutated[first_index] = first_base
            candidate = "".join(mutated)
            if candidate not in seen:
                seen.add(candidate)
                yield candidate

    if max_mismatches < 2:
        return

    length = len(kmer_list)
    for first_index in range(length - 1):
        original_first = kmer_list[first_index]
        for second_index in range(first_index + 1, length):
            original_second = kmer_list[second_index]
            for first_base in bases:
                if first_base == original_first:
                    continue
                for second_base in bases:
                    if second_base == original_second:
                        continue
                    mutated = kmer_list.copy()
                    mutated[first_index] = first_base
                    mutated[second_index] = second_base
                    candidate = "".join(mutated)
                    if candidate not in seen:
                        seen.add(candidate)
                        yield candidate


def _fuzzy_hits(
    *,
    query_kmer: str,
    diagnostics: dict[str, list[DiagnosticKmer]],
    max_mismatches: int,
) -> list[tuple[int, DiagnosticKmer]]:
    """Find fuzzy diagnostic k-mer hits using neighbourhood lookup."""
    if max_mismatches <= 0:
        return []
    hits: list[tuple[int, DiagnosticKmer]] = []
    for candidate in iter_mismatch_neighbourhood(
        kmer=query_kmer,
        max_mismatches=max_mismatches,
    ):
        panel_items = diagnostics.get(candidate, [])
        if not panel_items:
            continue
        distance = hamming_distance(left=query_kmer, right=candidate)
        if 0 < distance <= max_mismatches:
            hits.extend((distance, item) for item in panel_items)
    return hits


def screen_sequence_for_kmers(
    *,
    sequence_record: SequenceRecord,
    panel_index: dict[int, dict[str, list[DiagnosticKmer]]],
    sample_id: str,
    sequence_type: str,
    max_mismatches: int = 0,
    fuzzy_min_k: int = 71,
    exact_orientation_index: bool = False,
) -> list[KmerHit]:
    """Screen a sequence for diagnostic k-mers.

    Parameters
    ----------
    sequence_record : SequenceRecord
        Query sequence.
    panel_index : dict[int, dict[str, list[DiagnosticKmer]]]
        Diagnostic k-mer index.
    sample_id : str
        Sample identifier.
    sequence_type : str
        Query sequence type.
    max_mismatches : int, optional
        Maximum mismatches for optional fuzzy matching.
    fuzzy_min_k : int, optional
        Minimum k value eligible for fuzzy matching.
    exact_orientation_index : bool, optional
        If True, ``panel_index`` is assumed to contain both forward and
        reverse-complement diagnostic keys, so query k-mers are looked up in
        their raw forward orientation without canonicalisation. This is only
        valid for exact matching.

    Returns
    -------
    list[KmerHit]
        Query hits.
    """
    if max_mismatches < 0:
        raise ValueError("max_mismatches must be zero or greater")
    if max_mismatches > 2:
        raise ValueError("max_mismatches above two is not currently supported")

    if exact_orientation_index and max_mismatches > 0:
        raise ValueError("exact_orientation_index is only valid for exact matching")

    hits: list[KmerHit] = []
    for k, diagnostics in panel_index.items():
        if exact_orientation_index:
            kmer_iterator = iter_unambiguous_windows(
                sequence=sequence_record.sequence,
                k=k,
            )
        else:
            kmer_iterator = iter_kmers(sequence=sequence_record.sequence, k=k)

        for position, query_kmer in kmer_iterator:
            matched = _exact_hits(query_kmer=query_kmer, diagnostics=diagnostics)
            if max_mismatches > 0 and k >= fuzzy_min_k:
                matched.extend(
                    _fuzzy_hits(
                        query_kmer=query_kmer,
                        diagnostics=diagnostics,
                        max_mismatches=max_mismatches,
                    )
                )
            for mismatches, diagnostic in matched:
                hits.append(
                    KmerHit(
                        sample_id=sample_id,
                        sequence_id=sequence_record.identifier,
                        sequence_type=sequence_type,
                        k=k,
                        query_position=position,
                        matched_kmer=diagnostic.kmer,
                        query_kmer=query_kmer,
                        mismatches=mismatches,
                        panel_type=diagnostic.panel_type,
                        species_name=diagnostic.species_name,
                        clade=diagnostic.clade,
                        evidence_taxid=diagnostic.evidence_taxid,
                        evidence_name=diagnostic.evidence_name,
                        evidence_rank=diagnostic.evidence_rank,
                    )
                )
    return hits


def _chunks(*, records: Iterable[SequenceRecord], chunk_size: int) -> Iterator[list[SequenceRecord]]:
    """Yield fixed-size chunks from a record iterable.

    Parameters
    ----------
    records : iterable of SequenceRecord
        Records to chunk.
    chunk_size : int
        Number of records per chunk.

    Yields
    ------
    list[SequenceRecord]
        Record chunk.
    """
    chunk: list[SequenceRecord] = []
    for record in records:
        chunk.append(record)
        if len(chunk) >= chunk_size:
            yield chunk
            chunk = []
    if chunk:
        yield chunk


def _screen_chunk(
    *,
    records: list[SequenceRecord],
    panel_index: dict[int, dict[str, list[DiagnosticKmer]]],
    sample_id: str,
    sequence_type: str,
    max_mismatches: int,
    fuzzy_min_k: int,
    exact_orientation_index: bool,
) -> list[KmerHit]:
    """Screen one chunk of records."""
    hits: list[KmerHit] = []
    for record in records:
        hits.extend(
            screen_sequence_for_kmers(
                sequence_record=record,
                panel_index=panel_index,
                sample_id=sample_id,
                sequence_type=sequence_type,
                max_mismatches=max_mismatches,
                fuzzy_min_k=fuzzy_min_k,
                exact_orientation_index=exact_orientation_index,
            )
        )
    return hits


def _init_worker(
    panel_index: dict[int, dict[str, list[DiagnosticKmer]]],
    sample_id: str,
    sequence_type: str,
    max_mismatches: int,
    fuzzy_min_k: int,
    exact_orientation_index: bool,
) -> None:
    """Initialise worker process globals for screening."""
    global _GLOBAL_PANEL_INDEX
    global _GLOBAL_SAMPLE_ID
    global _GLOBAL_SEQUENCE_TYPE
    global _GLOBAL_MAX_MISMATCHES
    global _GLOBAL_FUZZY_MIN_K
    global _GLOBAL_EXACT_ORIENTATION_INDEX
    _GLOBAL_PANEL_INDEX = panel_index
    _GLOBAL_SAMPLE_ID = sample_id
    _GLOBAL_SEQUENCE_TYPE = sequence_type
    _GLOBAL_MAX_MISMATCHES = max_mismatches
    _GLOBAL_FUZZY_MIN_K = fuzzy_min_k
    _GLOBAL_EXACT_ORIENTATION_INDEX = exact_orientation_index


def _screen_chunk_worker(records: list[SequenceRecord]) -> list[KmerHit]:
    """Screen one record chunk inside a worker process."""
    if _GLOBAL_PANEL_INDEX is None:
        raise RuntimeError("Worker panel index has not been initialised")
    return _screen_chunk(
        records=records,
        panel_index=_GLOBAL_PANEL_INDEX,
        sample_id=_GLOBAL_SAMPLE_ID,
        sequence_type=_GLOBAL_SEQUENCE_TYPE,
        max_mismatches=_GLOBAL_MAX_MISMATCHES,
        fuzzy_min_k=_GLOBAL_FUZZY_MIN_K,
        exact_orientation_index=_GLOBAL_EXACT_ORIENTATION_INDEX,
    )


def screen_records_for_species_kmers(
    *,
    records: Iterable[SequenceRecord],
    panel_index: dict[int, dict[str, list[DiagnosticKmer]]],
    sample_id: str,
    sequence_type: str,
    max_mismatches: int = 0,
    fuzzy_min_k: int = 71,
    threads: int = 1,
    chunk_size: int = 1000,
    max_pending_chunks: int | None = None,
    logger: logging.Logger | None = None,
) -> list[KmerHit]:
    """Screen query records for species or clade diagnostic k-mers.

    Parameters
    ----------
    records : iterable of SequenceRecord
        Query records.
    panel_index : dict[int, dict[str, list[DiagnosticKmer]]]
        Diagnostic k-mer index.
    sample_id : str
        Sample identifier.
    sequence_type : str
        Query sequence type.
    max_mismatches : int, optional
        Maximum mismatches for fuzzy matching.
    fuzzy_min_k : int, optional
        Minimum k value eligible for fuzzy matching.
    threads : int, optional
        Number of worker processes.
    chunk_size : int, optional
        Number of records submitted per worker task.
    max_pending_chunks : int or None, optional
        Maximum number of submitted chunks kept in the worker queue. If omitted,
        twice the worker count is used. This avoids materialising all chunks in
        memory for large FASTQ/FASTA files.
    logger : logging.Logger | None, optional
        Logger for progress messages.

    Returns
    -------
    list[KmerHit]
        Query hits.
    """
    if threads <= 0:
        raise ValueError("threads must be a positive integer")
    if chunk_size <= 0:
        raise ValueError("chunk_size must be a positive integer")

    use_exact_orientation_index = max_mismatches == 0
    screening_panel_index = (
        build_orientation_aware_exact_index(panel_index=panel_index)
        if use_exact_orientation_index
        else panel_index
    )

    if logger:
        logger.info(
            "Screening records using %d worker(s), chunk_size=%d, max_mismatches=%d",
            threads,
            chunk_size,
            max_mismatches,
        )
        if use_exact_orientation_index:
            n_keys = sum(len(kmer_map) for kmer_map in screening_panel_index.values())
            logger.info(
                "Using orientation-aware exact screening index with %d lookup keys",
                n_keys,
            )
        else:
            logger.info(
                "Using canonical/fuzzy screening path because max_mismatches=%d",
                max_mismatches,
            )

    if threads == 1:
        all_hits: list[KmerHit] = []
        n_records = 0
        for chunk in _chunks(records=records, chunk_size=chunk_size):
            n_records += len(chunk)
            all_hits.extend(
                _screen_chunk(
                    records=chunk,
                    panel_index=screening_panel_index,
                    sample_id=sample_id,
                    sequence_type=sequence_type,
                    max_mismatches=max_mismatches,
                    fuzzy_min_k=fuzzy_min_k,
                    exact_orientation_index=use_exact_orientation_index,
                )
            )
        if logger:
            logger.info("Screened %d records and retained %d hits", n_records, len(all_hits))
        return all_hits

    all_hits: list[KmerHit] = []
    n_records = 0
    submitted = 0
    completed = 0
    max_pending = max(threads * 2, 1) if max_pending_chunks is None else max_pending_chunks
    if max_pending <= 0:
        raise ValueError("max_pending_chunks must be a positive integer")

    def submit_chunk(executor: ThreadPoolExecutor, chunk: list[SequenceRecord]):
        return executor.submit(
            _screen_chunk,
            records=chunk,
            panel_index=screening_panel_index,
            sample_id=sample_id,
            sequence_type=sequence_type,
            max_mismatches=max_mismatches,
            fuzzy_min_k=fuzzy_min_k,
            exact_orientation_index=use_exact_orientation_index,
        )

    with ThreadPoolExecutor(max_workers=threads) as executor:
        pending = set()
        for chunk in _chunks(records=records, chunk_size=chunk_size):
            n_records += len(chunk)
            pending.add(submit_chunk(executor, chunk))
            submitted += 1

            if len(pending) >= max_pending:
                done = next(as_completed(pending))
                pending.remove(done)
                all_hits.extend(done.result())
                completed += 1
                if logger and (completed == 1 or completed % 10 == 0):
                    logger.info(
                        "Completed %d chunks; submitted %d chunks; retained %d hits so far",
                        completed,
                        submitted,
                        len(all_hits),
                    )

        for future in as_completed(pending):
            all_hits.extend(future.result())
            completed += 1
            if logger and (completed == 1 or completed % 10 == 0 or completed == submitted):
                logger.info(
                    "Completed %d/%d chunks; retained %d hits so far",
                    completed,
                    submitted,
                    len(all_hits),
                )

    if logger:
        logger.info("Screened %d records and retained %d hits", n_records, len(all_hits))
    return all_hits



def screen_file_for_panel_index(
    *,
    input_path: str | Path,
    panel_index: dict[int, dict[str, list[DiagnosticKmer]]],
    sample_id: str,
    input_format: str,
    max_mismatches: int = 0,
    fuzzy_min_k: int = 71,
    threads: int = 1,
    chunk_size: int = 1000,
    max_pending_chunks: int | None = None,
    profile_records: list[dict[str, object]] | None = None,
    logger: logging.Logger | None = None,
) -> list[KmerHit]:
    """Screen a FASTA or FASTQ file against a pre-loaded panel index.

    Parameters
    ----------
    input_path : str or pathlib.Path
        Query FASTA/FASTQ path.
    panel_index : dict[int, dict[str, list[DiagnosticKmer]]]
        Pre-loaded diagnostic panel index. This is used by hierarchical
        two-pass screening to avoid writing temporary merged panel files and
        to avoid repeatedly loading many module panels.
    sample_id : str
        Sample identifier.
    input_format : str
        Input format, either ``fasta`` or ``fastq``.
    max_mismatches : int, optional
        Maximum mismatches for fuzzy matching.
    fuzzy_min_k : int, optional
        Minimum k value eligible for fuzzy matching.
    threads : int, optional
        Number of worker threads.
    chunk_size : int, optional
        Number of records submitted per worker task.
    max_pending_chunks : int or None, optional
        Maximum number of submitted chunks kept in the worker queue.
    profile_records : list[dict[str, object]] or None, optional
        Optional profile sink.
    logger : logging.Logger or None, optional
        Logger for progress messages.

    Returns
    -------
    list[KmerHit]
        Query hits.

    Raises
    ------
    ValueError
        If ``input_format`` is not ``fasta`` or ``fastq``.
    """
    if not panel_index:
        if logger:
            logger.info("Pre-loaded panel index is empty; returning no hits")
        return []

    n_panel_kmers = sum(len(kmer_map) for kmer_map in panel_index.values())
    if logger:
        logger.info(
            "Screening with pre-loaded panel index containing %d k values and "
            "%d unique panel k-mer keys",
            len(panel_index),
            n_panel_kmers,
        )

    parse_start = time.perf_counter()
    if input_format == "fasta":
        records = read_fasta_records(fasta_path=input_path)
        sequence_type = "contig"
    elif input_format == "fastq":
        records = read_fastq_records(fastq_path=input_path)
        sequence_type = "read"
    else:
        raise ValueError("input_format must be either fasta or fastq")
    parse_seconds = time.perf_counter() - parse_start
    if profile_records is not None:
        profile_records.append({
            "stage": "prepare_input_iterator",
            "seconds": f"{parse_seconds:.6f}",
            "detail": input_format,
        })

    screen_start = time.perf_counter()
    hits = screen_records_for_species_kmers(
        records=records,
        panel_index=panel_index,
        sample_id=sample_id,
        sequence_type=sequence_type,
        max_mismatches=max_mismatches,
        fuzzy_min_k=fuzzy_min_k,
        threads=threads,
        chunk_size=chunk_size,
        max_pending_chunks=max_pending_chunks,
        logger=logger,
    )
    screen_seconds = time.perf_counter() - screen_start
    if profile_records is not None:
        profile_records.append({
            "stage": "screen_records",
            "seconds": f"{screen_seconds:.6f}",
            "detail": f"hits={len(hits)}",
        })
    return hits

def screen_file_for_species_kmers(
    *,
    input_path: str | Path,
    panel_path: str | Path,
    sample_id: str,
    input_format: str,
    max_mismatches: int = 0,
    fuzzy_min_k: int = 71,
    threads: int = 1,
    chunk_size: int = 1000,
    max_pending_chunks: int | None = None,
    panel_cache_path: str | Path | None = None,
    use_panel_cache: bool = False,
    write_panel_cache: bool = False,
    profile_records: list[dict[str, object]] | None = None,
    logger: logging.Logger | None = None,
) -> list[KmerHit]:
    """Screen a FASTA or FASTQ file against a diagnostic panel.

    Parameters
    ----------
    input_path : str or pathlib.Path
        Query FASTA/FASTQ path.
    panel_path : str or pathlib.Path
        Diagnostic panel TSV path.
    sample_id : str
        Sample identifier.
    input_format : str
        Input format, either fasta or fastq.
    max_mismatches : int, optional
        Maximum mismatches for fuzzy matching.
    fuzzy_min_k : int, optional
        Minimum k value eligible for fuzzy matching.
    threads : int, optional
        Number of worker processes.
    chunk_size : int, optional
        Number of records submitted per worker task.
    max_pending_chunks : int or None, optional
        Maximum number of submitted chunks kept in the worker queue. If omitted,
        twice the worker count is used. This avoids materialising all chunks in
        memory for large FASTQ/FASTA files.
    logger : logging.Logger | None, optional
        Logger for progress messages.

    Returns
    -------
    list[KmerHit]
        Query hits.
    """
    if logger:
        logger.info("Loading panel: %s", panel_path)
        logger.info("Use panel cache: %s", use_panel_cache)
        logger.info("Write panel cache: %s", write_panel_cache)
    panel_start = time.perf_counter()
    panel_index, panel_source = load_panel_with_cache(
        panel_path=panel_path,
        cache_path=panel_cache_path,
        use_cache=use_panel_cache,
        write_cache=write_panel_cache,
    )
    panel_seconds = time.perf_counter() - panel_start
    if profile_records is not None:
        profile_records.append({
            "stage": "load_panel",
            "seconds": f"{panel_seconds:.6f}",
            "detail": panel_source,
        })
    if logger:
        n_panel_kmers = sum(len(kmer_map) for kmer_map in panel_index.values())
        logger.info(
            "Loaded panel from %s with %d k values and %d unique panel k-mer keys in %.3fs",
            panel_source,
            len(panel_index),
            n_panel_kmers,
            panel_seconds,
        )

    parse_start = time.perf_counter()
    if input_format == "fasta":
        records = read_fasta_records(fasta_path=input_path)
        sequence_type = "contig"
    elif input_format == "fastq":
        records = read_fastq_records(fastq_path=input_path)
        sequence_type = "read"
    else:
        raise ValueError("input_format must be either fasta or fastq")
    parse_seconds = time.perf_counter() - parse_start
    if profile_records is not None:
        profile_records.append({
            "stage": "prepare_input_iterator",
            "seconds": f"{parse_seconds:.6f}",
            "detail": input_format,
        })

    screen_start = time.perf_counter()
    hits = screen_records_for_species_kmers(
        records=records,
        panel_index=panel_index,
        sample_id=sample_id,
        sequence_type=sequence_type,
        max_mismatches=max_mismatches,
        fuzzy_min_k=fuzzy_min_k,
        threads=threads,
        chunk_size=chunk_size,
        max_pending_chunks=max_pending_chunks,
        logger=logger,
    )
    screen_seconds = time.perf_counter() - screen_start
    if profile_records is not None:
        profile_records.append({
            "stage": "screen_records",
            "seconds": f"{screen_seconds:.6f}",
            "detail": f"hits={len(hits)}",
        })
    return hits
