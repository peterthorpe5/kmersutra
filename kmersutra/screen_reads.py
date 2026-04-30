"""Screen reads or assemblies against a KmerSutra panel."""

from __future__ import annotations

import logging
from collections.abc import Iterable, Iterator
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path

from kmersutra.build_panel import DiagnosticKmer, load_panel
from kmersutra.fasta import SequenceRecord, read_fasta_records, read_fastq_records
from kmersutra.kmers import VALID_BASES, hamming_distance, iter_kmers

_GLOBAL_PANEL_INDEX: dict[int, dict[str, list[DiagnosticKmer]]] | None = None
_GLOBAL_SAMPLE_ID = ""
_GLOBAL_SEQUENCE_TYPE = ""
_GLOBAL_MAX_MISMATCHES = 0
_GLOBAL_FUZZY_MIN_K = 71


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
        }


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

    Returns
    -------
    list[KmerHit]
        Query hits.
    """
    if max_mismatches < 0:
        raise ValueError("max_mismatches must be zero or greater")
    if max_mismatches > 2:
        raise ValueError("max_mismatches above two is not currently supported")

    hits: list[KmerHit] = []
    for k, diagnostics in panel_index.items():
        for position, query_kmer in iter_kmers(sequence=sequence_record.sequence, k=k):
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
            )
        )
    return hits


def _init_worker(
    panel_index: dict[int, dict[str, list[DiagnosticKmer]]],
    sample_id: str,
    sequence_type: str,
    max_mismatches: int,
    fuzzy_min_k: int,
) -> None:
    """Initialise worker process globals for screening."""
    global _GLOBAL_PANEL_INDEX
    global _GLOBAL_SAMPLE_ID
    global _GLOBAL_SEQUENCE_TYPE
    global _GLOBAL_MAX_MISMATCHES
    global _GLOBAL_FUZZY_MIN_K
    _GLOBAL_PANEL_INDEX = panel_index
    _GLOBAL_SAMPLE_ID = sample_id
    _GLOBAL_SEQUENCE_TYPE = sequence_type
    _GLOBAL_MAX_MISMATCHES = max_mismatches
    _GLOBAL_FUZZY_MIN_K = fuzzy_min_k


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

    if logger:
        logger.info(
            "Screening records using %d worker(s), chunk_size=%d, max_mismatches=%d",
            threads,
            chunk_size,
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
                    panel_index=panel_index,
                    sample_id=sample_id,
                    sequence_type=sequence_type,
                    max_mismatches=max_mismatches,
                    fuzzy_min_k=fuzzy_min_k,
                )
            )
        if logger:
            logger.info("Screened %d records and retained %d hits", n_records, len(all_hits))
        return all_hits

    chunks = list(_chunks(records=records, chunk_size=chunk_size))
    if logger:
        logger.info("Prepared %d screening chunks", len(chunks))

    all_hits: list[KmerHit] = []
    n_records = sum(len(chunk) for chunk in chunks)
    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = [
            executor.submit(
                _screen_chunk,
                records=chunk,
                panel_index=panel_index,
                sample_id=sample_id,
                sequence_type=sequence_type,
                max_mismatches=max_mismatches,
                fuzzy_min_k=fuzzy_min_k,
            )
            for chunk in chunks
        ]
        for index, future in enumerate(as_completed(futures), start=1):
            chunk_hits = future.result()
            all_hits.extend(chunk_hits)
            if logger and (index == 1 or index % 10 == 0 or index == len(futures)):
                logger.info(
                    "Completed %d/%d chunks; retained %d hits so far",
                    index,
                    len(futures),
                    len(all_hits),
                )

    if logger:
        logger.info("Screened %d records and retained %d hits", n_records, len(all_hits))
    return all_hits


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
    logger : logging.Logger | None, optional
        Logger for progress messages.

    Returns
    -------
    list[KmerHit]
        Query hits.
    """
    if logger:
        logger.info("Loading panel: %s", panel_path)
    panel_index = load_panel(panel_path=panel_path)
    if logger:
        n_panel_kmers = sum(len(kmer_map) for kmer_map in panel_index.values())
        logger.info(
            "Loaded panel with %d k values and %d unique panel k-mer keys",
            len(panel_index),
            n_panel_kmers,
        )

    if input_format == "fasta":
        records = read_fasta_records(fasta_path=input_path)
        sequence_type = "contig"
    elif input_format == "fastq":
        records = read_fastq_records(fastq_path=input_path)
        sequence_type = "read"
    else:
        raise ValueError("input_format must be either fasta or fastq")

    return screen_records_for_species_kmers(
        records=records,
        panel_index=panel_index,
        sample_id=sample_id,
        sequence_type=sequence_type,
        max_mismatches=max_mismatches,
        fuzzy_min_k=fuzzy_min_k,
        threads=threads,
        chunk_size=chunk_size,
        logger=logger,
    )
