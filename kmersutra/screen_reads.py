"""Screen reads or assemblies against a KmerSutra panel."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

from kmersutra.build_panel import DiagnosticKmer, load_panel
from kmersutra.fasta import SequenceRecord, read_fasta_records, read_fastq_records
from kmersutra.kmers import hamming_distance, iter_kmers


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
    query_position: int,
    diagnostics: dict[str, list[DiagnosticKmer]],
) -> list[tuple[int, DiagnosticKmer]]:
    """Find exact diagnostic k-mer hits.

    Parameters
    ----------
    query_kmer : str
        Query k-mer.
    query_position : int
        Query position.
    diagnostics : dict[str, list[DiagnosticKmer]]
        Diagnostic lookup for one k.

    Returns
    -------
    list[tuple[int, DiagnosticKmer]]
        Mismatch count and diagnostic record pairs.
    """
    return [(0, item) for item in diagnostics.get(query_kmer, [])]


def _fuzzy_hits(
    *,
    query_kmer: str,
    diagnostics: dict[str, list[DiagnosticKmer]],
    max_mismatches: int,
) -> list[tuple[int, DiagnosticKmer]]:
    """Find fuzzy diagnostic k-mer hits by Hamming distance.

    Parameters
    ----------
    query_kmer : str
        Query k-mer.
    diagnostics : dict[str, list[DiagnosticKmer]]
        Diagnostic lookup for one k.
    max_mismatches : int
        Maximum allowed mismatches.

    Returns
    -------
    list[tuple[int, DiagnosticKmer]]
        Mismatch count and diagnostic record pairs.
    """
    if max_mismatches <= 0:
        return []
    hits: list[tuple[int, DiagnosticKmer]] = []
    for panel_kmer, panel_items in diagnostics.items():
        distance = hamming_distance(left=query_kmer, right=panel_kmer)
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
    hits: list[KmerHit] = []
    for k, diagnostics in panel_index.items():
        for position, query_kmer in iter_kmers(sequence=sequence_record.sequence, k=k):
            matched = _exact_hits(
                query_kmer=query_kmer,
                query_position=position,
                diagnostics=diagnostics,
            )
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


def screen_records_for_species_kmers(
    *,
    records: Iterable[SequenceRecord],
    panel_index: dict[int, dict[str, list[DiagnosticKmer]]],
    sample_id: str,
    sequence_type: str,
    max_mismatches: int = 0,
    fuzzy_min_k: int = 71,
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

    Returns
    -------
    list[KmerHit]
        Query hits.
    """
    all_hits: list[KmerHit] = []
    for record in records:
        all_hits.extend(
            screen_sequence_for_kmers(
                sequence_record=record,
                panel_index=panel_index,
                sample_id=sample_id,
                sequence_type=sequence_type,
                max_mismatches=max_mismatches,
                fuzzy_min_k=fuzzy_min_k,
            )
        )
    return all_hits


def screen_file_for_species_kmers(
    *,
    input_path: str | Path,
    panel_path: str | Path,
    sample_id: str,
    input_format: str,
    max_mismatches: int = 0,
    fuzzy_min_k: int = 71,
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

    Returns
    -------
    list[KmerHit]
        Query hits.
    """
    panel_index = load_panel(panel_path=panel_path)
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
    )
