"""FASTA and FASTQ parsing utilities."""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path

from kmersutra.io import open_text


@dataclass(frozen=True)
class SequenceRecord:
    """A simple sequence record.

    Attributes
    ----------
    identifier : str
        Record identifier, taken from the first token of the header.
    description : str
        Full header line without the leading marker.
    sequence : str
        Sequence string.
    """

    identifier: str
    description: str
    sequence: str


def normalise_sequence(sequence: str) -> str:
    """Normalise a nucleotide sequence for k-mer processing.

    Parameters
    ----------
    sequence : str
        Input nucleotide sequence.

    Returns
    -------
    str
        Uppercase sequence with whitespace removed and U converted to T.
    """
    return "".join(sequence.upper().replace("U", "T").split())


def read_fasta_records(*, fasta_path: str | Path) -> Iterator[SequenceRecord]:
    """Yield records from a FASTA or FASTA.GZ file.

    Parameters
    ----------
    fasta_path : str or pathlib.Path
        Path to the FASTA file.

    Yields
    ------
    SequenceRecord
        Parsed sequence record.
    """
    header: str | None = None
    sequence_parts: list[str] = []

    with open_text(fasta_path, "rt") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    identifier = header.split()[0]
                    yield SequenceRecord(
                        identifier=identifier,
                        description=header,
                        sequence=normalise_sequence("".join(sequence_parts)),
                    )
                header = line[1:].strip()
                sequence_parts = []
            else:
                sequence_parts.append(line)

    if header is not None:
        identifier = header.split()[0]
        yield SequenceRecord(
            identifier=identifier,
            description=header,
            sequence=normalise_sequence("".join(sequence_parts)),
        )


def read_fastq_records(*, fastq_path: str | Path) -> Iterator[SequenceRecord]:
    """Yield records from a FASTQ or FASTQ.GZ file.

    Parameters
    ----------
    fastq_path : str or pathlib.Path
        Path to FASTQ file.

    Yields
    ------
    SequenceRecord
        Parsed sequence record without quality scores.
    """
    with open_text(fastq_path, "rt") as handle:
        while True:
            header = handle.readline().rstrip("\n")
            if not header:
                break
            sequence = handle.readline().rstrip("\n")
            plus = handle.readline().rstrip("\n")
            quality = handle.readline().rstrip("\n")
            if not header.startswith("@") or not plus.startswith("+"):
                raise ValueError(f"Malformed FASTQ record in {fastq_path}")
            if len(quality) != len(sequence):
                raise ValueError(f"FASTQ sequence and quality length differ in {fastq_path}")
            description = header[1:].strip()
            identifier = description.split()[0]
            yield SequenceRecord(
                identifier=identifier,
                description=description,
                sequence=normalise_sequence(sequence),
            )
