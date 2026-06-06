"""FASTA and FASTQ parsing utilities."""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path

from kmersutra.io import open_text_reader


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


def normalise_sequence(sequence: str, mask_lowercase: bool = False) -> str:
    """Normalise a nucleotide sequence for k-mer processing.

    Parameters
    ----------
    sequence : str
        Input nucleotide sequence.
    mask_lowercase : bool, optional
        If true, lowercase ``a/c/g/t/u`` bases are converted to ``N`` before
        uppercasing. This allows repeat-masked FASTA regions to be excluded
        by the usual ambiguous-base k-mer filter without changing default
        behaviour.

    Returns
    -------
    str
        Uppercase sequence with whitespace removed and U converted to T.
    """
    normalised: list[str] = []
    for character in sequence:
        if character.isspace():
            continue
        if mask_lowercase and character in {"a", "c", "g", "t", "u"}:
            normalised.append("N")
            continue
        base = character.upper()
        normalised.append("T" if base == "U" else base)
    return "".join(normalised)


def read_fasta_records(
    *,
    fasta_path: str | Path,
    mask_lowercase: bool = False,
    decompressor: str = "python",
) -> Iterator[SequenceRecord]:
    """Yield records from a FASTA or FASTA.GZ file.

    Parameters
    ----------
    fasta_path : str or pathlib.Path
        Path to the FASTA file.
    mask_lowercase : bool, optional
        If true, lowercase repeat-masked bases are converted to ``N`` before
        sequence records are yielded.
    decompressor : str, optional
        Text decompressor for gzip inputs. ``python`` preserves historical
        behaviour; ``pigz`` uses external pigz; ``auto`` uses pigz when
        available and falls back to Python gzip.

    Yields
    ------
    SequenceRecord
        Parsed sequence record.
    """
    header: str | None = None
    sequence_parts: list[str] = []

    with open_text_reader(fasta_path, decompressor=decompressor) as handle:
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
                        sequence=normalise_sequence(
                "".join(sequence_parts),
                mask_lowercase=mask_lowercase,
            ),
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
            sequence=normalise_sequence(
                "".join(sequence_parts),
                mask_lowercase=mask_lowercase,
            ),
        )


def read_fastq_records(
    *,
    fastq_path: str | Path,
    decompressor: str = "python",
) -> Iterator[SequenceRecord]:
    """Yield records from a FASTQ or FASTQ.GZ file.

    Parameters
    ----------
    fastq_path : str or pathlib.Path
        Path to FASTQ file.
    decompressor : str, optional
        Text decompressor for gzip inputs. ``python`` preserves historical
        behaviour; ``pigz`` uses external pigz; ``auto`` uses pigz when
        available and falls back to Python gzip.

    Yields
    ------
    SequenceRecord
        Parsed sequence record without quality scores.
    """
    with open_text_reader(fastq_path, decompressor=decompressor) as handle:
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
                raise ValueError(
                    f"FASTQ sequence and quality length differ in {fastq_path}"
                )
            description = header[1:].strip()
            identifier = description.split()[0]
            yield SequenceRecord(
                identifier=identifier,
                description=description,
                sequence=normalise_sequence(sequence),
            )
