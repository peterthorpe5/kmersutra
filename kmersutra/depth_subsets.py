"""Deterministic one-pass depth subsets for FASTA and FASTQ benchmarks."""

from __future__ import annotations

import gzip
import hashlib
import logging
from contextlib import ExitStack
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO

from kmersutra.io import open_text_reader
from kmersutra.table_io import write_records_table


@dataclass(frozen=True)
class RawSequenceRecord:
    """A sequence record retaining its original quality string.

    Attributes:
        identifier: First token from the header.
        description: Complete header without the leading marker.
        sequence: Sequence text.
        quality: FASTQ quality text, or ``None`` for FASTA.
    """

    identifier: str
    description: str
    sequence: str
    quality: str | None


@dataclass
class SubsetWriter:
    """State for one deterministic subset output."""

    fraction: float
    seed: int
    sample_id: str
    final_path: Path
    temporary_path: Path
    handle: TextIO
    n_selected: int = 0


def validate_fractions(fractions: list[float]) -> list[float]:
    """Validate and sort unique sampling fractions.

    Args:
        fractions: Requested fractions.

    Returns:
        Sorted unique fractions.

    Raises:
        ValueError: If a fraction is outside (0, 1].
    """
    output = sorted(set(float(value) for value in fractions))
    if not output:
        raise ValueError("At least one fraction is required")
    invalid = [value for value in output if not 0.0 < value <= 1.0]
    if invalid:
        raise ValueError(
            "Fractions must be greater than zero and at most one: "
            + ", ".join(str(value) for value in invalid)
        )
    return output


def validate_seeds(seeds: list[int]) -> list[int]:
    """Validate and retain unique deterministic seeds.

    Args:
        seeds: Requested integer seeds.

    Returns:
        Unique seeds in first-seen order.

    Raises:
        ValueError: If no seeds are supplied.
    """
    output = list(dict.fromkeys(int(value) for value in seeds))
    if not output:
        raise ValueError("At least one seed is required")
    return output


def fraction_label(fraction: float) -> str:
    """Return a filesystem-safe fraction label.

    Args:
        fraction: Fraction in (0, 1].

    Returns:
        Stable label such as ``p010000`` for 1%.
    """
    return f"p{round(fraction * 1_000_000):06d}"


def deterministic_unit_interval(
    *,
    identifier: str,
    record_index: int,
    seed: int,
) -> float:
    """Map a record identity deterministically into [0, 1).

    Args:
        identifier: Sequence record identifier.
        record_index: Zero-based record index.
        seed: Sampling seed.

    Returns:
        Deterministic pseudo-random value in [0, 1).

    Raises:
        ValueError: If ``record_index`` is negative.
    """
    if record_index < 0:
        raise ValueError("record_index must be zero or greater")
    payload = f"{seed}\0{record_index}\0{identifier}".encode("utf-8")
    digest = hashlib.blake2b(payload, digest_size=8).digest()
    integer = int.from_bytes(digest, byteorder="big", signed=False)
    return integer / float(1 << 64)


def iter_fastq_raw(
    *,
    input_path: str | Path,
    decompressor: str,
):
    """Yield validated FASTQ records while retaining qualities.

    Args:
        input_path: FASTQ or FASTQ.GZ input.
        decompressor: Decompressor mode passed to KmerSutra I/O.

    Yields:
        Validated raw FASTQ records.

    Raises:
        ValueError: If a record is malformed or truncated.
    """
    with open_text_reader(input_path, decompressor=decompressor) as handle:
        record_number = 0
        while True:
            header = handle.readline()
            if not header:
                break
            sequence = handle.readline()
            plus = handle.readline()
            quality = handle.readline()
            record_number += 1
            if not sequence or not plus or not quality:
                raise ValueError(
                    f"Truncated FASTQ record {record_number} in {input_path}"
                )
            header = header.rstrip("\r\n")
            sequence = sequence.rstrip("\r\n")
            plus = plus.rstrip("\r\n")
            quality = quality.rstrip("\r\n")
            if not header.startswith("@") or not plus.startswith("+"):
                raise ValueError(
                    f"Malformed FASTQ record {record_number} in {input_path}"
                )
            if len(sequence) != len(quality):
                raise ValueError(
                    f"Sequence and quality lengths differ in FASTQ record "
                    f"{record_number} in {input_path}"
                )
            description = header[1:].strip()
            identifier = description.split()[0] if description else ""
            if not identifier:
                raise ValueError(
                    f"Blank FASTQ identifier at record {record_number} in {input_path}"
                )
            yield RawSequenceRecord(
                identifier=identifier,
                description=description,
                sequence=sequence,
                quality=quality,
            )


def iter_fasta_raw(
    *,
    input_path: str | Path,
    decompressor: str,
):
    """Yield validated FASTA records.

    Args:
        input_path: FASTA or FASTA.GZ input.
        decompressor: Decompressor mode passed to KmerSutra I/O.

    Yields:
        Validated FASTA records.

    Raises:
        ValueError: If sequence text occurs before the first header.
    """
    description: str | None = None
    sequence_parts: list[str] = []
    with open_text_reader(input_path, decompressor=decompressor) as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\r\n")
            if not line:
                continue
            if line.startswith(">"):
                if description is not None:
                    identifier = description.split()[0]
                    yield RawSequenceRecord(
                        identifier=identifier,
                        description=description,
                        sequence="".join(sequence_parts),
                        quality=None,
                    )
                description = line[1:].strip()
                if not description:
                    raise ValueError(
                        f"Blank FASTA header at line {line_number} in {input_path}"
                    )
                sequence_parts = []
            else:
                if description is None:
                    raise ValueError(
                        f"FASTA sequence before first header at line "
                        f"{line_number} in {input_path}"
                    )
                sequence_parts.append(line)
    if description is not None:
        yield RawSequenceRecord(
            identifier=description.split()[0],
            description=description,
            sequence="".join(sequence_parts),
            quality=None,
        )


def _open_output(path: Path, *, compress: bool):
    """Open a temporary sequence output."""
    if compress:
        return gzip.open(path, "wt", encoding="utf-8", newline="\n")
    return path.open("w", encoding="utf-8", newline="\n")


def _write_record(
    *,
    handle: TextIO,
    record: RawSequenceRecord,
    input_format: str,
) -> None:
    """Write one FASTA or FASTQ record."""
    if input_format == "fastq":
        if record.quality is None:
            raise ValueError("FASTQ output requires a quality string")
        handle.write(
            f"@{record.description}\n{record.sequence}\n+\n{record.quality}\n"
        )
    else:
        handle.write(f">{record.description}\n{record.sequence}\n")


def create_depth_subsets(
    *,
    input_path: str | Path,
    input_format: str,
    output_dir: str | Path,
    fractions: list[float],
    seeds: list[int],
    sample_prefix: str,
    compress: bool = True,
    decompressor: str = "auto",
    manifest_path: str | Path | None = None,
    progress_interval: int = 1_000_000,
    logger: logging.Logger | None = None,
) -> list[dict[str, object]]:
    """Create all depth/seed subsets in one input pass.

    Selection is deterministic for the same input order, identifiers, seed and
    KmerSutra version. Fractions for one seed are nested because they use the
    same hash value.

    Args:
        input_path: FASTA/FASTQ input.
        input_format: ``fasta`` or ``fastq``.
        output_dir: Output directory.
        fractions: Fractions in (0, 1].
        seeds: Integer seeds.
        sample_prefix: Prefix for generated sample identifiers.
        compress: Write gzip-compressed outputs.
        decompressor: Input gzip decompressor mode.
        manifest_path: Optional explicit manifest path.
        progress_interval: Record interval for progress logging.
        logger: Optional logger.

    Returns:
        Subset manifest records.

    Raises:
        FileNotFoundError: If the input is missing or empty.
        ValueError: If configuration values are invalid.
    """
    source = Path(input_path).expanduser().resolve()
    if not source.is_file() or source.stat().st_size <= 0:
        raise FileNotFoundError(f"Sequence input is missing or empty: {source}")
    if input_format not in {"fastq", "fasta"}:
        raise ValueError("input_format must be 'fastq' or 'fasta'")
    if progress_interval <= 0:
        raise ValueError("progress_interval must be positive")
    resolved_fractions = validate_fractions(fractions)
    resolved_seeds = validate_seeds(seeds)
    destination = Path(output_dir).expanduser().resolve()
    destination.mkdir(parents=True, exist_ok=True)
    extension = f".{input_format}" + (".gz" if compress else "")
    temporary_paths: list[Path] = []
    writers: list[SubsetWriter] = []
    n_total = 0
    iterator = iter_fastq_raw if input_format == "fastq" else iter_fasta_raw

    try:
        with ExitStack() as stack:
            for seed in resolved_seeds:
                for fraction in resolved_fractions:
                    label = fraction_label(fraction)
                    sample_id = f"{sample_prefix}_depth_{label}_seed_{seed}"
                    final_path = destination / f"{sample_id}{extension}"
                    temporary_path = destination / f".{final_path.name}.tmp"
                    temporary_paths.append(temporary_path)
                    handle = stack.enter_context(
                        _open_output(temporary_path, compress=compress)
                    )
                    writers.append(
                        SubsetWriter(
                            fraction=fraction,
                            seed=seed,
                            sample_id=sample_id,
                            final_path=final_path,
                            temporary_path=temporary_path,
                            handle=handle,
                        )
                    )

            writers_by_seed: dict[int, list[SubsetWriter]] = {
                seed: [
                    writer for writer in writers if writer.seed == seed
                ]
                for seed in resolved_seeds
            }
            for record_index, record in enumerate(
                iterator(input_path=source, decompressor=decompressor)
            ):
                n_total += 1
                for seed, seed_writers in writers_by_seed.items():
                    draw = deterministic_unit_interval(
                        identifier=record.identifier,
                        record_index=record_index,
                        seed=seed,
                    )
                    for writer in seed_writers:
                        if draw < writer.fraction:
                            _write_record(
                                handle=writer.handle,
                                record=record,
                                input_format=input_format,
                            )
                            writer.n_selected += 1
                if logger and n_total % progress_interval == 0:
                    logger.info("Processed %d input records for depth subsets", n_total)

        if n_total == 0:
            raise ValueError(f"Sequence input contains no records: {source}")
        for writer in writers:
            writer.temporary_path.replace(writer.final_path)
    except Exception:
        for path in temporary_paths:
            path.unlink(missing_ok=True)
        raise

    manifest = [
        {
            "sample_id": writer.sample_id,
            "fraction": f"{writer.fraction:.8f}",
            "seed": writer.seed,
            "input_format": input_format,
            "reads_path": str(writer.final_path),
            "n_input_records": n_total,
            "n_selected_records": writer.n_selected,
            "observed_fraction": f"{writer.n_selected / n_total:.8f}",
            "selection_method": "blake2b_identifier_record_index_v1",
        }
        for writer in writers
    ]
    output_manifest = (
        Path(manifest_path).expanduser().resolve()
        if manifest_path is not None
        else destination / "depth_subset_manifest.tsv"
    )
    write_records_table(
        records=manifest,
        output_path=output_manifest,
        fieldnames=list(manifest[0]),
        logger=logger,
    )
    if logger:
        logger.info(
            "Created %d deterministic subsets from %d input records",
            len(manifest),
            n_total,
        )
    return manifest
