"""Input and output helpers for KmerSutra."""

from __future__ import annotations

import gzip
import io
import json
import logging
import shutil
import subprocess
from contextlib import contextmanager
from pathlib import Path
from typing import Iterable, Iterator, Mapping, TextIO

VALID_DECOMPRESSORS = {"python", "pigz", "auto"}
_PIGZ_SIGPIPE_CODES = {-13, 141}


def open_text(path: str | Path, mode: str = "rt") -> TextIO:
    """Open plain-text or gzip-compressed files.

    Parameters
    ----------
    path : str or pathlib.Path
        Input or output path. Files ending in ``.gz`` are opened with gzip.
    mode : str, optional
        Text mode, such as ``rt`` or ``wt``.

    Returns
    -------
    TextIO
        Open file handle.
    """
    file_path = Path(path)
    if "b" in mode:
        raise ValueError("open_text only supports text modes")
    if file_path.suffix == ".gz":
        return gzip.open(file_path, mode)  # type: ignore[return-value]
    return file_path.open(mode, encoding="utf-8")


def resolve_decompressor(
    *,
    path: str | Path,
    decompressor: str = "python",
    logger: logging.Logger | None = None,
) -> str:
    """Resolve the text decompressor to use for an input path.

    Parameters
    ----------
    path : str or pathlib.Path
        Input path that will be opened for reading.
    decompressor : str, optional
        One of ``python``, ``pigz`` or ``auto``. ``python`` preserves the
        historical behaviour. ``pigz`` requires the external ``pigz`` command.
        ``auto`` uses ``pigz`` for gzip inputs when available, otherwise Python
        gzip.
    logger : logging.Logger or None, optional
        Optional logger for fallback messages.

    Returns
    -------
    str
        Resolved decompressor: ``plain``, ``python`` or ``pigz``.

    Raises
    ------
    ValueError
        If ``decompressor`` is not recognised.
    RuntimeError
        If ``pigz`` was requested explicitly but is not available.
    """
    mode = str(decompressor).strip().lower()
    if mode not in VALID_DECOMPRESSORS:
        raise ValueError(
            "decompressor must be one of: "
            f"{', '.join(sorted(VALID_DECOMPRESSORS))}"
        )

    file_path = Path(path)
    if file_path.suffix != ".gz":
        return "plain"

    pigz_path = shutil.which("pigz")
    if mode == "python":
        return "python"
    if mode == "pigz":
        if pigz_path is None:
            raise RuntimeError(
                "The pigz decompressor was requested, but pigz was not found "
                "on PATH. Use --decompressor python or install pigz."
            )
        return "pigz"
    if pigz_path is not None:
        return "pigz"
    if logger:
        logger.info("pigz not found on PATH; falling back to Python gzip")
    return "python"


@contextmanager
def open_text_reader(
    path: str | Path,
    *,
    decompressor: str = "python",
    logger: logging.Logger | None = None,
) -> Iterator[TextIO]:
    """Open a text input stream with optional pigz decompression.

    Parameters
    ----------
    path : str or pathlib.Path
        Input file path. Plain files are opened normally. Gzip-compressed files
        are opened with either Python gzip or external ``pigz -dc``.
    decompressor : str, optional
        One of ``python``, ``pigz`` or ``auto``. The default preserves previous
        package behaviour.
    logger : logging.Logger or None, optional
        Optional logger used to report the selected decompression path.

    Yields
    ------
    TextIO
        Text-mode input handle.
    """
    file_path = Path(path)
    resolved = resolve_decompressor(
        path=file_path,
        decompressor=decompressor,
        logger=logger,
    )
    if resolved == "plain":
        if logger:
            logger.info("Reading plain text input: %s", file_path)
        with file_path.open("rt", encoding="utf-8") as handle:
            yield handle
        return
    if resolved == "python":
        if logger:
            logger.info("Reading gzip input with Python gzip: %s", file_path)
        with gzip.open(file_path, "rt") as handle:
            yield handle  # type: ignore[misc]
        return

    if logger:
        logger.info("Reading gzip input with pigz: %s", file_path)
    process = subprocess.Popen(
        ["pigz", "-dc", str(file_path)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    if process.stdout is None:
        raise RuntimeError("pigz did not provide a stdout stream")
    text_handle = io.TextIOWrapper(process.stdout, encoding="utf-8")
    raised_in_context = False
    try:
        yield text_handle
    except Exception:
        raised_in_context = True
        raise
    finally:
        try:
            text_handle.close()
        finally:
            stderr_text = ""
            if process.stderr is not None:
                stderr_text = process.stderr.read().decode(
                    "utf-8",
                    errors="replace",
                ).strip()
            return_code = process.wait()
            if (
                return_code != 0
                and not raised_in_context
                and return_code not in _PIGZ_SIGPIPE_CODES
            ):
                detail = f": {stderr_text}" if stderr_text else ""
                raise RuntimeError(
                    f"pigz failed while reading {file_path} "
                    f"with exit code {return_code}{detail}"
                )


def write_tsv(
    *,
    records: Iterable[Mapping[str, object]],
    output_path: str | Path,
    fieldnames: list[str] | None = None,
) -> None:
    """Write dictionaries to a tab-separated file.

    Parameters
    ----------
    records : iterable of mappings
        Records to write.
    output_path : str or pathlib.Path
        Output path. May end in ``.gz``.
    fieldnames : list[str] | None, optional
        Optional explicit column order. If omitted, columns are inferred from
        the first record.
    """
    record_list = list(records)
    if not record_list and fieldnames is None:
        raise ValueError("Cannot infer TSV columns from zero records")
    columns = fieldnames or list(record_list[0].keys())

    with open_text(output_path, "wt") as handle:
        handle.write("\t".join(columns) + "\n")
        for record in record_list:
            values = [str(record.get(column, "")) for column in columns]
            handle.write("\t".join(values) + "\n")


def read_tsv(*, input_path: str | Path) -> list[dict[str, str]]:
    """Read a tab-separated file into dictionaries.

    Parameters
    ----------
    input_path : str or pathlib.Path
        Path to a TSV or TSV.GZ file.

    Returns
    -------
    list[dict[str, str]]
        Parsed records.
    """
    with open_text(input_path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        if not header or header == [""]:
            return []
        records: list[dict[str, str]] = []
        for line in handle:
            values = line.rstrip("\n").split("\t")
            if len(values) < len(header):
                values = values + ([""] * (len(header) - len(values)))
            elif len(values) > len(header):
                values = values[: len(header)]
            records.append(dict(zip(header, values)))
    return records


def write_json(*, data: object, output_path: str | Path) -> None:
    """Write JSON with stable indentation.

    Parameters
    ----------
    data : object
        JSON-serialisable object.
    output_path : str or pathlib.Path
        Output path.
    """
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def write_tsv_stream(
    *,
    records: Iterable[Mapping[str, object]],
    output_path: str | Path,
    fieldnames: list[str],
) -> int:
    """Write records to a tab-separated file without materialising them.

    Parameters
    ----------
    records : iterable of mappings
        Records to write.
    output_path : str or pathlib.Path
        Output path. May end in ``.gz``.
    fieldnames : list[str]
        Explicit output column order.

    Returns
    -------
    int
        Number of data records written.
    """
    n_records = 0
    with open_text(output_path, "wt") as handle:
        handle.write("\t".join(fieldnames) + "\n")
        for record in records:
            values = [str(record.get(column, "")) for column in fieldnames]
            handle.write("\t".join(values) + "\n")
            n_records += 1
    return n_records
