"""Input and output helpers for KmerSutra."""

from __future__ import annotations

import gzip
import json
from pathlib import Path
from typing import Iterable, Mapping, TextIO


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
            records.append(dict(zip(header, values, strict=False)))
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
