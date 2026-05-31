"""Generic tabular I/O helpers for KmerSutra result tables.

The package writes portable tab-separated text by default, but some KmerSutra
outputs can become large enough that compressed TSV or Parquet storage is more
appropriate. This module centralises format detection and column validation so
callers do not need to duplicate suffix handling.
"""

from __future__ import annotations

import logging
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path

from kmersutra.io import read_tsv, write_tsv
from kmersutra.table_parquet import read_records_parquet, write_records_parquet

SUPPORTED_TABLE_FORMATS = {"tsv", "tsv.gz", "parquet"}


def infer_table_format(*, path: str | Path) -> str:
    """Infer the table format from a file suffix.

    Parameters
    ----------
    path : str or pathlib.Path
        Input or output path.

    Returns
    -------
    str
        One of ``tsv``, ``tsv.gz`` or ``parquet``.

    Raises
    ------
    ValueError
        If the suffix is unsupported.
    """
    file_path = Path(path)
    suffixes = [suffix.lower() for suffix in file_path.suffixes]
    if suffixes[-2:] == [".tsv", ".gz"]:
        return "tsv.gz"
    if suffixes and suffixes[-1] == ".tsv":
        return "tsv"
    if suffixes and suffixes[-1] == ".parquet":
        return "parquet"
    raise ValueError(
        "Unsupported table format for path "
        f"{file_path}. Supported suffixes are .tsv, .tsv.gz and .parquet."
    )


def validate_required_columns(
    *,
    records: Sequence[Mapping[str, object]],
    required_columns: Sequence[str] | None = None,
    source_path: str | Path | None = None,
) -> None:
    """Validate that records contain required columns.

    Parameters
    ----------
    records : sequence of mappings
        Records to validate.
    required_columns : sequence of str or None, optional
        Required column names. If omitted, no validation is performed.
    source_path : str or pathlib.Path or None, optional
        Optional path included in error messages.

    Raises
    ------
    ValueError
        If one or more required columns are absent.
    """
    required = list(required_columns or [])
    if not required:
        return
    observed = set(records[0].keys()) if records else set()
    missing = [column for column in required if column not in observed]
    if missing:
        location = f" in {source_path}" if source_path is not None else ""
        raise ValueError(
            "Missing required table column(s)" + location + ": " + ", ".join(missing)
        )


def read_records_table(
    *,
    input_path: str | Path,
    required_columns: Sequence[str] | None = None,
    logger: logging.Logger | None = None,
) -> list[dict[str, str]]:
    """Read a KmerSutra records table from TSV, TSV.GZ or Parquet.

    Parameters
    ----------
    input_path : str or pathlib.Path
        Input table path. Supported suffixes are ``.tsv``, ``.tsv.gz`` and
        ``.parquet``.
    required_columns : sequence of str or None, optional
        Columns that must be present in the parsed table.
    logger : logging.Logger or None, optional
        Optional logger used for shape and format messages.

    Returns
    -------
    list[dict[str, str]]
        Parsed records with string values.
    """
    table_format = infer_table_format(path=input_path)
    if table_format == "parquet":
        records = read_records_parquet(input_path=input_path)
    else:
        records = read_tsv(input_path=input_path)
    validate_required_columns(
        records=records,
        required_columns=required_columns,
        source_path=input_path,
    )
    if logger:
        logger.info(
            "Read %d records from %s table %s",
            len(records),
            table_format,
            input_path,
        )
    return records


def write_records_table(
    *,
    records: Iterable[Mapping[str, object]],
    output_path: str | Path,
    fieldnames: Sequence[str] | None = None,
    logger: logging.Logger | None = None,
) -> int:
    """Write a KmerSutra records table as TSV, TSV.GZ or Parquet.

    Parameters
    ----------
    records : iterable of mappings
        Records to write.
    output_path : str or pathlib.Path
        Output table path. Supported suffixes are ``.tsv``, ``.tsv.gz`` and
        ``.parquet``.
    fieldnames : sequence of str or None, optional
        Optional output column order.
    logger : logging.Logger or None, optional
        Optional logger used for shape and format messages.

    Returns
    -------
    int
        Number of records written.
    """
    table_format = infer_table_format(path=output_path)
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    record_list = [dict(record) for record in records]
    columns = list(fieldnames) if fieldnames is not None else None
    if table_format == "parquet":
        n_records = write_records_parquet(
            records=record_list,
            output_path=output,
            fieldnames=columns,
        )
    else:
        write_tsv(records=record_list, output_path=output, fieldnames=columns)
        n_records = len(record_list)
    if logger:
        logger.info(
            "Wrote %d records to %s table %s",
            n_records,
            table_format,
            output,
        )
    return n_records
