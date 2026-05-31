"""Optional Parquet helpers for generic KmerSutra result tables."""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path

from kmersutra.parquet_modules import OptionalParquetDependencyError, require_pyarrow


def infer_fieldnames(*, records: Sequence[Mapping[str, object]]) -> list[str]:
    """Infer stable field names from generic record dictionaries.

    Parameters
    ----------
    records : sequence of mappings
        Records to inspect.

    Returns
    -------
    list[str]
        First-seen union of field names.
    """
    fieldnames: list[str] = []
    seen: set[str] = set()
    for record in records:
        for key in record:
            if key not in seen:
                seen.add(str(key))
                fieldnames.append(str(key))
    return fieldnames


def records_to_string_columns(
    *,
    records: Iterable[Mapping[str, object]],
    fieldnames: Sequence[str] | None = None,
) -> tuple[list[str], dict[str, list[str]]]:
    """Convert records to string columns suitable for Arrow serialisation.

    Parameters
    ----------
    records : iterable of mappings
        Input records.
    fieldnames : sequence of str or None, optional
        Explicit column order. If omitted, fields are inferred.

    Returns
    -------
    tuple[list[str], dict[str, list[str]]]
        Column order and values.
    """
    record_list = [dict(record) for record in records]
    columns = list(fieldnames) if fieldnames is not None else infer_fieldnames(records=record_list)
    values: dict[str, list[str]] = {column: [] for column in columns}
    for record in record_list:
        for column in columns:
            value = record.get(column, "")
            values[column].append("" if value is None else str(value))
    return columns, values


def write_records_parquet(
    *,
    records: Iterable[Mapping[str, object]],
    output_path: str | Path,
    fieldnames: Sequence[str] | None = None,
) -> int:
    """Write generic records to Parquet using optional ``pyarrow``.

    Parameters
    ----------
    records : iterable of mappings
        Records to write.
    output_path : str or pathlib.Path
        Parquet output path.
    fieldnames : sequence of str or None, optional
        Optional explicit column order.

    Returns
    -------
    int
        Number of records written.

    Raises
    ------
    OptionalParquetDependencyError
        If ``pyarrow`` is not installed.
    """
    pa, pq = require_pyarrow()
    record_list = [dict(record) for record in records]
    columns, values = records_to_string_columns(
        records=record_list,
        fieldnames=fieldnames,
    )
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    table = pa.table({column: values[column] for column in columns})
    pq.write_table(table, output)
    return len(record_list)


def read_records_parquet(*, input_path: str | Path) -> list[dict[str, str]]:
    """Read generic KmerSutra records from Parquet.

    Parameters
    ----------
    input_path : str or pathlib.Path
        Parquet input path.

    Returns
    -------
    list[dict[str, str]]
        Records with values converted to strings.

    Raises
    ------
    FileNotFoundError
        If the path is missing or empty.
    OptionalParquetDependencyError
        If ``pyarrow`` is not installed.
    """
    _pa, pq = require_pyarrow()
    path = Path(input_path)
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"Parquet table missing or empty: {path}")
    table = pq.read_table(path)
    data = table.to_pydict()
    rows = table.num_rows
    return [
        {
            column: "" if data[column][row_index] is None else str(data[column][row_index])
            for column in table.column_names
        }
        for row_index in range(rows)
    ]
