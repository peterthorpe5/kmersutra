"""Optional Parquet-backed diagnostic panel helpers.

The primary screening path needs an in-memory Python lookup index, regardless
of whether the source panel is stored as TSV.GZ or Parquet. Parquet therefore
mainly improves repeated loading, module storage and columnar interchange. The
helpers in this module keep ``pyarrow`` optional and provide clear failures when
Parquet is requested without the optional dependency.
"""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path
from typing import Any

from kmersutra.parquet_modules import OptionalParquetDependencyError, require_pyarrow

PANEL_FIELDNAMES = [
    "kmer",
    "k",
    "panel_type",
    "species_name",
    "clade",
    "source_genomes",
    "source_contigs",
    "example_position",
    "evidence_taxid",
    "evidence_name",
    "evidence_rank",
    "lineage_taxids",
    "source_taxids",
]


def is_parquet_path(*, path: str | Path) -> bool:
    """Return True if ``path`` looks like a Parquet panel path.

    Parameters
    ----------
    path : str or pathlib.Path
        File path to inspect.

    Returns
    -------
    bool
        True when the suffix is ``.parquet``.
    """
    return Path(path).suffix.lower() == ".parquet"


def pyarrow_available() -> bool:
    """Return whether the optional ``pyarrow`` dependency is importable."""
    try:
        require_pyarrow()
    except OptionalParquetDependencyError:
        return False
    return True


def records_to_arrow_table(*, records: Iterable[dict[str, object]]) -> Any:
    """Convert panel records to a pyarrow table with stable string columns.

    Parameters
    ----------
    records : iterable of dict[str, object]
        Panel records to serialise.

    Returns
    -------
    pyarrow.Table
        Arrow table with the KmerSutra panel schema.
    """
    pa, _ = require_pyarrow()
    columns: dict[str, list[str]] = {column: [] for column in PANEL_FIELDNAMES}
    for record in records:
        for column in PANEL_FIELDNAMES:
            columns[column].append(str(record.get(column, "")))
    return pa.table(columns)


def write_panel_parquet(
    *,
    records: Iterable[dict[str, object]],
    output_path: str | Path,
) -> int:
    """Write diagnostic panel records to Parquet.

    Parameters
    ----------
    records : iterable of dict[str, object]
        Records with KmerSutra panel fields.
    output_path : str or pathlib.Path
        Parquet output path.

    Returns
    -------
    int
        Number of records written.

    Raises
    ------
    OptionalParquetDependencyError
        If ``pyarrow`` is not installed.
    """
    _, pq = require_pyarrow()
    record_list = [
        {column: str(record.get(column, "")) for column in PANEL_FIELDNAMES}
        for record in records
    ]
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    table = records_to_arrow_table(records=record_list)
    pq.write_table(table, output)
    return len(record_list)


def read_panel_parquet(*, input_path: str | Path) -> list[dict[str, str]]:
    """Read KmerSutra diagnostic panel records from Parquet.

    Parameters
    ----------
    input_path : str or pathlib.Path
        Parquet file to read.

    Returns
    -------
    list[dict[str, str]]
        Panel records as strings.

    Raises
    ------
    FileNotFoundError
        If the Parquet path is missing or empty.
    ValueError
        If required panel columns are absent.
    OptionalParquetDependencyError
        If ``pyarrow`` is not installed.
    """
    _, pq = require_pyarrow()
    path = Path(input_path)
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"Parquet panel file is missing or empty: {path}")
    table = pq.read_table(path)
    names = set(table.column_names)
    missing = [column for column in PANEL_FIELDNAMES if column not in names]
    if missing:
        raise ValueError(
            "Parquet panel is missing required columns: " + ", ".join(missing)
        )
    data = table.select(PANEL_FIELDNAMES).to_pydict()
    n_rows = table.num_rows
    return [
        {column: str(data[column][row] if data[column][row] is not None else "") for column in PANEL_FIELDNAMES}
        for row in range(n_rows)
    ]


def derive_panel_parquet_path(*, panel_path: str | Path) -> Path:
    """Return the default Parquet companion path for a TSV panel."""
    path = Path(panel_path)
    name = path.name
    if name.endswith(".tsv.gz"):
        name = name[:-7]
    elif name.endswith(".tsv"):
        name = name[:-4]
    return path.with_name(f"{name}.parquet")


def select_screen_panel_path(
    *,
    tsv_path: str | Path,
    parquet_path: str | Path | None,
    panel_storage_format: str,
    prefer_parquet_if_available: bool = True,
) -> Path:
    """Choose which panel path should be placed in screen manifests.

    Parameters
    ----------
    tsv_path : str or pathlib.Path
        TSV.GZ fallback panel path.
    parquet_path : str or pathlib.Path or None
        Optional Parquet companion path.
    panel_storage_format : str
        One of ``tsv``, ``parquet``, ``both`` or ``auto``.
    prefer_parquet_if_available : bool, optional
        If True, ``auto`` and ``both`` prefer an existing Parquet path.

    Returns
    -------
    pathlib.Path
        Selected panel path.

    Raises
    ------
    ValueError
        If the storage format is unsupported or a required Parquet path is absent.
    """
    fmt = str(panel_storage_format).lower()
    tsv = Path(tsv_path)
    parquet = Path(parquet_path) if parquet_path else None
    if fmt == "tsv":
        return tsv
    if fmt == "parquet":
        if parquet is None or not parquet.exists():
            raise ValueError("Parquet panel was requested but no Parquet path exists")
        return parquet
    if fmt in {"auto", "both"}:
        if prefer_parquet_if_available and parquet is not None and parquet.exists():
            return parquet
        return tsv
    raise ValueError(f"Unsupported panel storage format: {panel_storage_format}")
