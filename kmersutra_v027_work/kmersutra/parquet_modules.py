"""Optional Parquet module storage and global revalidation helpers.

KmerSutra can build taxonomic modules separately, for example a Plasmodium
module, a host/background module and a bacterial module. To make those modules
mergeable without retaining very large SQLite databases, this module provides an
optional Parquet/Arrow storage layer for the global k-mer source index.

The Parquet dependency is intentionally optional. Importing KmerSutra does not
require ``pyarrow``. Commands that request Parquet export/import raise a clear
error when the optional dependency is not installed.
"""

from __future__ import annotations

import logging
import sqlite3
from collections import defaultdict
from collections.abc import Iterable, Iterator, Sequence
from pathlib import Path
from typing import Any

from kmersutra.global_candidate_evidence import (
    initialise_global_candidate_database,
)
from kmersutra.io import write_json
from kmersutra.target_evidence import _connect

GLOBAL_KMER_COLUMNS = [
    "k",
    "kmer",
    "species_names",
    "genome_ids",
    "contig_ids",
    "taxids",
    "clades",
    "roles",
    "example_position",
]

RETAINED_KMER_COLUMNS = [
    "k",
    "kmer",
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

BUILD_EVENT_COLUMNS = [
    "event_order",
    "timestamp_epoch",
    "stage",
    "detail",
    "n_records",
]


class OptionalParquetDependencyError(ImportError):
    """Raised when Parquet functionality is requested without pyarrow."""


def require_pyarrow() -> tuple[Any, Any]:
    """Import pyarrow lazily and return the table/parquet modules.

    Returns
    -------
    tuple[Any, Any]
        ``pyarrow`` and ``pyarrow.parquet`` modules.

    Raises
    ------
    OptionalParquetDependencyError
        If ``pyarrow`` is not installed.
    """
    try:
        import pyarrow as pa  # type: ignore[import-not-found]
        import pyarrow.parquet as pq  # type: ignore[import-not-found]
    except ImportError as exc:
        raise OptionalParquetDependencyError(
            "Parquet module storage requires the optional pyarrow dependency. "
            "Install KmerSutra with the parquet extra, for example: "
            "pip install '.[parquet]', or install pyarrow in the active environment."
        ) from exc
    return pa, pq


def split_semicolon(value: object) -> set[str]:
    """Split a semicolon-separated field into cleaned unique tokens.

    Parameters
    ----------
    value : object
        Value to split.

    Returns
    -------
    set[str]
        Cleaned non-empty values.
    """
    if value is None:
        return set()
    return {
        item.strip()
        for item in str(value).split(";")
        if item is not None and item.strip()
    }


def join_semicolon(values: Iterable[object]) -> str:
    """Join cleaned values as a deterministic semicolon-separated string.

    Parameters
    ----------
    values : iterable of object
        Values to join.

    Returns
    -------
    str
        Sorted semicolon-separated string.
    """
    clean_values = {
        str(value).strip()
        for value in values
        if value is not None and str(value).strip()
    }
    return ";".join(sorted(clean_values))


def normalise_global_kmer_record(record: dict[str, object]) -> dict[str, object]:
    """Normalise one global source-index record.

    Parameters
    ----------
    record : dict[str, object]
        Input record, typically from SQLite or Parquet.

    Returns
    -------
    dict[str, object]
        Record with the global_kmers schema and normalised types.
    """
    output = {column: record.get(column, "") for column in GLOBAL_KMER_COLUMNS}
    output["k"] = int(output["k"] or 0)
    output["kmer"] = str(output["kmer"] or "")
    output["example_position"] = int(output["example_position"] or 0)
    for column in [
        "species_names",
        "genome_ids",
        "contig_ids",
        "taxids",
        "clades",
        "roles",
    ]:
        output[column] = join_semicolon(split_semicolon(output[column]))
    return output


def merge_global_kmer_record_group(
    records: Iterable[dict[str, object]],
) -> dict[str, object]:
    """Merge records that describe the same ``(k, kmer)`` source key.

    Parameters
    ----------
    records : iterable of dict[str, object]
        Records sharing the same ``k`` and ``kmer``.

    Returns
    -------
    dict[str, object]
        Merged global source-index record.

    Raises
    ------
    ValueError
        If no records are supplied or records do not share the same key.
    """
    normalised = [normalise_global_kmer_record(record) for record in records]
    if not normalised:
        raise ValueError("Cannot merge an empty global k-mer record group")

    key = (normalised[0]["k"], normalised[0]["kmer"])
    for record in normalised:
        if (record["k"], record["kmer"]) != key:
            raise ValueError("All records in a group must share k and kmer")

    merged: dict[str, object] = {
        "k": int(key[0]),
        "kmer": str(key[1]),
        "example_position": min(int(record["example_position"]) for record in normalised),
    }
    for column in [
        "species_names",
        "genome_ids",
        "contig_ids",
        "taxids",
        "clades",
        "roles",
    ]:
        values: set[str] = set()
        for record in normalised:
            values.update(split_semicolon(record[column]))
        merged[column] = join_semicolon(values)
    return merged


def merge_global_kmer_records(
    records: Iterable[dict[str, object]],
) -> list[dict[str, object]]:
    """Merge duplicate global source-index records in memory.

    This helper is intended for tests and small module manifests. Large module
    merges should use SQLite-backed upserts via
    :func:`import_global_kmer_parquets_to_sqlite`.

    Parameters
    ----------
    records : iterable of dict[str, object]
        Global source-index records.

    Returns
    -------
    list[dict[str, object]]
        Merged records sorted by ``k`` and ``kmer``.
    """
    groups: dict[tuple[int, str], list[dict[str, object]]] = defaultdict(list)
    for record in records:
        normalised = normalise_global_kmer_record(record)
        groups[(int(normalised["k"]), str(normalised["kmer"]))].append(normalised)
    return [
        merge_global_kmer_record_group(records=groups[key])
        for key in sorted(groups)
    ]


def upsert_global_kmer_record(
    *,
    connection: sqlite3.Connection,
    record: dict[str, object],
) -> None:
    """Insert or merge one global source-index record into SQLite.

    Parameters
    ----------
    connection : sqlite3.Connection
        Open SQLite connection using the global candidate schema.
    record : dict[str, object]
        Global source-index record.
    """
    clean = normalise_global_kmer_record(record)
    connection.execute(
        """
        INSERT INTO global_kmers(
            k, kmer, species_names, genome_ids, contig_ids, taxids, clades,
            roles, example_position
        )
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
        ON CONFLICT(k, kmer) DO UPDATE SET
            species_names = CASE
                WHEN excluded.species_names = '' THEN species_names
                WHEN species_names = '' THEN excluded.species_names
                ELSE species_names || ';' || excluded.species_names
            END,
            genome_ids = CASE
                WHEN excluded.genome_ids = '' THEN genome_ids
                WHEN genome_ids = '' THEN excluded.genome_ids
                ELSE genome_ids || ';' || excluded.genome_ids
            END,
            contig_ids = CASE
                WHEN excluded.contig_ids = '' THEN contig_ids
                WHEN contig_ids = '' THEN excluded.contig_ids
                ELSE contig_ids || ';' || excluded.contig_ids
            END,
            taxids = CASE
                WHEN excluded.taxids = '' THEN taxids
                WHEN taxids = '' THEN excluded.taxids
                ELSE taxids || ';' || excluded.taxids
            END,
            clades = CASE
                WHEN excluded.clades = '' THEN clades
                WHEN clades = '' THEN excluded.clades
                ELSE clades || ';' || excluded.clades
            END,
            roles = CASE
                WHEN excluded.roles = '' THEN roles
                WHEN roles = '' THEN excluded.roles
                ELSE roles || ';' || excluded.roles
            END,
            example_position = MIN(example_position, excluded.example_position)
        """,
        (
            int(clean["k"]),
            str(clean["kmer"]),
            str(clean["species_names"]),
            str(clean["genome_ids"]),
            str(clean["contig_ids"]),
            str(clean["taxids"]),
            str(clean["clades"]),
            str(clean["roles"]),
            int(clean["example_position"]),
        ),
    )


def deduplicate_global_kmer_sqlite(*, sqlite_path: str | Path) -> None:
    """Deduplicate semicolon fields after module upserts.

    Parameters
    ----------
    sqlite_path : str or pathlib.Path
        SQLite database containing the merged global source index.
    """
    connection = _connect(sqlite_path=sqlite_path)
    connection.row_factory = sqlite3.Row
    try:
        cursor = connection.execute(
            """
            SELECT k, kmer, species_names, genome_ids, contig_ids, taxids,
                   clades, roles, example_position
            FROM global_kmers
            """
        )
        rows = cursor.fetchall()
        for row in rows:
            clean = normalise_global_kmer_record(dict(row))
            connection.execute(
                """
                UPDATE global_kmers
                SET species_names = ?, genome_ids = ?, contig_ids = ?,
                    taxids = ?, clades = ?, roles = ?, example_position = ?
                WHERE k = ? AND kmer = ?
                """,
                (
                    clean["species_names"],
                    clean["genome_ids"],
                    clean["contig_ids"],
                    clean["taxids"],
                    clean["clades"],
                    clean["roles"],
                    int(clean["example_position"]),
                    int(clean["k"]),
                    str(clean["kmer"]),
                ),
            )
        connection.commit()
    finally:
        connection.close()


def iter_sqlite_table_records(
    *,
    sqlite_path: str | Path,
    table_name: str,
    columns: Sequence[str],
    batch_size: int = 50000,
) -> Iterator[list[dict[str, object]]]:
    """Yield SQLite table rows in batches.

    Parameters
    ----------
    sqlite_path : str or pathlib.Path
        SQLite database path.
    table_name : str
        Table to read. Must be a trusted internal table name.
    columns : sequence of str
        Columns to select.
    batch_size : int, optional
        Number of rows per batch.

    Yields
    ------
    list[dict[str, object]]
        Batch of records.
    """
    if batch_size <= 0:
        raise ValueError("batch_size must be positive")
    allowed = {"global_kmers", "retained_kmers", "build_events"}
    if table_name not in allowed:
        raise ValueError(f"Unsupported table for export: {table_name}")

    connection = _connect(sqlite_path=sqlite_path)
    connection.row_factory = sqlite3.Row
    try:
        query = f"SELECT {', '.join(columns)} FROM {table_name}"
        cursor = connection.execute(query)
        while True:
            rows = cursor.fetchmany(batch_size)
            if not rows:
                break
            yield [dict(row) for row in rows]
    finally:
        connection.close()


def write_records_to_parquet(
    *,
    records: Iterable[dict[str, object]],
    output_path: str | Path,
    batch_size: int = 50000,
) -> int:
    """Write records to a Parquet file.

    Parameters
    ----------
    records : iterable of dict[str, object]
        Records to write.
    output_path : str or pathlib.Path
        Output Parquet path.
    batch_size : int, optional
        Number of records to buffer per Arrow table.

    Returns
    -------
    int
        Number of records written.
    """
    if batch_size <= 0:
        raise ValueError("batch_size must be positive")
    pa, pq = require_pyarrow()
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    writer = None
    total = 0
    buffer: list[dict[str, object]] = []
    try:
        for record in records:
            buffer.append(dict(record))
            if len(buffer) >= batch_size:
                table = pa.Table.from_pylist(buffer)
                if writer is None:
                    writer = pq.ParquetWriter(str(output), table.schema)
                writer.write_table(table)
                total += len(buffer)
                buffer.clear()
        if buffer:
            table = pa.Table.from_pylist(buffer)
            if writer is None:
                writer = pq.ParquetWriter(str(output), table.schema)
            writer.write_table(table)
            total += len(buffer)
        if writer is None:
            table = pa.Table.from_pylist([])
            pq.write_table(table, str(output))
    finally:
        if writer is not None:
            writer.close()
    return total


def export_sqlite_table_to_parquet(
    *,
    sqlite_path: str | Path,
    table_name: str,
    output_path: str | Path,
    batch_size: int = 50000,
) -> int:
    """Export a supported SQLite table to Parquet.

    Parameters
    ----------
    sqlite_path : str or pathlib.Path
        SQLite database path.
    table_name : str
        Supported table name.
    output_path : str or pathlib.Path
        Output Parquet path.
    batch_size : int, optional
        SQLite/Arrow batch size.

    Returns
    -------
    int
        Number of exported records.
    """
    columns_by_table = {
        "global_kmers": GLOBAL_KMER_COLUMNS,
        "retained_kmers": RETAINED_KMER_COLUMNS,
        "build_events": BUILD_EVENT_COLUMNS,
    }
    if table_name not in columns_by_table:
        raise ValueError(f"Unsupported table for export: {table_name}")

    def iter_records() -> Iterator[dict[str, object]]:
        for batch in iter_sqlite_table_records(
            sqlite_path=sqlite_path,
            table_name=table_name,
            columns=columns_by_table[table_name],
            batch_size=batch_size,
        ):
            yield from batch

    return write_records_to_parquet(
        records=iter_records(),
        output_path=output_path,
        batch_size=batch_size,
    )


def export_global_candidate_module(
    *,
    sqlite_path: str | Path,
    module_dir: str | Path,
    module_name: str,
    metadata: dict[str, object] | None = None,
    batch_size: int = 50000,
    logger: logging.Logger | None = None,
) -> dict[str, object]:
    """Export global candidate SQLite tables as a shareable Parquet module.

    Parameters
    ----------
    sqlite_path : str or pathlib.Path
        SQLite database created by ``--global_candidate_evidence``.
    module_dir : str or pathlib.Path
        Output module directory.
    module_name : str
        Human-readable module name.
    metadata : dict[str, object] or None, optional
        Additional metadata to write to ``module_metadata.json``.
    batch_size : int, optional
        Export batch size.
    logger : logging.Logger or None, optional
        Logger.

    Returns
    -------
    dict[str, object]
        Module metadata including exported row counts and output paths.
    """
    require_pyarrow()
    output_dir = Path(module_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    outputs = {
        "global_kmers": output_dir / "global_kmers.parquet",
        "retained_kmers": output_dir / "retained_kmers.parquet",
        "build_events": output_dir / "build_events.parquet",
    }
    counts: dict[str, int] = {}
    for table_name, output_path in outputs.items():
        if logger:
            logger.info("Exporting %s to %s", table_name, output_path)
        counts[table_name] = export_sqlite_table_to_parquet(
            sqlite_path=sqlite_path,
            table_name=table_name,
            output_path=output_path,
            batch_size=batch_size,
        )

    module_metadata: dict[str, object] = {
        "module_name": module_name,
        "sqlite_path": str(sqlite_path),
        "global_kmers_parquet": str(outputs["global_kmers"]),
        "retained_kmers_parquet": str(outputs["retained_kmers"]),
        "build_events_parquet": str(outputs["build_events"]),
        "n_global_kmers": counts["global_kmers"],
        "n_retained_kmers": counts["retained_kmers"],
        "n_build_events": counts["build_events"],
    }
    if metadata:
        module_metadata.update(metadata)
    write_json(data=module_metadata, output_path=output_dir / "module_metadata.json")
    return module_metadata


def read_parquet_records(
    *,
    parquet_path: str | Path,
    batch_size: int = 50000,
) -> Iterator[dict[str, object]]:
    """Yield records from a Parquet file in batches.

    Parameters
    ----------
    parquet_path : str or pathlib.Path
        Input Parquet path.
    batch_size : int, optional
        Number of rows per Arrow batch.

    Yields
    ------
    dict[str, object]
        Parquet row as a dictionary.
    """
    if batch_size <= 0:
        raise ValueError("batch_size must be positive")
    _pa, pq = require_pyarrow()
    path = Path(parquet_path)
    parquet_file = pq.ParquetFile(str(path))
    for batch in parquet_file.iter_batches(batch_size=batch_size):
        table = batch.to_pydict()
        keys = list(table)
        n_rows = len(next(iter(table.values()))) if table else 0
        for row_index in range(n_rows):
            yield {key: table[key][row_index] for key in keys}


def resolve_global_kmer_parquet_paths(
    *,
    module_dirs: Sequence[str | Path] | None = None,
    global_kmer_parquets: Sequence[str | Path] | None = None,
) -> list[Path]:
    """Resolve module directories and explicit Parquet paths.

    Parameters
    ----------
    module_dirs : sequence of str or pathlib.Path or None, optional
        Module directories containing ``global_kmers.parquet``.
    global_kmer_parquets : sequence of str or pathlib.Path or None, optional
        Explicit global-kmer Parquet files.

    Returns
    -------
    list[pathlib.Path]
        Existing Parquet paths.

    Raises
    ------
    FileNotFoundError
        If any requested path is missing.
    ValueError
        If no inputs are provided.
    """
    paths: list[Path] = []
    for module_dir in module_dirs or []:
        path = Path(module_dir) / "global_kmers.parquet"
        if not path.is_file():
            raise FileNotFoundError(f"Missing module global_kmers.parquet: {path}")
        paths.append(path)
    for parquet in global_kmer_parquets or []:
        path = Path(parquet)
        if not path.is_file():
            raise FileNotFoundError(f"Missing global-kmer Parquet file: {path}")
        paths.append(path)
    if not paths:
        raise ValueError("At least one module directory or global-kmer Parquet path is required")
    return paths


def import_global_kmer_parquets_to_sqlite(
    *,
    parquet_paths: Sequence[str | Path],
    sqlite_path: str | Path,
    batch_size: int = 50000,
    logger: logging.Logger | None = None,
) -> list[dict[str, object]]:
    """Import and merge global-kmer Parquet modules into SQLite.

    Parameters
    ----------
    parquet_paths : sequence of str or pathlib.Path
        Parquet files containing ``global_kmers`` source-index records.
    sqlite_path : str or pathlib.Path
        Output SQLite database path for the merged source index.
    batch_size : int, optional
        Import batch size.
    logger : logging.Logger or None, optional
        Logger.

    Returns
    -------
    list[dict[str, object]]
        Per-module import summary records.
    """
    if batch_size <= 0:
        raise ValueError("batch_size must be positive")
    require_pyarrow()
    output = Path(sqlite_path)
    if output.exists():
        output.unlink()
    initialise_global_candidate_database(sqlite_path=output)
    connection = _connect(sqlite_path=output)
    summaries: list[dict[str, object]] = []
    try:
        for parquet_path in parquet_paths:
            path = Path(parquet_path)
            if logger:
                logger.info("Importing global-kmer module table: %s", path)
            n_records = 0
            committed_at = 0
            for record in read_parquet_records(parquet_path=path, batch_size=batch_size):
                upsert_global_kmer_record(connection=connection, record=record)
                n_records += 1
                if n_records - committed_at >= batch_size:
                    connection.commit()
                    committed_at = n_records
            connection.commit()
            summaries.append(
                {
                    "module_parquet": str(path),
                    "imported_records": n_records,
                }
            )
            if logger:
                logger.info("Imported %d module rows from %s", n_records, path)
        connection.commit()
    finally:
        connection.close()
    deduplicate_global_kmer_sqlite(sqlite_path=output)
    return summaries
