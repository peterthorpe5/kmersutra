"""Global SQLite-backed all-candidate evidence building for KmerSutra.

This module implements the scalable successor to the v0.13 all-candidate
builder. The v0.13 builder was biologically correct, but it validated each
candidate species by rescanning every other genome. That made the build behave
like an all-versus-all loop and became too slow for realistic panels.

The global builder supports both exhaustive and bounded candidate-universe
construction. Exhaustive modes read each genome once per requested k value and
store source metadata for every distinct ``(k, kmer)`` key. The scalable
``candidate_universe`` mode first samples genome-spread candidate markers from
each genome/bin, then rescans references only to annotate conflicts for that
bounded candidate universe. Both routes assign k-mers to the most specific
supported taxonomic evidence level using observed source taxids.
"""

from __future__ import annotations

import logging
import sqlite3
import time
from collections import defaultdict
from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from pathlib import Path

from kmersutra.build_panel import DiagnosticKmer
from kmersutra.config import GenomeConfig
from kmersutra.fasta import read_fasta_records
from kmersutra.kmers import iter_kmers
from kmersutra.target_evidence import _connect
from kmersutra.taxonomy import CORE_RANK_ORDER, TaxonomyDatabase

DEFAULT_EXCLUDED_REPORTABLE_ROLES = {
    "exclude",
    "host",
    "host_or_background",
    "background_host",
    "host_background",
    "background",
    "environmental_background",
}

VALID_GLOBAL_SOURCE_INDEX_MODES = {
    "source_rows",
    "aggregated",
    "candidate_universe",
}


@dataclass(frozen=True)
class GlobalCandidateEvidenceBuildResult:
    """Result metadata from a global all-candidate evidence build.

    Attributes
    ----------
    sqlite_path : pathlib.Path
        SQLite database containing the global source index and retained
        diagnostic evidence.
    collection_summary : list[dict[str, object]]
        Per-genome collection summaries and final evidence summaries.
    build_summary : list[dict[str, object]]
        Coarse build summary records suitable for TSV output.
    panel_summary : list[dict[str, object]]
        Diagnostic k-mer counts grouped by evidence bucket.
    """

    sqlite_path: Path
    collection_summary: list[dict[str, object]]
    build_summary: list[dict[str, object]]
    panel_summary: list[dict[str, object]]


def _normalise_roles(roles: Iterable[str] | None) -> set[str]:
    """Return a cleaned set of role labels.

    Parameters
    ----------
    roles : iterable of str or None
        Role names supplied by a caller.

    Returns
    -------
    set[str]
        Non-empty role labels.
    """
    if roles is None:
        return set()
    return {str(role) for role in roles if str(role)}


def _normalise_semicolon_set(value: str | None) -> set[str]:
    """Convert a semicolon-separated value into a set.

    Parameters
    ----------
    value : str or None
        Existing semicolon-separated value.

    Returns
    -------
    set[str]
        Non-empty items.
    """
    if not value:
        return set()
    return {item for item in str(value).split(";") if item}


def _merge_semicolon_values(existing: str, new_value: str) -> str:
    """Merge a value into a semicolon-separated string.

    Parameters
    ----------
    existing : str
        Existing semicolon-separated values.
    new_value : str
        New value to add.

    Returns
    -------
    str
        Deterministically sorted semicolon-separated values.
    """
    values = _normalise_semicolon_set(existing)
    if new_value:
        values.add(str(new_value))
    return ";".join(sorted(values))


def _join_semicolon_from_csv(value: str | None) -> str:
    """Convert SQLite comma group-concat output into semicolon-separated text.

    Parameters
    ----------
    value : str or None
        Comma-separated SQLite ``GROUP_CONCAT`` value.

    Returns
    -------
    str
        Semicolon-separated unique values.
    """
    if not value:
        return ""
    return ";".join(sorted({item for item in str(value).split(",") if item}))


def initialise_global_candidate_database(*, sqlite_path: str | Path) -> None:
    """Create the SQLite schema for global all-candidate evidence.

    Parameters
    ----------
    sqlite_path : str or pathlib.Path
        SQLite database path to create.
    """
    connection = _connect(sqlite_path=sqlite_path)
    try:
        connection.execute(
            """
            CREATE TABLE IF NOT EXISTS global_kmers (
                k INTEGER NOT NULL,
                kmer TEXT NOT NULL,
                species_names TEXT NOT NULL DEFAULT '',
                genome_ids TEXT NOT NULL DEFAULT '',
                contig_ids TEXT NOT NULL DEFAULT '',
                taxids TEXT NOT NULL DEFAULT '',
                clades TEXT NOT NULL DEFAULT '',
                roles TEXT NOT NULL DEFAULT '',
                example_position INTEGER NOT NULL DEFAULT 0,
                PRIMARY KEY (k, kmer)
            )
            """
        )
        connection.execute(
            """
            CREATE TABLE IF NOT EXISTS retained_kmers (
                k INTEGER NOT NULL,
                kmer TEXT NOT NULL,
                panel_type TEXT NOT NULL,
                species_name TEXT NOT NULL,
                clade TEXT NOT NULL,
                source_genomes TEXT NOT NULL,
                source_contigs TEXT NOT NULL,
                example_position INTEGER NOT NULL,
                evidence_taxid TEXT NOT NULL,
                evidence_name TEXT NOT NULL,
                evidence_rank TEXT NOT NULL,
                lineage_taxids TEXT NOT NULL,
                source_taxids TEXT NOT NULL,
                PRIMARY KEY (k, kmer, evidence_taxid, evidence_rank)
            )
            """
        )
        connection.execute(
            """
            CREATE TABLE IF NOT EXISTS global_kmer_sources (
                k INTEGER NOT NULL,
                kmer TEXT NOT NULL,
                species_name TEXT NOT NULL DEFAULT '',
                genome_id TEXT NOT NULL DEFAULT '',
                contig_id TEXT NOT NULL DEFAULT '',
                taxid TEXT NOT NULL DEFAULT '',
                clade TEXT NOT NULL DEFAULT '',
                role TEXT NOT NULL DEFAULT '',
                example_position INTEGER NOT NULL DEFAULT 0,
                PRIMARY KEY (k, kmer, genome_id)
            )
            """
        )
        connection.execute(
            """
            CREATE TABLE IF NOT EXISTS candidate_kmers (
                k INTEGER NOT NULL,
                kmer TEXT NOT NULL,
                first_genome_id TEXT NOT NULL DEFAULT '',
                first_contig_id TEXT NOT NULL DEFAULT '',
                first_position INTEGER NOT NULL DEFAULT 0,
                PRIMARY KEY (k, kmer)
            )
            """
        )
        connection.execute(
            """
            CREATE TABLE IF NOT EXISTS build_events (
                event_order INTEGER PRIMARY KEY AUTOINCREMENT,
                timestamp_epoch REAL NOT NULL,
                stage TEXT NOT NULL,
                detail TEXT NOT NULL,
                n_records INTEGER NOT NULL DEFAULT 0
            )
            """
        )
        connection.execute(
            """
            CREATE INDEX IF NOT EXISTS retained_kmers_summary_idx
            ON retained_kmers(evidence_rank, evidence_name, species_name, k)
            """
        )
        connection.execute(
            """
            CREATE INDEX IF NOT EXISTS global_kmer_sources_key_idx
            ON global_kmer_sources(k, kmer)
            """
        )
        connection.execute(
            """
            CREATE INDEX IF NOT EXISTS candidate_kmers_k_idx
            ON candidate_kmers(k)
            """
        )
        connection.commit()
    finally:
        connection.close()


def _record_event(
    *,
    connection: sqlite3.Connection,
    stage: str,
    detail: str,
    n_records: int = 0,
) -> None:
    """Record one build event.

    Parameters
    ----------
    connection : sqlite3.Connection
        Open SQLite connection.
    stage : str
        Build stage.
    detail : str
        Human-readable detail.
    n_records : int, optional
        Relevant record count.
    """
    connection.execute(
        """
        INSERT INTO build_events(timestamp_epoch, stage, detail, n_records)
        VALUES (?, ?, ?, ?)
        """,
        (time.time(), stage, detail, int(n_records)),
    )


GLOBAL_KMER_UPSERT_SQL = """
    INSERT INTO global_kmers(
        k, kmer, species_names, genome_ids, contig_ids, taxids, clades,
        roles, example_position
    )
    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
    ON CONFLICT(k, kmer) DO UPDATE SET
        species_names = CASE
            WHEN instr(';' || species_names || ';', ';' || excluded.species_names || ';') > 0
            THEN species_names
            WHEN species_names = '' THEN excluded.species_names
            ELSE species_names || ';' || excluded.species_names
        END,
        genome_ids = CASE
            WHEN instr(';' || genome_ids || ';', ';' || excluded.genome_ids || ';') > 0
            THEN genome_ids
            WHEN genome_ids = '' THEN excluded.genome_ids
            ELSE genome_ids || ';' || excluded.genome_ids
        END,
        contig_ids = CASE
            WHEN instr(';' || contig_ids || ';', ';' || excluded.contig_ids || ';') > 0
            THEN contig_ids
            WHEN contig_ids = '' THEN excluded.contig_ids
            ELSE contig_ids || ';' || excluded.contig_ids
        END,
        taxids = CASE
            WHEN excluded.taxids = '' THEN taxids
            WHEN instr(';' || taxids || ';', ';' || excluded.taxids || ';') > 0
            THEN taxids
            WHEN taxids = '' THEN excluded.taxids
            ELSE taxids || ';' || excluded.taxids
        END,
        clades = CASE
            WHEN excluded.clades = '' THEN clades
            WHEN instr(';' || clades || ';', ';' || excluded.clades || ';') > 0
            THEN clades
            WHEN clades = '' THEN excluded.clades
            ELSE clades || ';' || excluded.clades
        END,
        roles = CASE
            WHEN excluded.roles = '' THEN roles
            WHEN instr(';' || roles || ';', ';' || excluded.roles || ';') > 0
            THEN roles
            WHEN roles = '' THEN excluded.roles
            ELSE roles || ';' || excluded.roles
        END
"""

GLOBAL_KMER_SOURCE_INSERT_SQL = """
    INSERT OR IGNORE INTO global_kmer_sources(
        k, kmer, species_name, genome_id, contig_id, taxid, clade, role,
        example_position
    )
    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
"""

CANDIDATE_KMER_INSERT_SQL = """
    INSERT OR IGNORE INTO candidate_kmers(
        k, kmer, first_genome_id, first_contig_id, first_position
    )
    VALUES (?, ?, ?, ?, ?)
"""


def configure_fast_global_sqlite(*, connection: sqlite3.Connection) -> None:
    """Apply pragmatic SQLite settings for rebuildable TMPDIR builds.

    Parameters
    ----------
    connection : sqlite3.Connection
        Open SQLite connection used during global database construction.

    Notes
    -----
    These settings are intended for intermediate build databases that can be
    regenerated if a job fails. Final KmerSutra outputs are still written after
    the build completes. The settings reduce fsync and journal overhead, which
    is usually the dominant cost when inserting millions of k-mer observations.
    """
    connection.execute("PRAGMA journal_mode=OFF")
    connection.execute("PRAGMA synchronous=OFF")
    connection.execute("PRAGMA temp_store=MEMORY")
    connection.execute("PRAGMA locking_mode=EXCLUSIVE")
    connection.execute("PRAGMA cache_size=-1000000")


def _global_row_tuple(
    *,
    k: int,
    kmer: str,
    genome_config: GenomeConfig,
    contig_id: str,
    position: int,
) -> tuple[object, ...]:
    """Return an aggregated global-kmer upsert tuple.

    Parameters
    ----------
    k : int
        K-mer length.
    kmer : str
        Canonical k-mer sequence.
    genome_config : GenomeConfig
        Source genome metadata.
    contig_id : str
        Source contig identifier.
    position : int
        Example source position.

    Returns
    -------
    tuple[object, ...]
        SQLite parameter tuple for ``GLOBAL_KMER_UPSERT_SQL``.
    """
    return (
        int(k),
        kmer,
        genome_config.species_name,
        genome_config.genome_id,
        contig_id,
        genome_config.taxid,
        genome_config.clade,
        genome_config.role,
        int(position),
    )


def _source_row_tuple(
    *,
    k: int,
    kmer: str,
    genome_config: GenomeConfig,
    contig_id: str,
    position: int,
) -> tuple[object, ...]:
    """Return a normalised source-row tuple for one genome/k-mer pair.

    Parameters
    ----------
    k : int
        K-mer length.
    kmer : str
        Canonical k-mer sequence.
    genome_config : GenomeConfig
        Source genome metadata.
    contig_id : str
        Source contig identifier.
    position : int
        Example source position.

    Returns
    -------
    tuple[object, ...]
        SQLite parameter tuple for ``GLOBAL_KMER_SOURCE_INSERT_SQL``.
    """
    return (
        int(k),
        kmer,
        genome_config.species_name,
        genome_config.genome_id,
        contig_id,
        genome_config.taxid,
        genome_config.clade,
        genome_config.role,
        int(position),
    )


def _flush_global_rows(
    *,
    connection: sqlite3.Connection,
    rows: list[tuple[object, ...]],
    source_index_mode: str,
) -> int:
    """Flush buffered global source rows to SQLite.

    Parameters
    ----------
    connection : sqlite3.Connection
        Open SQLite connection.
    rows : list[tuple[object, ...]]
        Buffered SQLite parameter tuples.
    source_index_mode : str
        Either ``source_rows`` or ``aggregated``.

    Returns
    -------
    int
        Number of buffered rows flushed.
    """
    if not rows:
        return 0
    if source_index_mode == "source_rows":
        connection.executemany(GLOBAL_KMER_SOURCE_INSERT_SQL, rows)
    elif source_index_mode == "aggregated":
        connection.executemany(GLOBAL_KMER_UPSERT_SQL, rows)
    else:
        raise ValueError(
            "source_index_mode must be one of: "
            + ", ".join(sorted(VALID_GLOBAL_SOURCE_INDEX_MODES))
        )
    flushed = len(rows)
    rows.clear()
    return flushed


def materialise_global_kmers_from_sources(
    *,
    sqlite_path: str | Path,
    logger: logging.Logger | None = None,
) -> dict[str, int]:
    """Materialise the aggregated ``global_kmers`` table from source rows.

    Parameters
    ----------
    sqlite_path : str or pathlib.Path
        SQLite database containing ``global_kmer_sources``.
    source_index_mode : str, optional
        Source-index implementation. ``source_rows`` stores one source row per
        genome/k-mer and materialises aggregated evidence after collection.
        ``aggregated`` preserves the legacy direct-upsert behaviour.
    progress_interval : int, optional
        Attempted k-mer interval for progress logging during genome indexing.
    logger : logging.Logger or None, optional
        Logger for progress messages.

    Returns
    -------
    dict[str, int]
        Counts for source rows and distinct global k-mer keys.
    """
    connection = _connect(sqlite_path=sqlite_path)
    try:
        configure_fast_global_sqlite(connection=connection)
        source_rows = int(
            connection.execute("SELECT COUNT(*) FROM global_kmer_sources").fetchone()[0]
        )
        if logger:
            logger.info(
                "Materialising aggregated global_kmers from %d source row(s)",
                source_rows,
            )
        # Rebuild the materialised table from scratch. This is faster than
        # inserting into the primary-keyed table created for the legacy
        # aggregated mode because SQLite does not have to maintain the global
        # k-mer key index during the large grouped INSERT. The assignment phase
        # only reads this table, so a normal post-build index is sufficient.
        connection.execute("DROP TABLE IF EXISTS global_kmers")
        connection.execute(
            """
            CREATE TABLE global_kmers (
                k INTEGER NOT NULL,
                kmer TEXT NOT NULL,
                species_names TEXT NOT NULL DEFAULT '',
                genome_ids TEXT NOT NULL DEFAULT '',
                contig_ids TEXT NOT NULL DEFAULT '',
                taxids TEXT NOT NULL DEFAULT '',
                clades TEXT NOT NULL DEFAULT '',
                roles TEXT NOT NULL DEFAULT '',
                example_position INTEGER NOT NULL DEFAULT 0
            )
            """
        )
        connection.execute(
            """
            INSERT INTO global_kmers(
                k, kmer, species_names, genome_ids, contig_ids, taxids, clades,
                roles, example_position
            )
            SELECT
                k,
                kmer,
                replace(group_concat(DISTINCT species_name), ',', ';'),
                replace(group_concat(DISTINCT genome_id), ',', ';'),
                replace(group_concat(DISTINCT contig_id), ',', ';'),
                replace(group_concat(DISTINCT taxid), ',', ';'),
                replace(group_concat(DISTINCT clade), ',', ';'),
                replace(group_concat(DISTINCT role), ',', ';'),
                min(example_position)
            FROM global_kmer_sources
            GROUP BY k, kmer
            """
        )
        connection.execute(
            """
            CREATE INDEX IF NOT EXISTS global_kmers_assignment_idx
            ON global_kmers(k, kmer)
            """
        )
        global_rows = int(
            connection.execute("SELECT COUNT(*) FROM global_kmers").fetchone()[0]
        )
        _record_event(
            connection=connection,
            stage="materialise_global_sources",
            detail="source_rows_to_global_kmers",
            n_records=global_rows,
        )
        connection.commit()
        if logger:
            logger.info(
                "Materialised %d distinct global k-mer key(s) from source rows",
                global_rows,
            )
    finally:
        connection.close()
    return {"source_rows": source_rows, "global_kmers": global_rows}


def _upsert_global_kmer(
    *,
    connection: sqlite3.Connection,
    k: int,
    kmer: str,
    genome_config: GenomeConfig,
    contig_id: str,
    position: int,
) -> None:
    """Insert or merge one observed k-mer into the global source index.

    Parameters
    ----------
    connection : sqlite3.Connection
        Open SQLite connection.
    k : int
        K-mer size.
    kmer : str
        Canonical k-mer.
    genome_config : GenomeConfig
        Source genome metadata.
    contig_id : str
        Source contig identifier.
    position : int
        Example position in the source contig.
    """
    connection.execute(
        """
        INSERT INTO global_kmers(
            k, kmer, species_names, genome_ids, contig_ids, taxids, clades,
            roles, example_position
        )
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
        ON CONFLICT(k, kmer) DO UPDATE SET
            species_names = CASE
                WHEN instr(';' || species_names || ';', ';' || excluded.species_names || ';') > 0
                THEN species_names
                WHEN species_names = '' THEN excluded.species_names
                ELSE species_names || ';' || excluded.species_names
            END,
            genome_ids = CASE
                WHEN instr(';' || genome_ids || ';', ';' || excluded.genome_ids || ';') > 0
                THEN genome_ids
                WHEN genome_ids = '' THEN excluded.genome_ids
                ELSE genome_ids || ';' || excluded.genome_ids
            END,
            contig_ids = CASE
                WHEN instr(';' || contig_ids || ';', ';' || excluded.contig_ids || ';') > 0
                THEN contig_ids
                WHEN contig_ids = '' THEN excluded.contig_ids
                ELSE contig_ids || ';' || excluded.contig_ids
            END,
            taxids = CASE
                WHEN excluded.taxids = '' THEN taxids
                WHEN instr(';' || taxids || ';', ';' || excluded.taxids || ';') > 0
                THEN taxids
                WHEN taxids = '' THEN excluded.taxids
                ELSE taxids || ';' || excluded.taxids
            END,
            clades = CASE
                WHEN excluded.clades = '' THEN clades
                WHEN instr(';' || clades || ';', ';' || excluded.clades || ';') > 0
                THEN clades
                WHEN clades = '' THEN excluded.clades
                ELSE clades || ';' || excluded.clades
            END,
            roles = CASE
                WHEN excluded.roles = '' THEN roles
                WHEN instr(';' || roles || ';', ';' || excluded.roles || ';') > 0
                THEN roles
                WHEN roles = '' THEN excluded.roles
                ELSE roles || ';' || excluded.roles
            END
        """,
        (
            int(k),
            kmer,
            genome_config.species_name,
            genome_config.genome_id,
            contig_id,
            genome_config.taxid,
            genome_config.clade,
            genome_config.role,
            int(position),
        ),
    )



def _candidate_row_tuple(
    *,
    k: int,
    kmer: str,
    genome_config: GenomeConfig,
    contig_id: str,
    position: int,
) -> tuple[object, ...]:
    """Return a candidate-universe row for one sampled k-mer.

    Parameters
    ----------
    k : int
        K-mer length.
    kmer : str
        Canonical k-mer sequence.
    genome_config : GenomeConfig
        Source genome metadata.
    contig_id : str
        Source contig identifier.
    position : int
        Zero-based position in the source contig.

    Returns
    -------
    tuple[object, ...]
        SQLite parameter tuple for ``CANDIDATE_KMER_INSERT_SQL``.
    """
    return (
        int(k),
        kmer,
        genome_config.genome_id,
        contig_id,
        int(position),
    )


def _flush_candidate_rows(
    *,
    connection: sqlite3.Connection,
    candidate_rows: list[tuple[object, ...]],
    source_rows: list[tuple[object, ...]],
) -> dict[str, int]:
    """Flush candidate and source rows used by candidate-universe mode.

    Parameters
    ----------
    connection : sqlite3.Connection
        Open SQLite connection.
    candidate_rows : list[tuple[object, ...]]
        Buffered candidate-kmer rows.
    source_rows : list[tuple[object, ...]]
        Buffered source rows for the same sampled candidates or hits.

    Returns
    -------
    dict[str, int]
        Counts of candidate and source rows flushed from Python buffers.
    """
    n_candidate = len(candidate_rows)
    n_source = len(source_rows)
    if candidate_rows:
        connection.executemany(CANDIDATE_KMER_INSERT_SQL, candidate_rows)
        candidate_rows.clear()
    if source_rows:
        connection.executemany(GLOBAL_KMER_SOURCE_INSERT_SQL, source_rows)
        source_rows.clear()
    return {"candidate_rows": n_candidate, "source_rows": n_source}


def _load_candidate_kmers_by_k(
    *,
    connection: sqlite3.Connection,
) -> dict[int, set[str]]:
    """Load the bounded candidate universe into memory grouped by k.

    Parameters
    ----------
    connection : sqlite3.Connection
        Open SQLite connection containing ``candidate_kmers``.

    Returns
    -------
    dict[int, set[str]]
        Candidate k-mers grouped by k value.
    """
    candidates: dict[int, set[str]] = defaultdict(set)
    cursor = connection.execute("SELECT k, kmer FROM candidate_kmers ORDER BY k")
    for k_value, kmer in cursor:
        candidates[int(k_value)].add(str(kmer))
    return dict(candidates)


def collect_candidate_universe_sqlite(
    *,
    genome_configs: Iterable[GenomeConfig],
    k_values: list[int],
    sqlite_path: str | Path,
    batch_size: int = 50000,
    genome_bin_size: int = 10000,
    max_per_genome_bin: int = 10,
    progress_interval: int = 1000000,
    logger: logging.Logger | None = None,
) -> list[dict[str, object]]:
    """Collect a bounded genome-spread candidate k-mer universe.

    This is the scalable alternative to storing every sliding-window k-mer.
    For each genome, contig, k value and genomic bin, at most
    ``max_per_genome_bin`` candidate k-mers are retained. Candidate source rows
    are written immediately so every sampled candidate has at least its origin
    genome represented before cross-genome conflict annotation.

    Parameters
    ----------
    genome_configs : iterable of GenomeConfig
        Genome records to sample.
    k_values : list[int]
        K-mer sizes to sample.
    sqlite_path : str or pathlib.Path
        SQLite database path.
    batch_size : int, optional
        SQLite flush interval for sampled candidate rows.
    genome_bin_size : int, optional
        Reference bases per candidate-sampling bin.
    max_per_genome_bin : int, optional
        Maximum sampled candidates per genome/contig/k/bin.
    progress_interval : int, optional
        Attempted k-mer interval for logging.
    logger : logging.Logger or None, optional
        Logger for progress messages.

    Returns
    -------
    list[dict[str, object]]
        Per-genome candidate-sampling summaries.
    """
    if batch_size <= 0:
        raise ValueError("batch_size must be positive")
    if genome_bin_size <= 0:
        raise ValueError("genome_bin_size must be positive")
    if max_per_genome_bin <= 0:
        raise ValueError("max_per_genome_bin must be positive")
    if progress_interval <= 0:
        raise ValueError("progress_interval must be positive")

    configs = list(genome_configs)
    initialise_global_candidate_database(sqlite_path=sqlite_path)
    summaries: list[dict[str, object]] = []
    connection = _connect(sqlite_path=sqlite_path)
    try:
        configure_fast_global_sqlite(connection=connection)
        for genome_index, config in enumerate(configs, start=1):
            if logger:
                logger.info(
                    "Sampling genome-spread candidate k-mers from genome %d/%d: %s (%s; role=%s)",
                    genome_index,
                    len(configs),
                    config.genome_id,
                    config.species_name,
                    config.role,
                )
            attempted = 0
            selected = 0
            last_logged_at = 0
            candidate_buffer: list[tuple[object, ...]] = []
            source_buffer: list[tuple[object, ...]] = []
            bin_counts: dict[tuple[int, str, int], int] = defaultdict(int)
            for record in read_fasta_records(fasta_path=config.genome_fasta):
                for k in k_values:
                    for position, kmer in iter_kmers(sequence=record.sequence, k=k):
                        attempted += 1
                        bin_id = int(position) // genome_bin_size
                        bin_key = (int(k), record.identifier, bin_id)
                        if bin_counts[bin_key] < max_per_genome_bin:
                            bin_counts[bin_key] += 1
                            selected += 1
                            candidate_buffer.append(
                                _candidate_row_tuple(
                                    k=k,
                                    kmer=kmer,
                                    genome_config=config,
                                    contig_id=record.identifier,
                                    position=position,
                                )
                            )
                            source_buffer.append(
                                _source_row_tuple(
                                    k=k,
                                    kmer=kmer,
                                    genome_config=config,
                                    contig_id=record.identifier,
                                    position=position,
                                )
                            )
                        if len(candidate_buffer) + len(source_buffer) >= batch_size:
                            _flush_candidate_rows(
                                connection=connection,
                                candidate_rows=candidate_buffer,
                                source_rows=source_buffer,
                            )
                            connection.commit()
                        if logger and attempted - last_logged_at >= progress_interval:
                            logger.info(
                                "Sampled %d candidate k-mer(s) from %d attempted observations in %s",
                                selected,
                                attempted,
                                config.genome_id,
                            )
                            last_logged_at = attempted
            _flush_candidate_rows(
                connection=connection,
                candidate_rows=candidate_buffer,
                source_rows=source_buffer,
            )
            connection.commit()
            _record_event(
                connection=connection,
                stage="sample_candidate_universe",
                detail=f"{config.genome_id}:{config.species_name}",
                n_records=selected,
            )
            connection.commit()
            summaries.append(
                {
                    "stage": "sample_candidate_universe",
                    "genome_id": config.genome_id,
                    "species_name": config.species_name,
                    "role": config.role,
                    "taxid": config.taxid,
                    "k_values": ";".join(str(value) for value in k_values),
                    "source_index_mode": "candidate_universe",
                    "attempted_kmers": attempted,
                    "sampled_candidate_kmers": selected,
                    "genome_bin_size": genome_bin_size,
                    "max_per_genome_bin": max_per_genome_bin,
                }
            )
    finally:
        connection.close()
    return summaries


def annotate_candidate_universe_sources_sqlite(
    *,
    genome_configs: Iterable[GenomeConfig],
    k_values: list[int],
    sqlite_path: str | Path,
    batch_size: int = 50000,
    progress_interval: int = 1000000,
    logger: logging.Logger | None = None,
) -> list[dict[str, object]]:
    """Annotate candidate-universe k-mers against all source genomes.

    Each genome is scanned, but only k-mers present in the bounded candidate
    universe are written to ``global_kmer_sources``. This preserves global
    near-neighbour/outgroup validation while avoiding SQLite writes for every
    non-candidate sliding-window observation.

    Parameters
    ----------
    genome_configs : iterable of GenomeConfig
        Genome records to scan for candidate hits.
    k_values : list[int]
        K-mer sizes to scan.
    sqlite_path : str or pathlib.Path
        SQLite database path.
    batch_size : int, optional
        SQLite flush interval for source-hit rows.
    progress_interval : int, optional
        Attempted k-mer interval for logging.
    logger : logging.Logger or None, optional
        Logger for progress messages.

    Returns
    -------
    list[dict[str, object]]
        Per-genome candidate-hit annotation summaries.
    """
    if batch_size <= 0:
        raise ValueError("batch_size must be positive")
    if progress_interval <= 0:
        raise ValueError("progress_interval must be positive")

    configs = list(genome_configs)
    connection = _connect(sqlite_path=sqlite_path)
    try:
        configure_fast_global_sqlite(connection=connection)
        candidates_by_k = _load_candidate_kmers_by_k(connection=connection)
        candidate_count = sum(len(values) for values in candidates_by_k.values())
        if logger:
            logger.info(
                "Loaded %d distinct candidate k-mer(s) across %d k value(s) for conflict annotation",
                candidate_count,
                len(candidates_by_k),
            )
    finally:
        connection.close()

    summaries: list[dict[str, object]] = []
    connection = _connect(sqlite_path=sqlite_path)
    try:
        configure_fast_global_sqlite(connection=connection)
        for genome_index, config in enumerate(configs, start=1):
            if logger:
                logger.info(
                    "Annotating candidate hits in genome %d/%d: %s (%s; role=%s)",
                    genome_index,
                    len(configs),
                    config.genome_id,
                    config.species_name,
                    config.role,
                )
            attempted = 0
            matched = 0
            last_logged_at = 0
            source_buffer: list[tuple[object, ...]] = []
            for record in read_fasta_records(fasta_path=config.genome_fasta):
                for k in k_values:
                    candidate_set = candidates_by_k.get(int(k))
                    if not candidate_set:
                        continue
                    for position, kmer in iter_kmers(sequence=record.sequence, k=k):
                        attempted += 1
                        if kmer not in candidate_set:
                            if logger and attempted - last_logged_at >= progress_interval:
                                logger.info(
                                    "Scanned %d observations and found %d candidate hit(s) in %s",
                                    attempted,
                                    matched,
                                    config.genome_id,
                                )
                                last_logged_at = attempted
                            continue
                        matched += 1
                        source_buffer.append(
                            _source_row_tuple(
                                k=k,
                                kmer=kmer,
                                genome_config=config,
                                contig_id=record.identifier,
                                position=position,
                            )
                        )
                        if len(source_buffer) >= batch_size:
                            connection.executemany(
                                GLOBAL_KMER_SOURCE_INSERT_SQL,
                                source_buffer,
                            )
                            source_buffer.clear()
                            connection.commit()
                        if logger and attempted - last_logged_at >= progress_interval:
                            logger.info(
                                "Scanned %d observations and found %d candidate hit(s) in %s",
                                attempted,
                                matched,
                                config.genome_id,
                            )
                            last_logged_at = attempted
            if source_buffer:
                connection.executemany(GLOBAL_KMER_SOURCE_INSERT_SQL, source_buffer)
                source_buffer.clear()
            connection.commit()
            _record_event(
                connection=connection,
                stage="annotate_candidate_universe",
                detail=f"{config.genome_id}:{config.species_name}",
                n_records=matched,
            )
            connection.commit()
            summaries.append(
                {
                    "stage": "annotate_candidate_universe",
                    "genome_id": config.genome_id,
                    "species_name": config.species_name,
                    "role": config.role,
                    "taxid": config.taxid,
                    "k_values": ";".join(str(value) for value in k_values),
                    "source_index_mode": "candidate_universe",
                    "attempted_kmers": attempted,
                    "candidate_hits": matched,
                }
            )
    finally:
        connection.close()
    return summaries

def collect_global_kmer_sources_sqlite(
    *,
    genome_configs: Iterable[GenomeConfig],
    k_values: list[int],
    sqlite_path: str | Path,
    batch_size: int = 50000,
    source_index_mode: str = "source_rows",
    progress_interval: int = 1000000,
    genome_bin_size: int = 10000,
    max_per_genome_bin: int = 10,
    logger: logging.Logger | None = None,
) -> list[dict[str, object]]:
    """Collect source metadata for all genomes in one global SQLite index.

    Parameters
    ----------
    genome_configs : iterable of GenomeConfig
        Genome records to index.
    k_values : list[int]
        K-mer sizes to collect.
    sqlite_path : str or pathlib.Path
        SQLite database path.
    batch_size : int, optional
        Commit interval in attempted k-mer observations.
    logger : logging.Logger or None, optional
        Logger for progress messages.

    Returns
    -------
    list[dict[str, object]]
        Per-genome collection summary records.
    """
    if batch_size <= 0:
        raise ValueError("batch_size must be positive")
    if progress_interval <= 0:
        raise ValueError("progress_interval must be positive")
    if genome_bin_size <= 0:
        raise ValueError("genome_bin_size must be positive")
    if max_per_genome_bin <= 0:
        raise ValueError("max_per_genome_bin must be positive")
    if source_index_mode not in VALID_GLOBAL_SOURCE_INDEX_MODES:
        raise ValueError(
            "source_index_mode must be one of: "
            + ", ".join(sorted(VALID_GLOBAL_SOURCE_INDEX_MODES))
        )

    configs = list(genome_configs)
    if source_index_mode == "candidate_universe":
        summaries = collect_candidate_universe_sqlite(
            genome_configs=configs,
            k_values=k_values,
            sqlite_path=sqlite_path,
            batch_size=batch_size,
            genome_bin_size=genome_bin_size,
            max_per_genome_bin=max_per_genome_bin,
            progress_interval=progress_interval,
            logger=logger,
        )
        summaries.extend(
            annotate_candidate_universe_sources_sqlite(
                genome_configs=configs,
                k_values=k_values,
                sqlite_path=sqlite_path,
                batch_size=batch_size,
                progress_interval=progress_interval,
                logger=logger,
            )
        )
        materialise_summary = materialise_global_kmers_from_sources(
            sqlite_path=sqlite_path,
            logger=logger,
        )
        summaries.append(
            {
                "stage": "materialise_global_sources",
                "genome_id": "",
                "species_name": "",
                "role": "",
                "taxid": "",
                "k_values": ";".join(str(value) for value in k_values),
                "source_index_mode": source_index_mode,
                "attempted_kmers": materialise_summary["source_rows"],
                "distinct_global_kmers": materialise_summary["global_kmers"],
            }
        )
        return summaries

    summaries: list[dict[str, object]] = []
    initialise_global_candidate_database(sqlite_path=sqlite_path)
    connection = _connect(sqlite_path=sqlite_path)
    try:
        configure_fast_global_sqlite(connection=connection)
        for genome_index, config in enumerate(configs, start=1):
            if logger:
                logger.info(
                    "Indexing genome %d/%d once for global evidence: %s (%s; role=%s)",
                    genome_index,
                    len(configs),
                    config.genome_id,
                    config.species_name,
                    config.role,
                )
            attempted = 0
            committed_at = 0
            last_logged_at = 0
            buffer: list[tuple[object, ...]] = []
            for record in read_fasta_records(fasta_path=config.genome_fasta):
                for k in k_values:
                    for position, kmer in iter_kmers(sequence=record.sequence, k=k):
                        if source_index_mode == "source_rows":
                            buffer.append(
                                _source_row_tuple(
                                    k=k,
                                    kmer=kmer,
                                    genome_config=config,
                                    contig_id=record.identifier,
                                    position=position,
                                )
                            )
                        else:
                            buffer.append(
                                _global_row_tuple(
                                    k=k,
                                    kmer=kmer,
                                    genome_config=config,
                                    contig_id=record.identifier,
                                    position=position,
                                )
                            )
                        attempted += 1
                        if len(buffer) >= batch_size:
                            _flush_global_rows(
                                connection=connection,
                                rows=buffer,
                                source_index_mode=source_index_mode,
                            )
                            connection.commit()
                            committed_at = attempted
                        if logger and attempted - last_logged_at >= progress_interval:
                            logger.info(
                                "Indexed %d attempted k-mer observations from %s",
                                attempted,
                                config.genome_id,
                            )
                            last_logged_at = attempted
            _flush_global_rows(
                connection=connection,
                rows=buffer,
                source_index_mode=source_index_mode,
            )
            connection.commit()
            _record_event(
                connection=connection,
                stage="collect_global_sources",
                detail=f"{config.genome_id}:{config.species_name}",
                n_records=attempted,
            )
            connection.commit()
            summaries.append(
                {
                    "stage": "collect_global_sources",
                    "genome_id": config.genome_id,
                    "species_name": config.species_name,
                    "role": config.role,
                    "taxid": config.taxid,
                    "k_values": ";".join(str(value) for value in k_values),
                    "source_index_mode": source_index_mode,
                    "attempted_kmers": attempted,
                }
            )
    finally:
        connection.close()
    if source_index_mode == "source_rows":
        materialise_summary = materialise_global_kmers_from_sources(
            sqlite_path=sqlite_path,
            logger=logger,
        )
        summaries.append(
            {
                "stage": "materialise_global_sources",
                "genome_id": "",
                "species_name": "",
                "role": "",
                "taxid": "",
                "k_values": ";".join(str(value) for value in k_values),
                "source_index_mode": source_index_mode,
                "attempted_kmers": materialise_summary["source_rows"],
                "distinct_global_kmers": materialise_summary["global_kmers"],
            }
        )
    return summaries


def _has_reportable_source(
    *,
    roles: set[str],
    role_filter: set[str],
    excluded_roles: set[str],
) -> bool:
    """Return whether an evidence group has at least one reportable source.

    Parameters
    ----------
    roles : set[str]
        Roles observed for a k-mer.
    role_filter : set[str]
        Optional whitelist of candidate/reportable roles.
    excluded_roles : set[str]
        Roles that should not be reportable candidates.

    Returns
    -------
    bool
        True when the k-mer has a reportable source.
    """
    if role_filter:
        return bool(roles & role_filter)
    return bool(roles - excluded_roles)


def _insert_retained_diagnostic(
    *,
    connection: sqlite3.Connection,
    diagnostic: DiagnosticKmer,
) -> bool:
    """Insert a retained diagnostic k-mer.

    Parameters
    ----------
    connection : sqlite3.Connection
        Open SQLite connection.
    diagnostic : DiagnosticKmer
        Diagnostic k-mer record.

    Returns
    -------
    bool
        True if a new row was inserted.
    """
    cursor = connection.execute(
        """
        INSERT OR IGNORE INTO retained_kmers(
            k, kmer, panel_type, species_name, clade, source_genomes,
            source_contigs, example_position, evidence_taxid, evidence_name,
            evidence_rank, lineage_taxids, source_taxids
        )
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            int(diagnostic.k),
            diagnostic.kmer,
            diagnostic.panel_type,
            diagnostic.species_name,
            diagnostic.clade,
            diagnostic.source_genomes,
            diagnostic.source_contigs,
            int(diagnostic.example_position),
            diagnostic.evidence_taxid,
            diagnostic.evidence_name,
            diagnostic.evidence_rank,
            diagnostic.lineage_taxids,
            diagnostic.source_taxids,
        ),
    )
    return cursor.rowcount == 1


def _diagnostic_from_global_row(
    *,
    row: sqlite3.Row,
    evidence_taxid: str,
    evidence_name: str,
    evidence_rank: str,
    lineage_taxids: str,
) -> DiagnosticKmer:
    """Create a DiagnosticKmer from one global source-index row.

    Parameters
    ----------
    row : sqlite3.Row
        Row from ``global_kmers``.
    evidence_taxid : str
        Taxid represented by the diagnostic k-mer.
    evidence_name : str
        Taxonomic name represented by the diagnostic k-mer.
    evidence_rank : str
        Evidence rank.
    lineage_taxids : str
        Semicolon-separated lineage taxids for evidence_taxid.

    Returns
    -------
    DiagnosticKmer
        Diagnostic k-mer record.
    """
    species_names = _normalise_semicolon_set(row["species_names"])
    clades = _normalise_semicolon_set(row["clades"])
    species_name = (
        next(iter(species_names))
        if evidence_rank == "species" and len(species_names) == 1
        else ""
    )
    clade = next(iter(clades)) if len(clades) == 1 else evidence_name
    panel_type = "species_unique" if evidence_rank == "species" else f"{evidence_rank}_core"
    return DiagnosticKmer(
        kmer=row["kmer"],
        k=int(row["k"]),
        panel_type=panel_type,
        species_name=species_name,
        clade=clade,
        source_genomes=row["genome_ids"],
        source_contigs=row["contig_ids"],
        example_position=int(row["example_position"]),
        evidence_taxid=evidence_taxid,
        evidence_name=evidence_name,
        evidence_rank=evidence_rank,
        lineage_taxids=lineage_taxids,
        source_taxids=row["taxids"],
    )


def _single_semicolon_value(value: str | None) -> str:
    """Return a value only when a semicolon field has exactly one item.

    Parameters
    ----------
    value : str or None
        Semicolon-separated field.

    Returns
    -------
    str
        The only non-empty value, or an empty string if zero or multiple values
        are present.
    """
    if not value:
        return ""
    text = str(value)
    if ";" not in text:
        return text if text else ""
    first = ""
    for item in text.split(";"):
        if not item:
            continue
        if first:
            return ""
        first = item
    return first


def _diagnostic_tuple_from_global_row(
    *,
    row: sqlite3.Row,
    evidence_taxid: str,
    evidence_name: str,
    evidence_rank: str,
    lineage_taxids: str,
) -> tuple[object, ...]:
    """Create a retained-kmer SQLite tuple from one global source row.

    Parameters
    ----------
    row : sqlite3.Row
        Row from ``global_kmers``.
    evidence_taxid : str
        Taxid represented by the diagnostic k-mer.
    evidence_name : str
        Taxonomic name represented by the diagnostic k-mer.
    evidence_rank : str
        Evidence rank.
    lineage_taxids : str
        Semicolon-separated lineage taxids for evidence_taxid.

    Returns
    -------
    tuple[object, ...]
        SQLite parameter tuple for ``retained_kmers`` insertion.

    Notes
    -----
    This helper intentionally avoids constructing a ``DiagnosticKmer`` for the
    evidence-assignment hot path. Large builds may evaluate millions of rows,
    so avoiding repeated dataclass construction and repeated set parsing gives
    a meaningful speed-up without changing the retained evidence semantics.
    """
    species_name = ""
    if evidence_rank == "species":
        species_name = _single_semicolon_value(row["species_names"])
    clade = _single_semicolon_value(row["clades"]) or evidence_name
    panel_type = "species_unique" if evidence_rank == "species" else f"{evidence_rank}_core"
    return (
        int(row["k"]),
        row["kmer"],
        panel_type,
        species_name,
        clade,
        row["genome_ids"],
        row["contig_ids"],
        int(row["example_position"]),
        evidence_taxid,
        evidence_name,
        evidence_rank,
        lineage_taxids,
        row["taxids"],
    )


RETAINED_KMER_INSERT_SQL = """
    INSERT OR IGNORE INTO retained_kmers(
        k, kmer, panel_type, species_name, clade, source_genomes,
        source_contigs, example_position, evidence_taxid, evidence_name,
        evidence_rank, lineage_taxids, source_taxids
    )
    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
"""


def _flush_retained_rows(
    *,
    connection: sqlite3.Connection,
    rows: list[tuple[object, ...]],
) -> int:
    """Flush buffered retained-kmer rows to SQLite.

    Parameters
    ----------
    connection : sqlite3.Connection
        Open SQLite connection.
    rows : list[tuple[object, ...]]
        Buffered retained-kmer rows.

    Returns
    -------
    int
        Number of newly inserted rows.
    """
    if not rows:
        return 0
    before = connection.total_changes
    connection.executemany(RETAINED_KMER_INSERT_SQL, rows)
    inserted = int(connection.total_changes - before)
    rows.clear()
    return inserted


def _cached_taxid_tuple(
    *,
    taxid_text: str,
    taxonomy_db: TaxonomyDatabase,
    cache: dict[str, tuple[str, ...]],
) -> tuple[str, ...]:
    """Return a cached normalised taxid tuple from a semicolon field.

    Parameters
    ----------
    taxid_text : str
        Semicolon-separated source taxids.
    taxonomy_db : TaxonomyDatabase
        Taxonomy database used for taxid normalisation.
    cache : dict[str, tuple[str, ...]]
        Cache keyed by the raw semicolon field.

    Returns
    -------
    tuple[str, ...]
        Sorted non-empty normalised taxids.
    """
    cached = cache.get(taxid_text)
    if cached is not None:
        return cached
    values = {
        taxonomy_db.normalise_taxid(item)
        for item in str(taxid_text or "").split(";")
        if item
    }
    result = tuple(sorted(taxid for taxid in values if taxid))
    cache[taxid_text] = result
    return result


def _cached_reportable_roles(
    *,
    role_text: str,
    role_filter: set[str],
    excluded_roles: set[str],
    cache: dict[str, bool],
) -> bool:
    """Return whether a semicolon role field has reportable evidence.

    Parameters
    ----------
    role_text : str
        Semicolon-separated roles.
    role_filter : set[str]
        Optional whitelist of candidate/reportable roles.
    excluded_roles : set[str]
        Roles that should not be reportable candidates.
    cache : dict[str, bool]
        Cache keyed by the raw role field.

    Returns
    -------
    bool
        True when the role field has at least one reportable source.
    """
    cached = cache.get(role_text)
    if cached is not None:
        return cached
    roles = {item for item in str(role_text or "").split(";") if item}
    result = _has_reportable_source(
        roles=roles,
        role_filter=role_filter,
        excluded_roles=excluded_roles,
    )
    cache[role_text] = result
    return result


def _cached_evidence_assignment(
    *,
    source_taxids: tuple[str, ...],
    taxonomy_db: TaxonomyDatabase,
    ranks: list[str],
    target_taxid: str,
    cache: dict[tuple[str, ...], tuple[str, str, str, str, str]],
) -> tuple[str, str, str, str, str]:
    """Return cached taxonomic evidence assignment for a source-taxid tuple.

    Parameters
    ----------
    source_taxids : tuple[str, ...]
        Sorted normalised source taxids.
    taxonomy_db : TaxonomyDatabase
        Taxonomy database.
    ranks : list[str]
        Preferred evidence ranks.
    target_taxid : str
        Optional subtree restriction.
    cache : dict[tuple[str, ...], tuple[str, str, str, str, str]]
        Cache keyed by source taxid tuple. The first returned field is a status
        label: ``ok``, ``unranked`` or ``outside_target``.

    Returns
    -------
    tuple[str, str, str, str, str]
        Status, evidence taxid, evidence name, evidence rank and semicolon
        lineage. Non-``ok`` statuses return empty evidence fields.
    """
    cached = cache.get(source_taxids)
    if cached is not None:
        return cached
    evidence_node = taxonomy_db.best_named_ancestor(
        taxids=source_taxids,
        preferred_ranks=ranks,
    )
    if evidence_node is None or evidence_node.rank not in ranks:
        result = ("unranked", "", "", "", "")
        cache[source_taxids] = result
        return result
    if target_taxid and not taxonomy_db.is_descendant(
        taxid=evidence_node.taxid,
        ancestor_taxid=target_taxid,
    ):
        result = ("outside_target", "", "", "", "")
        cache[source_taxids] = result
        return result
    result = (
        "ok",
        evidence_node.taxid,
        evidence_node.name,
        evidence_node.rank,
        ";".join(taxonomy_db.get_lineage(evidence_node.taxid)),
    )
    cache[source_taxids] = result
    return result

def assign_global_candidate_evidence_sqlite(
    *,
    sqlite_path: str | Path,
    taxonomy_db: TaxonomyDatabase,
    preferred_ranks: list[str] | None = None,
    target_taxid: str = "",
    candidate_roles: Iterable[str] | None = None,
    excluded_candidate_roles: Iterable[str] | None = None,
    batch_size: int = 50000,
    max_per_evidence_per_k: int | None = None,
    logger: logging.Logger | None = None,
) -> list[dict[str, object]]:
    """Assign taxonomic evidence from a global source index.

    Parameters
    ----------
    sqlite_path : str or pathlib.Path
        SQLite database containing ``global_kmers``.
    taxonomy_db : TaxonomyDatabase
        Parsed taxonomy database.
    preferred_ranks : list[str] or None, optional
        Evidence ranks to retain.
    target_taxid : str, optional
        Optional taxid subtree restriction.
    candidate_roles : iterable of str or None, optional
        Optional whitelist of reportable source roles.
    excluded_candidate_roles : iterable of str or None, optional
        Roles excluded from reportable evidence when no whitelist is supplied.
    batch_size : int, optional
        Commit interval for retained evidence.
    max_per_evidence_per_k : int or None, optional
        Optional cap per evidence taxon/rank/k bucket.
    logger : logging.Logger or None, optional
        Logger for progress messages.

    Returns
    -------
    list[dict[str, object]]
        Evidence assignment summary records.
    """
    if batch_size <= 0:
        raise ValueError("batch_size must be positive")
    if max_per_evidence_per_k is not None and max_per_evidence_per_k <= 0:
        raise ValueError("max_per_evidence_per_k must be positive")

    ranks = preferred_ranks or CORE_RANK_ORDER
    target_taxid = taxonomy_db.normalise_taxid(target_taxid)
    role_filter = _normalise_roles(candidate_roles)
    excluded_roles = _normalise_roles(excluded_candidate_roles) or DEFAULT_EXCLUDED_REPORTABLE_ROLES

    retained_counts: dict[tuple[str, str, str, str, int], int] = defaultdict(int)
    considered = 0
    retained = 0
    skipped_no_taxid = 0
    skipped_unranked = 0
    skipped_outside_target = 0
    skipped_non_reportable = 0
    skipped_by_limit = 0
    duplicates = 0

    taxid_tuple_cache: dict[str, tuple[str, ...]] = {}
    role_reportable_cache: dict[str, bool] = {}
    evidence_cache: dict[tuple[str, ...], tuple[str, str, str, str, str]] = {}

    connection = _connect(sqlite_path=sqlite_path)
    connection.row_factory = sqlite3.Row
    retained_buffer: list[tuple[object, ...]] = []
    try:
        configure_fast_global_sqlite(connection=connection)
        cursor = connection.execute(
            """
            SELECT
                k, kmer, species_names, genome_ids, contig_ids, taxids, clades,
                roles, example_position
            FROM global_kmers
            ORDER BY k, kmer
            """
        )
        while True:
            rows = cursor.fetchmany(batch_size)
            if not rows:
                break
            for row in rows:
                considered += 1
                source_taxids = _cached_taxid_tuple(
                    taxid_text=row["taxids"],
                    taxonomy_db=taxonomy_db,
                    cache=taxid_tuple_cache,
                )
                if not source_taxids:
                    skipped_no_taxid += 1
                    continue
                if not _cached_reportable_roles(
                    role_text=row["roles"],
                    role_filter=role_filter,
                    excluded_roles=excluded_roles,
                    cache=role_reportable_cache,
                ):
                    skipped_non_reportable += 1
                    continue
                (
                    evidence_status,
                    evidence_taxid,
                    evidence_name,
                    evidence_rank,
                    lineage_taxids,
                ) = _cached_evidence_assignment(
                    source_taxids=source_taxids,
                    taxonomy_db=taxonomy_db,
                    ranks=ranks,
                    target_taxid=target_taxid,
                    cache=evidence_cache,
                )
                if evidence_status == "unranked":
                    skipped_unranked += 1
                    continue
                if evidence_status == "outside_target":
                    skipped_outside_target += 1
                    continue

                diagnostic_tuple = _diagnostic_tuple_from_global_row(
                    row=row,
                    evidence_taxid=evidence_taxid,
                    evidence_name=evidence_name,
                    evidence_rank=evidence_rank,
                    lineage_taxids=lineage_taxids,
                )
                retention_key = (
                    diagnostic_tuple[2],
                    diagnostic_tuple[8],
                    diagnostic_tuple[3],
                    diagnostic_tuple[4],
                    diagnostic_tuple[0],
                )
                if max_per_evidence_per_k is not None:
                    if retained_counts[retention_key] >= max_per_evidence_per_k:
                        skipped_by_limit += 1
                        continue
                retained_counts[retention_key] += 1
                retained_buffer.append(diagnostic_tuple)
                if len(retained_buffer) >= batch_size:
                    n_buffered = len(retained_buffer)
                    inserted = _flush_retained_rows(
                        connection=connection,
                        rows=retained_buffer,
                    )
                    retained += inserted
                    duplicates += n_buffered - inserted
            if retained_buffer:
                n_buffered = len(retained_buffer)
                inserted = _flush_retained_rows(
                    connection=connection,
                    rows=retained_buffer,
                )
                retained += inserted
                duplicates += n_buffered - inserted
            connection.commit()
            if logger and (considered <= batch_size or considered % (batch_size * 20) == 0):
                logger.info(
                    "Assigned evidence for %d global k-mer keys; retained=%d; "
                    "skipped_by_limit=%d; taxid_cache=%d; evidence_cache=%d",
                    considered,
                    retained,
                    skipped_by_limit,
                    len(taxid_tuple_cache),
                    len(evidence_cache),
                )
        if retained_buffer:
            n_buffered = len(retained_buffer)
            inserted = _flush_retained_rows(
                connection=connection,
                rows=retained_buffer,
            )
            retained += inserted
            duplicates += n_buffered - inserted
            connection.commit()
        _record_event(
            connection=connection,
            stage="assign_global_evidence",
            detail="global all-candidate evidence assignment",
            n_records=retained,
        )
        connection.commit()
    finally:
        connection.close()

    return [
        {"stage": "assign_global_evidence", "metric": "global_kmers_considered", "value": considered},
        {"stage": "assign_global_evidence", "metric": "diagnostics_retained", "value": retained},
        {"stage": "assign_global_evidence", "metric": "duplicates_ignored", "value": duplicates},
        {"stage": "assign_global_evidence", "metric": "skipped_no_taxid", "value": skipped_no_taxid},
        {"stage": "assign_global_evidence", "metric": "skipped_unranked", "value": skipped_unranked},
        {"stage": "assign_global_evidence", "metric": "skipped_outside_target", "value": skipped_outside_target},
        {"stage": "assign_global_evidence", "metric": "skipped_non_reportable", "value": skipped_non_reportable},
        {"stage": "assign_global_evidence", "metric": "skipped_by_limit", "value": skipped_by_limit},
    ]


def iter_retained_global_candidate_diagnostics(
    *,
    sqlite_path: str | Path,
    chunk_size: int = 10000,
) -> Iterator[DiagnosticKmer]:
    """Yield retained diagnostics from a global all-candidate database.

    Parameters
    ----------
    sqlite_path : str or pathlib.Path
        SQLite database path.
    chunk_size : int, optional
        Number of rows fetched per SQLite batch.

    Yields
    ------
    DiagnosticKmer
        Retained diagnostic k-mer record.
    """
    if chunk_size <= 0:
        raise ValueError("chunk_size must be positive")
    connection = _connect(sqlite_path=sqlite_path)
    connection.row_factory = sqlite3.Row
    try:
        cursor = connection.execute(
            """
            SELECT
                k, kmer, panel_type, species_name, clade, source_genomes,
                source_contigs, example_position, evidence_taxid,
                evidence_name, evidence_rank, lineage_taxids, source_taxids
            FROM retained_kmers
            ORDER BY k, evidence_rank, evidence_name, kmer
            """
        )
        while True:
            rows = cursor.fetchmany(chunk_size)
            if not rows:
                break
            for row in rows:
                yield DiagnosticKmer(
                    kmer=row["kmer"],
                    k=int(row["k"]),
                    panel_type=row["panel_type"],
                    species_name=row["species_name"],
                    clade=row["clade"],
                    source_genomes=row["source_genomes"],
                    source_contigs=row["source_contigs"],
                    example_position=int(row["example_position"]),
                    evidence_taxid=row["evidence_taxid"],
                    evidence_name=row["evidence_name"],
                    evidence_rank=row["evidence_rank"],
                    lineage_taxids=row["lineage_taxids"],
                    source_taxids=row["source_taxids"],
                )
    finally:
        connection.close()


def summarise_retained_global_candidate_evidence(
    *,
    sqlite_path: str | Path,
) -> list[dict[str, object]]:
    """Summarise retained global all-candidate diagnostics.

    Parameters
    ----------
    sqlite_path : str or pathlib.Path
        SQLite database path.

    Returns
    -------
    list[dict[str, object]]
        Summary rows grouped by evidence bucket.
    """
    connection = _connect(sqlite_path=sqlite_path)
    try:
        rows = connection.execute(
            """
            SELECT
                panel_type,
                species_name,
                clade,
                evidence_taxid,
                evidence_rank,
                k,
                COUNT(*) AS diagnostic_kmers
            FROM retained_kmers
            GROUP BY
                panel_type, species_name, clade, evidence_taxid,
                evidence_rank, k
            ORDER BY evidence_rank, clade, species_name, k
            """
        ).fetchall()
    finally:
        connection.close()
    return [
        {
            "panel_type": row[0],
            "species_name": row[1],
            "clade": row[2],
            "evidence_taxid": row[3],
            "evidence_rank": row[4],
            "k": int(row[5]),
            "diagnostic_kmers": int(row[6]),
        }
        for row in rows
    ]


def build_global_candidate_evidence_sqlite(
    *,
    genome_configs: list[GenomeConfig],
    k_values: list[int],
    sqlite_path: str | Path,
    taxonomy_db: TaxonomyDatabase,
    target_taxid: str = "",
    preferred_ranks: list[str] | None = None,
    candidate_roles: Iterable[str] | None = None,
    excluded_candidate_roles: Iterable[str] | None = None,
    batch_size: int = 50000,
    max_per_evidence_per_k: int | None = None,
    source_index_mode: str = "source_rows",
    progress_interval: int = 1000000,
    genome_bin_size: int = 10000,
    max_per_genome_bin: int = 10,
    logger: logging.Logger | None = None,
) -> GlobalCandidateEvidenceBuildResult:
    """Build a global query-agnostic evidence database.

    Parameters
    ----------
    genome_configs : list of GenomeConfig
        Genome records to index once.
    k_values : list[int]
        K-mer sizes.
    sqlite_path : str or pathlib.Path
        SQLite path used for the global index and retained evidence.
    taxonomy_db : TaxonomyDatabase
        Parsed taxonomy database. Required for global evidence assignment.
    target_taxid : str, optional
        Optional subtree restriction. Leave blank for broad panels.
    preferred_ranks : list[str] or None, optional
        Evidence ranks to retain.
    candidate_roles : iterable of str or None, optional
        Optional whitelist of reportable roles.
    excluded_candidate_roles : iterable of str or None, optional
        Roles excluded from reportable candidate status when no whitelist is set.
    batch_size : int, optional
        SQLite commit and fetch interval.
    max_per_evidence_per_k : int or None, optional
        Optional cap per evidence bucket and k.
    source_index_mode : str, optional
        Source-index implementation. ``source_rows`` is the default faster mode.
        ``aggregated`` preserves the legacy direct-upsert implementation.
    progress_interval : int, optional
        Attempted k-mer interval for progress logging while indexing genomes.
    logger : logging.Logger or None, optional
        Logger.

    Returns
    -------
    GlobalCandidateEvidenceBuildResult
        Build metadata.
    """
    if taxonomy_db is None:
        raise ValueError("taxonomy_db is required for global candidate evidence")
    if batch_size <= 0:
        raise ValueError("batch_size must be positive")
    if max_per_evidence_per_k is not None and max_per_evidence_per_k <= 0:
        raise ValueError("max_per_evidence_per_k must be positive")
    if progress_interval <= 0:
        raise ValueError("progress_interval must be positive")
    if source_index_mode not in VALID_GLOBAL_SOURCE_INDEX_MODES:
        raise ValueError(
            "source_index_mode must be one of: "
            + ", ".join(sorted(VALID_GLOBAL_SOURCE_INDEX_MODES))
        )

    db_path = Path(sqlite_path)
    if db_path.exists():
        db_path.unlink()
    initialise_global_candidate_database(sqlite_path=db_path)

    if logger:
        logger.info(
            "Starting global all-candidate evidence build with %d genome(s), %d k value(s)",
            len(genome_configs),
            len(k_values),
        )
        logger.info(
            "Each genome will be indexed once before taxonomic evidence assignment"
        )
        logger.info("Global source-index mode: %s", source_index_mode)

    collection_summary = collect_global_kmer_sources_sqlite(
        genome_configs=genome_configs,
        k_values=k_values,
        sqlite_path=db_path,
        batch_size=batch_size,
        source_index_mode=source_index_mode,
        progress_interval=progress_interval,
        genome_bin_size=genome_bin_size,
        max_per_genome_bin=max_per_genome_bin,
        logger=logger,
    )
    assignment_summary = assign_global_candidate_evidence_sqlite(
        sqlite_path=db_path,
        taxonomy_db=taxonomy_db,
        preferred_ranks=preferred_ranks,
        target_taxid=target_taxid,
        candidate_roles=candidate_roles,
        excluded_candidate_roles=excluded_candidate_roles,
        batch_size=batch_size,
        max_per_evidence_per_k=max_per_evidence_per_k,
        logger=logger,
    )
    panel_summary = summarise_retained_global_candidate_evidence(sqlite_path=db_path)

    connection = _connect(sqlite_path=db_path)
    try:
        global_kmer_count = connection.execute(
            "SELECT COUNT(*) FROM global_kmers"
        ).fetchone()[0]
        retained_count = connection.execute(
            "SELECT COUNT(*) FROM retained_kmers"
        ).fetchone()[0]
    finally:
        connection.close()

    build_summary = [
        {"summary_name": "genome_records", "summary_value": len(genome_configs)},
        {"summary_name": "k_values", "summary_value": ";".join(map(str, k_values))},
        {"summary_name": "global_source_index_mode", "summary_value": source_index_mode},
        {"summary_name": "global_distinct_kmer_keys", "summary_value": int(global_kmer_count)},
        {"summary_name": "diagnostics_retained", "summary_value": int(retained_count)},
        {"summary_name": "sqlite_path", "summary_value": str(db_path)},
    ]
    build_summary.extend(
        {
            "summary_name": f"assignment::{row['metric']}",
            "summary_value": row["value"],
        }
        for row in assignment_summary
    )
    collection_summary.extend(assignment_summary)
    return GlobalCandidateEvidenceBuildResult(
        sqlite_path=db_path,
        collection_summary=collection_summary,
        build_summary=build_summary,
        panel_summary=panel_summary,
    )
