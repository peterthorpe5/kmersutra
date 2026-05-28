"""Low-memory target-aware panel building for KmerSutra.

This module provides an SQLite-backed build path that is designed for larger
panels where the global in-memory compact dictionary is too expensive. It keeps
candidate k-mers from target genomes on disk and streams all non-target genomes
against those candidates. This is most appropriate when the immediate question
is whether named target species retain species-level evidence after adding
near-neighbours and outgroups.
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
from kmersutra.taxonomy import CORE_RANK_ORDER, TaxonomyDatabase


@dataclass(frozen=True)
class TargetEvidenceBuildResult:
    """Result metadata from an SQLite-backed target-evidence build.

    Attributes
    ----------
    sqlite_path : pathlib.Path
        Path to the SQLite database containing target candidate k-mers.
    collection_summary : list[dict[str, object]]
        Per-genome k-mer collection and comparison summaries.
    build_summary : list[dict[str, object]]
        Coarse build-stage summary records.
    """

    sqlite_path: Path
    collection_summary: list[dict[str, object]]
    build_summary: list[dict[str, object]]


def _normalise_semicolon_set(value: str) -> set[str]:
    """Convert a semicolon-separated string to a set.

    Parameters
    ----------
    value : str
        Semicolon-separated value.

    Returns
    -------
    set[str]
        Non-empty values.
    """
    if not value:
        return set()
    return {item for item in value.split(";") if item}


def _merge_semicolon_value(value: str, new_value: str) -> str:
    """Merge a value into a semicolon-separated set string.

    Parameters
    ----------
    value : str
        Existing semicolon-separated values.
    new_value : str
        New value to add.

    Returns
    -------
    str
        Sorted semicolon-separated values.
    """
    values = _normalise_semicolon_set(value)
    if new_value:
        values.add(str(new_value))
    return ";".join(sorted(values))


def _connect(sqlite_path: str | Path) -> sqlite3.Connection:
    """Open an SQLite connection with pragmatic build settings.

    Parameters
    ----------
    sqlite_path : str or pathlib.Path
        SQLite database path.

    Returns
    -------
    sqlite3.Connection
        Open connection.
    """
    path = Path(sqlite_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    connection = sqlite3.connect(str(path))
    connection.execute("PRAGMA journal_mode=WAL")
    connection.execute("PRAGMA synchronous=NORMAL")
    connection.execute("PRAGMA temp_store=FILE")
    connection.execute("PRAGMA cache_size=-200000")
    return connection


def initialise_target_evidence_database(*, sqlite_path: str | Path) -> None:
    """Create the SQLite schema for target-aware evidence building.

    Parameters
    ----------
    sqlite_path : str or pathlib.Path
        SQLite database path.
    """
    connection = _connect(sqlite_path=sqlite_path)
    try:
        connection.execute(
            """
            CREATE TABLE IF NOT EXISTS target_kmers (
                k INTEGER NOT NULL,
                kmer TEXT NOT NULL,
                target_species TEXT NOT NULL DEFAULT '',
                target_genomes TEXT NOT NULL DEFAULT '',
                target_contigs TEXT NOT NULL DEFAULT '',
                target_taxids TEXT NOT NULL DEFAULT '',
                target_clades TEXT NOT NULL DEFAULT '',
                example_position INTEGER NOT NULL DEFAULT 0,
                non_target_species TEXT NOT NULL DEFAULT '',
                non_target_genomes TEXT NOT NULL DEFAULT '',
                non_target_taxids TEXT NOT NULL DEFAULT '',
                non_target_clades TEXT NOT NULL DEFAULT '',
                non_target_roles TEXT NOT NULL DEFAULT '',
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
    """Record a build event in the SQLite database.

    Parameters
    ----------
    connection : sqlite3.Connection
        Open SQLite connection.
    stage : str
        Build stage name.
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


def _upsert_target_row(
    *,
    connection: sqlite3.Connection,
    k: int,
    kmer: str,
    genome_config: GenomeConfig,
    contig_id: str,
    position: int,
) -> None:
    """Insert or update one target candidate k-mer row.

    Parameters
    ----------
    connection : sqlite3.Connection
        Open SQLite connection.
    k : int
        K-mer length.
    kmer : str
        Canonical k-mer.
    genome_config : GenomeConfig
        Source genome metadata.
    contig_id : str
        Source contig identifier.
    position : int
        Example k-mer position.
    """
    connection.execute(
        """
        INSERT INTO target_kmers(
            k, kmer, target_species, target_genomes, target_contigs,
            target_taxids, target_clades, example_position
        )
        VALUES (?, ?, ?, ?, ?, ?, ?, ?)
        ON CONFLICT(k, kmer) DO UPDATE SET
            target_species = CASE
                WHEN instr(';' || target_species || ';', ';' || excluded.target_species || ';') > 0
                THEN target_species
                WHEN target_species = '' THEN excluded.target_species
                ELSE target_species || ';' || excluded.target_species
            END,
            target_genomes = CASE
                WHEN instr(';' || target_genomes || ';', ';' || excluded.target_genomes || ';') > 0
                THEN target_genomes
                WHEN target_genomes = '' THEN excluded.target_genomes
                ELSE target_genomes || ';' || excluded.target_genomes
            END,
            target_contigs = CASE
                WHEN instr(';' || target_contigs || ';', ';' || excluded.target_contigs || ';') > 0
                THEN target_contigs
                WHEN target_contigs = '' THEN excluded.target_contigs
                ELSE target_contigs || ';' || excluded.target_contigs
            END,
            target_taxids = CASE
                WHEN excluded.target_taxids = '' THEN target_taxids
                WHEN instr(';' || target_taxids || ';', ';' || excluded.target_taxids || ';') > 0
                THEN target_taxids
                WHEN target_taxids = '' THEN excluded.target_taxids
                ELSE target_taxids || ';' || excluded.target_taxids
            END,
            target_clades = CASE
                WHEN instr(';' || target_clades || ';', ';' || excluded.target_clades || ';') > 0
                THEN target_clades
                WHEN target_clades = '' THEN excluded.target_clades
                ELSE target_clades || ';' || excluded.target_clades
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
            int(position),
        ),
    )


def collect_target_candidates_sqlite(
    *,
    genome_configs: Iterable[GenomeConfig],
    k_values: Iterable[int],
    sqlite_path: str | Path,
    batch_size: int = 50000,
    logger: logging.Logger | None = None,
) -> list[dict[str, object]]:
    """Collect candidate k-mers from target genomes into SQLite.

    Parameters
    ----------
    genome_configs : iterable of GenomeConfig
        Target genome configuration records.
    k_values : iterable of int
        K-mer lengths to collect.
    sqlite_path : str or pathlib.Path
        SQLite database path.
    batch_size : int, optional
        Commit interval in processed k-mers.
    logger : logging.Logger | None, optional
        Optional logger.

    Returns
    -------
    list[dict[str, object]]
        Per-genome collection summary records.
    """
    if batch_size <= 0:
        raise ValueError("batch_size must be positive")
    k_tuple = tuple(k_values)
    summaries: list[dict[str, object]] = []
    connection = _connect(sqlite_path=sqlite_path)
    try:
        for genome_index, genome_config in enumerate(genome_configs, start=1):
            if logger:
                logger.info(
                    "Collecting target candidates from genome %d: %s (%s)",
                    genome_index,
                    genome_config.genome_id,
                    genome_config.species_name,
                )
            contigs = 0
            total_bases = 0
            total_observations = 0
            per_k_counts = {k: 0 for k in k_tuple}
            since_commit = 0
            for fasta_record in read_fasta_records(fasta_path=genome_config.genome_fasta):
                contigs += 1
                total_bases += len(fasta_record.sequence)
                for k in k_tuple:
                    for position, kmer in iter_kmers(sequence=fasta_record.sequence, k=k):
                        _upsert_target_row(
                            connection=connection,
                            k=k,
                            kmer=kmer,
                            genome_config=genome_config,
                            contig_id=fasta_record.identifier,
                            position=position,
                        )
                        total_observations += 1
                        per_k_counts[k] += 1
                        since_commit += 1
                        if since_commit >= batch_size:
                            connection.commit()
                            since_commit = 0
            connection.commit()
            summary = {
                "genome_id": genome_config.genome_id,
                "species_name": genome_config.species_name,
                "taxid": genome_config.taxid,
                "role": genome_config.role,
                "clade": genome_config.clade,
                "genome_fasta": str(genome_config.genome_fasta),
                "collection_mode": "target_candidate_insert",
                "contigs": contigs,
                "total_bases": total_bases,
                "total_observations": total_observations,
                **{f"observations_k{k}": per_k_counts[k] for k in k_tuple},
            }
            summaries.append(summary)
            _record_event(
                connection=connection,
                stage="collect_target_candidates",
                detail=genome_config.genome_id,
                n_records=total_observations,
            )
            connection.commit()
    finally:
        connection.close()
    return summaries


def _update_non_target_row(
    *,
    connection: sqlite3.Connection,
    k: int,
    kmer: str,
    genome_config: GenomeConfig,
) -> None:
    """Mark a target candidate as also observed in a non-target genome.

    Parameters
    ----------
    connection : sqlite3.Connection
        Open SQLite connection.
    k : int
        K-mer length.
    kmer : str
        Canonical k-mer.
    genome_config : GenomeConfig
        Non-target source genome metadata.
    """
    connection.execute(
        """
        UPDATE target_kmers SET
            non_target_species = CASE
                WHEN instr(';' || non_target_species || ';', ';' || ? || ';') > 0
                THEN non_target_species
                WHEN non_target_species = '' THEN ?
                ELSE non_target_species || ';' || ?
            END,
            non_target_genomes = CASE
                WHEN instr(';' || non_target_genomes || ';', ';' || ? || ';') > 0
                THEN non_target_genomes
                WHEN non_target_genomes = '' THEN ?
                ELSE non_target_genomes || ';' || ?
            END,
            non_target_taxids = CASE
                WHEN ? = '' THEN non_target_taxids
                WHEN instr(';' || non_target_taxids || ';', ';' || ? || ';') > 0
                THEN non_target_taxids
                WHEN non_target_taxids = '' THEN ?
                ELSE non_target_taxids || ';' || ?
            END,
            non_target_clades = CASE
                WHEN instr(';' || non_target_clades || ';', ';' || ? || ';') > 0
                THEN non_target_clades
                WHEN non_target_clades = '' THEN ?
                ELSE non_target_clades || ';' || ?
            END,
            non_target_roles = CASE
                WHEN instr(';' || non_target_roles || ';', ';' || ? || ';') > 0
                THEN non_target_roles
                WHEN non_target_roles = '' THEN ?
                ELSE non_target_roles || ';' || ?
            END
        WHERE k = ? AND kmer = ?
        """,
        (
            genome_config.species_name,
            genome_config.species_name,
            genome_config.species_name,
            genome_config.genome_id,
            genome_config.genome_id,
            genome_config.genome_id,
            genome_config.taxid,
            genome_config.taxid,
            genome_config.taxid,
            genome_config.taxid,
            genome_config.clade,
            genome_config.clade,
            genome_config.clade,
            genome_config.role,
            genome_config.role,
            genome_config.role,
            int(k),
            kmer,
        ),
    )


def scan_non_target_matches_sqlite(
    *,
    genome_configs: Iterable[GenomeConfig],
    k_values: Iterable[int],
    sqlite_path: str | Path,
    batch_size: int = 50000,
    logger: logging.Logger | None = None,
) -> list[dict[str, object]]:
    """Stream non-target genomes and mark overlaps with target candidates.

    Parameters
    ----------
    genome_configs : iterable of GenomeConfig
        Non-target genome records.
    k_values : iterable of int
        K-mer lengths to scan.
    sqlite_path : str or pathlib.Path
        SQLite database path containing target candidates.
    batch_size : int, optional
        Commit interval in processed k-mers.
    logger : logging.Logger | None, optional
        Optional logger.

    Returns
    -------
    list[dict[str, object]]
        Per-genome scan summary records.
    """
    if batch_size <= 0:
        raise ValueError("batch_size must be positive")
    k_tuple = tuple(k_values)
    summaries: list[dict[str, object]] = []
    connection = _connect(sqlite_path=sqlite_path)
    try:
        for genome_index, genome_config in enumerate(genome_configs, start=1):
            if logger:
                logger.info(
                    "Scanning non-target genome %d: %s (%s; role=%s)",
                    genome_index,
                    genome_config.genome_id,
                    genome_config.species_name,
                    genome_config.role,
                )
            contigs = 0
            total_bases = 0
            total_observations = 0
            candidate_hits = 0
            per_k_counts = {k: 0 for k in k_tuple}
            since_commit = 0
            for fasta_record in read_fasta_records(fasta_path=genome_config.genome_fasta):
                contigs += 1
                total_bases += len(fasta_record.sequence)
                for k in k_tuple:
                    for _, kmer in iter_kmers(sequence=fasta_record.sequence, k=k):
                        before = connection.total_changes
                        _update_non_target_row(
                            connection=connection,
                            k=k,
                            kmer=kmer,
                            genome_config=genome_config,
                        )
                        if connection.total_changes > before:
                            candidate_hits += 1
                        total_observations += 1
                        per_k_counts[k] += 1
                        since_commit += 1
                        if since_commit >= batch_size:
                            connection.commit()
                            since_commit = 0
            connection.commit()
            summary = {
                "genome_id": genome_config.genome_id,
                "species_name": genome_config.species_name,
                "taxid": genome_config.taxid,
                "role": genome_config.role,
                "clade": genome_config.clade,
                "genome_fasta": str(genome_config.genome_fasta),
                "collection_mode": "non_target_overlap_scan",
                "contigs": contigs,
                "total_bases": total_bases,
                "total_observations": total_observations,
                "candidate_hits": candidate_hits,
                **{f"observations_k{k}": per_k_counts[k] for k in k_tuple},
            }
            summaries.append(summary)
            _record_event(
                connection=connection,
                stage="scan_non_target_matches",
                detail=genome_config.genome_id,
                n_records=total_observations,
            )
            connection.commit()
    finally:
        connection.close()
    return summaries


def _row_to_diagnostic(
    *,
    row: sqlite3.Row,
    taxonomy_db: TaxonomyDatabase | None,
    target_taxid: str = "",
    preferred_ranks: list[str] | None = None,
) -> DiagnosticKmer | None:
    """Convert an SQLite row into a diagnostic k-mer if retained.

    Parameters
    ----------
    row : sqlite3.Row
        Row from ``target_kmers``.
    taxonomy_db : TaxonomyDatabase | None
        Optional NCBI taxonomy database.
    target_taxid : str, optional
        Taxid root to retain.
    preferred_ranks : list[str] | None, optional
        Retained evidence ranks.

    Returns
    -------
    DiagnosticKmer | None
        Diagnostic record, or None if the k-mer is not retained.
    """
    target_species = _normalise_semicolon_set(row["target_species"])
    target_taxids = _normalise_semicolon_set(row["target_taxids"])
    non_target_taxids = _normalise_semicolon_set(row["non_target_taxids"])
    all_taxids = target_taxids | non_target_taxids
    clades = _normalise_semicolon_set(row["target_clades"]) | _normalise_semicolon_set(
        row["non_target_clades"]
    )

    if taxonomy_db is not None:
        ranks = preferred_ranks or CORE_RANK_ORDER
        normalised_taxids = {taxonomy_db.normalise_taxid(taxid) for taxid in all_taxids}
        normalised_taxids = {taxid for taxid in normalised_taxids if taxid}
        if not normalised_taxids:
            return None
        evidence_node = taxonomy_db.best_named_ancestor(
            taxids=normalised_taxids,
            preferred_ranks=ranks,
        )
        if evidence_node is None:
            return None
        normalised_target_taxid = taxonomy_db.normalise_taxid(target_taxid)
        if normalised_target_taxid and not taxonomy_db.is_descendant(
            taxid=evidence_node.taxid,
            ancestor_taxid=normalised_target_taxid,
        ):
            return None
        if evidence_node.rank not in ranks:
            return None
        evidence_rank = evidence_node.rank
        evidence_name = evidence_node.name
        evidence_taxid = evidence_node.taxid
        lineage_taxids = ";".join(taxonomy_db.get_lineage(evidence_node.taxid))
    else:
        if row["non_target_taxids"]:
            return None
        evidence_rank = "species" if len(target_species) == 1 else "clade"
        evidence_name = next(iter(sorted(target_species))) if len(target_species) == 1 else "shared_target"
        evidence_taxid = next(iter(sorted(target_taxids))) if len(target_taxids) == 1 else ""
        lineage_taxids = ""

    species_name = ""
    if evidence_rank == "species" and len(target_species) == 1:
        species_name = next(iter(sorted(target_species)))
    panel_type = "species_unique" if evidence_rank == "species" else f"{evidence_rank}_core"

    return DiagnosticKmer(
        kmer=row["kmer"],
        k=int(row["k"]),
        panel_type=panel_type,
        species_name=species_name,
        clade=next(iter(sorted(clades))) if len(clades) == 1 else evidence_name,
        source_genomes=";".join(
            sorted(
                _normalise_semicolon_set(row["target_genomes"])
                | _normalise_semicolon_set(row["non_target_genomes"])
            )
        ),
        source_contigs=row["target_contigs"],
        example_position=int(row["example_position"]),
        evidence_taxid=evidence_taxid,
        evidence_name=evidence_name,
        evidence_rank=evidence_rank,
        lineage_taxids=lineage_taxids,
        source_taxids=";".join(sorted(all_taxids)),
    )


def iter_target_evidence_diagnostics(
    *,
    sqlite_path: str | Path,
    taxonomy_db: TaxonomyDatabase | None = None,
    target_taxid: str = "",
    preferred_ranks: list[str] | None = None,
    chunk_size: int = 10000,
) -> Iterator[DiagnosticKmer]:
    """Yield retained diagnostic k-mers from the SQLite candidate store.

    Parameters
    ----------
    sqlite_path : str or pathlib.Path
        SQLite candidate database.
    taxonomy_db : TaxonomyDatabase | None, optional
        Optional NCBI taxonomy database.
    target_taxid : str, optional
        Taxid root to retain.
    preferred_ranks : list[str] | None, optional
        Retained evidence ranks.
    chunk_size : int, optional
        Number of SQLite rows fetched per batch.

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
                k, kmer, target_species, target_genomes, target_contigs,
                target_taxids, target_clades, example_position,
                non_target_species, non_target_genomes, non_target_taxids,
                non_target_clades, non_target_roles
            FROM target_kmers
            ORDER BY k, kmer
            """
        )
        while True:
            rows = cursor.fetchmany(chunk_size)
            if not rows:
                break
            for row in rows:
                diagnostic = _row_to_diagnostic(
                    row=row,
                    taxonomy_db=taxonomy_db,
                    target_taxid=target_taxid,
                    preferred_ranks=preferred_ranks,
                )
                if diagnostic is not None:
                    yield diagnostic
    finally:
        connection.close()


def summarise_diagnostics_stream(
    *,
    diagnostics: Iterable[DiagnosticKmer],
    max_per_evidence_per_k: int | None = None,
) -> tuple[list[DiagnosticKmer], list[dict[str, object]]]:
    """Collect and summarise diagnostics from an iterator.

    This helper is intended for tests and small builds. The CLI uses a more
    direct streaming writer to avoid storing all diagnostics in memory.

    Parameters
    ----------
    diagnostics : iterable of DiagnosticKmer
        Diagnostic k-mer records.
    max_per_evidence_per_k : int | None, optional
        Optional retention limit per evidence bucket.

    Returns
    -------
    tuple[list[DiagnosticKmer], list[dict[str, object]]]
        Retained diagnostics and summary rows.
    """
    counts: dict[tuple[str, str, str, str, int], int] = defaultdict(int)
    summary_counts: dict[tuple[str, str, str, str, str, int], int] = defaultdict(int)
    retained: list[DiagnosticKmer] = []
    for diagnostic in diagnostics:
        retention_key = (
            diagnostic.panel_type,
            diagnostic.evidence_taxid,
            diagnostic.species_name,
            diagnostic.clade,
            diagnostic.k,
        )
        if max_per_evidence_per_k is not None:
            if max_per_evidence_per_k <= 0:
                raise ValueError("max_per_evidence_per_k must be positive")
            if counts[retention_key] >= max_per_evidence_per_k:
                continue
        counts[retention_key] += 1
        retained.append(diagnostic)
        summary_key = (
            diagnostic.panel_type,
            diagnostic.species_name,
            diagnostic.clade,
            diagnostic.evidence_taxid,
            diagnostic.evidence_rank,
            diagnostic.k,
        )
        summary_counts[summary_key] += 1
    summary = [
        {
            "panel_type": panel_type,
            "species_name": species_name,
            "clade": clade,
            "evidence_taxid": evidence_taxid,
            "evidence_rank": evidence_rank,
            "k": k,
            "diagnostic_kmers": count,
        }
        for (panel_type, species_name, clade, evidence_taxid, evidence_rank, k), count
        in sorted(summary_counts.items())
    ]
    return retained, summary


def build_target_evidence_sqlite(
    *,
    genome_configs: list[GenomeConfig],
    k_values: list[int],
    sqlite_path: str | Path,
    batch_size: int = 50000,
    logger: logging.Logger | None = None,
) -> TargetEvidenceBuildResult:
    """Build an SQLite target-candidate database from genome configs.

    Parameters
    ----------
    genome_configs : list of GenomeConfig
        Genome records. Records with ``is_target`` are inserted as candidates;
        all others are streamed as filters/downgrade evidence.
    k_values : list of int
        K-mer lengths.
    sqlite_path : str or pathlib.Path
        SQLite output path.
    batch_size : int, optional
        Commit interval.
    logger : logging.Logger | None, optional
        Optional logger.

    Returns
    -------
    TargetEvidenceBuildResult
        Build metadata and collection summaries.
    """
    target_configs = [config for config in genome_configs if config.is_target]
    non_target_configs = [config for config in genome_configs if not config.is_target]
    if not target_configs:
        raise ValueError("At least one genome must have role target_species")
    if logger:
        logger.info(
            "Starting SQLite target-evidence build with %d target and %d non-target genome(s)",
            len(target_configs),
            len(non_target_configs),
        )
        logger.info("SQLite candidate database: %s", sqlite_path)

    sqlite_file = Path(sqlite_path)
    if sqlite_file.exists():
        sqlite_file.unlink()
    initialise_target_evidence_database(sqlite_path=sqlite_file)
    target_summary = collect_target_candidates_sqlite(
        genome_configs=target_configs,
        k_values=k_values,
        sqlite_path=sqlite_file,
        batch_size=batch_size,
        logger=logger,
    )
    non_target_summary = scan_non_target_matches_sqlite(
        genome_configs=non_target_configs,
        k_values=k_values,
        sqlite_path=sqlite_file,
        batch_size=batch_size,
        logger=logger,
    )
    connection = _connect(sqlite_path=sqlite_file)
    try:
        total_candidates = connection.execute(
            "SELECT COUNT(*) FROM target_kmers"
        ).fetchone()[0]
        shared_candidates = connection.execute(
            "SELECT COUNT(*) FROM target_kmers WHERE non_target_taxids != ''"
        ).fetchone()[0]
        build_summary = [
            {
                "summary_name": "target_genomes",
                "summary_value": len(target_configs),
            },
            {
                "summary_name": "non_target_genomes",
                "summary_value": len(non_target_configs),
            },
            {
                "summary_name": "target_candidate_kmers",
                "summary_value": int(total_candidates),
            },
            {
                "summary_name": "target_candidates_seen_in_non_targets",
                "summary_value": int(shared_candidates),
            },
        ]
        _record_event(
            connection=connection,
            stage="complete_sqlite_target_evidence_build",
            detail=str(sqlite_file),
            n_records=int(total_candidates),
        )
        connection.commit()
    finally:
        connection.close()
    return TargetEvidenceBuildResult(
        sqlite_path=sqlite_file,
        collection_summary=[*target_summary, *non_target_summary],
        build_summary=build_summary,
    )
