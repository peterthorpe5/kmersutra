"""Low-memory all-candidate panel building for KmerSutra.

This module extends the SQLite-backed target-evidence build into a query-
agnostic build mode. Instead of requiring the user to name the expected target
species in advance, it iterates over each eligible candidate species, treats
that species as the temporary candidate, streams all other genomes as filters,
and retains validated evidence for every candidate taxon.

The design is deliberately slower than the target-only build, but it avoids the
large global in-memory dictionary used by the compact builder and therefore
supports realistic unknown-sample panels where the true species is not known in
advance.
"""

from __future__ import annotations

import logging
import sqlite3
from collections import defaultdict
from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from pathlib import Path

from kmersutra.build_panel import DiagnosticKmer
from kmersutra.config import GenomeConfig
from kmersutra.target_evidence import (
    _connect,
    collect_target_candidates_sqlite,
    initialise_target_evidence_database,
    iter_target_evidence_diagnostics,
    scan_non_target_matches_sqlite,
)
from kmersutra.taxonomy import CORE_RANK_ORDER, TaxonomyDatabase

DEFAULT_EXCLUDED_CANDIDATE_ROLES = {
    "exclude",
    "host",
    "host_or_background",
    "background_host",
    "host_background",
    "background",
    "environmental_background",
}


@dataclass(frozen=True)
class AllCandidateEvidenceBuildResult:
    """Result metadata from an all-candidate evidence build.

    Attributes
    ----------
    retained_sqlite_path : pathlib.Path
        SQLite database containing retained diagnostic evidence.
    collection_summary : list[dict[str, object]]
        Per-candidate and per-filter collection summaries.
    build_summary : list[dict[str, object]]
        Coarse build summary records.
    panel_summary : list[dict[str, object]]
        Diagnostic k-mer counts grouped by evidence bucket.
    """

    retained_sqlite_path: Path
    collection_summary: list[dict[str, object]]
    build_summary: list[dict[str, object]]
    panel_summary: list[dict[str, object]]


def _retention_key(item: DiagnosticKmer) -> tuple[str, str, str, str, int]:
    """Return the quota key for one diagnostic k-mer.

    Parameters
    ----------
    item : DiagnosticKmer
        Diagnostic k-mer record.

    Returns
    -------
    tuple[str, str, str, str, int]
        Retention key based on evidence type, taxon and k value.
    """
    return (
        item.panel_type,
        item.evidence_taxid,
        item.species_name,
        item.clade,
        item.k,
    )


def _normalise_roles(roles: Iterable[str] | None) -> set[str]:
    """Normalise optional role strings.

    Parameters
    ----------
    roles : iterable of str or None
        Role names.

    Returns
    -------
    set[str]
        Non-empty role names.
    """
    if roles is None:
        return set()
    return {str(role) for role in roles if str(role)}


def select_candidate_genomes(
    *,
    genome_configs: Iterable[GenomeConfig],
    candidate_roles: Iterable[str] | None = None,
    excluded_roles: Iterable[str] | None = None,
) -> list[GenomeConfig]:
    """Select genomes that should become reportable candidates.

    Parameters
    ----------
    genome_configs : iterable of GenomeConfig
        Genome configuration records.
    candidate_roles : iterable of str or None, optional
        If provided, only genomes with these roles are reportable candidates.
    excluded_roles : iterable of str or None, optional
        Roles to exclude when candidate_roles is not provided.

    Returns
    -------
    list[GenomeConfig]
        Candidate genome records.
    """
    configs = list(genome_configs)
    role_filter = _normalise_roles(candidate_roles)
    excluded = _normalise_roles(excluded_roles) or DEFAULT_EXCLUDED_CANDIDATE_ROLES
    if role_filter:
        return [config for config in configs if config.role in role_filter]
    return [config for config in configs if config.role not in excluded]


def group_candidate_genomes_by_species(
    *,
    genome_configs: Iterable[GenomeConfig],
) -> list[tuple[str, list[GenomeConfig]]]:
    """Group candidate genome records by species name.

    Parameters
    ----------
    genome_configs : iterable of GenomeConfig
        Candidate genome records.

    Returns
    -------
    list[tuple[str, list[GenomeConfig]]]
        Sorted species-to-genome groups.
    """
    grouped: dict[str, list[GenomeConfig]] = defaultdict(list)
    for config in genome_configs:
        grouped[config.species_name].append(config)
    return [(species, grouped[species]) for species in sorted(grouped)]


def initialise_retained_evidence_database(*, sqlite_path: str | Path) -> None:
    """Create the SQLite schema for retained all-candidate evidence.

    Parameters
    ----------
    sqlite_path : str or pathlib.Path
        SQLite database path.
    """
    connection = _connect(sqlite_path=sqlite_path)
    try:
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
                candidate_species TEXT NOT NULL,
                PRIMARY KEY (k, kmer, evidence_taxid, evidence_rank)
            )
            """
        )
        connection.execute(
            """
            CREATE INDEX IF NOT EXISTS retained_kmers_evidence_idx
            ON retained_kmers(evidence_rank, evidence_name, species_name, k)
            """
        )
        connection.commit()
    finally:
        connection.close()


def _insert_retained_diagnostic(
    *,
    connection: sqlite3.Connection,
    diagnostic: DiagnosticKmer,
    candidate_species: str,
) -> bool:
    """Insert one retained diagnostic k-mer.

    Parameters
    ----------
    connection : sqlite3.Connection
        Open retained-evidence database connection.
    diagnostic : DiagnosticKmer
        Diagnostic k-mer record.
    candidate_species : str
        Candidate species group that produced the diagnostic.

    Returns
    -------
    bool
        True if a new row was inserted, otherwise False for duplicate evidence.
    """
    cursor = connection.execute(
        """
        INSERT OR IGNORE INTO retained_kmers(
            k, kmer, panel_type, species_name, clade, source_genomes,
            source_contigs, example_position, evidence_taxid, evidence_name,
            evidence_rank, lineage_taxids, source_taxids, candidate_species
        )
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
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
            candidate_species,
        ),
    )
    return cursor.rowcount == 1


def iter_retained_all_candidate_diagnostics(
    *,
    sqlite_path: str | Path,
    chunk_size: int = 10000,
) -> Iterator[DiagnosticKmer]:
    """Yield retained all-candidate diagnostics from SQLite.

    Parameters
    ----------
    sqlite_path : str or pathlib.Path
        Retained evidence SQLite path.
    chunk_size : int, optional
        Number of rows to fetch per SQLite batch.

    Yields
    ------
    DiagnosticKmer
        Retained diagnostic k-mer records.
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


def summarise_retained_all_candidate_evidence(
    *,
    sqlite_path: str | Path,
) -> list[dict[str, object]]:
    """Summarise retained all-candidate evidence records.

    Parameters
    ----------
    sqlite_path : str or pathlib.Path
        Retained evidence SQLite path.

    Returns
    -------
    list[dict[str, object]]
        Evidence-count summary rows.
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
            ORDER BY evidence_rank, evidence_name, species_name, k
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


def build_all_candidate_evidence_sqlite(
    *,
    genome_configs: list[GenomeConfig],
    k_values: list[int],
    retained_sqlite_path: str | Path,
    work_sqlite_path: str | Path,
    taxonomy_db: TaxonomyDatabase | None = None,
    target_taxid: str = "",
    preferred_ranks: list[str] | None = None,
    candidate_roles: Iterable[str] | None = None,
    excluded_candidate_roles: Iterable[str] | None = None,
    batch_size: int = 50000,
    max_per_evidence_per_k: int | None = None,
    logger: logging.Logger | None = None,
) -> AllCandidateEvidenceBuildResult:
    """Build a low-memory query-agnostic evidence database.

    Parameters
    ----------
    genome_configs : list of GenomeConfig
        Genome records used for candidate and filter evidence.
    k_values : list of int
        K-mer lengths.
    retained_sqlite_path : str or pathlib.Path
        SQLite path for retained diagnostic evidence.
    work_sqlite_path : str or pathlib.Path
        Temporary SQLite path reused for each candidate species.
    taxonomy_db : TaxonomyDatabase or None, optional
        Optional taxonomy database used for evidence-level assignment.
    target_taxid : str, optional
        Optional taxonomic subtree restriction. Leave empty for true
        query-agnostic broad panels including outgroups.
    preferred_ranks : list[str] or None, optional
        Evidence ranks to retain.
    candidate_roles : iterable of str or None, optional
        Optional role whitelist for reportable candidates.
    excluded_candidate_roles : iterable of str or None, optional
        Roles excluded from candidate status when no whitelist is supplied.
    batch_size : int, optional
        SQLite commit interval.
    max_per_evidence_per_k : int or None, optional
        Optional cap per evidence bucket and k value.
    logger : logging.Logger or None, optional
        Optional logger.

    Returns
    -------
    AllCandidateEvidenceBuildResult
        Build metadata and retained evidence summaries.
    """
    if batch_size <= 0:
        raise ValueError("batch_size must be positive")
    if max_per_evidence_per_k is not None and max_per_evidence_per_k <= 0:
        raise ValueError("max_per_evidence_per_k must be positive")

    candidate_configs = select_candidate_genomes(
        genome_configs=genome_configs,
        candidate_roles=candidate_roles,
        excluded_roles=excluded_candidate_roles,
    )
    candidate_groups = group_candidate_genomes_by_species(
        genome_configs=candidate_configs,
    )
    if not candidate_groups:
        raise ValueError("No candidate genomes were selected for all-candidate build")

    retained_path = Path(retained_sqlite_path)
    work_path = Path(work_sqlite_path)
    if retained_path.exists():
        retained_path.unlink()
    if work_path.exists():
        work_path.unlink()
    initialise_retained_evidence_database(sqlite_path=retained_path)

    ranks = preferred_ranks or CORE_RANK_ORDER
    retained_counts: dict[tuple[str, str, str, str, int], int] = defaultdict(int)
    collection_summary: list[dict[str, object]] = []
    candidate_build_rows: list[dict[str, object]] = []
    total_considered = 0
    total_retained = 0
    total_duplicates = 0
    total_skipped_by_limit = 0

    retained_connection = _connect(sqlite_path=retained_path)
    try:
        for group_index, (candidate_species, group_configs) in enumerate(
            candidate_groups,
            start=1,
        ):
            filter_configs = [
                config
                for config in genome_configs
                if config.species_name != candidate_species
                and config.role not in {"exclude"}
            ]
            if logger:
                logger.info(
                    "Building all-candidate evidence for species %d/%d: %s "
                    "(%d candidate genome(s), %d filter genome(s))",
                    group_index,
                    len(candidate_groups),
                    candidate_species,
                    len(group_configs),
                    len(filter_configs),
                )
            if work_path.exists():
                work_path.unlink()
            initialise_target_evidence_database(sqlite_path=work_path)
            candidate_summary = collect_target_candidates_sqlite(
                genome_configs=group_configs,
                k_values=k_values,
                sqlite_path=work_path,
                batch_size=batch_size,
                logger=logger,
            )
            filter_summary = scan_non_target_matches_sqlite(
                genome_configs=filter_configs,
                k_values=k_values,
                sqlite_path=work_path,
                batch_size=batch_size,
                logger=logger,
            )
            for record in [*candidate_summary, *filter_summary]:
                record = dict(record)
                record["candidate_species_group"] = candidate_species
                collection_summary.append(record)

            considered = 0
            retained = 0
            duplicates = 0
            skipped_by_limit = 0
            for diagnostic in iter_target_evidence_diagnostics(
                sqlite_path=work_path,
                taxonomy_db=taxonomy_db,
                target_taxid=target_taxid,
                preferred_ranks=ranks,
            ):
                considered += 1
                retention_key = _retention_key(diagnostic)
                if max_per_evidence_per_k is not None:
                    if retained_counts[retention_key] >= max_per_evidence_per_k:
                        skipped_by_limit += 1
                        continue
                inserted = _insert_retained_diagnostic(
                    connection=retained_connection,
                    diagnostic=diagnostic,
                    candidate_species=candidate_species,
                )
                if inserted:
                    retained_counts[retention_key] += 1
                    retained += 1
                    if retained % batch_size == 0:
                        retained_connection.commit()
                else:
                    duplicates += 1
            retained_connection.commit()
            candidate_build_rows.append(
                {
                    "candidate_species": candidate_species,
                    "candidate_genomes": len(group_configs),
                    "filter_genomes": len(filter_configs),
                    "diagnostics_considered": considered,
                    "diagnostics_retained": retained,
                    "diagnostic_duplicates_ignored": duplicates,
                    "diagnostics_skipped_by_limit": skipped_by_limit,
                }
            )
            total_considered += considered
            total_retained += retained
            total_duplicates += duplicates
            total_skipped_by_limit += skipped_by_limit
            if logger:
                logger.info(
                    "Finished candidate species %s: considered=%d, "
                    "retained=%d, duplicates=%d, skipped_by_limit=%d",
                    candidate_species,
                    considered,
                    retained,
                    duplicates,
                    skipped_by_limit,
                )
    finally:
        retained_connection.close()
        work_path.unlink(missing_ok=True)

    panel_summary = summarise_retained_all_candidate_evidence(
        sqlite_path=retained_path,
    )
    build_summary = [
        {
            "summary_name": "genome_records",
            "summary_value": len(genome_configs),
        },
        {
            "summary_name": "candidate_genomes",
            "summary_value": len(candidate_configs),
        },
        {
            "summary_name": "candidate_species_groups",
            "summary_value": len(candidate_groups),
        },
        {
            "summary_name": "diagnostics_considered",
            "summary_value": total_considered,
        },
        {
            "summary_name": "diagnostics_retained",
            "summary_value": total_retained,
        },
        {
            "summary_name": "diagnostic_duplicates_ignored",
            "summary_value": total_duplicates,
        },
        {
            "summary_name": "diagnostics_skipped_by_limit",
            "summary_value": total_skipped_by_limit,
        },
        {
            "summary_name": "retained_sqlite_path",
            "summary_value": str(retained_path),
        },
    ]
    build_summary.extend(
        {
            "summary_name": f"candidate::{row['candidate_species']}::retained",
            "summary_value": row["diagnostics_retained"],
        }
        for row in candidate_build_rows
    )
    collection_summary.extend(candidate_build_rows)
    return AllCandidateEvidenceBuildResult(
        retained_sqlite_path=retained_path,
        collection_summary=collection_summary,
        build_summary=build_summary,
        panel_summary=panel_summary,
    )
