"""Hierarchical cascade screening for KmerSutra modules.

This module implements a sample-level cascade screen. A small broad-rank gate
panel is used to decide which taxonomic modules should be activated. Activated
modules are then screened with their more detailed genus/species panels. The
current implementation deliberately reuses the existing exact/fuzzy screening
engine so that hierarchical mode changes panel selection and reporting, not the
underlying k-mer matching semantics.
"""

from __future__ import annotations

import logging
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

from kmersutra.io import read_tsv
from kmersutra.screen_reads import KmerHit, screen_file_for_species_kmers


MODULE_MANIFEST_FIELDNAMES = [
    "module_id",
    "module_name",
    "rank",
    "parent_module_id",
    "gate_panel_path",
    "module_panel_path",
    "min_gate_unique_kmers",
    "min_gate_positive_sequences",
    "min_gate_k_values",
    "min_gate_best_k",
]

MODULE_ACTIVATION_FIELDNAMES = [
    "sample_id",
    "module_id",
    "module_name",
    "rank",
    "parent_module_id",
    "gate_panel_path",
    "module_panel_path",
    "gate_n_hits",
    "gate_unique_kmers",
    "gate_positive_sequences",
    "gate_k_values",
    "gate_best_k",
    "activated",
    "activation_reason",
]


@dataclass(frozen=True)
class ModuleDefinition:
    """Definition of a screenable KmerSutra taxonomic module.

    Parameters
    ----------
    module_id : str
        Stable module identifier, for example ``apicomplexa`` or
        ``plasmodium``.
    module_name : str
        Human-readable module name.
    rank : str
        Taxonomic or organisational rank represented by the module.
    parent_module_id : str
        Parent module identifier. Empty means root-level module.
    gate_panel_path : str
        Optional broad-rank gate panel used to activate this module.
    module_panel_path : str
        Detailed panel screened if the module is activated.
    min_gate_unique_kmers : int
        Minimum unique gate k-mers required for activation.
    min_gate_positive_sequences : int
        Minimum independent gate-positive reads/contigs required.
    min_gate_k_values : int
        Minimum number of positive k values required.
    min_gate_best_k : int
        Minimum longest positive k value required.
    """

    module_id: str
    module_name: str
    rank: str
    parent_module_id: str
    gate_panel_path: str
    module_panel_path: str
    min_gate_unique_kmers: int = 1
    min_gate_positive_sequences: int = 1
    min_gate_k_values: int = 1
    min_gate_best_k: int = 0

    @classmethod
    def from_record(
        cls,
        *,
        record: dict[str, str],
        manifest_dir: Path,
    ) -> "ModuleDefinition":
        """Create a module definition from a manifest record.

        Parameters
        ----------
        record : dict[str, str]
            Manifest row.
        manifest_dir : pathlib.Path
            Directory used to resolve relative panel paths.

        Returns
        -------
        ModuleDefinition
            Parsed module definition.

        Raises
        ------
        ValueError
            If required identifiers or panel paths are missing.
        """
        module_id = record.get("module_id", "").strip()
        if not module_id:
            raise ValueError("Module manifest row is missing module_id")

        gate_panel = _resolve_optional_path(
            value=record.get("gate_panel_path", ""),
            base_dir=manifest_dir,
        )
        module_panel = _resolve_optional_path(
            value=record.get("module_panel_path", ""),
            base_dir=manifest_dir,
        )
        if not gate_panel and not module_panel:
            raise ValueError(
                f"Module {module_id} must define gate_panel_path, module_panel_path, or both"
            )

        return cls(
            module_id=module_id,
            module_name=record.get("module_name", module_id).strip() or module_id,
            rank=record.get("rank", "module").strip() or "module",
            parent_module_id=record.get("parent_module_id", "").strip(),
            gate_panel_path=gate_panel,
            module_panel_path=module_panel,
            min_gate_unique_kmers=_parse_positive_int(
                value=record.get("min_gate_unique_kmers", ""),
                default=1,
                name="min_gate_unique_kmers",
            ),
            min_gate_positive_sequences=_parse_positive_int(
                value=record.get("min_gate_positive_sequences", ""),
                default=1,
                name="min_gate_positive_sequences",
            ),
            min_gate_k_values=_parse_positive_int(
                value=record.get("min_gate_k_values", ""),
                default=1,
                name="min_gate_k_values",
            ),
            min_gate_best_k=_parse_non_negative_int(
                value=record.get("min_gate_best_k", ""),
                default=0,
                name="min_gate_best_k",
            ),
        )


def _parse_positive_int(*, value: object, default: int, name: str) -> int:
    """Parse a positive integer with a default."""
    parsed = default if value in (None, "") else int(str(value))
    if parsed <= 0:
        raise ValueError(f"{name} must be positive")
    return parsed


def _parse_non_negative_int(*, value: object, default: int, name: str) -> int:
    """Parse a non-negative integer with a default."""
    parsed = default if value in (None, "") else int(str(value))
    if parsed < 0:
        raise ValueError(f"{name} must be zero or greater")
    return parsed


def _resolve_optional_path(*, value: str | None, base_dir: Path) -> str:
    """Resolve an optional manifest path."""
    if value is None or not str(value).strip():
        return ""
    path = Path(str(value).strip())
    if not path.is_absolute():
        path = base_dir / path
    return str(path)


def load_module_manifest(*, manifest_path: str | Path) -> list[ModuleDefinition]:
    """Load a hierarchical module manifest.

    Parameters
    ----------
    manifest_path : str or pathlib.Path
        Module manifest TSV.

    Returns
    -------
    list[ModuleDefinition]
        Parsed modules in manifest order.

    Raises
    ------
    FileNotFoundError
        If the manifest or referenced panels are missing.
    ValueError
        If the manifest is empty, malformed, or contains duplicate module IDs.
    """
    path = Path(manifest_path)
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"Module manifest is missing or empty: {path}")

    records = read_tsv(input_path=path)
    if not records:
        raise ValueError(f"Module manifest contains no modules: {path}")

    modules = [
        ModuleDefinition.from_record(record=record, manifest_dir=path.parent)
        for record in records
    ]
    module_ids = [module.module_id for module in modules]
    duplicates = sorted({module_id for module_id in module_ids if module_ids.count(module_id) > 1})
    if duplicates:
        raise ValueError(f"Duplicate module_id values in manifest: {', '.join(duplicates)}")

    known_ids = set(module_ids)
    for module in modules:
        if module.parent_module_id and module.parent_module_id not in known_ids:
            raise ValueError(
                f"Module {module.module_id} references unknown parent_module_id "
                f"{module.parent_module_id}"
            )
        for panel_path in [module.gate_panel_path, module.module_panel_path]:
            if panel_path and not Path(panel_path).is_file():
                raise FileNotFoundError(
                    f"Panel path for module {module.module_id} is missing: {panel_path}"
                )
    return modules


def order_modules_by_parentage(*, modules: list[ModuleDefinition]) -> list[ModuleDefinition]:
    """Return modules sorted so parents are screened before children.

    Parameters
    ----------
    modules : list[ModuleDefinition]
        Parsed module definitions.

    Returns
    -------
    list[ModuleDefinition]
        Modules sorted by parent depth and then module identifier.

    Raises
    ------
    ValueError
        If a cyclic parent relationship is detected.
    """
    by_id = {module.module_id: module for module in modules}
    visiting: set[str] = set()
    depth_cache: dict[str, int] = {}

    def depth(module_id: str) -> int:
        if module_id in depth_cache:
            return depth_cache[module_id]
        if module_id in visiting:
            raise ValueError("Cyclic module parent relationship detected")
        visiting.add(module_id)
        module = by_id[module_id]
        if not module.parent_module_id:
            value = 0
        else:
            value = depth(module.parent_module_id) + 1
        visiting.remove(module_id)
        depth_cache[module_id] = value
        return value

    return sorted(modules, key=lambda module: (depth(module.module_id), module.module_id))


def summarise_gate_hits(*, hits: Iterable[KmerHit]) -> dict[str, int]:
    """Summarise gate hits for module activation.

    Parameters
    ----------
    hits : iterable of KmerHit
        Hits observed against a gate panel.

    Returns
    -------
    dict[str, int]
        Hit, unique-marker, independent-sequence, k-value and best-k summary.
    """
    hit_list = list(hits)
    unique_kmers = {hit.matched_kmer for hit in hit_list}
    positive_sequences = {hit.sequence_id for hit in hit_list}
    k_values = {int(hit.k) for hit in hit_list}
    return {
        "gate_n_hits": len(hit_list),
        "gate_unique_kmers": len(unique_kmers),
        "gate_positive_sequences": len(positive_sequences),
        "gate_k_values": len(k_values),
        "gate_best_k": max(k_values) if k_values else 0,
    }


def gate_summary_passes(*, module: ModuleDefinition, summary: dict[str, int]) -> bool:
    """Return whether a gate summary activates a module.

    Parameters
    ----------
    module : ModuleDefinition
        Module thresholds.
    summary : dict[str, int]
        Gate hit summary from :func:`summarise_gate_hits`.

    Returns
    -------
    bool
        True when all module activation thresholds are met.
    """
    return (
        summary["gate_unique_kmers"] >= module.min_gate_unique_kmers
        and summary["gate_positive_sequences"] >= module.min_gate_positive_sequences
        and summary["gate_k_values"] >= module.min_gate_k_values
        and summary["gate_best_k"] >= module.min_gate_best_k
    )


def unique_hits(*, hits: Iterable[KmerHit]) -> list[KmerHit]:
    """Deduplicate hits produced by overlapping gate and module panels.

    Parameters
    ----------
    hits : iterable of KmerHit
        Candidate hit records.

    Returns
    -------
    list[KmerHit]
        First-seen unique hits in deterministic insertion order.
    """
    seen: set[tuple[object, ...]] = set()
    output: list[KmerHit] = []
    for hit in hits:
        key = (
            hit.sample_id,
            hit.sequence_id,
            hit.sequence_type,
            hit.k,
            hit.query_position,
            hit.matched_kmer,
            hit.query_kmer,
            hit.mismatches,
            hit.panel_type,
            hit.species_name,
            hit.clade,
            hit.evidence_taxid,
            hit.evidence_name,
            hit.evidence_rank,
        )
        if key in seen:
            continue
        seen.add(key)
        output.append(hit)
    return output


@dataclass(frozen=True)
class HierarchicalScreenResult:
    """Result from hierarchical screening.

    Attributes
    ----------
    hits : list[KmerHit]
        Combined gate and activated-module hits.
    activation_records : list[dict[str, object]]
        Per-module activation summary records.
    """

    hits: list[KmerHit]
    activation_records: list[dict[str, object]]


def _module_cache_path(
    *,
    panel_path: str,
    panel_cache_path: str | Path | None,
) -> str | Path | None:
    """Resolve a safe cache path for a module panel."""
    if panel_cache_path is None:
        return None
    base = Path(panel_cache_path)
    safe_name = Path(panel_path).name.replace(".gz", "").replace(".tsv", "")
    return base.with_name(f"{base.stem}.{safe_name}{base.suffix}")


def _screen_panel(
    *,
    input_path: str | Path,
    panel_path: str,
    sample_id: str,
    input_format: str,
    max_mismatches: int,
    fuzzy_min_k: int,
    threads: int,
    chunk_size: int,
    max_pending_chunks: int | None,
    panel_cache_path: str | Path | None,
    use_panel_cache: bool,
    write_panel_cache: bool,
    profile_records: list[dict[str, object]] | None,
    logger: logging.Logger | None,
) -> list[KmerHit]:
    """Screen one panel with cache handling suitable for module mode."""
    return screen_file_for_species_kmers(
        input_path=input_path,
        panel_path=panel_path,
        sample_id=sample_id,
        input_format=input_format,
        max_mismatches=max_mismatches,
        fuzzy_min_k=fuzzy_min_k,
        threads=threads,
        chunk_size=chunk_size,
        max_pending_chunks=max_pending_chunks,
        panel_cache_path=_module_cache_path(
            panel_path=panel_path,
            panel_cache_path=panel_cache_path,
        ),
        use_panel_cache=use_panel_cache,
        write_panel_cache=write_panel_cache,
        profile_records=profile_records,
        logger=logger,
    )


def screen_file_hierarchical(
    *,
    input_path: str | Path,
    module_manifest_path: str | Path,
    sample_id: str,
    input_format: str,
    max_mismatches: int = 0,
    fuzzy_min_k: int = 71,
    threads: int = 1,
    chunk_size: int = 1000,
    max_pending_chunks: int | None = None,
    panel_cache_path: str | Path | None = None,
    use_panel_cache: bool = False,
    write_panel_cache: bool = False,
    hierarchical_fail_open: bool = False,
    profile_records: list[dict[str, object]] | None = None,
    logger: logging.Logger | None = None,
) -> HierarchicalScreenResult:
    """Run sample-level hierarchical cascade screening.

    Parameters
    ----------
    input_path : str or pathlib.Path
        Query FASTA/FASTQ path.
    module_manifest_path : str or pathlib.Path
        Hierarchical module manifest TSV.
    sample_id : str
        Sample identifier.
    input_format : str
        Input format, either ``fastq`` or ``fasta``.
    max_mismatches : int, optional
        Maximum mismatches for fuzzy matching.
    fuzzy_min_k : int, optional
        Minimum k value eligible for fuzzy matching.
    threads : int, optional
        Worker count for each panel screen.
    chunk_size : int, optional
        Records per chunk.
    max_pending_chunks : int or None, optional
        Maximum queued chunks.
    panel_cache_path : str or pathlib.Path or None, optional
        Optional base cache path. In module mode, panel-specific suffixes are
        derived from this path so modules do not overwrite each other's caches.
    use_panel_cache : bool, optional
        Use current panel caches when available.
    write_panel_cache : bool, optional
        Write panel caches after TSV loading.
    hierarchical_fail_open : bool, optional
        If no module passes its gate but weak gate evidence is observed,
        activate detailed module panels in fail-open mode. This preserves
        sensitivity for unresolved/novel signals while recording the reason.
    profile_records : list[dict[str, object]] or None, optional
        Optional profile sink.
    logger : logging.Logger or None, optional
        Logger.

    Returns
    -------
    HierarchicalScreenResult
        Combined hits and activation records.
    """
    modules = load_module_manifest(manifest_path=module_manifest_path)
    if logger:
        logger.info("Loaded %d hierarchical module definitions", len(modules))

    modules = order_modules_by_parentage(modules=modules)
    all_hits: list[KmerHit] = []
    activation_records: list[dict[str, object]] = []
    active_modules: set[str] = set()
    weak_gate_modules: set[str] = set()

    for module in modules:
        parent_active = not module.parent_module_id or module.parent_module_id in active_modules
        if not parent_active:
            summary = {
                "gate_n_hits": 0,
                "gate_unique_kmers": 0,
                "gate_positive_sequences": 0,
                "gate_k_values": 0,
                "gate_best_k": 0,
            }
            activation_records.append(_activation_record(
                sample_id=sample_id,
                module=module,
                summary=summary,
                activated=False,
                reason="parent_not_active",
            ))
            continue

        gate_hits: list[KmerHit] = []
        if module.gate_panel_path:
            if logger:
                logger.info(
                    "Screening gate panel for module %s: %s",
                    module.module_id,
                    module.gate_panel_path,
                )
            gate_hits = _screen_panel(
                input_path=input_path,
                panel_path=module.gate_panel_path,
                sample_id=sample_id,
                input_format=input_format,
                max_mismatches=max_mismatches,
                fuzzy_min_k=fuzzy_min_k,
                threads=threads,
                chunk_size=chunk_size,
                max_pending_chunks=max_pending_chunks,
                panel_cache_path=panel_cache_path,
                use_panel_cache=use_panel_cache,
                write_panel_cache=write_panel_cache,
                profile_records=profile_records,
                logger=logger,
            )
            all_hits.extend(gate_hits)
            summary = summarise_gate_hits(hits=gate_hits)
            passes = gate_summary_passes(module=module, summary=summary)
            if gate_hits and not passes:
                weak_gate_modules.add(module.module_id)
        else:
            summary = {
                "gate_n_hits": 0,
                "gate_unique_kmers": 0,
                "gate_positive_sequences": 0,
                "gate_k_values": 0,
                "gate_best_k": 0,
            }
            passes = True

        if passes:
            active_modules.add(module.module_id)
            activation_records.append(_activation_record(
                sample_id=sample_id,
                module=module,
                summary=summary,
                activated=True,
                reason="gate_passed" if module.gate_panel_path else "root_without_gate",
            ))
        else:
            activation_records.append(_activation_record(
                sample_id=sample_id,
                module=module,
                summary=summary,
                activated=False,
                reason="gate_below_threshold",
            ))

    if hierarchical_fail_open and not active_modules and weak_gate_modules:
        if logger:
            logger.info(
                "Hierarchical fail-open activated because weak gate evidence was observed"
            )
        for module in modules:
            if module.module_panel_path:
                active_modules.add(module.module_id)
                activation_records.append(_activation_record(
                    sample_id=sample_id,
                    module=module,
                    summary={
                        "gate_n_hits": 0,
                        "gate_unique_kmers": 0,
                        "gate_positive_sequences": 0,
                        "gate_k_values": 0,
                        "gate_best_k": 0,
                    },
                    activated=True,
                    reason="fail_open_weak_gate_signal",
                ))

    screened_module_panels: set[str] = set()
    for module in modules:
        if module.module_id not in active_modules or not module.module_panel_path:
            continue
        if module.module_panel_path in screened_module_panels:
            continue
        screened_module_panels.add(module.module_panel_path)
        if logger:
            logger.info(
                "Screening activated module %s panel: %s",
                module.module_id,
                module.module_panel_path,
            )
        module_hits = _screen_panel(
            input_path=input_path,
            panel_path=module.module_panel_path,
            sample_id=sample_id,
            input_format=input_format,
            max_mismatches=max_mismatches,
            fuzzy_min_k=fuzzy_min_k,
            threads=threads,
            chunk_size=chunk_size,
            max_pending_chunks=max_pending_chunks,
            panel_cache_path=panel_cache_path,
            use_panel_cache=use_panel_cache,
            write_panel_cache=write_panel_cache,
            profile_records=profile_records,
            logger=logger,
        )
        all_hits.extend(module_hits)

    return HierarchicalScreenResult(
        hits=unique_hits(hits=all_hits),
        activation_records=activation_records,
    )


def _activation_record(
    *,
    sample_id: str,
    module: ModuleDefinition,
    summary: dict[str, int],
    activated: bool,
    reason: str,
) -> dict[str, object]:
    """Build one module activation record."""
    return {
        "sample_id": sample_id,
        "module_id": module.module_id,
        "module_name": module.module_name,
        "rank": module.rank,
        "parent_module_id": module.parent_module_id,
        "gate_panel_path": module.gate_panel_path,
        "module_panel_path": module.module_panel_path,
        "gate_n_hits": summary["gate_n_hits"],
        "gate_unique_kmers": summary["gate_unique_kmers"],
        "gate_positive_sequences": summary["gate_positive_sequences"],
        "gate_k_values": summary["gate_k_values"],
        "gate_best_k": summary["gate_best_k"],
        "activated": str(bool(activated)),
        "activation_reason": reason,
    }
