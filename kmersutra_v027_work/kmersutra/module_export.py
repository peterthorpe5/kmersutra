"""Automatic hierarchical module export for KmerSutra panels.

The scalable screening engine expects a module manifest plus gate/detail panel
files. This module derives those files from a completed flat KmerSutra panel so
large databases can be screened in hierarchical mode without manual splitting.
The export is intentionally conservative: every marker written to a module is
copied from the already validated panel, and module generation never changes
marker evidence ranks.
"""

from __future__ import annotations

import logging
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from kmersutra.hierarchical import MODULE_MANIFEST_FIELDNAMES
from kmersutra.io import open_text, write_tsv

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

MODULE_EXPORT_SUMMARY_FIELDNAMES = [
    "module_id",
    "module_name",
    "rank",
    "parent_module_id",
    "gate_records",
    "module_records",
    "gate_panel_path",
    "module_panel_path",
]

DEFAULT_GATE_RANKS = {"phylum", "class", "order", "family", "genus"}
_SPECIES_NAME_RE = re.compile(r"^[A-Z][A-Za-z0-9_.-]+\s+.+")
_SAFE_ID_RE = re.compile(r"[^A-Za-z0-9_.-]+")


@dataclass(frozen=True)
class ModuleExportConfig:
    """Configuration for automatic hierarchical module export.

    Parameters
    ----------
    gate_ranks : set[str]
        Evidence ranks that may be used as broad gate markers.
    min_gate_unique_kmers : int
        Default module activation threshold for unique gate k-mers.
    min_gate_positive_sequences : int
        Default module activation threshold for positive reads/contigs.
    min_gate_k_values : int
        Default module activation threshold for distinct positive k values.
    min_gate_best_k : int
        Default module activation threshold for the longest positive k value.
    allow_species_gate_fallback : bool
        If True, modules with no broad-rank gate records receive a small
        species-marker gate copied from their detail panel. This preserves
        activation for taxa without genus/core markers, while the summary makes
        the fallback explicit via gate counts.
    max_gate_records_per_module_per_k : int
        Optional cap for gate panels. Use 0 for no cap.
    """

    gate_ranks: set[str]
    min_gate_unique_kmers: int = 1
    min_gate_positive_sequences: int = 1
    min_gate_k_values: int = 1
    min_gate_best_k: int = 0
    allow_species_gate_fallback: bool = True
    max_gate_records_per_module_per_k: int = 0

    def validate(self) -> None:
        """Validate the export configuration.

        Raises
        ------
        ValueError
            If numeric thresholds are invalid.
        """
        if self.min_gate_unique_kmers <= 0:
            raise ValueError("min_gate_unique_kmers must be positive")
        if self.min_gate_positive_sequences <= 0:
            raise ValueError("min_gate_positive_sequences must be positive")
        if self.min_gate_k_values <= 0:
            raise ValueError("min_gate_k_values must be positive")
        if self.min_gate_best_k < 0:
            raise ValueError("min_gate_best_k must be zero or greater")
        if self.max_gate_records_per_module_per_k < 0:
            raise ValueError("max_gate_records_per_module_per_k must be zero or greater")


@dataclass(frozen=True)
class ModuleExportResult:
    """Paths and counts from an automatic module export."""

    manifest_path: Path
    summary_path: Path
    n_modules: int
    n_gate_panels: int
    n_module_panels: int


def safe_module_id(*, value: str, prefix: str) -> str:
    """Return a stable file/module-safe identifier.

    Parameters
    ----------
    value : str
        Human-readable taxon or clade label.
    prefix : str
        Prefix identifying module type.

    Returns
    -------
    str
        Lower-case safe identifier.
    """
    cleaned = _SAFE_ID_RE.sub("_", str(value).strip().lower()).strip("_")
    if not cleaned:
        cleaned = "unclassified"
    return f"{prefix}_{cleaned}"


def infer_genus_from_record(*, record: dict[str, str]) -> str:
    """Infer the genus label represented by a panel record.

    Parameters
    ----------
    record : dict[str, str]
        Panel record.

    Returns
    -------
    str
        Genus-like label, or an empty string if none can be inferred.
    """
    evidence_rank = record.get("evidence_rank", "").strip().lower()
    evidence_name = record.get("evidence_name", "").strip()
    species_name = record.get("species_name", "").strip()
    clade = record.get("clade", "").strip()

    if evidence_rank == "genus" and evidence_name:
        return evidence_name
    if species_name and _SPECIES_NAME_RE.match(species_name):
        return species_name.split()[0]
    if clade and " " not in clade:
        return clade
    return ""


def read_panel_records(*, panel_path: str | Path) -> list[dict[str, str]]:
    """Read a KmerSutra panel into dictionaries.

    Parameters
    ----------
    panel_path : str or pathlib.Path
        Panel TSV or TSV.GZ path.

    Returns
    -------
    list[dict[str, str]]
        Panel records.

    Raises
    ------
    FileNotFoundError
        If the panel path is missing or empty.
    ValueError
        If required columns are missing.
    """
    path = Path(panel_path)
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"Panel file is missing or empty: {path}")
    with open_text(path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        missing = [column for column in PANEL_FIELDNAMES if column not in header]
        if missing:
            raise ValueError(f"Panel is missing required columns: {', '.join(missing)}")
        records: list[dict[str, str]] = []
        for line_number, line in enumerate(handle, start=2):
            if not line.strip():
                continue
            values = line.rstrip("\n").split("\t")
            if len(values) < len(header):
                values += [""] * (len(header) - len(values))
            elif len(values) > len(header):
                values = values[: len(header)]
            record = dict(zip(header, values))
            if not record.get("kmer") or not record.get("k"):
                raise ValueError(f"Malformed panel record at line {line_number}: missing kmer/k")
            records.append(record)
    if not records:
        raise ValueError(f"Panel contains no diagnostic records: {path}")
    return records


def _write_panel_subset(
    *,
    records: Iterable[dict[str, str]],
    output_path: Path,
) -> int:
    """Write a panel subset and return the number of records written."""
    n_records = 0
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open_text(output_path, "wt") as handle:
        handle.write("\t".join(PANEL_FIELDNAMES) + "\n")
        for record in records:
            handle.write(
                "\t".join(str(record.get(column, "")) for column in PANEL_FIELDNAMES)
                + "\n"
            )
            n_records += 1
    return n_records


def _limit_gate_records(
    *,
    records: list[dict[str, str]],
    max_per_module_per_k: int,
) -> list[dict[str, str]]:
    """Apply a deterministic per-k cap to gate records."""
    if max_per_module_per_k <= 0:
        return records
    retained: list[dict[str, str]] = []
    counts: dict[str, int] = defaultdict(int)
    for record in records:
        k_value = str(record.get("k", ""))
        if counts[k_value] >= max_per_module_per_k:
            continue
        counts[k_value] += 1
        retained.append(record)
    return retained


def _records_by_clade(
    *,
    records: Iterable[dict[str, str]],
) -> dict[str, list[dict[str, str]]]:
    """Group panel records by clade label."""
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for record in records:
        clade = record.get("clade", "").strip() or "unclassified"
        grouped[clade].append(record)
    return dict(grouped)


def _records_by_genus(
    *,
    records: Iterable[dict[str, str]],
) -> dict[str, list[dict[str, str]]]:
    """Group panel records by inferred genus label."""
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for record in records:
        genus = infer_genus_from_record(record=record)
        if genus:
            grouped[genus].append(record)
    return dict(grouped)


def export_hierarchical_modules_from_panel(
    *,
    panel_path: str | Path,
    out_dir: str | Path,
    config: ModuleExportConfig,
    logger: logging.Logger | None = None,
) -> ModuleExportResult:
    """Export gate/detail module panels and a manifest from one flat panel.

    Parameters
    ----------
    panel_path : str or pathlib.Path
        Existing flat KmerSutra panel.
    out_dir : str or pathlib.Path
        Output directory for module panels and manifest.
    config : ModuleExportConfig
        Module export settings.
    logger : logging.Logger or None, optional
        Logger.

    Returns
    -------
    ModuleExportResult
        Manifest and summary paths plus counts.
    """
    config.validate()
    records = read_panel_records(panel_path=panel_path)
    output_dir = Path(out_dir)
    gates_dir = output_dir / "gates"
    modules_dir = output_dir / "modules"
    manifest_path = output_dir / "kmersutra_module_manifest.tsv"
    summary_path = output_dir / "module_export_summary.tsv"
    output_dir.mkdir(parents=True, exist_ok=True)

    if logger:
        logger.info(
            "Exporting hierarchical modules from %s into %s (%d records)",
            panel_path,
            output_dir,
            len(records),
        )

    manifest_records: list[dict[str, object]] = []
    summary_records: list[dict[str, object]] = []
    n_gate_panels = 0
    n_module_panels = 0

    root_id = ""

    clade_parent_ids: dict[str, str] = {}
    for clade, clade_records in sorted(_records_by_clade(records=records).items()):
        broad_records = [
            record
            for record in clade_records
            if record.get("evidence_rank", "").strip().lower() in config.gate_ranks
        ]
        module_records = broad_records or clade_records
        module_id = safe_module_id(value=clade, prefix="clade")
        clade_parent_ids[clade] = module_id
        gate_path = gates_dir / f"{module_id}.gate.tsv.gz"
        module_path = modules_dir / f"{module_id}.module.tsv.gz"
        gate_records = _limit_gate_records(
            records=broad_records or clade_records,
            max_per_module_per_k=config.max_gate_records_per_module_per_k,
        )
        n_gate = _write_panel_subset(records=gate_records, output_path=gate_path)
        n_module = _write_panel_subset(records=module_records, output_path=module_path)
        n_gate_panels += 1
        n_module_panels += 1
        manifest_records.append(
            _manifest_record(
                module_id=module_id,
                module_name=clade,
                rank="clade",
                parent_module_id=root_id,
                gate_path=gate_path,
                module_path=module_path,
                manifest_dir=output_dir,
                config=config,
            )
        )
        summary_records.append(
            _summary_record(
                module_id=module_id,
                module_name=clade,
                rank="clade",
                parent_module_id=root_id,
                gate_records=n_gate,
                module_records=n_module,
                gate_path=gate_path,
                module_path=module_path,
                manifest_dir=output_dir,
            )
        )

    for genus, genus_records in sorted(_records_by_genus(records=records).items()):
        genus_gate_records = [
            record
            for record in genus_records
            if record.get("evidence_rank", "").strip().lower() in config.gate_ranks
            and record.get("evidence_name", "").strip() == genus
        ]
        if not genus_gate_records and config.allow_species_gate_fallback:
            genus_gate_records = [
                record
                for record in genus_records
                if record.get("evidence_rank", "").strip().lower() == "species"
            ]
        if not genus_gate_records:
            if logger:
                logger.warning(
                    "Skipping genus module %s because no gate records were available",
                    genus,
                )
            continue

        parent_clade = _dominant_clade(records=genus_records)
        parent_module_id = clade_parent_ids.get(parent_clade, root_id)
        module_id = safe_module_id(value=genus, prefix="genus")
        gate_path = gates_dir / f"{module_id}.gate.tsv.gz"
        module_path = modules_dir / f"{module_id}.module.tsv.gz"
        gate_records = _limit_gate_records(
            records=genus_gate_records,
            max_per_module_per_k=config.max_gate_records_per_module_per_k,
        )
        n_gate = _write_panel_subset(records=gate_records, output_path=gate_path)
        n_module = _write_panel_subset(records=genus_records, output_path=module_path)
        n_gate_panels += 1
        n_module_panels += 1
        manifest_records.append(
            _manifest_record(
                module_id=module_id,
                module_name=genus,
                rank="genus",
                parent_module_id=parent_module_id,
                gate_path=gate_path,
                module_path=module_path,
                manifest_dir=output_dir,
                config=config,
            )
        )
        summary_records.append(
            _summary_record(
                module_id=module_id,
                module_name=genus,
                rank="genus",
                parent_module_id=parent_module_id,
                gate_records=n_gate,
                module_records=n_module,
                gate_path=gate_path,
                module_path=module_path,
                manifest_dir=output_dir,
            )
        )

    write_tsv(
        records=manifest_records,
        output_path=manifest_path,
        fieldnames=MODULE_MANIFEST_FIELDNAMES,
    )
    write_tsv(
        records=summary_records,
        output_path=summary_path,
        fieldnames=MODULE_EXPORT_SUMMARY_FIELDNAMES,
    )
    if logger:
        logger.info(
            "Hierarchical module manifest written to %s with %d module(s)",
            manifest_path,
            len(manifest_records),
        )
    return ModuleExportResult(
        manifest_path=manifest_path,
        summary_path=summary_path,
        n_modules=len(manifest_records),
        n_gate_panels=n_gate_panels,
        n_module_panels=n_module_panels,
    )


def _dominant_clade(*, records: Iterable[dict[str, str]]) -> str:
    """Return the most frequent non-empty clade among records."""
    counts: dict[str, int] = defaultdict(int)
    for record in records:
        clade = record.get("clade", "").strip() or "unclassified"
        counts[clade] += 1
    if not counts:
        return "unclassified"
    return sorted(counts.items(), key=lambda item: (-item[1], item[0]))[0][0]


def _relative_path(*, path: Path, base_dir: Path) -> str:
    """Return a POSIX-style path relative to a manifest directory."""
    return path.relative_to(base_dir).as_posix()


def _manifest_record(
    *,
    module_id: str,
    module_name: str,
    rank: str,
    parent_module_id: str,
    gate_path: Path,
    module_path: Path,
    manifest_dir: Path,
    config: ModuleExportConfig,
) -> dict[str, object]:
    """Create one module manifest record."""
    return {
        "module_id": module_id,
        "module_name": module_name,
        "rank": rank,
        "parent_module_id": parent_module_id,
        "gate_panel_path": _relative_path(path=gate_path, base_dir=manifest_dir),
        "module_panel_path": _relative_path(path=module_path, base_dir=manifest_dir),
        "min_gate_unique_kmers": config.min_gate_unique_kmers,
        "min_gate_positive_sequences": config.min_gate_positive_sequences,
        "min_gate_k_values": config.min_gate_k_values,
        "min_gate_best_k": config.min_gate_best_k,
    }


def _summary_record(
    *,
    module_id: str,
    module_name: str,
    rank: str,
    parent_module_id: str,
    gate_records: int,
    module_records: int,
    gate_path: Path,
    module_path: Path,
    manifest_dir: Path,
) -> dict[str, object]:
    """Create one module export summary record."""
    return {
        "module_id": module_id,
        "module_name": module_name,
        "rank": rank,
        "parent_module_id": parent_module_id,
        "gate_records": gate_records,
        "module_records": module_records,
        "gate_panel_path": _relative_path(path=gate_path, base_dir=manifest_dir),
        "module_panel_path": _relative_path(path=module_path, base_dir=manifest_dir),
    }
