"""Prepare and validate leakage-controlled ATCC reference panels.

The public ATCC MSA-1003 benchmark requires reference genomes for all expected
species, but must not use the exact sequences listed in the benchmark truth
manifest.  This module creates a reproducible NCBI download plan, combines the
downloaded held-out genomes with an established background genome collection,
and evaluates the reference audit without relying on shell-specific AWK code.
"""

from __future__ import annotations

import logging
import re
from collections import Counter
from collections.abc import Iterable, Mapping
from pathlib import Path

from kmersutra.config import VALID_ROLES
from kmersutra.mock_community import (
    MockCommunityTruth,
    load_truth_manifest,
    normalise_accession,
    normalise_species_key,
)
from kmersutra.table_io import read_records_table, write_records_table
from kmersutra.taxonomy import TaxonomyDatabase

LOGGER = logging.getLogger("kmersutra.atcc_reference")

TAXID_PLAN_FIELDS = [
    "taxid",
    "role",
    "clade",
    "group_label",
    "max_assemblies",
    "best_per_species",
    "min_total_length",
    "max_total_length",
    "min_scaffold_n50",
    "min_contig_n50",
]

FINAL_CONFIG_FIELDS = [
    "genome_fasta",
    "species_name",
    "strain_name",
    "taxid",
    "assembly_accession",
    "role",
    "clade",
    "source",
    "query_taxid",
    "assembly_level",
    "scaffold_n50",
    "contig_n50",
]

COVERAGE_FIELDS = [
    "species_name",
    "ncbi_taxid",
    "reference_count",
    "panel_status",
]

GATE_SUMMARY_FIELDS = ["metric", "value"]


def _expected_truth_records(
    *, truth_manifest: str | Path
) -> list[MockCommunityTruth]:
    """Load expected truth records with usable NCBI taxids.

    Args:
        truth_manifest: Mock-community truth manifest.

    Returns:
        Expected truth records.

    Raises:
        ValueError: If an expected organism has no NCBI taxid.
    """
    records = [
        record
        for record in load_truth_manifest(manifest_path=truth_manifest)
        if record.expected
    ]
    missing_taxids = [record.species_name for record in records if not record.ncbi_taxid]
    if missing_taxids:
        raise ValueError(
            "Expected truth organisms lack NCBI taxids: "
            + ", ".join(sorted(missing_taxids))
        )
    return records


def _safe_label(value: str) -> str:
    """Return a stable label containing only portable filename characters.

    Args:
        value: Free-text label.

    Returns:
        Sanitised label.
    """
    label = re.sub(r"[^A-Za-z0-9._-]+", "_", value.strip())
    return re.sub(r"_+", "_", label).strip("_") or "unknown"


def _genus_node(
    *, taxonomy: TaxonomyDatabase, taxid: str, species_name: str
):
    """Resolve the genus node for one species taxid.

    Args:
        taxonomy: Loaded NCBI taxonomy database.
        taxid: Species NCBI taxid.
        species_name: Species name used in diagnostics.

    Returns:
        The genus-level taxonomy node.

    Raises:
        ValueError: If the taxid is unavailable or has no genus ancestor.
    """
    normalised_taxid = taxonomy.normalise_taxid(taxid)
    if not normalised_taxid or taxonomy.get_node(normalised_taxid) is None:
        raise ValueError(
            f"Truth taxid is absent from the taxonomy database: "
            f"{species_name} ({taxid})"
        )
    for node in reversed(taxonomy.get_ranked_lineage(normalised_taxid)):
        if node.rank == "genus":
            return node
    raise ValueError(
        f"No genus ancestor was found for {species_name} ({normalised_taxid})"
    )


def create_taxid_plan(
    *,
    truth_manifest: str | Path,
    taxonomy: TaxonomyDatabase,
    output_path: str | Path,
    target_assemblies: int = 5,
    near_neighbour_assemblies_per_genus: int = 25,
    logger: logging.Logger | None = None,
) -> list[dict[str, object]]:
    """Create target-species and genus-neighbour NCBI download plans.

    Each expected species is queried directly to obtain held-out same-species
    references.  Each distinct parent genus is also queried to provide absent
    near-neighbour species.  Target-species rows are later removed from the
    genus-derived neighbour set to avoid duplicate or conflicting roles.

    Args:
        truth_manifest: ATCC truth manifest.
        taxonomy: Loaded NCBI taxonomy database.
        output_path: Destination TSV path.
        target_assemblies: Maximum alternative assemblies requested per target.
        near_neighbour_assemblies_per_genus: Maximum assemblies requested per
            distinct target genus.
        logger: Optional logger.

    Returns:
        Written taxid-plan records.

    Raises:
        ValueError: If an assembly limit is not positive.
    """
    if target_assemblies <= 0:
        raise ValueError("target_assemblies must be positive")
    if near_neighbour_assemblies_per_genus <= 0:
        raise ValueError("near_neighbour_assemblies_per_genus must be positive")

    truth_records = _expected_truth_records(truth_manifest=truth_manifest)
    rows: list[dict[str, object]] = []
    genera: dict[str, str] = {}
    for record in sorted(truth_records, key=lambda item: item.species_name):
        species_taxid = taxonomy.normalise_taxid(record.ncbi_taxid)
        genus = _genus_node(
            taxonomy=taxonomy,
            taxid=species_taxid,
            species_name=record.species_name,
        )
        genera[genus.taxid] = genus.name
        rows.append(
            {
                "taxid": species_taxid,
                "role": "target_species",
                "clade": "ATCC_MSA1003",
                "group_label": f"target_{_safe_label(record.species_name)}",
                "max_assemblies": target_assemblies,
                "best_per_species": target_assemblies,
                "min_total_length": "",
                "max_total_length": "",
                "min_scaffold_n50": "",
                "min_contig_n50": "",
            }
        )

    for genus_taxid, genus_name in sorted(
        genera.items(), key=lambda item: (item[1], item[0])
    ):
        rows.append(
            {
                "taxid": genus_taxid,
                "role": "near_neighbour",
                "clade": "ATCC_target_genera",
                "group_label": f"genus_{_safe_label(genus_name)}",
                "max_assemblies": near_neighbour_assemblies_per_genus,
                "best_per_species": 1,
                "min_total_length": "",
                "max_total_length": "",
                "min_scaffold_n50": "",
                "min_contig_n50": "",
            }
        )

    write_records_table(
        records=rows,
        output_path=output_path,
        fieldnames=TAXID_PLAN_FIELDS,
        logger=logger,
    )
    if logger:
        logger.info(
            "Created ATCC taxid plan with %d target species and %d genera",
            len(truth_records),
            len(genera),
        )
    return rows


def _resolve_fasta_path(*, config_path: Path, value: object) -> Path:
    """Resolve and validate a FASTA path from a genome configuration.

    Args:
        config_path: Source genome-configuration path.
        value: FASTA path value.

    Returns:
        Absolute FASTA path.

    Raises:
        FileNotFoundError: If the FASTA is missing or empty.
    """
    text = str(value or "").strip()
    if not text:
        raise FileNotFoundError(f"Blank genome_fasta in {config_path}")
    fasta_path = Path(text).expanduser()
    if not fasta_path.is_absolute():
        fasta_path = config_path.parent / fasta_path
    fasta_path = fasta_path.resolve()
    if not fasta_path.is_file() or fasta_path.stat().st_size <= 0:
        raise FileNotFoundError(f"Genome FASTA is missing or empty: {fasta_path}")
    return fasta_path


def _normalise_config_row(
    *, row: Mapping[str, object], config_path: Path, source_label: str
) -> dict[str, str]:
    """Normalise one genome-configuration row.

    Args:
        row: Source row.
        config_path: Source configuration path.
        source_label: Provenance label used when the row lacks ``source``.

    Returns:
        Normalised row using :data:`FINAL_CONFIG_FIELDS`.

    Raises:
        ValueError: If the row uses an unsupported role.
    """
    role = str(row.get("role", "")).strip() or "downloaded"
    if role not in VALID_ROLES:
        raise ValueError(f"Unsupported genome role in {config_path}: {role}")
    normalised = {
        field: str(row.get(field, "") or "").strip()
        for field in FINAL_CONFIG_FIELDS
    }
    normalised["genome_fasta"] = str(
        _resolve_fasta_path(config_path=config_path, value=row.get("genome_fasta", ""))
    )
    normalised["role"] = role
    normalised["source"] = normalised["source"] or source_label
    return normalised


def _truth_lookup(
    *, truth_records: Iterable[MockCommunityTruth], taxonomy: TaxonomyDatabase
) -> tuple[dict[str, MockCommunityTruth], dict[str, MockCommunityTruth]]:
    """Create taxid and accepted-name truth lookup tables.

    Args:
        truth_records: Expected truth records.
        taxonomy: Loaded NCBI taxonomy database.

    Returns:
        Taxid and normalised-name lookup dictionaries.
    """
    by_taxid: dict[str, MockCommunityTruth] = {}
    by_name: dict[str, MockCommunityTruth] = {}
    for record in truth_records:
        by_taxid[taxonomy.normalise_taxid(record.ncbi_taxid)] = record
        for accepted_name in record.accepted_species_names:
            by_name[accepted_name] = record
    return by_taxid, by_name


def _matching_truth(
    *,
    row: Mapping[str, object],
    taxonomy: TaxonomyDatabase,
    truth_by_taxid: Mapping[str, MockCommunityTruth],
    truth_by_name: Mapping[str, MockCommunityTruth],
) -> MockCommunityTruth | None:
    """Return the expected truth species represented by a config row.

    Args:
        row: Genome-configuration row.
        taxonomy: Loaded NCBI taxonomy database.
        truth_by_taxid: Truth records keyed by current taxid.
        truth_by_name: Truth records keyed by accepted species name.

    Returns:
        Matching truth record, or ``None``.
    """
    row_taxid = taxonomy.normalise_taxid(str(row.get("taxid", "")))
    if row_taxid and row_taxid in truth_by_taxid:
        return truth_by_taxid[row_taxid]
    species_key = normalise_species_key(row.get("species_name", ""))
    return truth_by_name.get(species_key)


def _deduplication_key(row: Mapping[str, object]) -> tuple[str, str]:
    """Return an accession-first stable genome deduplication key.

    Args:
        row: Normalised genome-configuration row.

    Returns:
        Key type and value.
    """
    accession = normalise_accession(row.get("assembly_accession", ""))
    if accession:
        return "accession", accession
    return "fasta", str(row.get("genome_fasta", ""))


def finalise_reference_config(
    *,
    downloaded_config: str | Path,
    background_config: str | Path,
    truth_manifest: str | Path,
    taxonomy: TaxonomyDatabase,
    output_config: str | Path,
    coverage_table: str | Path,
    minimum_target_references: int = 1,
    logger: logging.Logger | None = None,
) -> tuple[list[dict[str, str]], list[dict[str, object]]]:
    """Build the final held-out ATCC-plus-background genome configuration.

    Target species downloaded through a genus query are discarded from the
    near-neighbour set.  Expected species are also removed from the background
    collection so only the explicitly downloaded, accession-excluded copies can
    represent ATCC targets.  Duplicate assemblies are resolved by precedence:
    target, near neighbour, then background.

    Args:
        downloaded_config: Downloader-produced genome configuration.
        background_config: Established v4 genome configuration.
        truth_manifest: ATCC truth manifest.
        taxonomy: Loaded NCBI taxonomy database.
        output_config: Destination final genome configuration.
        coverage_table: Destination target-coverage table.
        minimum_target_references: Minimum non-leaking references per target.
        logger: Optional logger.

    Returns:
        Final genome rows and coverage rows.

    Raises:
        ValueError: If target coverage is incomplete or role assignments are
            inconsistent.
    """
    if minimum_target_references <= 0:
        raise ValueError("minimum_target_references must be positive")

    truth_records = _expected_truth_records(truth_manifest=truth_manifest)
    truth_by_taxid, truth_by_name = _truth_lookup(
        truth_records=truth_records,
        taxonomy=taxonomy,
    )
    downloaded_path = Path(downloaded_config).expanduser().resolve()
    background_path = Path(background_config).expanduser().resolve()
    downloaded_rows = read_records_table(
        input_path=downloaded_path,
        required_columns=[
            "genome_fasta",
            "species_name",
            "taxid",
            "assembly_accession",
            "role",
        ],
        logger=logger,
    )
    background_rows = read_records_table(
        input_path=background_path,
        required_columns=[
            "genome_fasta",
            "species_name",
            "taxid",
            "assembly_accession",
            "role",
        ],
        logger=logger,
    )

    candidates: list[tuple[int, dict[str, str]]] = []
    for source_row in downloaded_rows:
        row = _normalise_config_row(
            row=source_row,
            config_path=downloaded_path,
            source_label="NCBI_ATCC_heldout",
        )
        truth = _matching_truth(
            row=row,
            taxonomy=taxonomy,
            truth_by_taxid=truth_by_taxid,
            truth_by_name=truth_by_name,
        )
        if row["role"] == "target_species":
            if truth is None:
                raise ValueError(
                    "A target_species download row is not an expected ATCC species: "
                    f"{row['species_name']} ({row['taxid']})"
                )
            candidates.append((0, row))
        elif row["role"] in {"near_neighbour", "near_neighbor"}:
            if truth is not None:
                continue
            row["role"] = "near_neighbour"
            candidates.append((1, row))
        else:
            raise ValueError(
                "Downloaded ATCC rows must use target_species or near_neighbour "
                f"roles, observed {row['role']!r}"
            )

    excluded_background_targets = 0
    for source_row in background_rows:
        if str(source_row.get("role", "")).strip() == "exclude":
            continue
        row = _normalise_config_row(
            row=source_row,
            config_path=background_path,
            source_label="established_v4_background",
        )
        truth = _matching_truth(
            row=row,
            taxonomy=taxonomy,
            truth_by_taxid=truth_by_taxid,
            truth_by_name=truth_by_name,
        )
        if truth is not None:
            excluded_background_targets += 1
            continue
        candidates.append((2, row))

    selected: dict[tuple[str, str], tuple[int, dict[str, str]]] = {}
    for priority, row in sorted(
        candidates,
        key=lambda item: (
            item[0],
            item[1]["species_name"],
            item[1]["assembly_accession"],
            item[1]["genome_fasta"],
        ),
    ):
        key = _deduplication_key(row)
        if key not in selected:
            selected[key] = (priority, row)

    final_rows = [
        item[1]
        for item in sorted(
            selected.values(),
            key=lambda item: (
                item[0],
                item[1]["species_name"],
                item[1]["assembly_accession"],
                item[1]["genome_fasta"],
            ),
        )
    ]
    reference_counts: Counter[str] = Counter()
    for row in final_rows:
        if row["role"] != "target_species":
            continue
        truth = _matching_truth(
            row=row,
            taxonomy=taxonomy,
            truth_by_taxid=truth_by_taxid,
            truth_by_name=truth_by_name,
        )
        if truth is not None:
            reference_counts[truth.species_name] += 1

    coverage_rows: list[dict[str, object]] = []
    missing: list[str] = []
    for truth in sorted(truth_records, key=lambda item: item.species_name):
        count = reference_counts[truth.species_name]
        status = "represented" if count >= minimum_target_references else "missing"
        coverage_rows.append(
            {
                "species_name": truth.species_name,
                "ncbi_taxid": taxonomy.normalise_taxid(truth.ncbi_taxid),
                "reference_count": count,
                "panel_status": status,
            }
        )
        if status == "missing":
            missing.append(truth.species_name)

    write_records_table(
        records=final_rows,
        output_path=output_config,
        fieldnames=FINAL_CONFIG_FIELDS,
        logger=logger,
    )
    write_records_table(
        records=coverage_rows,
        output_path=coverage_table,
        fieldnames=COVERAGE_FIELDS,
        logger=logger,
    )
    if logger:
        logger.info(
            "Final configuration contains %d genomes; removed %d ATCC-like "
            "background row(s)",
            len(final_rows),
            excluded_background_targets,
        )
    if missing:
        raise ValueError(
            "Held-out reference coverage is incomplete: " + ", ".join(missing)
        )
    return final_rows, coverage_rows


def evaluate_reference_gate(
    *,
    audit_table: str | Path,
    truth_manifest: str | Path,
    coverage_table: str | Path,
    summary_table: str | Path,
    minimum_target_references: int = 1,
    logger: logging.Logger | None = None,
) -> bool:
    """Evaluate the final ATCC reference gate from the leakage audit.

    Args:
        audit_table: Output from ``kmersutra-audit-reference-panel``.
        truth_manifest: ATCC truth manifest.
        coverage_table: Destination per-species coverage table.
        summary_table: Destination gate summary table.
        minimum_target_references: Minimum references required per target.
        logger: Optional logger.

    Returns:
        ``True`` when the reference gate passes.

    Raises:
        ValueError: If ``minimum_target_references`` is not positive.
    """
    if minimum_target_references <= 0:
        raise ValueError("minimum_target_references must be positive")
    truth_records = _expected_truth_records(truth_manifest=truth_manifest)
    audit_rows = read_records_table(
        input_path=audit_table,
        required_columns=[
            "matched_truth_species",
            "same_species_reference",
            "fasta_status",
            "audit_status",
        ],
        logger=logger,
    )
    canonical_names = {record.species_name for record in truth_records}
    counts: Counter[str] = Counter(
        str(row.get("matched_truth_species", "")).strip()
        for row in audit_rows
        if str(row.get("matched_truth_species", "")).strip() in canonical_names
        and str(row.get("audit_status", "")).strip() == "same_species_reference"
    )
    coverage_rows: list[dict[str, object]] = []
    for record in sorted(truth_records, key=lambda item: item.species_name):
        count = counts[record.species_name]
        coverage_rows.append(
            {
                "species_name": record.species_name,
                "ncbi_taxid": record.ncbi_taxid,
                "reference_count": count,
                "panel_status": (
                    "represented"
                    if count >= minimum_target_references
                    else "missing"
                ),
            }
        )

    leakage_count = sum(
        str(row.get("audit_status", "")).strip() == "leakage"
        for row in audit_rows
    )
    incomplete_count = sum(
        str(row.get("fasta_status", "")).strip() != "present"
        for row in audit_rows
    )
    represented_count = sum(
        row["panel_status"] == "represented" for row in coverage_rows
    )
    expected_count = len(truth_records)
    missing_count = expected_count - represented_count
    minimum_observed = min(
        (int(row["reference_count"]) for row in coverage_rows),
        default=0,
    )
    passed = (
        represented_count == expected_count
        and missing_count == 0
        and leakage_count == 0
        and incomplete_count == 0
        and minimum_observed >= minimum_target_references
    )
    summary_rows = [
        {"metric": "expected_species", "value": expected_count},
        {"metric": "represented_species", "value": represented_count},
        {"metric": "missing_species", "value": missing_count},
        {"metric": "leakage_rows", "value": leakage_count},
        {"metric": "incomplete_fasta_rows", "value": incomplete_count},
        {
            "metric": "minimum_references_per_species",
            "value": minimum_observed,
        },
        {
            "metric": "required_references_per_species",
            "value": minimum_target_references,
        },
        {"metric": "gate_status", "value": "PASS" if passed else "BLOCKED"},
    ]
    write_records_table(
        records=coverage_rows,
        output_path=coverage_table,
        fieldnames=COVERAGE_FIELDS,
        logger=logger,
    )
    write_records_table(
        records=summary_rows,
        output_path=summary_table,
        fieldnames=GATE_SUMMARY_FIELDS,
        logger=logger,
    )
    if logger:
        logger.info(
            "ATCC reference gate: %s (%d/%d represented, %d leakage, %d incomplete)",
            "PASS" if passed else "BLOCKED",
            represented_count,
            expected_count,
            leakage_count,
            incomplete_count,
        )
    return passed
