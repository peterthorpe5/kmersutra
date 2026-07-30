"""Reference-panel leakage auditing for independent benchmarks."""

from __future__ import annotations

import gzip
import logging
import re
from collections.abc import Iterable
from pathlib import Path

from kmersutra.mock_community import (
    MockCommunityTruth,
    load_truth_manifest,
    normalise_accession,
    normalise_species_key,
    sha256_file,
    truth_by_species,
)
from kmersutra.table_io import read_records_table, write_records_table

ACCESSION_PATTERN = re.compile(r"\b(?:GC[AF]_\d+|[A-Z]{1,4}_?[A-Z0-9]*\d+)(?:\.\d+)?\b")

AUDIT_FIELDNAMES = [
    "genome_fasta",
    "species_name",
    "assembly_accession",
    "normalised_assembly_accession",
    "matched_truth_species",
    "matched_truth_organism_id",
    "same_species_reference",
    "assembly_accession_leakage",
    "sequence_accession_leakage",
    "checksum_leakage",
    "matched_truth_accessions",
    "panel_fasta_sha256",
    "fasta_status",
    "audit_status",
    "audit_reason",
]


def open_text_auto(path: Path):
    """Open a plain or gzip-compressed text file.

    Args:
        path: Text file path.

    Returns:
        Readable text handle.
    """
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def fasta_header_accessions(path: Path) -> set[str]:
    """Collect normalised accession-like tokens from FASTA headers.

    Args:
        path: FASTA or FASTA.GZ file.

    Returns:
        Normalised accession tokens found in header lines.
    """
    accessions: set[str] = set()
    with open_text_auto(path) as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            for match in ACCESSION_PATTERN.findall(line):
                accession = normalise_accession(match)
                if accession:
                    accessions.add(accession)
    return accessions


def truth_accession_sets(
    truth_records: Iterable[MockCommunityTruth],
) -> tuple[set[str], set[str]]:
    """Collect truth assembly and sequence accession sets.

    Args:
        truth_records: Truth records.

    Returns:
        Assembly accession set and sequence accession set.
    """
    assembly_accessions: set[str] = set()
    sequence_accessions: set[str] = set()
    for record in truth_records:
        assembly_accessions.update(record.truth_assembly_accessions)
        sequence_accessions.update(record.truth_sequence_accessions)
    return assembly_accessions, sequence_accessions


def truth_checksums(
    *,
    truth_records: Iterable[MockCommunityTruth],
    logger: logging.Logger | None = None,
) -> dict[str, MockCommunityTruth]:
    """Resolve exact truth FASTA checksums.

    Args:
        truth_records: Truth records.
        logger: Optional logger.

    Returns:
        Mapping from SHA-256 digest to truth record.

    Raises:
        FileNotFoundError: If a declared truth FASTA is missing.
        ValueError: If a declared truth checksum does not match its FASTA.
    """
    checksums: dict[str, MockCommunityTruth] = {}
    for record in truth_records:
        declared = record.truth_fasta_sha256
        observed = ""
        if record.truth_fasta is not None:
            if not record.truth_fasta.is_file():
                raise FileNotFoundError(
                    f"Truth FASTA is missing for {record.species_name}: "
                    f"{record.truth_fasta}"
                )
            observed = sha256_file(record.truth_fasta)
            if declared and observed != declared:
                raise ValueError(
                    f"Truth FASTA checksum mismatch for {record.species_name}: "
                    f"declared={declared}, observed={observed}"
                )
        checksum = declared or observed
        if checksum:
            checksums[checksum] = record
            if logger:
                logger.debug(
                    "Registered truth checksum for %s: %s",
                    record.species_name,
                    checksum,
                )
    return checksums


def _resolve_fasta_path(*, genome_config: Path, value: object) -> Path | None:
    """Resolve a genome FASTA path relative to its configuration."""
    text = str(value or "").strip()
    if not text:
        return None
    path = Path(text).expanduser()
    if not path.is_absolute():
        path = genome_config.parent / path
    return path.resolve()


def audit_reference_panel(
    *,
    genome_config: str | Path,
    truth_manifest: str | Path,
    allow_missing_fastas: bool = False,
    logger: logging.Logger | None = None,
) -> list[dict[str, object]]:
    """Audit a KmerSutra genome configuration for benchmark leakage.

    Args:
        genome_config: KmerSutra genome configuration table.
        truth_manifest: Mock-community truth manifest.
        allow_missing_fastas: Keep metadata-only rows when FASTAs are missing.
        logger: Optional logger.

    Returns:
        One audit row per panel genome.

    Raises:
        ValueError: If the genome configuration is empty or invalid.
        FileNotFoundError: If a required panel FASTA is missing.
    """
    config_path = Path(genome_config).expanduser().resolve()
    rows = read_records_table(
        input_path=config_path,
        required_columns=["genome_fasta", "species_name", "assembly_accession"],
        logger=logger,
    )
    if not rows:
        raise ValueError(f"Genome configuration contains no rows: {config_path}")
    truth_records = load_truth_manifest(
        manifest_path=truth_manifest,
        logger=logger,
    )
    species_mapping = truth_by_species(truth_records=truth_records)
    assembly_truth, sequence_truth = truth_accession_sets(truth_records)
    checksum_truth = truth_checksums(
        truth_records=truth_records,
        logger=logger,
    )

    audit_rows: list[dict[str, object]] = []
    for row in rows:
        species_name = str(row.get("species_name", "")).strip()
        truth = species_mapping.get(normalise_species_key(species_name))
        assembly_accession = str(row.get("assembly_accession", "")).strip()
        normalised_assembly = normalise_accession(assembly_accession)
        assembly_leakage = bool(
            normalised_assembly and normalised_assembly in assembly_truth
        )
        fasta_path = _resolve_fasta_path(
            genome_config=config_path,
            value=row.get("genome_fasta", ""),
        )
        fasta_status = "not_declared"
        panel_checksum = ""
        matching_sequence_accessions: set[str] = set()
        checksum_leakage = False
        if fasta_path is not None:
            if not fasta_path.is_file() or fasta_path.stat().st_size <= 0:
                fasta_status = "missing_or_empty"
                if not allow_missing_fastas:
                    raise FileNotFoundError(
                        f"Panel FASTA is missing or empty: {fasta_path}"
                    )
            else:
                fasta_status = "present"
                header_accessions = fasta_header_accessions(fasta_path)
                matching_sequence_accessions = header_accessions.intersection(
                    sequence_truth
                )
                panel_checksum = sha256_file(fasta_path)
                checksum_leakage = panel_checksum in checksum_truth

        sequence_leakage = bool(matching_sequence_accessions)
        leakage_reasons: list[str] = []
        if assembly_leakage:
            leakage_reasons.append("truth_assembly_accession")
        if sequence_leakage:
            leakage_reasons.append("truth_sequence_accession")
        if checksum_leakage:
            leakage_reasons.append("truth_fasta_checksum")
        if leakage_reasons:
            status = "leakage"
            reason = ";".join(leakage_reasons)
        elif fasta_status == "missing_or_empty":
            status = "incomplete"
            reason = "panel_fasta_missing_or_empty"
        elif truth is not None:
            status = "same_species_reference"
            reason = "held_out_same_species_reference"
        else:
            status = "clear"
            reason = "no_exact_truth_identity"

        audit_rows.append(
            {
                "genome_fasta": str(fasta_path or ""),
                "species_name": species_name,
                "assembly_accession": assembly_accession,
                "normalised_assembly_accession": normalised_assembly,
                "matched_truth_species": truth.species_name if truth else "",
                "matched_truth_organism_id": truth.organism_id if truth else "",
                "same_species_reference": str(bool(truth)),
                "assembly_accession_leakage": str(assembly_leakage),
                "sequence_accession_leakage": str(sequence_leakage),
                "checksum_leakage": str(checksum_leakage),
                "matched_truth_accessions": ";".join(
                    sorted(matching_sequence_accessions)
                ),
                "panel_fasta_sha256": panel_checksum,
                "fasta_status": fasta_status,
                "audit_status": status,
                "audit_reason": reason,
            }
        )
    if logger:
        leakage_count = sum(
            row["audit_status"] == "leakage" for row in audit_rows
        )
        logger.info(
            "Reference audit completed: %d rows, %d leakage rows",
            len(audit_rows),
            leakage_count,
        )
    return audit_rows


def write_reference_panel_audit(
    *,
    genome_config: str | Path,
    truth_manifest: str | Path,
    output_table: str | Path,
    allow_missing_fastas: bool = False,
    fail_on_leakage: bool = False,
    logger: logging.Logger | None = None,
) -> list[dict[str, object]]:
    """Run and write a reference-panel leakage audit.

    Args:
        genome_config: KmerSutra genome configuration table.
        truth_manifest: Mock-community truth manifest.
        output_table: Audit output table.
        allow_missing_fastas: Permit metadata-only audit rows.
        fail_on_leakage: Raise when exact truth identity is found.
        logger: Optional logger.

    Returns:
        Written audit rows.

    Raises:
        ValueError: If ``fail_on_leakage`` is true and leakage is detected.
    """
    rows = audit_reference_panel(
        genome_config=genome_config,
        truth_manifest=truth_manifest,
        allow_missing_fastas=allow_missing_fastas,
        logger=logger,
    )
    write_records_table(
        records=rows,
        output_path=output_table,
        fieldnames=AUDIT_FIELDNAMES,
        logger=logger,
    )
    leaking = [
        row for row in rows if str(row.get("audit_status", "")) == "leakage"
    ]
    if fail_on_leakage and leaking:
        identifiers = [
            str(row.get("assembly_accession") or row.get("genome_fasta"))
            for row in leaking[:10]
        ]
        raise ValueError(
            "Reference-panel leakage detected in "
            f"{len(leaking)} row(s): {', '.join(identifiers)}"
        )
    return rows
