"""Manifest-driven truth labelling for mock-community benchmarks."""

from __future__ import annotations

import hashlib
import logging
from collections import Counter
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

from kmersutra.ai_calibration import build_call_feature_record
from kmersutra.table_io import read_records_table, write_records_table
from kmersutra.zymo_truth import (
    first_non_empty,
    infer_species_name,
    row_has_evidence,
    row_is_reportable,
)

TRUTH_REQUIRED_COLUMNS = {
    "benchmark_id",
    "species_name",
    "expected_abundance_fraction",
    "abundance_tier",
}


@dataclass(frozen=True)
class MockCommunityTruth:
    """One expected or deliberately absent organism in a benchmark.

    Attributes:
        benchmark_id: Stable benchmark identifier.
        organism_id: Stable organism identifier within the benchmark.
        species_name: Canonical species name used in outputs.
        accepted_species_names: Normalised canonical name and accepted aliases.
        expected: Whether the organism is expected to be present.
        expected_abundance_fraction: Expected fraction between zero and one.
        abundance_tier: Pre-specified abundance tier label.
        ncbi_taxid: Optional NCBI taxonomy identifier.
        truth_strain: Reference-material strain label.
        truth_assembly_accessions: Assembly accessions that must be excluded.
        truth_sequence_accessions: Sequence accessions that must be excluded.
        truth_fasta: Optional local truth FASTA used for exact checksum auditing.
        truth_fasta_sha256: Optional expected SHA-256 for ``truth_fasta``.
    """

    benchmark_id: str
    organism_id: str
    species_name: str
    accepted_species_names: frozenset[str]
    expected: bool
    expected_abundance_fraction: float
    abundance_tier: str
    ncbi_taxid: str
    truth_strain: str
    truth_assembly_accessions: frozenset[str]
    truth_sequence_accessions: frozenset[str]
    truth_fasta: Path | None
    truth_fasta_sha256: str


def normalise_species_key(value: object) -> str:
    """Normalise a species label for robust exact matching.

    Args:
        value: Species-like value.

    Returns:
        Case-folded, whitespace-normalised species label.
    """
    return " ".join(str(value or "").replace("_", " ").strip().casefold().split())


def parse_bool(value: object, *, default: bool = True) -> bool:
    """Parse a conservative boolean value.

    Args:
        value: Boolean-like value.
        default: Value returned for blank input.

    Returns:
        Parsed boolean.

    Raises:
        ValueError: If a non-blank value is not recognised.
    """
    text = str(value or "").strip().casefold()
    if not text:
        return default
    if text in {"1", "true", "yes", "y", "expected", "present"}:
        return True
    if text in {"0", "false", "no", "n", "absent", "negative"}:
        return False
    raise ValueError(f"Unsupported boolean value: {value!r}")


def parse_fraction(value: object, *, field_name: str) -> float:
    """Parse a fraction bounded between zero and one.

    Args:
        value: Numeric fraction.
        field_name: Field name used in diagnostics.

    Returns:
        Parsed floating-point fraction.

    Raises:
        ValueError: If the value is missing, non-numeric or outside [0, 1].
    """
    try:
        fraction = float(str(value).strip())
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{field_name} must be numeric: {value!r}") from exc
    if not 0.0 <= fraction <= 1.0:
        raise ValueError(f"{field_name} must be between zero and one: {fraction}")
    return fraction


def split_values(value: object) -> frozenset[str]:
    """Split semicolon-delimited manifest values.

    Args:
        value: Semicolon-delimited text.

    Returns:
        Non-empty stripped values.
    """
    return frozenset(
        item.strip() for item in str(value or "").split(";") if item.strip()
    )


def normalise_accession(value: object) -> str:
    """Return an upper-case accession without a version suffix.

    Args:
        value: Assembly or sequence accession.

    Returns:
        Normalised accession.
    """
    return str(value or "").strip().upper().split(".", 1)[0]


def sha256_file(path: Path, *, chunk_size: int = 1024 * 1024) -> str:
    """Calculate a file SHA-256 digest.

    Args:
        path: File to hash.
        chunk_size: Read size in bytes.

    Returns:
        Lower-case hexadecimal SHA-256 digest.

    Raises:
        ValueError: If ``chunk_size`` is not positive.
        FileNotFoundError: If the file does not exist.
    """
    if chunk_size <= 0:
        raise ValueError("chunk_size must be positive")
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def _resolve_truth_fasta(*, manifest_path: Path, value: object) -> Path | None:
    """Resolve an optional truth FASTA relative to its manifest."""
    text = str(value or "").strip()
    if not text:
        return None
    path = Path(text).expanduser()
    if not path.is_absolute():
        path = manifest_path.parent / path
    return path.resolve()


def load_truth_manifest(
    *,
    manifest_path: str | Path,
    require_expected: bool = True,
    logger: logging.Logger | None = None,
) -> list[MockCommunityTruth]:
    """Load and validate a mock-community truth manifest.

    Args:
        manifest_path: TSV, TSV.GZ or Parquet truth manifest.
        require_expected: Require at least one expected organism.
        logger: Optional logger.

    Returns:
        Validated truth records.

    Raises:
        ValueError: If the manifest is empty or internally inconsistent.
    """
    path = Path(manifest_path).expanduser().resolve()
    rows = read_records_table(input_path=path, logger=logger)
    if not rows:
        raise ValueError(f"Truth manifest contains no records: {path}")
    missing = TRUTH_REQUIRED_COLUMNS.difference(rows[0])
    if missing:
        raise ValueError(
            "Truth manifest is missing required columns: "
            + ", ".join(sorted(missing))
        )

    records: list[MockCommunityTruth] = []
    seen_organisms: set[tuple[str, str]] = set()
    for line_number, row in enumerate(rows, start=2):
        if not any(str(value or "").strip() for value in row.values()):
            continue
        benchmark_id = str(row.get("benchmark_id", "")).strip()
        species_name = str(row.get("species_name", "")).strip()
        organism_id = str(row.get("organism_id", "")).strip() or species_name
        abundance_tier = str(row.get("abundance_tier", "")).strip()
        if not benchmark_id:
            raise ValueError(f"Blank benchmark_id at manifest row {line_number}")
        if not species_name:
            raise ValueError(f"Blank species_name at manifest row {line_number}")
        if not abundance_tier:
            raise ValueError(f"Blank abundance_tier at manifest row {line_number}")
        key = (benchmark_id, organism_id)
        if key in seen_organisms:
            raise ValueError(f"Duplicate truth organism {key!r}")
        seen_organisms.add(key)

        accepted_names = {
            normalise_species_key(species_name),
            *(
                normalise_species_key(name)
                for name in split_values(row.get("accepted_species_names", ""))
            ),
        }
        accepted_names.discard("")
        truth_fasta = _resolve_truth_fasta(
            manifest_path=path,
            value=row.get("truth_fasta", ""),
        )
        truth_fasta_sha256 = str(row.get("truth_fasta_sha256", "")).strip().lower()
        if truth_fasta_sha256 and len(truth_fasta_sha256) != 64:
            raise ValueError(
                f"truth_fasta_sha256 must contain 64 hexadecimal characters "
                f"at row {line_number}"
            )
        records.append(
            MockCommunityTruth(
                benchmark_id=benchmark_id,
                organism_id=organism_id,
                species_name=species_name,
                accepted_species_names=frozenset(accepted_names),
                expected=parse_bool(row.get("expected", ""), default=True),
                expected_abundance_fraction=parse_fraction(
                    row.get("expected_abundance_fraction", ""),
                    field_name="expected_abundance_fraction",
                ),
                abundance_tier=abundance_tier,
                ncbi_taxid=str(row.get("ncbi_taxid", "")).strip(),
                truth_strain=str(row.get("truth_strain", "")).strip(),
                truth_assembly_accessions=frozenset(
                    normalise_accession(value)
                    for value in split_values(
                        row.get("truth_assembly_accessions", "")
                    )
                    if normalise_accession(value)
                ),
                truth_sequence_accessions=frozenset(
                    normalise_accession(value)
                    for value in split_values(
                        row.get("truth_sequence_accessions", "")
                    )
                    if normalise_accession(value)
                ),
                truth_fasta=truth_fasta,
                truth_fasta_sha256=truth_fasta_sha256,
            )
        )

    expected = [record for record in records if record.expected]
    if require_expected and not expected:
        raise ValueError("Truth manifest has no expected organisms")
    total_fraction = sum(record.expected_abundance_fraction for record in expected)
    if total_fraction > 1.000001:
        raise ValueError(
            "Expected abundance fractions exceed one: "
            f"{total_fraction:.8f}"
        )
    benchmark_ids = {record.benchmark_id for record in records}
    if len(benchmark_ids) != 1:
        raise ValueError(
            "One truth manifest must describe exactly one benchmark_id; found "
            + ", ".join(sorted(benchmark_ids))
        )
    if logger:
        logger.info(
            "Loaded %d truth records (%d expected) for %s",
            len(records),
            len(expected),
            next(iter(benchmark_ids)),
        )
    return records


def truth_by_species(
    *,
    truth_records: Iterable[MockCommunityTruth],
) -> dict[str, MockCommunityTruth]:
    """Map canonical and accepted species labels to truth records.

    Args:
        truth_records: Validated truth records.

    Returns:
        Mapping from normalised species label to truth record.

    Raises:
        ValueError: If the same accepted name maps to multiple truth records.
    """
    mapping: dict[str, MockCommunityTruth] = {}
    for record in truth_records:
        for species_key in record.accepted_species_names:
            previous = mapping.get(species_key)
            if previous is not None and previous.organism_id != record.organism_id:
                raise ValueError(
                    f"Accepted species name {species_key!r} maps to multiple organisms"
                )
            mapping[species_key] = record
    return mapping


def classify_mock_truth_category(
    *,
    record: dict[str, object],
    truth_mapping: dict[str, MockCommunityTruth],
) -> dict[str, object]:
    """Classify one KmerSutra call row against mock-community truth.

    Args:
        record: KmerSutra species-call row.
        truth_mapping: Mapping returned by :func:`truth_by_species`.

    Returns:
        Truth categories and benchmark metadata for the call row.
    """
    species_name = infer_species_name(record=record)
    truth = truth_mapping.get(normalise_species_key(species_name))
    has_evidence = row_has_evidence(record=record)
    is_reportable = row_is_reportable(record=record)
    is_expected = bool(truth and truth.expected)

    if is_expected and is_reportable:
        fine_category = "expected_species_reportable"
        coarse_label = "expected_target"
    elif is_expected and has_evidence:
        fine_category = "expected_species_below_threshold"
        coarse_label = "observed_below_threshold"
    elif is_expected:
        fine_category = "expected_species_not_detected"
        coarse_label = "not_detected"
    elif is_reportable:
        fine_category = "unexpected_species_reportable"
        coarse_label = "reportable_off_target_species"
    elif has_evidence:
        fine_category = "unexpected_species_below_threshold"
        coarse_label = "observed_below_threshold"
    else:
        fine_category = "not_detected"
        coarse_label = "not_detected"

    return {
        "mock_truth_category": fine_category,
        "ml_report_label": coarse_label,
        "benchmark_id": truth.benchmark_id if truth else "",
        "truth_organism_id": truth.organism_id if truth else "",
        "truth_species_name": truth.species_name if truth else "",
        "expected_abundance_fraction": (
            f"{truth.expected_abundance_fraction:.8f}" if truth else ""
        ),
        "abundance_tier": truth.abundance_tier if truth else "",
        "is_expected_species": str(is_expected),
        "has_mock_evidence": str(bool(has_evidence)),
        "is_mock_reportable": str(bool(is_reportable)),
    }


def build_mock_ai_feature_records(
    *,
    records: Iterable[dict[str, object]],
    truth_records: Iterable[MockCommunityTruth],
    logger: logging.Logger | None = None,
) -> list[dict[str, object]]:
    """Build AI-ready records with generic mock-community truth labels.

    Args:
        records: KmerSutra call rows.
        truth_records: Validated truth records.
        logger: Optional logger.

    Returns:
        AI-ready records with truth metadata.
    """
    truth_mapping = truth_by_species(truth_records=truth_records)
    output: list[dict[str, object]] = []
    for record in records:
        feature_record = build_call_feature_record(record=record)
        feature_record.update(
            classify_mock_truth_category(
                record=record,
                truth_mapping=truth_mapping,
            )
        )
        feature_record["public_call"] = first_non_empty(
            record=record,
            columns=["call", "call_status", "species_call", "report_call"],
        )
        feature_record["public_species_name"] = infer_species_name(record=record)
        feature_record["public_report_label"] = record.get("report_label", "")
        feature_record["public_reference_label"] = record.get("reference_label", "")
        output.append(feature_record)
    if logger:
        counts = Counter(str(row["mock_truth_category"]) for row in output)
        for category, count in sorted(counts.items()):
            logger.info("Mock truth category %s: %d", category, count)
    return output


def count_records_by_column(
    *,
    records: Iterable[dict[str, object]],
    column: str,
) -> list[dict[str, object]]:
    """Count records by one output column.

    Args:
        records: Input rows.
        column: Column to count.

    Returns:
        Sorted count records.
    """
    counts = Counter(str(row.get(column, "")) for row in records)
    return [
        {column: value, "n_records": count}
        for value, count in sorted(counts.items())
    ]


def write_mock_ai_feature_table(
    *,
    calls_table: str | Path,
    truth_manifest: str | Path,
    output_table: str | Path,
    category_counts_table: str | Path | None = None,
    coarse_label_counts_table: str | Path | None = None,
    logger: logging.Logger | None = None,
) -> list[dict[str, object]]:
    """Label a call table using a generic mock-community manifest.

    Args:
        calls_table: KmerSutra species-call table.
        truth_manifest: Mock-community truth manifest.
        output_table: Labelled AI feature table.
        category_counts_table: Optional fine-category count output.
        coarse_label_counts_table: Optional coarse-label count output.
        logger: Optional logger.

    Returns:
        Written feature records.

    Raises:
        ValueError: If no feature records are generated.
    """
    truth_records = load_truth_manifest(
        manifest_path=truth_manifest,
        logger=logger,
    )
    call_records = read_records_table(input_path=calls_table, logger=logger)
    feature_records = build_mock_ai_feature_records(
        records=call_records,
        truth_records=truth_records,
        logger=logger,
    )
    if not feature_records:
        raise ValueError("No mock-community feature records were generated")
    fieldnames = list(feature_records[0])
    write_records_table(
        records=feature_records,
        output_path=output_table,
        fieldnames=fieldnames,
        logger=logger,
    )
    if category_counts_table is not None:
        write_records_table(
            records=count_records_by_column(
                records=feature_records,
                column="mock_truth_category",
            ),
            output_path=category_counts_table,
            fieldnames=["mock_truth_category", "n_records"],
            logger=logger,
        )
    if coarse_label_counts_table is not None:
        write_records_table(
            records=count_records_by_column(
                records=feature_records,
                column="ml_report_label",
            ),
            output_path=coarse_label_counts_table,
            fieldnames=["ml_report_label", "n_records"],
            logger=logger,
        )
    return feature_records
