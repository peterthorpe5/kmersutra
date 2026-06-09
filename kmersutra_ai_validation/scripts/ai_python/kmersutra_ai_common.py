#!/usr/bin/env python3
"""Common helper functions for KmerSutra AI validation workflows.

The helpers in this module are deliberately lightweight and mostly independent
of the main KmerSutra package so that they can be unit-tested without requiring
large benchmark files or a full KmerSutra installation.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import logging
import sys
from collections import Counter
from pathlib import Path
from typing import Iterable, Sequence


DEFAULT_SAFE_FEATURE_COLUMNS = [
    "n_hits",
    "n_unique_kmers",
    "n_positive_sequences",
    "n_k_values_positive",
    "best_k",
    "n_exact_hits",
    "n_fuzzy_hits",
    "conflicting_unique_kmers",
    "conflict_ratio",
    "reportable_conflicting_unique_kmers",
    "reportable_conflict_ratio",
    "mixed_species_support_fraction",
    "confidence_score",
    "signal_confidence_score",
    "has_long_k_support",
    "has_multi_k_support",
    "exact_hit_fraction",
]

LEAKAGE_RISK_COLUMNS = {
    "spike_n",
    "spike_n_per_genome",
    "total_spike_n",
    "sample_id",
    "species_name",
    "benchmark_family",
    "panel",
    "replicate",
    "is_negative",
    "is_shuffled_control",
    "expected_target",
    "expected_species",
    "truth_label",
    "ml_report_label",
    "call",
    "call_status",
    "report_label",
    "reference_label",
}

EXPECTED_ZYMO_D6300_SPECIES = {
    "Pseudomonas aeruginosa",
    "Escherichia coli",
    "Salmonella enterica",
    "Limosilactobacillus fermentum",
    "Enterococcus faecalis",
    "Staphylococcus aureus",
    "Listeria monocytogenes",
    "Bacillus subtilis",
    "Saccharomyces cerevisiae",
    "Cryptococcus neoformans",
}


class ValidationError(RuntimeError):
    """Error raised for invalid validation inputs."""


def configure_logging(*, log_path: Path, verbose: bool) -> logging.Logger:
    """Configure console and file logging.

    Parameters
    ----------
    log_path : pathlib.Path
        Path to the log file to write.
    verbose : bool
        If true, write debug-level messages to the console.

    Returns
    -------
    logging.Logger
        Configured logger.
    """
    logger = logging.getLogger(log_path.stem)
    logger.setLevel(logging.DEBUG)
    logger.handlers.clear()
    logger.propagate = False

    formatter = logging.Formatter(
        "%(asctime)s %(levelname)s [%(name)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    stream_handler = logging.StreamHandler(sys.stderr)
    stream_handler.setLevel(logging.DEBUG if verbose else logging.INFO)
    stream_handler.setFormatter(formatter)
    logger.addHandler(stream_handler)

    log_path.parent.mkdir(parents=True, exist_ok=True)
    file_handler = logging.FileHandler(log_path)
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    return logger


def existing_file(path: str | Path, description: str) -> Path:
    """Validate that a file exists and is not empty.

    Parameters
    ----------
    path : str or pathlib.Path
        File path to validate.
    description : str
        Human-readable description for error messages.

    Returns
    -------
    pathlib.Path
        Validated path.

    Raises
    ------
    ValidationError
        If the path is missing or empty.
    """
    parsed = Path(path)
    if not parsed.exists():
        raise ValidationError(f"{description} does not exist: {parsed}")
    if not parsed.is_file():
        raise ValidationError(f"{description} is not a file: {parsed}")
    if parsed.stat().st_size == 0:
        raise ValidationError(f"{description} is empty: {parsed}")
    return parsed


def ensure_output_dir(path: str | Path) -> Path:
    """Create and return an output directory.

    Parameters
    ----------
    path : str or pathlib.Path
        Output directory path.

    Returns
    -------
    pathlib.Path
        Created directory.
    """
    parsed = Path(path)
    parsed.mkdir(parents=True, exist_ok=True)
    return parsed


def normalise_name(value: object) -> str:
    """Return a stripped string representation of a label or name.

    Parameters
    ----------
    value : object
        Input value.

    Returns
    -------
    str
        Normalised string.
    """
    return str(value or "").strip()


def normalise_status(value: object) -> str:
    """Return a lower-case status string.

    Parameters
    ----------
    value : object
        Input status value.

    Returns
    -------
    str
        Lower-case status string.
    """
    return normalise_name(value).lower().replace(" ", "_")


def safe_value(record: dict[str, object], column: str) -> str:
    """Return a non-empty string value for a split/group column.

    Parameters
    ----------
    record : dict[str, object]
        Input record.
    column : str
        Column name.

    Returns
    -------
    str
        Clean value, or ``missing``.
    """
    value = normalise_name(record.get(column, ""))
    return value if value else "missing"


def stable_group_hash(value: str) -> int:
    """Return a deterministic hash for grouped splitting.

    Parameters
    ----------
    value : str
        Value to hash.

    Returns
    -------
    int
        Stable integer hash.
    """
    digest = hashlib.sha256(value.encode("utf-8")).hexdigest()
    return int(digest[:12], 16)


def sanitise_name(value: str) -> str:
    """Return a filesystem-safe name fragment.

    Parameters
    ----------
    value : str
        Input value.

    Returns
    -------
    str
        Sanitised value.
    """
    cleaned = normalise_name(value)
    for old, new in [
        (" ", "_"),
        ("/", "_"),
        ("\\", "_"),
        ("(", ""),
        (")", ""),
        (":", "_"),
        (";", "_"),
        (",", "_"),
    ]:
        cleaned = cleaned.replace(old, new)
    return cleaned or "missing"


def label_counts(records: Iterable[dict[str, object]], label_column: str) -> Counter:
    """Count labels in records.

    Parameters
    ----------
    records : iterable of dict[str, object]
        Input records.
    label_column : str
        Label column name.

    Returns
    -------
    collections.Counter
        Label counts.
    """
    return Counter(
        normalise_name(row.get(label_column, ""))
        for row in records
        if normalise_name(row.get(label_column, ""))
    )


def label_counts_text(records: Iterable[dict[str, object]], label_column: str) -> str:
    """Return compact label-count text.

    Parameters
    ----------
    records : iterable of dict[str, object]
        Input records.
    label_column : str
        Label column name.

    Returns
    -------
    str
        Semicolon-delimited label counts.
    """
    counts = label_counts(records=records, label_column=label_column)
    return ";".join(f"{key}:{counts[key]}" for key in sorted(counts))


def count_records_by_column(
    *,
    records: Iterable[dict[str, object]],
    column: str,
) -> list[dict[str, object]]:
    """Count records by a column.

    Parameters
    ----------
    records : iterable of dict[str, object]
        Input records.
    column : str
        Column to count.

    Returns
    -------
    list[dict[str, object]]
        Count table records.
    """
    counts = Counter(normalise_name(row.get(column, "")) for row in records)
    return [
        {column: value, "n_records": count}
        for value, count in sorted(counts.items())
    ]


def write_json(data: object, output_path: str | Path) -> None:
    """Write JSON using stable formatting.

    Parameters
    ----------
    data : object
        JSON-serialisable object.
    output_path : str or pathlib.Path
        Output path.
    """
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True), encoding="utf-8")


def open_text_auto(path: str | Path, mode: str = "rt"):
    """Open plain or gzip-compressed text.

    Parameters
    ----------
    path : str or pathlib.Path
        Input path.
    mode : str, optional
        Text mode.

    Returns
    -------
    file-like
        Open handle.
    """
    parsed = Path(path)
    if parsed.suffix == ".gz":
        return gzip.open(parsed, mode)
    return parsed.open(mode, newline="")


def read_tsv_simple(path: str | Path) -> list[dict[str, str]]:
    """Read a small TSV or TSV.GZ file without pandas.

    Parameters
    ----------
    path : str or pathlib.Path
        Input table.

    Returns
    -------
    list[dict[str, str]]
        Parsed rows.
    """
    with open_text_auto(path, "rt") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv_simple(
    *,
    records: Iterable[dict[str, object]],
    output_path: str | Path,
    fieldnames: Sequence[str] | None = None,
) -> int:
    """Write a TSV or TSV.GZ file without pandas.

    Parameters
    ----------
    records : iterable of dict[str, object]
        Records to write.
    output_path : str or pathlib.Path
        Output path.
    fieldnames : sequence of str or None, optional
        Output columns. If omitted, the first row keys are used.

    Returns
    -------
    int
        Number of rows written.
    """
    rows = [dict(record) for record in records]
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)

    if fieldnames is None:
        fieldnames = list(rows[0].keys()) if rows else []

    with open_text_auto(path, "wt") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=list(fieldnames),
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    return len(rows)


def read_reference_label_map(
    *,
    reference_label_map: str | Path | None,
    logger: logging.Logger | None = None,
) -> tuple[set[str], set[str], dict[str, str]]:
    """Read expected Zymo reference labels/species from a map table.

    Parameters
    ----------
    reference_label_map : str, pathlib.Path or None
        Optional reference-label map.
    logger : logging.Logger or None, optional
        Logger.

    Returns
    -------
    tuple[set[str], set[str], dict[str, str]]
        Expected reference labels, expected species names and a reference-label
        to role mapping.
    """
    expected_reference_labels: set[str] = set()
    expected_species_names: set[str] = set(EXPECTED_ZYMO_D6300_SPECIES)
    reference_label_roles: dict[str, str] = {}

    if reference_label_map is None or not str(reference_label_map):
        if logger:
            logger.warning("No reference-label map supplied")
        return expected_reference_labels, expected_species_names, reference_label_roles

    path = Path(reference_label_map)
    if not path.exists() or path.stat().st_size == 0:
        if logger:
            logger.warning("Reference-label map missing or empty: %s", path)
        return expected_reference_labels, expected_species_names, reference_label_roles

    for row in read_tsv_simple(path):
        role = normalise_status(row.get("role"))
        reference_label = normalise_name(
            row.get("reference_label")
            or row.get("label")
            or row.get("species_name")
        )
        original_species = normalise_name(
            row.get("original_species_name")
            or row.get("species")
            or row.get("species_name")
        )
        species_name = normalise_name(row.get("species_name"))

        if reference_label:
            reference_label_roles[reference_label] = role
        if role == "target_species":
            for value in [reference_label, species_name]:
                if value:
                    expected_reference_labels.add(value)
            for value in [original_species, species_name]:
                if value:
                    expected_species_names.add(value)

    if logger:
        logger.info("Expected reference labels: %d", len(expected_reference_labels))
        logger.info("Expected species names: %d", len(expected_species_names))

    return expected_reference_labels, expected_species_names, reference_label_roles


def infer_public_truth_label(
    *,
    row: dict[str, object],
    expected_reference_labels: set[str],
    expected_species_names: set[str],
) -> str:
    """Infer a public Zymo species-level truth/evidence label.

    This function is intentionally conservative. It is for species/sample-level
    validation of a previously trained evidence calibrator, not for read-level
    truth labelling.

    Parameters
    ----------
    row : dict[str, object]
        KmerSutra species-detection row.
    expected_reference_labels : set[str]
        Expected reference labels from ``reference_label_map.tsv``.
    expected_species_names : set[str]
        Expected species names.

    Returns
    -------
    str
        Label compatible with the internal calibrator.
    """
    candidate_values = {
        normalise_name(row.get("species_name")),
        normalise_name(row.get("report_label")),
        normalise_name(row.get("reference_label")),
        normalise_name(row.get("original_species_name")),
    }
    candidate_values = {value for value in candidate_values if value}

    if candidate_values & expected_reference_labels:
        return "expected_target"
    if candidate_values & expected_species_names:
        return "expected_target"

    call = normalise_status(row.get("call")) or normalise_status(row.get("call_status"))
    decision = normalise_status(row.get("decision"))

    combined_status = " ".join([call, decision])
    if "present" in combined_status or "reportable" in combined_status:
        return "reportable_off_target_species"
    if "below" in combined_status or "neighbour" in combined_status:
        return "observed_below_threshold"
    if "lineage" in combined_status or "background" in combined_status:
        return "observed_below_threshold"

    n_unique = parse_float(row.get("n_unique_kmers", 0.0))
    n_sequences = parse_float(row.get("n_positive_sequences", 0.0))
    if n_unique > 0.0 or n_sequences > 0.0:
        return "observed_below_threshold"

    return "not_detected"


def parse_float(value: object) -> float:
    """Parse a float with defensive handling of missing values.

    Parameters
    ----------
    value : object
        Input value.

    Returns
    -------
    float
        Parsed float or zero.
    """
    if value is None:
        return 0.0
    text = str(value).strip()
    if not text or text.lower() in {"na", "nan", "none"}:
        return 0.0
    try:
        return float(text)
    except ValueError:
        return 0.0


def find_present_features(
    *,
    records: list[dict[str, object]],
    requested_features: Sequence[str],
) -> list[str]:
    """Return requested features present in the first record.

    Parameters
    ----------
    records : list of dict[str, object]
        Input records.
    requested_features : sequence of str
        Requested feature names.

    Returns
    -------
    list[str]
        Present feature names.
    """
    if not records:
        return []
    first = records[0]
    return [feature for feature in requested_features if feature in first]


def audit_feature_columns(features: Sequence[str]) -> list[dict[str, object]]:
    """Audit whether feature columns are likely safe or leakage-prone.

    Parameters
    ----------
    features : sequence of str
        Feature names.

    Returns
    -------
    list[dict[str, object]]
        Audit rows.
    """
    rows = []
    for feature in features:
        leakage_risk = feature in LEAKAGE_RISK_COLUMNS
        rows.append(
            {
                "feature": feature,
                "is_leakage_risk": "yes" if leakage_risk else "no",
                "recommended_for_training": "no" if leakage_risk else "yes",
            }
        )
    return rows


def add_common_parser_arguments(parser: argparse.ArgumentParser) -> None:
    """Add common verbosity argument to a parser.

    Parameters
    ----------
    parser : argparse.ArgumentParser
        Parser to modify.
    """
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Write verbose progress messages to stderr and the log file.",
    )
