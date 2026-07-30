#!/usr/bin/env python3
"""Download NCBI genomes and create KmerSutra-ready metadata.

This script replaces older pyani-style genome download helpers for the
KmerSutra workflow. It searches NCBI Assembly for one or more taxonomy IDs,
downloads genome FASTA files from GenBank or RefSeq FTP locations, organises
assemblies into taxon/species/accession folders, and writes TSV metadata files
including a KmerSutra genome configuration table.

The script intentionally writes tab-separated files rather than comma-separated
files so that outputs are consistent with the wider KmerSutra project.
"""

from __future__ import annotations

import argparse
import gzip
import http.client
import json
import logging
import re
import shutil
import sys
import time
from xml.parsers.expat import ExpatError
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

Entrez: object | None = None

LOGGER = logging.getLogger("download_ncbi_genomes_for_kmersutra")

ASSEMBLY_LEVEL_RANK = {
    "complete genome": 4,
    "chromosome": 3,
    "scaffold": 2,
    "contig": 1,
    "": 0,
}

DEFAULT_FORMATS = ("genomic_fna",)
FORMAT_SUFFIXES = {
    "genomic_fna": "genomic.fna.gz",
    "genomic_gff": "genomic.gff.gz",
    "genomic_gbff": "genomic.gbff.gz",
    "protein_faa": "protein.faa.gz",
    "assembly_report": "assembly_report.txt",
    "assembly_stats": "assembly_stats.txt",
}

TRANSIENT_ENTREZ_ERRORS = (
    HTTPError,
    URLError,
    TimeoutError,
    OSError,
    http.client.IncompleteRead,
    ExpatError,
    ValueError,
)

KMERSUTRA_CONFIG_HEADER = [
    "genome_fasta",
    "species_name",
    "strain_name",
    "taxid",
    "role",
    "clade",
    "assembly_accession",
    "query_taxid",
    "assembly_level",
    "scaffold_n50",
    "contig_n50",
]


@dataclass(frozen=True)
class TaxonPlan:
    """Describe how one queried taxon should be downloaded and labelled."""

    taxid: str
    role: str = "downloaded"
    clade: str = ""
    group_label: str = ""
    max_assemblies: int | None = None
    best_per_species: int | None = None
    min_total_length: int | None = None
    max_total_length: int | None = None
    min_scaffold_n50: int | None = None
    min_contig_n50: int | None = None


@dataclass(frozen=True)
class AssemblyRecord:
    """Represent selected NCBI Assembly metadata used by KmerSutra."""

    query_taxid: str
    assembly_uid: str
    assembly_accession: str
    assembly_name: str
    organism_name: str
    species_name: str
    species_taxid: str
    taxid: str
    strain_name: str
    assembly_level: str
    refseq_category: str
    scaffold_n50: int
    contig_n50: int
    total_length: int
    ftp_path_refseq: str
    ftp_path_genbank: str
    selected_source: str
    selected_ftp_path: str
    role: str
    clade: str
    group_label: str


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Download genomes from NCBI Assembly by taxid and write "
            "KmerSutra-ready genome configuration files."
        )
    )
    parser.add_argument(
        "--taxids",
        nargs="+",
        default=None,
        help="One or more NCBI taxonomy IDs to search below.",
    )
    parser.add_argument(
        "--taxid_plan",
        default=None,
        help=(
            "Optional TSV with columns: taxid, role, clade, group_label, "
            "max_assemblies, best_per_species, min_total_length, "
            "max_total_length, min_scaffold_n50, min_contig_n50. Values in "
            "this table can override command-line defaults per taxid."
        ),
    )
    parser.add_argument(
        "--custom_genome_config",
        default=None,
        help=(
            "Optional user-supplied KmerSutra genome_config TSV to append to "
            "downloaded NCBI assemblies. This allows local/custom genomes to be "
            "included without going through Entrez."
        ),
    )
    parser.add_argument(
        "--exclude_accessions_table",
        default=None,
        help=(
            "Optional TSV, TSV.GZ or Parquet table of truth accessions that "
            "must not enter the reference panel. Recognised columns include "
            "assembly_accession, truth_assembly_accessions, sequence_accession, "
            "truth_sequence_accessions and truth_reference_accession. Multiple "
            "values in one cell must be separated by semicolons."
        ),
    )
    parser.add_argument(
        "--allow_missing_custom_fastas",
        action="store_true",
        help=(
            "Allow custom genome-config rows whose genome_fasta paths are missing "
            "or empty. By default these are rejected defensively."
        ),
    )
    parser.add_argument(
        "--out_dir",
        required=True,
        help="Output directory for organised genome downloads and metadata.",
    )
    parser.add_argument(
        "--email",
        default=None,
        help=(
            "Email address supplied to NCBI Entrez. Required when downloading "
            "from NCBI, but not needed when only normalising a custom genome config."
        ),
    )
    parser.add_argument(
        "--api_key",
        default=None,
        help="Optional NCBI API key for higher request-rate limits.",
    )
    parser.add_argument(
        "--source",
        choices=("prefer_refseq", "prefer_genbank", "refseq", "genbank"),
        default="prefer_refseq",
        help="Assembly source preference. Default: prefer_refseq.",
    )
    parser.add_argument(
        "--formats",
        nargs="+",
        default=list(DEFAULT_FORMATS),
        choices=sorted(FORMAT_SUFFIXES),
        help="NCBI assembly file formats to download. Default: genomic_fna.",
    )
    parser.add_argument(
        "--assembly_levels",
        nargs="+",
        default=None,
        help=(
            "Optional assembly-level filter. Matching is case-insensitive "
            "and whitespace-normalised, e.g. 'complete genome', "
            "'Complete Genome', 'chromosome' or 'Scaffold'. If omitted, "
            "all levels are kept."
        ),
    )
    parser.add_argument(
        "--min_total_length",
        type=int,
        default=None,
        help=(
            "Optional minimum assembly length retained across all taxids. "
            "Per-taxid values in --taxid_plan override this."
        ),
    )
    parser.add_argument(
        "--max_total_length",
        type=int,
        default=None,
        help=(
            "Optional maximum assembly length retained across all taxids. "
            "Per-taxid values in --taxid_plan override this."
        ),
    )
    parser.add_argument(
        "--min_scaffold_n50",
        type=int,
        default=None,
        help=(
            "Optional minimum scaffold N50 retained across all taxids. "
            "Per-taxid values in --taxid_plan override this."
        ),
    )
    parser.add_argument(
        "--min_contig_n50",
        type=int,
        default=None,
        help=(
            "Optional minimum contig N50 retained across all taxids. "
            "Per-taxid values in --taxid_plan override this."
        ),
    )
    parser.add_argument(
        "--max_assemblies_per_taxid",
        type=int,
        default=None,
        help="Optional maximum assemblies retained per queried taxid.",
    )
    parser.add_argument(
        "--best_per_species",
        type=int,
        default=None,
        help=(
            "Optionally retain the best N assemblies per species within "
            "each queried taxid, ranked by assembly level, scaffold N50, "
            "contig N50, and total length."
        ),
    )
    parser.add_argument(
        "--default_role",
        default="downloaded",
        help="Default KmerSutra role used when no taxid plan is supplied.",
    )
    parser.add_argument(
        "--default_clade",
        default="",
        help="Default clade label used when no taxid plan is supplied.",
    )
    parser.add_argument(
        "--include_unplaced",
        action="store_true",
        help="Include assemblies without a downloadable FTP path.",
    )
    parser.add_argument(
        "--metadata_only",
        action="store_true",
        help="Only write metadata and KmerSutra config paths; do not download files.",
    )
    parser.add_argument(
        "--decompress",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Decompress downloaded genomic FASTA archives. Default: true.",
    )
    parser.add_argument(
        "--delete_archives",
        action="store_true",
        help="Delete .gz archives after successful decompression.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing downloaded files.",
    )
    parser.add_argument(
        "--sleep_seconds",
        type=float,
        default=0.34,
        help="Delay between NCBI Entrez requests. Default: 0.34 seconds.",
    )
    parser.add_argument(
        "--retries",
        type=int,
        default=5,
        help="Number of retries for Entrez and file download requests.",
    )
    parser.add_argument(
        "--batch_size",
        type=int,
        default=200,
        help=(
            "Number of assembly records retrieved per Entrez summary batch. "
            "Default: 200, chosen to reduce truncated chunked Entrez responses."
        ),
    )
    parser.add_argument(
        "--max_entrez_records_per_taxid",
        type=int,
        default=None,
        help=(
            "Optional cap on assembly UIDs fetched from Entrez for each queried "
            "taxid before summary retrieval. This is a safety valve for very broad "
            "taxa that return many thousands of assemblies."
        ),
    )
    parser.add_argument(
        "--log_file",
        default=None,
        help="Optional path for a detailed log file. Defaults to out_dir/logs/download.log.",
    )
    parser.add_argument(
        "--dry_run",
        action="store_true",
        help="Search/filter and write planned metadata, but do not download files.",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print detailed logging to stderr.",
    )
    return parser.parse_args(argv)


def configure_logging(log_file: Path, verbose: bool = False) -> None:
    """Configure console and file logging."""
    for handler in list(LOGGER.handlers):
        LOGGER.removeHandler(handler)
        handler.close()
    LOGGER.setLevel(logging.DEBUG)
    log_file.parent.mkdir(parents=True, exist_ok=True)

    formatter = logging.Formatter(
        "%(asctime)s %(levelname)s [%(name)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    stream_handler = logging.StreamHandler(sys.stderr)
    stream_handler.setLevel(logging.INFO if verbose else logging.WARNING)
    stream_handler.setFormatter(formatter)
    LOGGER.addHandler(stream_handler)

    file_handler = logging.FileHandler(log_file, mode="w", encoding="utf-8")
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)
    LOGGER.addHandler(file_handler)


def read_tsv(path: Path) -> list[dict[str, str]]:
    """Read a tab-separated file into dictionaries."""
    rows: list[dict[str, str]] = []
    with path.open("r", encoding="utf-8") as handle:
        header_line = handle.readline().rstrip("\n")
        if not header_line:
            return rows
        header = header_line.split("\t")
        for line_number, line in enumerate(handle, start=2):
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            values = line.split("\t")
            if len(values) < len(header):
                values.extend([""] * (len(header) - len(values)))
            if len(values) > len(header):
                LOGGER.warning(
                    "Ignoring extra fields in %s at line %s", path, line_number
                )
                values = values[: len(header)]
            rows.append(dict(zip(header, values)))
    return rows


EXCLUSION_ACCESSION_COLUMNS = {
    "assembly_accession",
    "truth_assembly_accession",
    "truth_assembly_accessions",
    "sequence_accession",
    "truth_sequence_accession",
    "truth_sequence_accessions",
    "truth_reference_accession",
}


def normalise_accession(value: Any) -> str:
    """Return an upper-case accession without a version suffix.

    Args:
        value: Assembly or sequence accession.

    Returns:
        Normalised accession.
    """
    return str(value or "").strip().upper().split(".", 1)[0]


def load_excluded_accessions(path: str | Path | None) -> set[str]:
    """Load accession exclusions from a tabular manifest.

    Args:
        path: Optional TSV, TSV.GZ or Parquet table.

    Returns:
        Normalised accession exclusions.

    Raises:
        ValueError: If the table has no recognised accession columns.
    """
    if path is None:
        return set()
    from kmersutra.table_io import read_records_table

    rows = read_records_table(input_path=path)
    if not rows:
        raise ValueError(f"Accession exclusion table contains no rows: {path}")
    observed = EXCLUSION_ACCESSION_COLUMNS.intersection(rows[0])
    if not observed:
        raise ValueError(
            "Accession exclusion table has no recognised columns. Expected one "
            "of: " + ", ".join(sorted(EXCLUSION_ACCESSION_COLUMNS))
        )
    accessions: set[str] = set()
    for row in rows:
        for column in observed:
            for value in str(row.get(column, "") or "").split(";"):
                accession = normalise_accession(value)
                if accession:
                    accessions.add(accession)
    LOGGER.info("Loaded %d explicit truth-accession exclusions", len(accessions))
    return accessions


def exclude_records_by_accession(
    *,
    records: list[AssemblyRecord],
    audit_rows: list[dict[str, Any]],
    excluded_accessions: set[str],
) -> list[AssemblyRecord]:
    """Remove candidate assemblies listed in an exclusion manifest.

    Args:
        records: Candidate assemblies.
        audit_rows: Candidate audit rows to update in place.
        excluded_accessions: Normalised accessions that must be excluded.

    Returns:
        Records not explicitly excluded.
    """
    if not excluded_accessions:
        return records
    retained: list[AssemblyRecord] = []
    excluded: set[str] = set()
    for record in records:
        accession = normalise_accession(record.assembly_accession)
        if accession and accession in excluded_accessions:
            excluded.add(record.assembly_accession)
            continue
        retained.append(record)
    for row in audit_rows:
        if str(row.get("assembly_accession", "")) in excluded:
            row["filter_status"] = "excluded"
            row["filter_reason"] = "excluded_by_accession_table"
    if excluded:
        LOGGER.info(
            "Excluded %d candidate assemblies by explicit accession",
            len(excluded),
        )
    return retained


def fasta_excluded_accessions(
    *,
    fasta_path: str | Path,
    excluded_accessions: set[str],
) -> set[str]:
    """Return truth accessions found in FASTA headers.

    Args:
        fasta_path: Panel FASTA path.
        excluded_accessions: Normalised truth accessions.

    Returns:
        Matching accessions.
    """
    if not excluded_accessions:
        return set()
    from kmersutra.reference_audit import fasta_header_accessions

    return fasta_header_accessions(Path(fasta_path)).intersection(
        excluded_accessions
    )


def write_tsv(path: Path, rows: list[dict[str, Any]], header: list[str]) -> None:
    """Write dictionaries to a tab-separated file with a fixed header."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("\t".join(header) + "\n")
        for row in rows:
            values = [normalise_tsv_value(row.get(column, "")) for column in header]
            handle.write("\t".join(values) + "\n")


def normalise_tsv_value(value: Any) -> str:
    """Convert a value to TSV-safe text."""
    if value is None:
        return ""
    text = str(value)
    return text.replace("\t", " ").replace("\n", " ").replace("\r", " ")


def safe_name(value: str, fallback: str = "unknown") -> str:
    """Return a filesystem-safe label."""
    text = value.strip() or fallback
    text = re.sub(r"[^A-Za-z0-9._-]+", "_", text)
    text = re.sub(r"_+", "_", text)
    return text.strip("_") or fallback


def parse_int(value: Any, default: int = 0) -> int:
    """Parse an integer-like value safely."""
    try:
        if value in (None, ""):
            return default
        return int(float(str(value).replace(",", "")))
    except (TypeError, ValueError):
        return default


def build_taxon_plan(args: argparse.Namespace) -> list[TaxonPlan]:
    """Build a list of taxon download plans from CLI arguments."""
    plans: list[TaxonPlan] = []
    if args.taxid_plan:
        for row in read_tsv(Path(args.taxid_plan)):
            taxid = row.get("taxid", "").strip()
            if not taxid:
                continue
            plans.append(
                TaxonPlan(
                    taxid=taxid,
                    role=row.get("role", args.default_role).strip()
                    or args.default_role,
                    clade=row.get("clade", args.default_clade).strip()
                    or args.default_clade,
                    group_label=row.get("group_label", "").strip(),
                    max_assemblies=parse_optional_int(row.get("max_assemblies", "")),
                    best_per_species=parse_optional_int(row.get("best_per_species", "")),
                    min_total_length=parse_optional_int(row.get("min_total_length", "")),
                    max_total_length=parse_optional_int(row.get("max_total_length", "")),
                    min_scaffold_n50=parse_optional_int(row.get("min_scaffold_n50", "")),
                    min_contig_n50=parse_optional_int(row.get("min_contig_n50", "")),
                )
            )
    if args.taxids:
        for taxid in args.taxids:
            plans.append(
                TaxonPlan(
                    taxid=str(taxid),
                    role=args.default_role,
                    clade=args.default_clade,
                    group_label="",
                    max_assemblies=args.max_assemblies_per_taxid,
                    best_per_species=args.best_per_species,
                    min_total_length=args.min_total_length,
                    max_total_length=args.max_total_length,
                    min_scaffold_n50=args.min_scaffold_n50,
                    min_contig_n50=args.min_contig_n50,
                )
            )
    if not plans and not getattr(args, "custom_genome_config", None):
        raise ValueError(
            "At least one taxid must be supplied via --taxids/--taxid_plan "
            "or provide --custom_genome_config"
        )
    return plans


def parse_optional_int(value: Any) -> int | None:
    """Parse an optional integer from TSV configuration text."""
    if value in (None, ""):
        return None
    return parse_int(value, default=0)


def normalise_assembly_level(value: Any) -> str:
    """Return a case-insensitive normal form for an NCBI assembly level.

    Parameters
    ----------
    value : object
        Raw assembly level from NCBI, such as ``Complete Genome`` or
        ``Scaffold``.

    Returns
    -------
    str
        Lowercase, whitespace-normalised assembly level.
    """
    return " ".join(str(value or "").strip().lower().split())


def normalise_assembly_level_set(levels: list[str] | None) -> set[str]:
    """Normalise requested assembly levels for robust filtering."""
    return {
        normalise_assembly_level(level)
        for level in (levels or [])
        if normalise_assembly_level(level)
    }


def normalise_download_url(path_or_url: str, suffix: str | None = None) -> str:
    """Normalise an NCBI FTP/HTTPS assembly path to a download URL.

    Parameters
    ----------
    path_or_url : str
        NCBI assembly directory or file URL. Historical NCBI metadata fields
        are named ``ftp_path`` and may contain ``ftp://`` URLs.
    suffix : str | None, optional
        Filename suffix such as ``genomic.fna.gz``. When supplied, the final
        filename is derived from the assembly directory basename.

    Returns
    -------
    str
        HTTPS download URL, or an empty string for blank input.
    """
    raw = str(path_or_url or "").strip()
    if not raw:
        return ""
    if raw.startswith("ftp://"):
        raw = "https://" + raw[len("ftp://") :]
    if suffix is None:
        return raw
    base_url = raw.rstrip("/")
    base_name = base_url.split("/")[-1]
    return f"{base_url}/{base_name}_{suffix}"


def count_values(rows: list[dict[str, Any]], column: str) -> dict[str, int]:
    """Count string values from a list of dictionaries."""
    counts: dict[str, int] = {}
    for row in rows:
        key = str(row.get(column, "") or "")
        counts[key] = counts.get(key, 0) + 1
    return counts


def configure_entrez(email: str, api_key: str | None = None) -> None:
    """Configure Bio.Entrez identity metadata.

    Parameters
    ----------
    email : str
        Contact email address required by NCBI Entrez.
    api_key : str | None
        Optional NCBI API key.

    Raises
    ------
    ImportError
        Raised if BioPython is not installed.
    """
    global Entrez
    try:
        from Bio import Entrez as bio_entrez  # type: ignore
    except ImportError as exc:  # pragma: no cover - environment dependent
        raise ImportError(
            "BioPython is required for Entrez access. Install with: "
            "pip install biopython"
        ) from exc

    Entrez = bio_entrez
    Entrez.email = email
    Entrez.tool = "download_ncbi_genomes_for_kmersutra"
    if api_key:
        Entrez.api_key = api_key
    LOGGER.info("Configured Entrez with email: %s", email)


def entrez_retry(
    func: Any,
    *args: Any,
    retries: int = 5,
    sleep_seconds: float = 0.34,
    **kwargs: Any,
) -> Any:
    """Run an Entrez call with retry and polite sleep handling."""
    last_error: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            handle = func(*args, **kwargs)
            time.sleep(sleep_seconds)
            return handle
        except TRANSIENT_ENTREZ_ERRORS as exc:
            last_error = exc
            wait_time = sleep_seconds * attempt
            LOGGER.warning(
                "Entrez request failed on attempt %s/%s: %s",
                attempt,
                retries,
                exc,
            )
            time.sleep(wait_time)
    raise RuntimeError(f"Entrez request failed after {retries} attempts") from last_error



def entrez_read_retry(
    func: Any,
    *args: Any,
    retries: int = 5,
    sleep_seconds: float = 0.34,
    validate: bool = False,
    context: str = "Entrez request",
    **kwargs: Any,
) -> Any:
    """Run an Entrez request and parse the response with retries.

    Parameters
    ----------
    func : callable
        Entrez function such as ``esearch`` or ``esummary``.
    *args : object
        Positional arguments passed to ``func``.
    retries : int, optional
        Maximum number of attempts.
    sleep_seconds : float, optional
        Polite base delay between attempts. Backoff is linear by attempt.
    validate : bool, optional
        Value passed to ``Entrez.read``.
    context : str, optional
        Human-readable context included in warning/error messages, such as a
        taxid, query or summary batch range.
    **kwargs : object
        Keyword arguments passed to ``func``.

    Returns
    -------
    object
        Parsed Entrez record.

    Raises
    ------
    RuntimeError
        Raised after all attempts fail.

    Notes
    -----
    Biopython can raise ``http.client.IncompleteRead`` while parsing a
    truncated chunked Entrez response. Retrying only handle creation is not
    sufficient in that case because the failure occurs during ``Entrez.read``.
    This helper retries the complete request-plus-parse operation.
    """
    entrez = require_entrez()
    last_error: Exception | None = None
    retries = max(1, int(retries))
    for attempt in range(1, retries + 1):
        handle = None
        try:
            LOGGER.debug("Starting %s, attempt %s/%s", context, attempt, retries)
            handle = func(*args, **kwargs)
            record = entrez.read(handle, validate=validate)
            time.sleep(sleep_seconds)
            return record
        except TRANSIENT_ENTREZ_ERRORS as exc:
            last_error = exc
            wait_time = max(0.0, sleep_seconds) * attempt
            LOGGER.warning(
                "%s failed on attempt %s/%s: %s",
                context,
                attempt,
                retries,
                exc,
            )
            if attempt < retries:
                time.sleep(wait_time)
        finally:
            if handle is not None:
                try:
                    handle.close()
                except Exception:  # pragma: no cover - defensive close only
                    pass
    raise RuntimeError(
        f"{context} failed after {retries} attempt(s): {last_error!r}"
    ) from last_error

def require_entrez() -> object:
    """Return the configured Entrez module.

    Returns
    -------
    object
        Bio.Entrez module.

    Raises
    ------
    RuntimeError
        Raised if Entrez has not been configured.
    """
    if Entrez is None:
        raise RuntimeError("Entrez has not been configured. Call configure_entrez() first.")
    return Entrez


def search_assembly_uids(
    taxid: str,
    retries: int,
    sleep_seconds: float,
    max_records: int | None = None,
) -> list[str]:
    """Search NCBI Assembly for all assemblies under a taxonomy subtree.

    Parameters
    ----------
    taxid : str
        NCBI taxonomy ID to search using the organism-expansion query.
    retries : int
        Maximum Entrez retry attempts.
    sleep_seconds : float
        Polite delay and retry backoff base in seconds.
    max_records : int | None, optional
        Optional cap on the number of assembly UIDs returned for very broad
        taxid searches.

    Returns
    -------
    list[str]
        Assembly UIDs returned by Entrez.
    """
    query = f"txid{taxid}[Organism:exp]"
    LOGGER.info("Searching NCBI Assembly with query: %s", query)
    entrez = require_entrez()
    count_context = f"Entrez esearch count taxid={taxid} query={query!r}"
    count_record = entrez_read_retry(
        entrez.esearch,
        db="assembly",
        term=query,
        retmax=0,
        usehistory="y",
        retries=retries,
        sleep_seconds=sleep_seconds,
        context=count_context,
    )
    count = int(count_record.get("Count", 0))
    LOGGER.info("Taxid %s returned %s assembly records", taxid, count)
    if count == 0:
        return []
    retmax = count
    if max_records is not None and max_records > 0 and count > max_records:
        LOGGER.warning(
            "Taxid %s returned %s records; capping Entrez UID retrieval to %s",
            taxid,
            count,
            max_records,
        )
        retmax = max_records
    fetch_context = (
        f"Entrez esearch UIDs taxid={taxid} query={query!r} "
        f"retmax={retmax} count={count}"
    )
    uid_record = entrez_read_retry(
        entrez.esearch,
        db="assembly",
        term=query,
        retmax=retmax,
        retries=retries,
        sleep_seconds=sleep_seconds,
        context=fetch_context,
    )
    return list(uid_record.get("IdList", []))


def fetch_assembly_summaries(
    assembly_uids: list[str],
    batch_size: int,
    retries: int,
    sleep_seconds: float,
    taxid: str = "",
    query: str = "",
) -> list[dict[str, Any]]:
    """Fetch full NCBI Assembly summaries in retryable batches.

    Parameters
    ----------
    assembly_uids : list[str]
        Assembly UIDs to fetch.
    batch_size : int
        Number of UIDs per Entrez summary request.
    retries : int
        Maximum Entrez retry attempts per batch.
    sleep_seconds : float
        Polite delay and retry backoff base in seconds.
    taxid : str, optional
        Query taxid for log/error context.
    query : str, optional
        Original Entrez query for log/error context.

    Returns
    -------
    list[dict[str, Any]]
        NCBI Assembly document summaries.
    """
    if batch_size <= 0:
        raise ValueError("batch_size must be a positive integer")
    summaries: list[dict[str, Any]] = []
    entrez = require_entrez()
    for start in range(0, len(assembly_uids), batch_size):
        batch = assembly_uids[start : start + batch_size]
        end = start + len(batch)
        LOGGER.info("Fetching assembly summaries %s-%s", start + 1, end)
        context = (
            "Entrez esummary "
            f"taxid={taxid or 'unknown'} query={query!r} "
            f"retstart={start} retmax={len(batch)} "
            f"batch_range={start + 1}-{end}"
        )
        record = entrez_read_retry(
            entrez.esummary,
            db="assembly",
            id=",".join(batch),
            report="full",
            retries=retries,
            sleep_seconds=sleep_seconds,
            context=context,
        )
        document_set = record.get("DocumentSummarySet", {})
        summaries.extend(document_set.get("DocumentSummary", []))
    return summaries

def extract_strain_name(summary: dict[str, Any]) -> str:
    """Extract a strain or isolate label from an Assembly summary."""
    biosource = summary.get("Biosource", {})
    infraspecies = biosource.get("InfraspeciesList", []) if biosource else []
    for entry in infraspecies:
        subtype = str(entry.get("Sub_type", "")).lower()
        if subtype in {"strain", "isolate", "cultivar", "ecotype"}:
            return str(entry.get("Sub_value", ""))
    if infraspecies:
        return str(infraspecies[0].get("Sub_value", ""))
    return ""


def choose_ftp_path(summary: dict[str, Any], source: str) -> tuple[str, str]:
    """Choose the best FTP path from an Assembly summary."""
    refseq_path = str(summary.get("FtpPath_RefSeq", "") or "")
    genbank_path = str(summary.get("FtpPath_GenBank", "") or "")
    if source == "refseq":
        return "refseq", refseq_path
    if source == "genbank":
        return "genbank", genbank_path
    if source == "prefer_refseq" and refseq_path:
        return "refseq", refseq_path
    if source == "prefer_refseq":
        return "genbank", genbank_path
    if source == "prefer_genbank" and genbank_path:
        return "genbank", genbank_path
    return "refseq", refseq_path


def summary_to_record(
    summary: dict[str, Any],
    query_taxid: str,
    assembly_uid: str,
    source: str,
    plan: TaxonPlan,
) -> AssemblyRecord:
    """Convert an NCBI Assembly summary into an AssemblyRecord."""
    selected_source, selected_path = choose_ftp_path(summary=summary, source=source)
    return AssemblyRecord(
        query_taxid=query_taxid,
        assembly_uid=assembly_uid,
        assembly_accession=str(summary.get("AssemblyAccession", "")),
        assembly_name=str(summary.get("AssemblyName", "")),
        organism_name=str(summary.get("Organism", "")),
        species_name=str(summary.get("SpeciesName", "")),
        species_taxid=str(summary.get("SpeciesTaxid", "")),
        taxid=str(summary.get("Taxid", "") or summary.get("SpeciesTaxid", "")),
        strain_name=extract_strain_name(summary),
        assembly_level=str(summary.get("AssemblyLevel", "")),
        refseq_category=str(summary.get("RefSeq_category", "")),
        scaffold_n50=parse_int(summary.get("ScaffoldN50", 0)),
        contig_n50=parse_int(summary.get("ContigN50", 0)),
        total_length=parse_int(summary.get("SeqLength", 0)),
        ftp_path_refseq=str(summary.get("FtpPath_RefSeq", "") or ""),
        ftp_path_genbank=str(summary.get("FtpPath_GenBank", "") or ""),
        selected_source=selected_source,
        selected_ftp_path=selected_path,
        role=plan.role,
        clade=plan.clade or plan.group_label,
        group_label=plan.group_label,
    )



def record_to_candidate_audit_row(
    record: AssemblyRecord,
    filter_status: str,
    filter_reason: str,
) -> dict[str, Any]:
    """Convert an assembly record to a pre-filter audit row.

    Parameters
    ----------
    record : AssemblyRecord
        Candidate assembly record.
    filter_status : str
        ``retained`` or ``excluded``.
    filter_reason : str
        Explicit reason for the final candidate decision.

    Returns
    -------
    dict[str, Any]
        TSV-ready candidate audit row.
    """
    suffix = FORMAT_SUFFIXES["genomic_fna"]
    return {
        "taxid": record.query_taxid,
        "query": f"txid{record.query_taxid}[Organism:exp]",
        "assembly_uid": record.assembly_uid,
        "assembly_accession": record.assembly_accession,
        "organism_name": record.organism_name,
        "species_name": record.species_name,
        "species_taxid": record.species_taxid,
        "assembly_level": record.assembly_level,
        "normalised_assembly_level": normalise_assembly_level(record.assembly_level),
        "refseq_category": record.refseq_category,
        "assembly_status": record.assembly_level,
        "ftp_path_refseq": record.ftp_path_refseq,
        "ftp_path_genbank": record.ftp_path_genbank,
        "selected_source": record.selected_source,
        "selected_path": record.selected_ftp_path,
        "normalised_download_url": normalise_download_url(
            record.selected_ftp_path,
            suffix=suffix,
        ),
        "scaffold_n50": record.scaffold_n50,
        "contig_n50": record.contig_n50,
        "total_length": record.total_length,
        "filter_status": filter_status,
        "filter_reason": filter_reason,
        "role": record.role,
        "clade": record.clade,
        "group_label": record.group_label,
    }


def evaluate_record_filter_reason(
    record: AssemblyRecord,
    assembly_levels: list[str] | None = None,
    include_unplaced: bool = False,
    min_total_length: int | None = None,
    max_total_length: int | None = None,
    min_scaffold_n50: int | None = None,
    min_contig_n50: int | None = None,
) -> str:
    """Return the first basic filter reason for an assembly record.

    A return value of ``retained`` means the record passes basic filtering.
    Later selection steps may still exclude it by best-per-species or maximum
    assembly caps.
    """
    allowed_levels = normalise_assembly_level_set(assembly_levels)
    record_level = normalise_assembly_level(record.assembly_level)
    if allowed_levels and record_level not in allowed_levels:
        return "excluded_by_assembly_level"
    if not include_unplaced and not record.selected_ftp_path:
        return "excluded_missing_download_path"
    if min_total_length is not None and record.total_length < min_total_length:
        return "excluded_by_assembly_quality"
    if max_total_length is not None and record.total_length > max_total_length:
        return "excluded_by_assembly_quality"
    if min_scaffold_n50 is not None and record.scaffold_n50 < min_scaffold_n50:
        return "excluded_by_assembly_quality"
    if min_contig_n50 is not None and record.contig_n50 < min_contig_n50:
        return "excluded_by_assembly_quality"
    return "retained"


def filter_records_with_audit(
    records: list[AssemblyRecord],
    assembly_levels: list[str] | None = None,
    include_unplaced: bool = False,
    min_total_length: int | None = None,
    max_total_length: int | None = None,
    min_scaffold_n50: int | None = None,
    min_contig_n50: int | None = None,
) -> tuple[list[AssemblyRecord], list[dict[str, Any]]]:
    """Filter assembly records and build candidate-audit rows.

    Parameters are equivalent to :func:`filter_records`. The audit table is
    written before download so zero-retained taxids still explain why records
    were excluded.
    """
    retained: list[AssemblyRecord] = []
    audit_rows: list[dict[str, Any]] = []
    for record in records:
        reason = evaluate_record_filter_reason(
            record=record,
            assembly_levels=assembly_levels,
            include_unplaced=include_unplaced,
            min_total_length=min_total_length,
            max_total_length=max_total_length,
            min_scaffold_n50=min_scaffold_n50,
            min_contig_n50=min_contig_n50,
        )
        status = "retained" if reason == "retained" else "excluded"
        audit_rows.append(
            record_to_candidate_audit_row(
                record=record,
                filter_status=status,
                filter_reason=reason,
            )
        )
        if reason == "retained":
            retained.append(record)
    return retained, audit_rows


def filter_records(
    records: list[AssemblyRecord],
    assembly_levels: list[str] | None = None,
    include_unplaced: bool = False,
    min_total_length: int | None = None,
    max_total_length: int | None = None,
    min_scaffold_n50: int | None = None,
    min_contig_n50: int | None = None,
) -> list[AssemblyRecord]:
    """Filter assembly records by level, downloadability, and quality."""
    retained, _audit_rows = filter_records_with_audit(
        records=records,
        assembly_levels=assembly_levels,
        include_unplaced=include_unplaced,
        min_total_length=min_total_length,
        max_total_length=max_total_length,
        min_scaffold_n50=min_scaffold_n50,
        min_contig_n50=min_contig_n50,
    )
    return retained

def sort_records_by_quality(records: list[AssemblyRecord]) -> list[AssemblyRecord]:
    """Sort records from highest to lowest assembly quality."""
    return sorted(
        records,
        key=lambda rec: (
            ASSEMBLY_LEVEL_RANK.get(normalise_assembly_level(rec.assembly_level), 0),
            rec.scaffold_n50,
            rec.contig_n50,
            rec.total_length,
            rec.assembly_accession,
        ),
        reverse=True,
    )


def select_best_per_species(
    records: list[AssemblyRecord],
    best_per_species: int | None,
) -> list[AssemblyRecord]:
    """Retain the best N assemblies per species taxid."""
    if not best_per_species or best_per_species <= 0:
        return records
    grouped: dict[str, list[AssemblyRecord]] = {}
    for record in records:
        key = record.species_taxid or record.species_name or record.taxid
        grouped.setdefault(key, []).append(record)
    selected: list[AssemblyRecord] = []
    for species_records in grouped.values():
        selected.extend(sort_records_by_quality(species_records)[:best_per_species])
    return sort_records_by_quality(selected)


def limit_records(
    records: list[AssemblyRecord],
    max_assemblies: int | None,
) -> list[AssemblyRecord]:
    """Limit records after quality sorting."""
    sorted_records = sort_records_by_quality(records)
    if not max_assemblies or max_assemblies <= 0:
        return sorted_records
    return sorted_records[:max_assemblies]



def ftp_url_to_file_url(ftp_path: str, suffix: str) -> str:
    """Convert an NCBI assembly FTP directory to a downloadable HTTPS URL."""
    return normalise_download_url(ftp_path, suffix=suffix)

def assembly_output_dir(out_dir: Path, record: AssemblyRecord) -> Path:
    """Return the organised output directory for one assembly."""
    taxid_label = safe_name(f"taxid_{record.query_taxid}_{record.group_label or record.clade}")
    species_label = safe_name(record.species_name or record.organism_name)
    accession_label = safe_name(record.assembly_accession)
    return out_dir / "genomes" / taxid_label / species_label / accession_label


def download_file(
    url: str,
    out_path: Path,
    force: bool,
    retries: int,
    sleep_seconds: float,
) -> str:
    """Download a URL to a local path with retries."""
    if not url:
        return "missing_url"
    if out_path.exists() and out_path.stat().st_size > 0 and not force:
        LOGGER.info("Keeping existing file: %s", out_path)
        return "exists"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    headers = {"User-Agent": "KmerSutra genome downloader"}
    last_error: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            LOGGER.info("Downloading %s", url)
            request = Request(url, headers=headers)
            with urlopen(request, timeout=60) as response:
                with out_path.open("wb") as output_handle:
                    shutil.copyfileobj(response, output_handle, length=1024 * 1024)
            return "downloaded"
        except (HTTPError, URLError, TimeoutError, OSError) as exc:
            last_error = exc
            LOGGER.warning(
                "Download failed on attempt %s/%s for %s: %s",
                attempt,
                retries,
                url,
                exc,
            )
            time.sleep(sleep_seconds * attempt)
    LOGGER.error("Download failed after %s attempts: %s", retries, url)
    if out_path.exists() and out_path.stat().st_size == 0:
        out_path.unlink()
    return f"failed:{last_error}"


def decompress_gzip(gzip_path: Path, out_path: Path, force: bool = False) -> str:
    """Decompress a gzip file to a target path."""
    if out_path.exists() and out_path.stat().st_size > 0 and not force:
        return "exists"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(gzip_path, "rb") as input_handle:
        with out_path.open("wb") as output_handle:
            shutil.copyfileobj(input_handle, output_handle, length=1024 * 1024)
    return "decompressed"


def download_record_files(
    record: AssemblyRecord,
    out_dir: Path,
    formats: list[str],
    force: bool,
    retries: int,
    sleep_seconds: float,
    decompress: bool,
    delete_archives: bool,
    metadata_only: bool,
    dry_run: bool,
) -> dict[str, str]:
    """Download requested files for one assembly and return path metadata."""
    record_dir = assembly_output_dir(out_dir=out_dir, record=record)
    result: dict[str, str] = {
        "record_dir": str(record_dir),
        "genome_fasta": "",
        "genome_fasta_gz": "",
        "download_status": "not_requested",
    }
    for fmt in formats:
        suffix = FORMAT_SUFFIXES[fmt]
        url = ftp_url_to_file_url(record.selected_ftp_path, suffix=suffix)
        result[f"{fmt}_url"] = url
        if not url:
            result[f"{fmt}_status"] = "missing_url"
            continue
        base_name = url.rstrip("/").split("/")[-1]
        archive_path = record_dir / base_name
        result[f"{fmt}_path"] = str(archive_path)
        if metadata_only or dry_run:
            result[f"{fmt}_status"] = "planned"
            if fmt == "genomic_fna":
                result["genome_fasta_gz"] = str(archive_path)
                result["genome_fasta"] = str(archive_path.with_suffix(""))
            continue
        status = download_file(
            url=url,
            out_path=archive_path,
            force=force,
            retries=retries,
            sleep_seconds=sleep_seconds,
        )
        result[f"{fmt}_status"] = status
        if fmt == "genomic_fna":
            if status.startswith("failed"):
                result["download_status"] = status
                result["genome_fasta_gz"] = ""
                result["genome_fasta"] = ""
                continue
            result["download_status"] = status
            result["genome_fasta_gz"] = str(archive_path)
            fasta_path = archive_path.with_suffix("")
            if decompress and archive_path.exists() and archive_path.stat().st_size > 0:
                decompress_status = decompress_gzip(
                    gzip_path=archive_path,
                    out_path=fasta_path,
                    force=force,
                )
                result["genomic_fna_decompress_status"] = decompress_status
                result["genome_fasta"] = str(fasta_path)
                if delete_archives and fasta_path.exists():
                    archive_path.unlink()
            else:
                result["genome_fasta"] = str(archive_path)
    return result


def record_to_metadata_row(record: AssemblyRecord, file_info: dict[str, str]) -> dict[str, Any]:
    """Convert one assembly record and file-info dictionary to TSV metadata."""
    row = {
        "query_taxid": record.query_taxid,
        "assembly_uid": record.assembly_uid,
        "assembly_accession": record.assembly_accession,
        "assembly_name": record.assembly_name,
        "organism_name": record.organism_name,
        "species_name": record.species_name,
        "species_taxid": record.species_taxid,
        "taxid": record.taxid,
        "strain_name": record.strain_name,
        "assembly_level": record.assembly_level,
        "refseq_category": record.refseq_category,
        "scaffold_n50": record.scaffold_n50,
        "contig_n50": record.contig_n50,
        "total_length": record.total_length,
        "ftp_path_refseq": record.ftp_path_refseq,
        "ftp_path_genbank": record.ftp_path_genbank,
        "selected_source": record.selected_source,
        "selected_ftp_path": record.selected_ftp_path,
        "role": record.role,
        "clade": record.clade,
        "group_label": record.group_label,
    }
    row.update(file_info)
    return row


def metadata_to_kmersutra_config_rows(
    metadata_rows: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    """Create KmerSutra genome-config rows from download metadata."""
    rows: list[dict[str, Any]] = []
    for row in metadata_rows:
        fasta_path = row.get("genome_fasta", "")
        if not fasta_path:
            continue
        rows.append(
            {
                "genome_fasta": fasta_path,
                "species_name": row.get("species_name", "") or row.get("organism_name", ""),
                "strain_name": row.get("strain_name", "") or row.get("assembly_name", ""),
                "taxid": row.get("species_taxid", "") or row.get("taxid", ""),
                "role": row.get("role", "downloaded"),
                "clade": row.get("clade", "") or row.get("group_label", ""),
                "assembly_accession": row.get("assembly_accession", ""),
                "query_taxid": row.get("query_taxid", ""),
                "assembly_level": row.get("assembly_level", ""),
                "scaffold_n50": row.get("scaffold_n50", ""),
                "contig_n50": row.get("contig_n50", ""),
            }
        )
    return rows


def write_run_config(path: Path, args: argparse.Namespace, plans: list[TaxonPlan]) -> None:
    """Write a JSON run configuration for reproducibility."""
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "command": sys.argv,
        "arguments": vars(args),
        "taxon_plans": [plan.__dict__ for plan in plans],
        "created_at": time.strftime("%Y-%m-%d %H:%M:%S"),
    }
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def load_custom_genome_config_rows(
    path: Path,
    allow_missing_fastas: bool = False,
) -> list[dict[str, Any]]:
    """Load and validate user-supplied KmerSutra genome-config rows.

    Parameters
    ----------
    path : Path
        TSV file with KmerSutra genome-config columns.
    allow_missing_fastas : bool, optional
        If true, keep rows whose FASTA paths do not exist. This is intended
        only for dry-run/planning workflows.

    Returns
    -------
    list[dict[str, Any]]
        Normalised genome-config rows ready to append to downloaded rows.

    Raises
    ------
    ValueError
        Raised if required columns are missing or FASTA paths are unusable.
    """
    required = set(KMERSUTRA_CONFIG_HEADER[:4])
    rows = read_tsv(path)
    if not rows:
        LOGGER.warning("Custom genome config %s contains no rows", path)
        return []
    missing = required.difference(rows[0])
    if missing:
        raise ValueError(
            f"Custom genome config {path} is missing required columns: "
            f"{', '.join(sorted(missing))}"
        )
    normalised_rows: list[dict[str, Any]] = []
    for index, row in enumerate(rows, start=2):
        fasta = Path(str(row.get("genome_fasta", "") or "")).expanduser()
        if not str(fasta):
            raise ValueError(f"Custom genome config {path} line {index} has blank genome_fasta")
        if not allow_missing_fastas and (not fasta.exists() or fasta.stat().st_size <= 0):
            raise ValueError(
                f"Custom genome FASTA is missing or empty at line {index}: {fasta}"
            )
        normalised = {column: row.get(column, "") for column in KMERSUTRA_CONFIG_HEADER}
        normalised["genome_fasta"] = str(fasta)
        normalised.setdefault("role", "custom")
        normalised_rows.append(normalised)
    LOGGER.info("Loaded %s custom genome-config row(s) from %s", len(normalised_rows), path)
    return normalised_rows


def custom_rows_to_metadata_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Convert custom genome-config rows to metadata-like rows."""
    metadata_rows: list[dict[str, Any]] = []
    for row in rows:
        metadata_rows.append(
            {
                "query_taxid": row.get("query_taxid", row.get("taxid", "")),
                "assembly_uid": "custom",
                "assembly_accession": row.get("assembly_accession", "custom"),
                "assembly_name": row.get("strain_name", "custom"),
                "organism_name": row.get("species_name", ""),
                "species_name": row.get("species_name", ""),
                "species_taxid": row.get("taxid", ""),
                "taxid": row.get("taxid", ""),
                "strain_name": row.get("strain_name", ""),
                "assembly_level": row.get("assembly_level", "custom"),
                "refseq_category": "custom",
                "scaffold_n50": row.get("scaffold_n50", ""),
                "contig_n50": row.get("contig_n50", ""),
                "total_length": "",
                "selected_source": "custom",
                "selected_ftp_path": "",
                "role": row.get("role", "custom"),
                "clade": row.get("clade", ""),
                "group_label": "custom",
                "record_dir": str(Path(row.get("genome_fasta", "")).parent),
                "genome_fasta": row.get("genome_fasta", ""),
                "genome_fasta_gz": "",
                "genomic_fna_url": "custom",
                "genomic_fna_status": "custom",
                "genomic_fna_decompress_status": "not_required",
                "ftp_path_refseq": "",
                "ftp_path_genbank": "",
            }
        )
    return metadata_rows



def collect_records_for_plan(
    plan: TaxonPlan,
    args: argparse.Namespace,
    excluded_accessions: set[str] | None = None,
) -> tuple[list[AssemblyRecord], list[dict[str, Any]]]:
    """Search, fetch, convert, filter, and select records for one taxon plan.

    Returns both retained records and a full candidate audit. The audit is
    emitted even when no assemblies survive filtering, so broad NCBI searches
    and source/quality filters can be diagnosed without rerunning Entrez.
    """
    query = f"txid{plan.taxid}[Organism:exp]"
    assembly_uids = search_assembly_uids(
        taxid=plan.taxid,
        retries=args.retries,
        sleep_seconds=args.sleep_seconds,
        max_records=args.max_entrez_records_per_taxid,
    )
    if not assembly_uids:
        LOGGER.warning("No assemblies found for taxid %s", plan.taxid)
        return [], []
    summaries = fetch_assembly_summaries(
        assembly_uids=assembly_uids,
        batch_size=args.batch_size,
        retries=args.retries,
        sleep_seconds=args.sleep_seconds,
        taxid=plan.taxid,
        query=query,
    )
    records = [
        summary_to_record(
            summary=summary,
            query_taxid=plan.taxid,
            assembly_uid=assembly_uids[index] if index < len(assembly_uids) else "",
            source=args.source,
            plan=plan,
        )
        for index, summary in enumerate(summaries)
    ]
    LOGGER.info("Taxid %s converted %s assembly summary records", plan.taxid, len(records))
    records, audit_rows = filter_records_with_audit(
        records=records,
        assembly_levels=args.assembly_levels,
        include_unplaced=args.include_unplaced,
        min_total_length=plan.min_total_length or args.min_total_length,
        max_total_length=plan.max_total_length or args.max_total_length,
        min_scaffold_n50=plan.min_scaffold_n50 or args.min_scaffold_n50,
        min_contig_n50=plan.min_contig_n50 or args.min_contig_n50,
    )
    records = exclude_records_by_accession(
        records=records,
        audit_rows=audit_rows,
        excluded_accessions=excluded_accessions or set(),
    )
    best_n = plan.best_per_species or args.best_per_species
    records_after_best = select_best_per_species(records=records, best_per_species=best_n)
    best_selected = {record.assembly_accession for record in records_after_best}
    for row in audit_rows:
        if row.get("filter_reason") == "retained" and row.get("assembly_accession") not in best_selected:
            row["filter_status"] = "excluded"
            row["filter_reason"] = "excluded_by_best_per_species"
    max_n = plan.max_assemblies or args.max_assemblies_per_taxid
    final_records = limit_records(records=records_after_best, max_assemblies=max_n)
    final_selected = {record.assembly_accession for record in final_records}
    for row in audit_rows:
        if row.get("filter_reason") == "retained" and row.get("assembly_accession") not in final_selected:
            row["filter_status"] = "excluded"
            row["filter_reason"] = "excluded_by_max_assemblies"
    reason_counts = count_values(audit_rows, "filter_reason")
    LOGGER.info("Taxid %s retained %s assemblies", plan.taxid, len(final_records))
    LOGGER.info(
        "Taxid %s filter reason counts: %s",
        plan.taxid,
        ", ".join(f"{key}={value}" for key, value in sorted(reason_counts.items())),
    )
    return final_records, audit_rows


def main(argv: list[str] | None = None) -> int:
    """Run the NCBI genome downloader."""
    args = parse_args(argv)
    out_dir = Path(args.out_dir).resolve()
    log_file = Path(args.log_file) if args.log_file else out_dir / "logs" / "download.log"
    configure_logging(log_file=log_file, verbose=args.verbose)

    plans = build_taxon_plan(args)
    if plans:
        if not args.email:
            raise ValueError("--email is required when downloading/querying NCBI Entrez")
        configure_entrez(email=args.email, api_key=args.api_key)

    out_dir.mkdir(parents=True, exist_ok=True)
    write_run_config(path=out_dir / "run_config.json", args=args, plans=plans)
    excluded_accessions = load_excluded_accessions(args.exclude_accessions_table)

    all_records: list[AssemblyRecord] = []
    candidate_audit_rows: list[dict[str, Any]] = []
    for plan in plans:
        LOGGER.info("Processing taxid %s", plan.taxid)
        retained_records, audit_rows = collect_records_for_plan(
            plan=plan,
            args=args,
            excluded_accessions=excluded_accessions,
        )
        all_records.extend(retained_records)
        candidate_audit_rows.extend(audit_rows)

    candidate_audit_header = [
        "taxid",
        "query",
        "assembly_uid",
        "assembly_accession",
        "organism_name",
        "species_name",
        "species_taxid",
        "assembly_level",
        "normalised_assembly_level",
        "refseq_category",
        "assembly_status",
        "ftp_path_refseq",
        "ftp_path_genbank",
        "selected_source",
        "selected_path",
        "normalised_download_url",
        "scaffold_n50",
        "contig_n50",
        "total_length",
        "filter_status",
        "filter_reason",
        "role",
        "clade",
        "group_label",
    ]
    metadata_rows: list[dict[str, Any]] = []
    for record in all_records:
        file_info = download_record_files(
            record=record,
            out_dir=out_dir,
            formats=args.formats,
            force=args.force,
            retries=args.retries,
            sleep_seconds=args.sleep_seconds,
            decompress=args.decompress,
            delete_archives=args.delete_archives,
            metadata_only=args.metadata_only,
            dry_run=args.dry_run,
        )
        if str(file_info.get("genomic_fna_status", "")).startswith("failed"):
            for audit_row in candidate_audit_rows:
                if audit_row.get("assembly_accession") == record.assembly_accession:
                    audit_row["filter_status"] = "excluded"
                    audit_row["filter_reason"] = "download_failed"
        metadata_row = record_to_metadata_row(record=record, file_info=file_info)
        genome_fasta = str(metadata_row.get("genome_fasta", "") or "")
        matching_accessions: set[str] = set()
        if genome_fasta and Path(genome_fasta).is_file():
            matching_accessions = fasta_excluded_accessions(
                fasta_path=genome_fasta,
                excluded_accessions=excluded_accessions,
            )
        if matching_accessions:
            LOGGER.warning(
                "Excluding %s after FASTA-header truth-accession match: %s",
                record.assembly_accession,
                ", ".join(sorted(matching_accessions)),
            )
            for audit_row in candidate_audit_rows:
                if audit_row.get("assembly_accession") == record.assembly_accession:
                    audit_row["filter_status"] = "excluded"
                    audit_row["filter_reason"] = "excluded_by_sequence_accession_table"
            continue
        metadata_rows.append(metadata_row)

    custom_config_rows: list[dict[str, Any]] = []
    if args.custom_genome_config:
        custom_config_rows = load_custom_genome_config_rows(
            path=Path(args.custom_genome_config),
            allow_missing_fastas=args.allow_missing_custom_fastas,
        )
        metadata_rows.extend(custom_rows_to_metadata_rows(custom_config_rows))

    write_tsv(
        path=out_dir / "ncbi_assembly_candidate_audit.tsv",
        rows=candidate_audit_rows,
        header=candidate_audit_header,
    )

    metadata_header = [
        "query_taxid",
        "assembly_uid",
        "assembly_accession",
        "assembly_name",
        "organism_name",
        "species_name",
        "species_taxid",
        "taxid",
        "strain_name",
        "assembly_level",
        "refseq_category",
        "scaffold_n50",
        "contig_n50",
        "total_length",
        "selected_source",
        "selected_ftp_path",
        "role",
        "clade",
        "group_label",
        "record_dir",
        "genome_fasta",
        "genome_fasta_gz",
        "genomic_fna_url",
        "genomic_fna_status",
        "genomic_fna_decompress_status",
        "ftp_path_refseq",
        "ftp_path_genbank",
    ]
    write_tsv(
        path=out_dir / "ncbi_download_metadata.tsv",
        rows=metadata_rows,
        header=metadata_header,
    )

    kmersutra_rows = metadata_to_kmersutra_config_rows(metadata_rows)
    write_tsv(
        path=out_dir / "kmersutra_genome_config.tsv",
        rows=kmersutra_rows,
        header=KMERSUTRA_CONFIG_HEADER,
    )

    query_summary_rows = []
    for plan in plans:
        retained = [row for row in metadata_rows if row.get("query_taxid") == plan.taxid]
        audit_for_taxid = [
            row for row in candidate_audit_rows if row.get("taxid") == plan.taxid
        ]
        reason_counts = count_values(audit_for_taxid, "filter_reason")
        query_summary_rows.append(
            {
                "query_taxid": plan.taxid,
                "role": plan.role,
                "clade": plan.clade,
                "group_label": plan.group_label,
                "n_candidate_assemblies": len(audit_for_taxid),
                "n_retained_assemblies": len(retained),
                "n_species": len({row.get("species_taxid") for row in retained}),
                "filter_reason_counts": ";".join(
                    f"{key}={value}" for key, value in sorted(reason_counts.items())
                ),
            }
        )
    if custom_config_rows:
        query_summary_rows.append(
            {
                "query_taxid": "custom",
                "role": "custom",
                "clade": "custom",
                "group_label": "custom",
                "n_candidate_assemblies": len(custom_config_rows),
                "n_retained_assemblies": len(custom_config_rows),
                "n_species": len({row.get("taxid") for row in custom_config_rows}),
                "filter_reason_counts": f"custom={len(custom_config_rows)}",
            }
        )
    write_tsv(
        path=out_dir / "query_summary.tsv",
        rows=query_summary_rows,
        header=[
            "query_taxid",
            "role",
            "clade",
            "group_label",
            "n_candidate_assemblies",
            "n_retained_assemblies",
            "n_species",
            "filter_reason_counts",
        ],
    )

    LOGGER.info("Wrote candidate audit rows: %s", len(candidate_audit_rows))
    LOGGER.info("Wrote metadata rows: %s", len(metadata_rows))
    LOGGER.info("Wrote KmerSutra config rows: %s", len(kmersutra_rows))
    LOGGER.info("Done")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
