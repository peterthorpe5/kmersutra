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
import json
import logging
import re
import shutil
import sys
import time
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
        "--out_dir",
        required=True,
        help="Output directory for organised genome downloads and metadata.",
    )
    parser.add_argument(
        "--email",
        required=True,
        help="Email address supplied to NCBI Entrez.",
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
        choices=("complete genome", "chromosome", "scaffold", "contig"),
        help="Optional assembly-level filter. If omitted, all levels are kept.",
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
        default=500,
        help="Number of assembly records retrieved per Entrez summary batch.",
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
    LOGGER.handlers.clear()
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
    if not plans:
        raise ValueError("At least one taxid must be supplied via --taxids or --taxid_plan")
    return plans


def parse_optional_int(value: Any) -> int | None:
    """Parse an optional integer from TSV configuration text."""
    if value in (None, ""):
        return None
    return parse_int(value, default=0)


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
        except (HTTPError, URLError, TimeoutError, OSError) as exc:
            last_error = exc
            wait_time = sleep_seconds * attempt
            LOGGER.warning(
                "Entrez request failed on attempt %s/%s: %s", attempt, retries, exc
            )
            time.sleep(wait_time)
    raise RuntimeError(f"Entrez request failed after {retries} attempts") from last_error


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
) -> list[str]:
    """Search NCBI Assembly for all assemblies under a taxonomy subtree."""
    query = f"txid{taxid}[Organism:exp]"
    LOGGER.info("Searching NCBI Assembly with query: %s", query)
    entrez = require_entrez()
    handle = entrez_retry(
        entrez.esearch,
        db="assembly",
        term=query,
        retmax=0,
        usehistory="y",
        retries=retries,
        sleep_seconds=sleep_seconds,
    )
    record = entrez.read(handle, validate=False)
    count = int(record.get("Count", 0))
    LOGGER.info("Taxid %s returned %s assembly records", taxid, count)
    if count == 0:
        return []
    entrez = require_entrez()
    handle = entrez_retry(
        entrez.esearch,
        db="assembly",
        term=query,
        retmax=count,
        retries=retries,
        sleep_seconds=sleep_seconds,
    )
    record = entrez.read(handle, validate=False)
    return list(record.get("IdList", []))


def fetch_assembly_summaries(
    assembly_uids: list[str],
    batch_size: int,
    retries: int,
    sleep_seconds: float,
) -> list[dict[str, Any]]:
    """Fetch full NCBI Assembly summaries in batches."""
    summaries: list[dict[str, Any]] = []
    for start in range(0, len(assembly_uids), batch_size):
        batch = assembly_uids[start : start + batch_size]
        LOGGER.info("Fetching assembly summaries %s-%s", start + 1, start + len(batch))
        entrez = require_entrez()
        handle = entrez_retry(
            entrez.esummary,
            db="assembly",
            id=",".join(batch),
            report="full",
            retries=retries,
            sleep_seconds=sleep_seconds,
        )
        record = entrez.read(handle, validate=False)
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


def filter_records(
    records: list[AssemblyRecord],
    assembly_levels: list[str] | None = None,
    include_unplaced: bool = False,
    min_total_length: int | None = None,
    max_total_length: int | None = None,
    min_scaffold_n50: int | None = None,
    min_contig_n50: int | None = None,
) -> list[AssemblyRecord]:
    """Filter assembly records by level, downloadability, and quality.

    Parameters
    ----------
    records : list[AssemblyRecord]
        Candidate assembly records.
    assembly_levels : list[str] | None
        Optional retained assembly levels.
    include_unplaced : bool
        Whether to retain records without a selected FTP path.
    min_total_length : int | None
        Optional minimum total assembly length.
    max_total_length : int | None
        Optional maximum total assembly length.
    min_scaffold_n50 : int | None
        Optional minimum scaffold N50.
    min_contig_n50 : int | None
        Optional minimum contig N50.

    Returns
    -------
    list[AssemblyRecord]
        Records passing all filters.
    """
    filtered: list[AssemblyRecord] = []
    allowed_levels = {level.lower() for level in assembly_levels or []}
    for record in records:
        if allowed_levels and record.assembly_level.lower() not in allowed_levels:
            continue
        if not include_unplaced and not record.selected_ftp_path:
            continue
        if min_total_length is not None and record.total_length < min_total_length:
            continue
        if max_total_length is not None and record.total_length > max_total_length:
            continue
        if min_scaffold_n50 is not None and record.scaffold_n50 < min_scaffold_n50:
            continue
        if min_contig_n50 is not None and record.contig_n50 < min_contig_n50:
            continue
        filtered.append(record)
    return filtered


def sort_records_by_quality(records: list[AssemblyRecord]) -> list[AssemblyRecord]:
    """Sort records from highest to lowest assembly quality."""
    return sorted(
        records,
        key=lambda rec: (
            ASSEMBLY_LEVEL_RANK.get(rec.assembly_level.lower(), 0),
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
    """Convert an NCBI assembly FTP directory to a downloadable file URL."""
    if not ftp_path:
        return ""
    base_url = ftp_path.replace("ftp://", "https://")
    base_name = base_url.rstrip("/").split("/")[-1]
    return f"{base_url}/{base_name}_{suffix}"


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


def collect_records_for_plan(
    plan: TaxonPlan,
    args: argparse.Namespace,
) -> list[AssemblyRecord]:
    """Search, fetch, convert, filter, and select records for one taxon plan."""
    assembly_uids = search_assembly_uids(
        taxid=plan.taxid,
        retries=args.retries,
        sleep_seconds=args.sleep_seconds,
    )
    if not assembly_uids:
        LOGGER.warning("No assemblies found for taxid %s", plan.taxid)
        return []
    summaries = fetch_assembly_summaries(
        assembly_uids=assembly_uids,
        batch_size=args.batch_size,
        retries=args.retries,
        sleep_seconds=args.sleep_seconds,
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
    records = filter_records(
        records=records,
        assembly_levels=args.assembly_levels,
        include_unplaced=args.include_unplaced,
        min_total_length=plan.min_total_length or args.min_total_length,
        max_total_length=plan.max_total_length or args.max_total_length,
        min_scaffold_n50=plan.min_scaffold_n50 or args.min_scaffold_n50,
        min_contig_n50=plan.min_contig_n50 or args.min_contig_n50,
    )
    best_n = plan.best_per_species or args.best_per_species
    records = select_best_per_species(records=records, best_per_species=best_n)
    max_n = plan.max_assemblies or args.max_assemblies_per_taxid
    records = limit_records(records=records, max_assemblies=max_n)
    LOGGER.info("Taxid %s retained %s assemblies", plan.taxid, len(records))
    return records


def main(argv: list[str] | None = None) -> int:
    """Run the NCBI genome downloader."""
    args = parse_args(argv)
    out_dir = Path(args.out_dir).resolve()
    log_file = Path(args.log_file) if args.log_file else out_dir / "logs" / "download.log"
    configure_logging(log_file=log_file, verbose=args.verbose)
    configure_entrez(email=args.email, api_key=args.api_key)

    plans = build_taxon_plan(args)
    out_dir.mkdir(parents=True, exist_ok=True)
    write_run_config(path=out_dir / "run_config.json", args=args, plans=plans)

    all_records: list[AssemblyRecord] = []
    for plan in plans:
        LOGGER.info("Processing taxid %s", plan.taxid)
        all_records.extend(collect_records_for_plan(plan=plan, args=args))

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
        metadata_rows.append(record_to_metadata_row(record=record, file_info=file_info))

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
    kmersutra_header = [
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
    write_tsv(
        path=out_dir / "kmersutra_genome_config.tsv",
        rows=kmersutra_rows,
        header=kmersutra_header,
    )

    query_summary_rows = []
    for plan in plans:
        retained = [row for row in metadata_rows if row.get("query_taxid") == plan.taxid]
        query_summary_rows.append(
            {
                "query_taxid": plan.taxid,
                "role": plan.role,
                "clade": plan.clade,
                "group_label": plan.group_label,
                "n_retained_assemblies": len(retained),
                "n_species": len({row.get("species_taxid") for row in retained}),
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
            "n_retained_assemblies",
            "n_species",
        ],
    )

    LOGGER.info("Wrote metadata rows: %s", len(metadata_rows))
    LOGGER.info("Wrote KmerSutra config rows: %s", len(kmersutra_rows))
    LOGGER.info("Done")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
