#!/usr/bin/env python3
"""Download NCBI assemblies for one or more taxonomy IDs.

This script is a KmerSutra-facing replacement for older taxon download scripts.
It writes a genome configuration table that can be passed directly to
``build_clade_kmer_panel.py``.
"""

from __future__ import annotations

import argparse
import gzip
import shutil
import time
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import urlopen

from kmersutra.io import write_tsv
from kmersutra.logging_utils import configure_logging
from kmersutra.ncbi import AssemblyRecord, build_ncbi_output_paths


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed arguments.
    """
    parser = argparse.ArgumentParser(
        description="Download NCBI assemblies under taxon IDs for KmerSutra."
    )
    parser.add_argument("--taxids", nargs="+", required=True)
    parser.add_argument("--email", required=True)
    parser.add_argument("--out_dir", required=True)
    parser.add_argument("--role", default="target_clade_member")
    parser.add_argument("--clade", default="")
    parser.add_argument("--retmax", type=int, default=100000)
    parser.add_argument("--timeout", type=int, default=30)
    parser.add_argument("--sleep_seconds", type=float, default=0.34)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--verbose", action="store_true")
    return parser.parse_args()


def _require_biopython() -> object:
    """Import Bio.Entrez lazily.

    Returns
    -------
    object
        Bio.Entrez module.
    """
    try:
        from Bio import Entrez  # type: ignore
    except ImportError as exc:
        raise ImportError(
            "Biopython is required for NCBI downloads. Install biopython or "
            "prepare genome_config.tsv manually."
        ) from exc
    return Entrez


def _assembly_url(*, ftp_path: str) -> str:
    """Build the genomic FASTA URL from an NCBI FTP path.

    Parameters
    ----------
    ftp_path : str
        Assembly FTP path from NCBI esummary.

    Returns
    -------
    str
        URL to genomic FASTA gzip.
    """
    stem = ftp_path.rstrip("/").split("/")[-1]
    return f"{ftp_path}/{stem}_genomic.fna.gz".replace("ftp://", "https://")


def search_assembly_ids(*, entrez: object, taxid: str, retmax: int) -> list[str]:
    """Search NCBI Assembly for assemblies under a taxon.

    Parameters
    ----------
    entrez : object
        Bio.Entrez module.
    taxid : str
        NCBI taxonomy identifier.
    retmax : int
        Maximum records to return.

    Returns
    -------
    list[str]
        Assembly UIDs.
    """
    query = f"txid{taxid}[Organism:exp]"
    handle = entrez.esearch(db="assembly", term=query, retmax=retmax)
    record = entrez.read(handle)
    return list(record.get("IdList", []))


def fetch_assembly_summary(*, entrez: object, assembly_uid: str) -> dict[str, object]:
    """Fetch one NCBI Assembly summary.

    Parameters
    ----------
    entrez : object
        Bio.Entrez module.
    assembly_uid : str
        Assembly UID.

    Returns
    -------
    dict[str, object]
        Assembly summary record.
    """
    handle = entrez.esummary(db="assembly", id=assembly_uid, report="full")
    summary = entrez.read(handle)
    return dict(summary["DocumentSummarySet"]["DocumentSummary"][0])


def download_and_extract_fasta(
    *,
    url: str,
    output_path: Path,
    timeout: int,
    force: bool,
) -> None:
    """Download and extract a gzipped FASTA file.

    Parameters
    ----------
    url : str
        URL to ``genomic.fna.gz``.
    output_path : pathlib.Path
        Output uncompressed FASTA path.
    timeout : int
        URL timeout in seconds.
    force : bool
        If true, overwrite existing output.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if output_path.exists() and not force:
        return
    gz_path = output_path.with_suffix(output_path.suffix + ".gz")
    with urlopen(url, timeout=timeout) as response, gz_path.open("wb") as handle:
        shutil.copyfileobj(response, handle)
    with gzip.open(gz_path, "rb") as source, output_path.open("wb") as target:
        shutil.copyfileobj(source, target)


def main() -> None:
    """Run the NCBI taxon downloader."""
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = build_ncbi_output_paths(out_dir=out_dir)
    logger = configure_logging(log_file=paths["log"], verbose=args.verbose)
    entrez = _require_biopython()
    entrez.email = args.email
    entrez.tool = "KmerSutra"

    manifest: list[dict[str, object]] = []
    genome_config: list[dict[str, object]] = []
    failures: list[dict[str, object]] = []

    for taxid in args.taxids:
        logger.info("Searching NCBI Assembly for taxid %s", taxid)
        assembly_ids = search_assembly_ids(entrez=entrez, taxid=taxid, retmax=args.retmax)
        logger.info("Found %d assemblies for taxid %s", len(assembly_ids), taxid)
        for assembly_uid in assembly_ids:
            try:
                summary = fetch_assembly_summary(entrez=entrez, assembly_uid=assembly_uid)
                accession = str(summary.get("AssemblyAccession", ""))
                species_name = str(summary.get("SpeciesName", summary.get("Organism", "")))
                species_taxid = str(summary.get("SpeciesTaxid", ""))
                ftp_path = str(summary.get("FtpPath_RefSeq") or summary.get("FtpPath_GenBank") or "")
                if not ftp_path:
                    raise ValueError("No RefSeq or GenBank FTP path available")
                strain = ""
                biosource = summary.get("Biosource", {})
                if isinstance(biosource, dict):
                    infra = biosource.get("InfraspeciesList", [])
                    if infra:
                        strain = str(infra[0].get("Sub_value", ""))
                fasta_name = f"{accession}_{species_name.replace(' ', '_')}.fna"
                fasta_path = out_dir / "fastas" / fasta_name
                download_and_extract_fasta(
                    url=_assembly_url(ftp_path=ftp_path),
                    output_path=fasta_path,
                    timeout=args.timeout,
                    force=args.force,
                )
                record = AssemblyRecord(
                    query_taxid=taxid,
                    assembly_uid=assembly_uid,
                    assembly_accession=accession,
                    species_name=species_name,
                    strain_name=strain,
                    taxid=species_taxid,
                    fasta_path=str(fasta_path),
                    role=args.role,
                    clade=args.clade,
                    status="downloaded",
                )
                manifest.append(record.to_manifest_record())
                genome_config.append(record.to_genome_config_record())
                logger.info("Downloaded %s", accession)
            except (HTTPError, URLError, TimeoutError, ValueError, OSError) as exc:
                logger.warning("Failed assembly UID %s: %s", assembly_uid, exc)
                failures.append(
                    {
                        "query_taxid": taxid,
                        "assembly_uid": assembly_uid,
                        "reason": str(exc),
                    }
                )
            time.sleep(args.sleep_seconds)

    write_tsv(records=manifest, output_path=paths["download_manifest"], fieldnames=[
        "query_taxid", "assembly_uid", "assembly_accession", "species_name",
        "strain_name", "taxid", "fasta_path", "role", "clade", "status",
    ])
    write_tsv(records=genome_config, output_path=paths["genome_config"], fieldnames=[
        "genome_fasta", "species_name", "strain_name", "taxid",
        "assembly_accession", "role", "clade", "source",
    ])
    write_tsv(records=failures, output_path=paths["download_failures"], fieldnames=[
        "query_taxid", "assembly_uid", "reason",
    ])
    logger.info("Done")


if __name__ == "__main__":
    main()
