"""NCBI download helpers for KmerSutra.

This module provides a modernised interface intended to replace legacy taxon
assembly download scripts. Network-dependent functions are deliberately small
so they can be mocked in tests.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class AssemblyRecord:
    """Metadata for one downloaded or attempted NCBI assembly.

    Attributes
    ----------
    query_taxid : str
        Taxonomy identifier used for the search.
    assembly_uid : str
        NCBI assembly UID.
    assembly_accession : str
        Assembly accession.
    species_name : str
        Species name.
    strain_name : str
        Strain name if available.
    taxid : str
        Species taxid if available.
    fasta_path : str
        Local FASTA path.
    role : str
        Intended KmerSutra role.
    clade : str
        Clade label.
    status : str
        Download status.
    """

    query_taxid: str
    assembly_uid: str
    assembly_accession: str
    species_name: str
    strain_name: str
    taxid: str
    fasta_path: str
    role: str
    clade: str
    status: str

    def to_manifest_record(self) -> dict[str, object]:
        """Return a download manifest row.

        Returns
        -------
        dict[str, object]
            Manifest record.
        """
        return self.__dict__.copy()

    def to_genome_config_record(self) -> dict[str, object]:
        """Return a genome configuration row for KmerSutra.

        Returns
        -------
        dict[str, object]
            Genome configuration record.
        """
        return {
            "genome_fasta": self.fasta_path,
            "species_name": self.species_name,
            "strain_name": self.strain_name,
            "taxid": self.taxid,
            "assembly_accession": self.assembly_accession,
            "role": self.role,
            "clade": self.clade,
            "source": "NCBI",
        }


def build_ncbi_output_paths(*, out_dir: str | Path) -> dict[str, Path]:
    """Build standard NCBI downloader output paths.

    Parameters
    ----------
    out_dir : str or pathlib.Path
        Output directory.

    Returns
    -------
    dict[str, pathlib.Path]
        Named output paths.
    """
    root = Path(out_dir)
    return {
        "genome_config": root / "genome_config.tsv",
        "download_manifest": root / "download_manifest.tsv",
        "download_failures": root / "download_failures.tsv",
        "log": root / "download_ncbi_taxon_genomes.log",
    }
