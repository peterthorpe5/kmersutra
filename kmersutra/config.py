"""Configuration parsing for KmerSutra genome panels."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from kmersutra.io import read_tsv

VALID_ROLES = {
    "target_species",
    "target_clade_member",
    "target",
    "near_neighbour",
    "near_neighbor",
    "outgroup",
    "apicomplexan_outgroup",
    "distant_outgroup",
    "host",
    "host_or_background",
    "background_host",
    "host_background",
    "background",
    "background_pathogen",
    "environmental_background",
    "non_target",
    "downloaded",
    "exclude",
}
TARGET_ROLES = {"target_species", "target", "target_clade_member"}
NON_TARGET_ROLES = VALID_ROLES - TARGET_ROLES - {"exclude"}


@dataclass(frozen=True)
class GenomeConfig:
    """Genome metadata used to build a k-mer panel.

    Attributes
    ----------
    genome_fasta : pathlib.Path
        Path to genome FASTA.
    species_name : str
        Species or taxon label.
    strain_name : str
        Optional strain label.
    taxid : str
        Optional NCBI taxonomy identifier.
    role : str
        Biological role in the panel. Roles are used to distinguish target
        species from near-neighbours, outgroups, host/background records, and
        other non-target genomes.
    clade : str
        Optional clade/group label.
    assembly_accession : str
        Optional assembly accession.
    source : str
        Optional source database label.
    """

    genome_fasta: Path
    species_name: str
    strain_name: str = ""
    taxid: str = ""
    role: str = "target_species"
    clade: str = ""
    assembly_accession: str = ""
    source: str = ""

    @property
    def genome_id(self) -> str:
        """Return a stable genome identifier.

        Returns
        -------
        str
            Assembly accession if available, otherwise FASTA stem.
        """
        return self.assembly_accession or self.genome_fasta.stem.replace(".fna", "")

    @property
    def is_target(self) -> bool:
        """Return whether this genome represents a target species.

        Returns
        -------
        bool
            True for target roles.
        """
        return self.role in TARGET_ROLES


def load_genome_config(*, config_path: str | Path, require_target: bool = True) -> list[GenomeConfig]:
    """Load genome configuration records from a TSV file.

    Parameters
    ----------
    config_path : str or pathlib.Path
        Path to genome configuration TSV.
    require_target : bool, optional
        Require at least one record with a target role. Set to False for
        query-agnostic all-candidate panel builds.

    Returns
    -------
    list[GenomeConfig]
        Parsed genome records.
    """
    records = read_tsv(input_path=config_path)
    if not records:
        raise ValueError("Genome configuration contains no records")

    required = {"genome_fasta", "species_name", "role"}
    missing = required - set(records[0].keys())
    if missing:
        raise ValueError(
            "Genome configuration is missing required columns: "
            + ", ".join(sorted(missing))
        )

    configs: list[GenomeConfig] = []
    for record in records:
        role = record.get("role", "target_species")
        if role not in VALID_ROLES:
            raise ValueError(f"Unsupported genome role: {role}")
        if role == "exclude":
            continue
        configs.append(
            GenomeConfig(
                genome_fasta=Path(record["genome_fasta"]),
                species_name=record["species_name"],
                strain_name=record.get("strain_name", ""),
                taxid=record.get("taxid", ""),
                role=role,
                clade=record.get("clade", ""),
                assembly_accession=record.get("assembly_accession", ""),
                source=record.get("source", ""),
            )
        )
    if require_target and not any(config.is_target for config in configs):
        raise ValueError("At least one target species is required")
    return configs
