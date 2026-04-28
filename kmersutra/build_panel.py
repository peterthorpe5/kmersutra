"""Build clade-aware and species-aware diagnostic k-mer panels."""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

from kmersutra.config import GenomeConfig
from kmersutra.fasta import read_fasta_records
from kmersutra.kmers import iter_kmers


@dataclass(frozen=True)
class KmerObservation:
    """A k-mer occurrence in a reference genome.

    Attributes
    ----------
    kmer : str
        Canonical k-mer sequence.
    k : int
        K-mer length.
    species_name : str
        Source species label.
    genome_id : str
        Source genome identifier.
    contig_id : str
        Source contig identifier.
    position : int
        Zero-based start position on the source contig.
    role : str
        Source genome role.
    clade : str
        Source clade label.
    """

    kmer: str
    k: int
    species_name: str
    genome_id: str
    contig_id: str
    position: int
    role: str
    clade: str


@dataclass(frozen=True)
class DiagnosticKmer:
    """A diagnostic k-mer retained in a panel.

    Attributes
    ----------
    kmer : str
        Canonical k-mer sequence.
    k : int
        K-mer length.
    panel_type : str
        Diagnostic class, such as species_unique or clade_core.
    species_name : str
        Species label for species-level k-mers. Empty for clade-level k-mers.
    clade : str
        Clade/group label.
    source_genomes : str
        Semicolon-separated source genome identifiers.
    source_contigs : str
        Semicolon-separated source contig identifiers.
    example_position : int
        Example zero-based source position.
    """

    kmer: str
    k: int
    panel_type: str
    species_name: str
    clade: str
    source_genomes: str
    source_contigs: str
    example_position: int

    def to_record(self) -> dict[str, object]:
        """Convert the diagnostic k-mer to a serialisable record.

        Returns
        -------
        dict[str, object]
            Dictionary representation.
        """
        return {
            "kmer": self.kmer,
            "k": self.k,
            "panel_type": self.panel_type,
            "species_name": self.species_name,
            "clade": self.clade,
            "source_genomes": self.source_genomes,
            "source_contigs": self.source_contigs,
            "example_position": self.example_position,
        }


def collect_kmer_observations(
    *,
    genome_configs: Iterable[GenomeConfig],
    k_values: Iterable[int],
) -> list[KmerObservation]:
    """Collect reference k-mer observations from configured genomes.

    Parameters
    ----------
    genome_configs : iterable of GenomeConfig
        Genome metadata records.
    k_values : iterable of int
        K-mer lengths to process.

    Returns
    -------
    list[KmerObservation]
        Observed reference k-mers.
    """
    observations: list[KmerObservation] = []
    for genome_config in genome_configs:
        for fasta_record in read_fasta_records(fasta_path=genome_config.genome_fasta):
            for k in k_values:
                for position, kmer in iter_kmers(sequence=fasta_record.sequence, k=k):
                    observations.append(
                        KmerObservation(
                            kmer=kmer,
                            k=k,
                            species_name=genome_config.species_name,
                            genome_id=genome_config.genome_id,
                            contig_id=fasta_record.identifier,
                            position=position,
                            role=genome_config.role,
                            clade=genome_config.clade,
                        )
                    )
    return observations


def identify_species_unique_kmers(
    *,
    observations: Iterable[KmerObservation],
    target_species: set[str],
) -> list[DiagnosticKmer]:
    """Identify k-mers unique to each target species.

    Parameters
    ----------
    observations : iterable of KmerObservation
        Reference k-mer observations.
    target_species : set[str]
        Species labels treated as target species.

    Returns
    -------
    list[DiagnosticKmer]
        Species-unique diagnostic k-mers.
    """
    by_key: dict[tuple[int, str], list[KmerObservation]] = defaultdict(list)
    for observation in observations:
        by_key[(observation.k, observation.kmer)].append(observation)

    diagnostics: list[DiagnosticKmer] = []
    for (k, kmer), group in by_key.items():
        species = {item.species_name for item in group}
        target_overlap = species & target_species
        if len(species) == 1 and len(target_overlap) == 1:
            first = group[0]
            diagnostics.append(
                DiagnosticKmer(
                    kmer=kmer,
                    k=k,
                    panel_type="species_unique",
                    species_name=first.species_name,
                    clade=first.clade,
                    source_genomes=";".join(sorted({item.genome_id for item in group})),
                    source_contigs=";".join(sorted({item.contig_id for item in group})),
                    example_position=first.position,
                )
            )
    return diagnostics


def identify_clade_core_kmers(
    *,
    observations: Iterable[KmerObservation],
    target_clade: str,
) -> list[DiagnosticKmer]:
    """Identify k-mers present only within a named clade.

    Parameters
    ----------
    observations : iterable of KmerObservation
        Reference k-mer observations.
    target_clade : str
        Clade label defining the clade-level panel.

    Returns
    -------
    list[DiagnosticKmer]
        Clade-level diagnostic k-mers.
    """
    if not target_clade:
        return []

    by_key: dict[tuple[int, str], list[KmerObservation]] = defaultdict(list)
    for observation in observations:
        by_key[(observation.k, observation.kmer)].append(observation)

    diagnostics: list[DiagnosticKmer] = []
    for (k, kmer), group in by_key.items():
        clades = {item.clade for item in group}
        if clades == {target_clade}:
            first = group[0]
            species = {item.species_name for item in group}
            diagnostics.append(
                DiagnosticKmer(
                    kmer=kmer,
                    k=k,
                    panel_type="clade_core",
                    species_name="" if len(species) > 1 else first.species_name,
                    clade=target_clade,
                    source_genomes=";".join(sorted({item.genome_id for item in group})),
                    source_contigs=";".join(sorted({item.contig_id for item in group})),
                    example_position=first.position,
                )
            )
    return diagnostics


def thin_diagnostic_kmers(
    *,
    diagnostic_kmers: Iterable[DiagnosticKmer],
    max_per_species_per_k: int | None = None,
) -> list[DiagnosticKmer]:
    """Thin diagnostic k-mers to a maximum count per species and k.

    Parameters
    ----------
    diagnostic_kmers : iterable of DiagnosticKmer
        Diagnostic k-mers.
    max_per_species_per_k : int | None, optional
        Maximum number retained per panel type, species/clade and k.

    Returns
    -------
    list[DiagnosticKmer]
        Thinned diagnostic k-mers.
    """
    if max_per_species_per_k is None:
        return list(diagnostic_kmers)
    if max_per_species_per_k <= 0:
        raise ValueError("max_per_species_per_k must be positive")

    counts: dict[tuple[str, str, str, int], int] = defaultdict(int)
    retained: list[DiagnosticKmer] = []
    for item in sorted(
        diagnostic_kmers,
        key=lambda x: (x.panel_type, x.species_name, x.clade, x.k, x.kmer),
    ):
        key = (item.panel_type, item.species_name, item.clade, item.k)
        if counts[key] >= max_per_species_per_k:
            continue
        retained.append(item)
        counts[key] += 1
    return retained


def build_panel(
    *,
    genome_configs: list[GenomeConfig],
    k_values: list[int],
    target_clade: str = "",
    max_per_species_per_k: int | None = None,
) -> tuple[list[DiagnosticKmer], list[dict[str, object]]]:
    """Build species and optional clade diagnostic k-mer panels.

    Parameters
    ----------
    genome_configs : list[GenomeConfig]
        Genome metadata records.
    k_values : list[int]
        K-mer lengths.
    target_clade : str, optional
        Optional clade label for clade-core k-mers.
    max_per_species_per_k : int | None, optional
        Optional thinning limit.

    Returns
    -------
    tuple[list[DiagnosticKmer], list[dict[str, object]]]
        Diagnostic k-mers and summary records.
    """
    observations = collect_kmer_observations(
        genome_configs=genome_configs,
        k_values=k_values,
    )
    target_species = {config.species_name for config in genome_configs if config.is_target}
    species_unique = identify_species_unique_kmers(
        observations=observations,
        target_species=target_species,
    )
    clade_core = identify_clade_core_kmers(
        observations=observations,
        target_clade=target_clade,
    )
    diagnostics = thin_diagnostic_kmers(
        diagnostic_kmers=[*species_unique, *clade_core],
        max_per_species_per_k=max_per_species_per_k,
    )
    summary = summarise_panel(diagnostic_kmers=diagnostics)
    return diagnostics, summary


def summarise_panel(
    *,
    diagnostic_kmers: Iterable[DiagnosticKmer],
) -> list[dict[str, object]]:
    """Summarise diagnostic k-mer counts.

    Parameters
    ----------
    diagnostic_kmers : iterable of DiagnosticKmer
        Diagnostic k-mers.

    Returns
    -------
    list[dict[str, object]]
        Summary records by panel type, species/clade and k.
    """
    counts: dict[tuple[str, str, str, int], int] = defaultdict(int)
    for item in diagnostic_kmers:
        key = (item.panel_type, item.species_name, item.clade, item.k)
        counts[key] += 1
    return [
        {
            "panel_type": panel_type,
            "species_name": species_name,
            "clade": clade,
            "k": k,
            "diagnostic_kmers": count,
        }
        for (panel_type, species_name, clade, k), count in sorted(counts.items())
    ]


def load_panel(*, panel_path: str | Path) -> dict[int, dict[str, list[DiagnosticKmer]]]:
    """Load a diagnostic panel into a k-mer lookup index.

    Parameters
    ----------
    panel_path : str or pathlib.Path
        Panel TSV or TSV.GZ path.

    Returns
    -------
    dict[int, dict[str, list[DiagnosticKmer]]]
        Mapping of k to k-mer to diagnostic records.
    """
    from kmersutra.io import read_tsv

    index: dict[int, dict[str, list[DiagnosticKmer]]] = defaultdict(lambda: defaultdict(list))
    for record in read_tsv(input_path=panel_path):
        item = DiagnosticKmer(
            kmer=record["kmer"],
            k=int(record["k"]),
            panel_type=record["panel_type"],
            species_name=record.get("species_name", ""),
            clade=record.get("clade", ""),
            source_genomes=record.get("source_genomes", ""),
            source_contigs=record.get("source_contigs", ""),
            example_position=int(record.get("example_position", "0") or 0),
        )
        index[item.k][item.kmer].append(item)
    return index
