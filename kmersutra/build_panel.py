"""Build taxonomically aware diagnostic k-mer panels."""

from __future__ import annotations

import logging
from collections import defaultdict
from collections.abc import Iterable
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path

from kmersutra.config import GenomeConfig
from kmersutra.fasta import read_fasta_records
from kmersutra.kmers import iter_kmers
from kmersutra.taxonomy import CORE_RANK_ORDER, TaxonomyDatabase


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
    taxid : str
        Source NCBI taxid, when available.
    """

    kmer: str
    k: int
    species_name: str
    genome_id: str
    contig_id: str
    position: int
    role: str
    clade: str
    taxid: str = ""


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
        Evidence class, such as species_unique, genus_core or clade_core.
    species_name : str
        Species label for species-level k-mers. Empty for broader evidence.
    clade : str
        Clade/group label.
    source_genomes : str
        Semicolon-separated source genome identifiers.
    source_contigs : str
        Semicolon-separated source contig identifiers.
    example_position : int
        Example zero-based source position.
    evidence_taxid : str
        Taxid represented by this evidence level, when available.
    evidence_name : str
        Taxonomic name represented by this evidence level, when available.
    evidence_rank : str
        Taxonomic rank represented by this evidence level, when available.
    lineage_taxids : str
        Semicolon-separated lineage taxids for the evidence taxid.
    source_taxids : str
        Semicolon-separated source taxids in which this k-mer was observed.
    """

    kmer: str
    k: int
    panel_type: str
    species_name: str
    clade: str
    source_genomes: str
    source_contigs: str
    example_position: int
    evidence_taxid: str = ""
    evidence_name: str = ""
    evidence_rank: str = ""
    lineage_taxids: str = ""
    source_taxids: str = ""

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
            "evidence_taxid": self.evidence_taxid,
            "evidence_name": self.evidence_name,
            "evidence_rank": self.evidence_rank,
            "lineage_taxids": self.lineage_taxids,
            "source_taxids": self.source_taxids,
        }


def _collect_for_genome(
    genome_config: GenomeConfig,
    k_values: tuple[int, ...],
) -> tuple[list[KmerObservation], dict[str, object]]:
    """Collect k-mer observations for one genome.

    Parameters
    ----------
    genome_config : GenomeConfig
        Genome metadata record.
    k_values : tuple[int, ...]
        K-mer lengths to collect.

    Returns
    -------
    tuple[list[KmerObservation], dict[str, object]]
        Observations and a compact per-genome collection summary.
    """
    observations: list[KmerObservation] = []
    contigs = 0
    total_bases = 0
    per_k_counts = {k: 0 for k in k_values}

    for fasta_record in read_fasta_records(fasta_path=genome_config.genome_fasta):
        contigs += 1
        total_bases += len(fasta_record.sequence)
        for k in k_values:
            count = 0
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
                        taxid=genome_config.taxid,
                    )
                )
                count += 1
            per_k_counts[k] += count

    summary = {
        "genome_id": genome_config.genome_id,
        "species_name": genome_config.species_name,
        "taxid": genome_config.taxid,
        "role": genome_config.role,
        "clade": genome_config.clade,
        "genome_fasta": str(genome_config.genome_fasta),
        "contigs": contigs,
        "total_bases": total_bases,
        "total_observations": len(observations),
        **{f"observations_k{k}": per_k_counts[k] for k in k_values},
    }
    return observations, summary


def collect_kmer_observations(
    *,
    genome_configs: Iterable[GenomeConfig],
    k_values: Iterable[int],
    threads: int = 1,
    logger: logging.Logger | None = None,
) -> tuple[list[KmerObservation], list[dict[str, object]]]:
    """Collect reference k-mer observations from configured genomes.

    Parameters
    ----------
    genome_configs : iterable of GenomeConfig
        Genome metadata records.
    k_values : iterable of int
        K-mer lengths to process.
    threads : int, optional
        Number of worker processes for per-genome collection.
    logger : logging.Logger | None, optional
        Logger for progress messages.

    Returns
    -------
    tuple[list[KmerObservation], list[dict[str, object]]]
        Observed reference k-mers and per-genome collection summaries.
    """
    configs = list(genome_configs)
    k_tuple = tuple(k_values)
    if threads <= 0:
        raise ValueError("threads must be a positive integer")

    if logger:
        logger.info(
            "Collecting k-mers for %d genomes, %d k values, using %d worker(s)",
            len(configs),
            len(k_tuple),
            threads,
        )

    observations: list[KmerObservation] = []
    summaries: list[dict[str, object]] = []

    if threads == 1 or len(configs) <= 1:
        for index, genome_config in enumerate(configs, start=1):
            if logger:
                logger.info(
                    "Collecting genome %d/%d: %s (%s; taxid=%s)",
                    index,
                    len(configs),
                    genome_config.genome_id,
                    genome_config.species_name,
                    genome_config.taxid or "NA",
                )
            genome_observations, summary = _collect_for_genome(genome_config, k_tuple)
            observations.extend(genome_observations)
            summaries.append(summary)
            if logger:
                logger.info(
                    "Collected %d observations from %s",
                    len(genome_observations),
                    genome_config.genome_id,
                )
        return observations, summaries

    max_workers = min(threads, len(configs))
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_config = {
            executor.submit(_collect_for_genome, genome_config, k_tuple): genome_config
            for genome_config in configs
        }
        for future in as_completed(future_to_config):
            genome_config = future_to_config[future]
            genome_observations, summary = future.result()
            observations.extend(genome_observations)
            summaries.append(summary)
            if logger:
                logger.info(
                    "Collected %d observations from %s (%s; taxid=%s)",
                    len(genome_observations),
                    genome_config.genome_id,
                    genome_config.species_name,
                    genome_config.taxid or "NA",
                )

    summaries.sort(key=lambda record: str(record["genome_id"]))
    return observations, summaries


def _group_observations_by_kmer(
    observations: Iterable[KmerObservation],
) -> dict[tuple[int, str], list[KmerObservation]]:
    """Group observations by k value and k-mer sequence."""
    by_key: dict[tuple[int, str], list[KmerObservation]] = defaultdict(list)
    for observation in observations:
        by_key[(observation.k, observation.kmer)].append(observation)
    return by_key


def _diagnostic_from_group(
    *,
    k: int,
    kmer: str,
    group: list[KmerObservation],
    panel_type: str,
    species_name: str,
    clade: str,
    evidence_taxid: str = "",
    evidence_name: str = "",
    evidence_rank: str = "",
    lineage_taxids: str = "",
) -> DiagnosticKmer:
    """Create a diagnostic k-mer record from grouped observations."""
    first = group[0]
    return DiagnosticKmer(
        kmer=kmer,
        k=k,
        panel_type=panel_type,
        species_name=species_name,
        clade=clade,
        source_genomes=";".join(sorted({item.genome_id for item in group})),
        source_contigs=";".join(sorted({item.contig_id for item in group})),
        example_position=first.position,
        evidence_taxid=evidence_taxid,
        evidence_name=evidence_name,
        evidence_rank=evidence_rank,
        lineage_taxids=lineage_taxids,
        source_taxids=";".join(sorted({item.taxid for item in group if item.taxid})),
    )


def identify_species_unique_kmers(
    *,
    observations: Iterable[KmerObservation],
    target_species: set[str],
    logger: logging.Logger | None = None,
) -> list[DiagnosticKmer]:
    """Identify k-mers unique to each target species.

    Parameters
    ----------
    observations : iterable of KmerObservation
        Reference k-mer observations.
    target_species : set[str]
        Species labels treated as target species.
    logger : logging.Logger | None, optional
        Logger for progress messages.

    Returns
    -------
    list[DiagnosticKmer]
        Species-unique diagnostic k-mers.
    """
    by_key = _group_observations_by_kmer(observations)

    if logger:
        logger.info("Testing %d distinct reference k-mer keys for species uniqueness", len(by_key))

    diagnostics: list[DiagnosticKmer] = []
    for (k, kmer), group in by_key.items():
        species = {item.species_name for item in group}
        target_overlap = species & target_species
        if len(species) == 1 and len(target_overlap) == 1:
            first = group[0]
            diagnostics.append(
                _diagnostic_from_group(
                    k=k,
                    kmer=kmer,
                    group=group,
                    panel_type="species_unique",
                    species_name=first.species_name,
                    clade=first.clade,
                    evidence_taxid=first.taxid,
                    evidence_name=first.species_name,
                    evidence_rank="species",
                )
            )
    if logger:
        logger.info("Identified %d species-unique diagnostic k-mers", len(diagnostics))
    return diagnostics


def identify_clade_core_kmers(
    *,
    observations: Iterable[KmerObservation],
    target_clade: str,
    logger: logging.Logger | None = None,
) -> list[DiagnosticKmer]:
    """Identify k-mers present only within a named clade.

    Parameters
    ----------
    observations : iterable of KmerObservation
        Reference k-mer observations.
    target_clade : str
        Clade label defining the clade-level panel.
    logger : logging.Logger | None, optional
        Logger for progress messages.

    Returns
    -------
    list[DiagnosticKmer]
        Clade-level diagnostic k-mers.
    """
    if not target_clade:
        return []

    by_key = _group_observations_by_kmer(observations)

    if logger:
        logger.info("Testing %d distinct reference k-mer keys for clade specificity", len(by_key))

    diagnostics: list[DiagnosticKmer] = []
    for (k, kmer), group in by_key.items():
        clades = {item.clade for item in group}
        if clades == {target_clade}:
            species = {item.species_name for item in group}
            first = group[0]
            diagnostics.append(
                _diagnostic_from_group(
                    k=k,
                    kmer=kmer,
                    group=group,
                    panel_type="clade_core",
                    species_name="" if len(species) > 1 else first.species_name,
                    clade=target_clade,
                    evidence_name=target_clade,
                    evidence_rank="clade",
                )
            )
    if logger:
        logger.info("Identified %d clade-core diagnostic k-mers", len(diagnostics))
    return diagnostics


def identify_taxonomic_evidence_kmers(
    *,
    observations: Iterable[KmerObservation],
    taxonomy_db: TaxonomyDatabase,
    target_taxid: str = "",
    preferred_ranks: list[str] | None = None,
    logger: logging.Logger | None = None,
) -> list[DiagnosticKmer]:
    """Identify diagnostic k-mers at their best supported taxonomic level.

    Parameters
    ----------
    observations : iterable of KmerObservation
        Reference k-mer observations.
    taxonomy_db : TaxonomyDatabase
        Parsed NCBI taxonomy database.
    target_taxid : str, optional
        Optional taxid whose subtree should be retained.
    preferred_ranks : list[str] | None, optional
        Taxonomic ranks to retain as evidence levels.
    logger : logging.Logger | None, optional
        Logger for progress messages.

    Returns
    -------
    list[DiagnosticKmer]
        Taxonomically assigned diagnostic k-mers.
    """
    ranks = preferred_ranks or CORE_RANK_ORDER
    by_key = _group_observations_by_kmer(observations)
    target_taxid = taxonomy_db.normalise_taxid(target_taxid)

    if logger:
        logger.info(
            "Testing %d distinct reference k-mer keys for taxonomic evidence levels",
            len(by_key),
        )
        if target_taxid:
            logger.info(
                "Restricting taxonomic evidence to descendants of taxid %s (%s)",
                target_taxid,
                taxonomy_db.get_name(target_taxid) or "unnamed",
            )

    diagnostics: list[DiagnosticKmer] = []
    skipped_no_taxid = 0
    skipped_outside_target = 0
    skipped_unranked = 0

    for (k, kmer), group in by_key.items():
        source_taxids = {
            taxonomy_db.normalise_taxid(item.taxid) for item in group if item.taxid
        }
        source_taxids = {taxid for taxid in source_taxids if taxid}
        if not source_taxids:
            skipped_no_taxid += 1
            continue

        evidence_node = taxonomy_db.best_named_ancestor(
            taxids=source_taxids,
            preferred_ranks=ranks,
        )
        if evidence_node is None:
            skipped_unranked += 1
            continue
        if target_taxid and not taxonomy_db.is_descendant(
            taxid=evidence_node.taxid,
            ancestor_taxid=target_taxid,
        ):
            skipped_outside_target += 1
            continue
        if evidence_node.rank not in ranks:
            skipped_unranked += 1
            continue

        species = {item.species_name for item in group}
        clades = {item.clade for item in group if item.clade}
        species_name = next(iter(species)) if evidence_node.rank == "species" and len(species) == 1 else ""
        clade = next(iter(clades)) if len(clades) == 1 else evidence_node.name
        panel_type = "species_unique" if evidence_node.rank == "species" else f"{evidence_node.rank}_core"

        diagnostics.append(
            _diagnostic_from_group(
                k=k,
                kmer=kmer,
                group=group,
                panel_type=panel_type,
                species_name=species_name,
                clade=clade,
                evidence_taxid=evidence_node.taxid,
                evidence_name=evidence_node.name,
                evidence_rank=evidence_node.rank,
                lineage_taxids=";".join(taxonomy_db.get_lineage(evidence_node.taxid)),
            )
        )

    if logger:
        logger.info("Identified %d taxonomic evidence k-mers", len(diagnostics))
        logger.info("Skipped %d k-mers without usable taxids", skipped_no_taxid)
        logger.info("Skipped %d k-mers outside target taxid", skipped_outside_target)
        logger.info("Skipped %d k-mers without a retained evidence rank", skipped_unranked)
    return diagnostics


def thin_diagnostic_kmers(
    *,
    diagnostic_kmers: Iterable[DiagnosticKmer],
    max_per_species_per_k: int | None = None,
    logger: logging.Logger | None = None,
) -> list[DiagnosticKmer]:
    """Thin diagnostic k-mers to a maximum count per evidence level and k.

    Parameters
    ----------
    diagnostic_kmers : iterable of DiagnosticKmer
        Diagnostic k-mers.
    max_per_species_per_k : int | None, optional
        Maximum number retained per panel type, taxon/species/clade and k.
    logger : logging.Logger | None, optional
        Logger for progress messages.

    Returns
    -------
    list[DiagnosticKmer]
        Thinned diagnostic k-mers.
    """
    diagnostics = list(diagnostic_kmers)
    if max_per_species_per_k is None:
        return diagnostics
    if max_per_species_per_k <= 0:
        raise ValueError("max_per_species_per_k must be positive")

    counts: dict[tuple[str, str, str, str, int], int] = defaultdict(int)
    retained: list[DiagnosticKmer] = []
    for item in sorted(
        diagnostics,
        key=lambda x: (x.panel_type, x.evidence_taxid, x.species_name, x.clade, x.k, x.kmer),
    ):
        key = (item.panel_type, item.evidence_taxid, item.species_name, item.clade, item.k)
        if counts[key] >= max_per_species_per_k:
            continue
        retained.append(item)
        counts[key] += 1
    if logger:
        logger.info(
            "Thinned diagnostic k-mers from %d to %d using max_per_species_per_k=%d",
            len(diagnostics),
            len(retained),
            max_per_species_per_k,
        )
    return retained


def build_panel(
    *,
    genome_configs: list[GenomeConfig],
    k_values: list[int],
    target_clade: str = "",
    max_per_species_per_k: int | None = None,
    threads: int = 1,
    taxonomy_db: TaxonomyDatabase | None = None,
    target_taxid: str = "",
    preferred_ranks: list[str] | None = None,
    logger: logging.Logger | None = None,
) -> tuple[list[DiagnosticKmer], list[dict[str, object]], list[dict[str, object]]]:
    """Build diagnostic k-mer panels.

    Parameters
    ----------
    genome_configs : list[GenomeConfig]
        Genome metadata records.
    k_values : list[int]
        K-mer lengths.
    target_clade : str, optional
        Optional clade label for clade-core k-mers when no taxonomy is used.
    max_per_species_per_k : int | None, optional
        Optional thinning limit.
    threads : int, optional
        Number of worker processes used during reference k-mer collection.
    taxonomy_db : TaxonomyDatabase | None, optional
        Optional NCBI taxonomy database for taxonomic evidence assignment.
    target_taxid : str, optional
        Optional root taxid for retained evidence when taxonomy is used.
    preferred_ranks : list[str] | None, optional
        Retained evidence ranks when taxonomy is used.
    logger : logging.Logger | None, optional
        Logger for progress messages.

    Returns
    -------
    tuple[list[DiagnosticKmer], list[dict[str, object]], list[dict[str, object]]]
        Diagnostic k-mers, panel summary records and collection summaries.
    """
    observations, collection_summary = collect_kmer_observations(
        genome_configs=genome_configs,
        k_values=k_values,
        threads=threads,
        logger=logger,
    )
    if logger:
        logger.info("Collected %d total reference k-mer observations", len(observations))

    if taxonomy_db is not None:
        diagnostics = identify_taxonomic_evidence_kmers(
            observations=observations,
            taxonomy_db=taxonomy_db,
            target_taxid=target_taxid,
            preferred_ranks=preferred_ranks,
            logger=logger,
        )
    else:
        target_species = {config.species_name for config in genome_configs if config.is_target}
        if logger:
            logger.info("Target species: %s", "; ".join(sorted(target_species)))
        species_unique = identify_species_unique_kmers(
            observations=observations,
            target_species=target_species,
            logger=logger,
        )
        clade_core = identify_clade_core_kmers(
            observations=observations,
            target_clade=target_clade,
            logger=logger,
        )
        diagnostics = [*species_unique, *clade_core]

    thinned = thin_diagnostic_kmers(
        diagnostic_kmers=diagnostics,
        max_per_species_per_k=max_per_species_per_k,
        logger=logger,
    )
    summary = summarise_panel(diagnostic_kmers=thinned)
    return thinned, summary, collection_summary


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
        Summary records by panel type, taxonomic evidence and k.
    """
    counts: dict[tuple[str, str, str, str, str, int], int] = defaultdict(int)
    for item in diagnostic_kmers:
        key = (
            item.panel_type,
            item.species_name,
            item.clade,
            item.evidence_taxid,
            item.evidence_rank,
            item.k,
        )
        counts[key] += 1
    return [
        {
            "panel_type": panel_type,
            "species_name": species_name,
            "clade": clade,
            "evidence_taxid": evidence_taxid,
            "evidence_rank": evidence_rank,
            "k": k,
            "diagnostic_kmers": count,
        }
        for (panel_type, species_name, clade, evidence_taxid, evidence_rank, k), count
        in sorted(counts.items())
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
            evidence_taxid=record.get("evidence_taxid", ""),
            evidence_name=record.get("evidence_name", ""),
            evidence_rank=record.get("evidence_rank", ""),
            lineage_taxids=record.get("lineage_taxids", ""),
            source_taxids=record.get("source_taxids", ""),
        )
        index[item.k][item.kmer].append(item)
    return index
