"""NCBI taxonomy utilities for KmerSutra.

This module provides a small, dependency-light parser for the NCBI taxonomy
``taxdump`` files. It is intentionally focused on the pieces needed by
KmerSutra: lineage lookup, rank lookup, lowest-common-ancestor assignment,
and automatic download of the standard taxonomy dump when local files are
missing.
"""

from __future__ import annotations

import logging
import shutil
import urllib.request
import zipfile
from dataclasses import dataclass
from pathlib import Path
from tempfile import NamedTemporaryFile

NCBI_TAXONOMY_URL = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip"
REQUIRED_TAXDUMP_FILES = ("nodes.dmp", "names.dmp", "merged.dmp", "delnodes.dmp")
CORE_RANK_ORDER = [
    "species",
    "genus",
    "family",
    "order",
    "class",
    "phylum",
    "superkingdom",
]


@dataclass(frozen=True)
class TaxonomyNode:
    """A single node in the NCBI taxonomy tree.

    Attributes
    ----------
    taxid : str
        NCBI taxon identifier.
    parent_taxid : str
        Parent NCBI taxon identifier.
    rank : str
        NCBI taxonomic rank.
    name : str
        Scientific name, when available.
    """

    taxid: str
    parent_taxid: str
    rank: str
    name: str = ""


def _split_dmp_line(line: str) -> list[str]:
    """Split one NCBI ``.dmp`` line into stripped fields.

    Parameters
    ----------
    line : str
        Raw line from an NCBI taxonomy dump file.

    Returns
    -------
    list[str]
        Parsed fields with terminators removed.
    """
    return [field.strip() for field in line.rstrip("\n").split("|")[:-1]]


def download_taxdump(
    *,
    taxonomy_dir: str | Path,
    url: str = NCBI_TAXONOMY_URL,
    overwrite: bool = False,
    logger: logging.Logger | None = None,
) -> Path:
    """Download and extract the NCBI taxonomy dump.

    Parameters
    ----------
    taxonomy_dir : str or pathlib.Path
        Directory where taxonomy files should be written.
    url : str, optional
        URL for the NCBI ``taxdmp.zip`` file.
    overwrite : bool, optional
        Replace existing required files when true.
    logger : logging.Logger | None, optional
        Logger for progress messages.

    Returns
    -------
    pathlib.Path
        Directory containing the extracted taxonomy files.
    """
    output_dir = Path(taxonomy_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    existing = [output_dir / filename for filename in REQUIRED_TAXDUMP_FILES]
    if all(path.exists() for path in existing) and not overwrite:
        if logger:
            logger.info("NCBI taxonomy files already exist in %s", output_dir)
        return output_dir

    if logger:
        logger.info("Downloading NCBI taxonomy dump from %s", url)

    with NamedTemporaryFile(suffix=".zip", delete=False) as tmp_handle:
        tmp_zip = Path(tmp_handle.name)

    try:
        with urllib.request.urlopen(url) as response, tmp_zip.open("wb") as out_handle:
            shutil.copyfileobj(response, out_handle)

        if logger:
            logger.info("Extracting required taxonomy files to %s", output_dir)

        with zipfile.ZipFile(tmp_zip) as zip_handle:
            names = set(zip_handle.namelist())
            missing = [filename for filename in REQUIRED_TAXDUMP_FILES if filename not in names]
            if missing:
                raise ValueError(
                    "Downloaded taxonomy archive is missing required files: "
                    + ", ".join(missing)
                )
            for filename in REQUIRED_TAXDUMP_FILES:
                target_path = output_dir / filename
                if target_path.exists() and not overwrite:
                    continue
                with zip_handle.open(filename) as source, target_path.open("wb") as dest:
                    shutil.copyfileobj(source, dest)
    finally:
        tmp_zip.unlink(missing_ok=True)

    return output_dir


def ensure_taxdump_files(
    *,
    taxonomy_dir: str | Path,
    download_if_missing: bool = False,
    url: str = NCBI_TAXONOMY_URL,
    logger: logging.Logger | None = None,
) -> dict[str, Path]:
    """Ensure that required NCBI taxonomy files are available.

    Parameters
    ----------
    taxonomy_dir : str or pathlib.Path
        Directory expected to contain NCBI taxdump files.
    download_if_missing : bool, optional
        Download ``taxdmp.zip`` when one or more files are missing.
    url : str, optional
        URL for the taxonomy archive.
    logger : logging.Logger | None, optional
        Logger for progress messages.

    Returns
    -------
    dict[str, pathlib.Path]
        Mapping from required filename to local path.
    """
    root = Path(taxonomy_dir)
    paths = {filename: root / filename for filename in REQUIRED_TAXDUMP_FILES}
    missing = [filename for filename, path in paths.items() if not path.exists()]

    if missing and download_if_missing:
        if logger:
            logger.info(
                "Missing taxonomy files (%s); downloading taxdump",
                ", ".join(missing),
            )
        download_taxdump(taxonomy_dir=root, url=url, logger=logger)
        missing = [filename for filename, path in paths.items() if not path.exists()]

    if missing:
        raise FileNotFoundError(
            "Missing NCBI taxonomy files in "
            f"{root}: "
            + ", ".join(missing)
            + ". Use --download_taxonomy_if_missing to fetch taxdmp.zip."
        )

    return paths


class TaxonomyDatabase:
    """In-memory NCBI taxonomy database.

    The class stores only the fields needed by KmerSutra and is deliberately
    lightweight enough for command-line workflows.
    """

    def __init__(
        self,
        *,
        nodes: dict[str, TaxonomyNode],
        merged: dict[str, str] | None = None,
        deleted: set[str] | None = None,
    ) -> None:
        """Create a taxonomy database.

        Parameters
        ----------
        nodes : dict[str, TaxonomyNode]
            Taxonomy nodes keyed by taxid.
        merged : dict[str, str] | None, optional
            Mapping from obsolete taxids to replacement taxids.
        deleted : set[str] | None, optional
            Deleted taxids.
        """
        self.nodes = nodes
        self.merged = merged or {}
        self.deleted = deleted or set()

    @classmethod
    def from_taxdump(
        cls,
        *,
        taxonomy_dir: str | Path,
        download_if_missing: bool = False,
        url: str = NCBI_TAXONOMY_URL,
        logger: logging.Logger | None = None,
    ) -> "TaxonomyDatabase":
        """Load an NCBI taxonomy database from taxdump files.

        Parameters
        ----------
        taxonomy_dir : str or pathlib.Path
            Directory containing ``nodes.dmp``, ``names.dmp``, ``merged.dmp``
            and ``delnodes.dmp``.
        download_if_missing : bool, optional
            Download missing taxonomy files when true.
        url : str, optional
            URL for the taxonomy archive.
        logger : logging.Logger | None, optional
            Logger for progress messages.

        Returns
        -------
        TaxonomyDatabase
            Parsed taxonomy database.
        """
        paths = ensure_taxdump_files(
            taxonomy_dir=taxonomy_dir,
            download_if_missing=download_if_missing,
            url=url,
            logger=logger,
        )
        names = cls._read_scientific_names(paths["names.dmp"])
        nodes = cls._read_nodes(paths["nodes.dmp"], names)
        merged = cls._read_merged(paths["merged.dmp"])
        deleted = cls._read_deleted(paths["delnodes.dmp"])
        if logger:
            logger.info(
                "Loaded taxonomy database with %d nodes, %d merged taxids and %d deleted taxids",
                len(nodes),
                len(merged),
                len(deleted),
            )
        return cls(nodes=nodes, merged=merged, deleted=deleted)

    @staticmethod
    def _read_scientific_names(path: Path) -> dict[str, str]:
        """Read scientific names from ``names.dmp``."""
        names: dict[str, str] = {}
        with path.open("rt", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                fields = _split_dmp_line(line)
                if len(fields) < 4:
                    continue
                taxid, name_txt, _, name_class = fields[:4]
                if name_class == "scientific name":
                    names[taxid] = name_txt
        return names

    @staticmethod
    def _read_nodes(path: Path, names: dict[str, str]) -> dict[str, TaxonomyNode]:
        """Read taxonomy nodes from ``nodes.dmp``."""
        nodes: dict[str, TaxonomyNode] = {}
        with path.open("rt", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                fields = _split_dmp_line(line)
                if len(fields) < 3:
                    continue
                taxid, parent_taxid, rank = fields[:3]
                nodes[taxid] = TaxonomyNode(
                    taxid=taxid,
                    parent_taxid=parent_taxid,
                    rank=rank,
                    name=names.get(taxid, ""),
                )
        return nodes

    @staticmethod
    def _read_merged(path: Path) -> dict[str, str]:
        """Read obsolete-to-current taxid mappings from ``merged.dmp``."""
        merged: dict[str, str] = {}
        with path.open("rt", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                fields = _split_dmp_line(line)
                if len(fields) >= 2:
                    merged[fields[0]] = fields[1]
        return merged

    @staticmethod
    def _read_deleted(path: Path) -> set[str]:
        """Read deleted taxids from ``delnodes.dmp``."""
        deleted: set[str] = set()
        with path.open("rt", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                fields = _split_dmp_line(line)
                if fields:
                    deleted.add(fields[0])
        return deleted

    def normalise_taxid(self, taxid: str | int | None) -> str:
        """Return the current taxid for a possibly obsolete taxid.

        Parameters
        ----------
        taxid : str, int or None
            Input taxid.

        Returns
        -------
        str
            Current taxid, or an empty string if unavailable/deleted.
        """
        if taxid is None:
            return ""
        value = str(taxid).strip()
        if not value or value in self.deleted:
            return ""
        return self.merged.get(value, value)

    def get_node(self, taxid: str | int | None) -> TaxonomyNode | None:
        """Return a taxonomy node by taxid."""
        value = self.normalise_taxid(taxid)
        return self.nodes.get(value)

    def get_lineage(self, taxid: str | int | None) -> list[str]:
        """Return lineage taxids from root to the requested taxid.

        Parameters
        ----------
        taxid : str, int or None
            Taxid of interest.

        Returns
        -------
        list[str]
            Lineage taxids ordered root-to-leaf.
        """
        value = self.normalise_taxid(taxid)
        if not value or value not in self.nodes:
            return []

        lineage: list[str] = []
        seen: set[str] = set()
        current = value
        while current and current not in seen and current in self.nodes:
            seen.add(current)
            lineage.append(current)
            parent = self.nodes[current].parent_taxid
            if parent == current:
                break
            current = parent
        lineage.reverse()
        return lineage

    def get_ranked_lineage(self, taxid: str | int | None) -> list[TaxonomyNode]:
        """Return lineage nodes from root to requested taxid."""
        return [self.nodes[item] for item in self.get_lineage(taxid)]

    def get_name(self, taxid: str | int | None) -> str:
        """Return the scientific name for a taxid if available."""
        node = self.get_node(taxid)
        return node.name if node else ""

    def get_rank(self, taxid: str | int | None) -> str:
        """Return the rank for a taxid if available."""
        node = self.get_node(taxid)
        return node.rank if node else ""

    def is_descendant(self, *, taxid: str | int | None, ancestor_taxid: str | int | None) -> bool:
        """Return whether one taxid descends from another.

        Parameters
        ----------
        taxid : str, int or None
            Potential descendant taxid.
        ancestor_taxid : str, int or None
            Potential ancestor taxid.

        Returns
        -------
        bool
            True if ancestor taxid is in the descendant lineage.
        """
        ancestor = self.normalise_taxid(ancestor_taxid)
        return bool(ancestor and ancestor in self.get_lineage(taxid))

    def lowest_common_ancestor(self, taxids: list[str] | set[str] | tuple[str, ...]) -> str:
        """Return the lowest common ancestor for a set of taxids.

        Parameters
        ----------
        taxids : list, set or tuple of str
            Taxids to compare.

        Returns
        -------
        str
            Lowest common ancestor taxid, or empty string if unavailable.
        """
        lineages = [self.get_lineage(taxid) for taxid in taxids if self.get_lineage(taxid)]
        if not lineages:
            return ""

        common = lineages[0]
        for lineage in lineages[1:]:
            limit = min(len(common), len(lineage))
            new_common: list[str] = []
            for index in range(limit):
                if common[index] != lineage[index]:
                    break
                new_common.append(common[index])
            common = new_common
            if not common:
                return ""
        return common[-1] if common else ""

    def best_named_ancestor(
        self,
        *,
        taxids: list[str] | set[str] | tuple[str, ...],
        preferred_ranks: list[str] | None = None,
    ) -> TaxonomyNode | None:
        """Return the best ranked ancestor shared by taxids.

        Parameters
        ----------
        taxids : list, set or tuple of str
            Taxids to compare.
        preferred_ranks : list[str] | None, optional
            Ranks considered useful as evidence levels.

        Returns
        -------
        TaxonomyNode | None
            Shared ancestor node, preferably at a named/core rank.
        """
        lca = self.lowest_common_ancestor(taxids)
        if not lca:
            return None
        ranks = preferred_ranks or CORE_RANK_ORDER
        lineage = self.get_ranked_lineage(lca)
        for node in reversed(lineage):
            if node.rank in ranks:
                return node
        return self.get_node(lca)
