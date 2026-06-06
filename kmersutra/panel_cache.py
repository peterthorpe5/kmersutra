"""Panel index caching utilities for faster KmerSutra screening.

The cache stores the already parsed in-memory panel index used by screening.
It is deliberately validated against panel metadata before use so array tasks
cannot silently reuse an index built from a different panel.
"""

from __future__ import annotations

import hashlib
import pickle
from pathlib import Path
from typing import Any

from kmersutra import __version__
from kmersutra.build_panel import DiagnosticKmer, load_panel

CACHE_FORMAT_VERSION = 2
_HASH_CHUNK_SIZE = 1024 * 1024


def get_default_panel_cache_path(*, panel_path: str | Path) -> Path:
    """Return the default cache path for a panel file.

    Parameters
    ----------
    panel_path : str or pathlib.Path
        Source panel file path.

    Returns
    -------
    pathlib.Path
        Default cache path next to the panel.
    """
    return Path(f"{Path(panel_path)}.index.pkl")


def _sha256_file(*, path: str | Path) -> str:
    """Return the SHA-256 digest of a file.

    Parameters
    ----------
    path : str or pathlib.Path
        File to hash.

    Returns
    -------
    str
        Hexadecimal SHA-256 digest.
    """
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(_HASH_CHUNK_SIZE), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _panel_signature(*, panel_path: str | Path) -> dict[str, object]:
    """Return a panel signature for stale-cache detection.

    Parameters
    ----------
    panel_path : str or pathlib.Path
        Source panel file path.

    Returns
    -------
    dict[str, object]
        File path, size, modification time and SHA-256 digest.
    """
    panel_file = Path(panel_path)
    stat = panel_file.stat()
    return {
        "path_name": panel_file.name,
        "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
        "sha256": _sha256_file(path=panel_file),
    }


def _plain_panel_index(
    *,
    panel_index: dict[int, dict[str, list[DiagnosticKmer]]],
) -> dict[int, dict[str, list[DiagnosticKmer]]]:
    """Convert nested defaultdict objects to plain dictionaries.

    Parameters
    ----------
    panel_index : dict[int, dict[str, list[DiagnosticKmer]]]
        Panel index that may contain ``defaultdict`` instances.

    Returns
    -------
    dict[int, dict[str, list[DiagnosticKmer]]]
        Plain nested dictionary suitable for pickling.
    """
    return {int(k): dict(kmer_map) for k, kmer_map in panel_index.items()}


def _panel_index_metadata(
    *,
    panel_index: dict[int, dict[str, list[DiagnosticKmer]]],
) -> dict[str, object]:
    """Return compact metadata describing a parsed panel index.

    Parameters
    ----------
    panel_index : dict[int, dict[str, list[DiagnosticKmer]]]
        Parsed panel index.

    Returns
    -------
    dict[str, object]
        K values, number of unique k-mer keys and number of taxa/species.
    """
    species_names: set[str] = set()
    evidence_taxids: set[str] = set()
    n_diagnostic_rows = 0
    for kmer_map in panel_index.values():
        for diagnostics in kmer_map.values():
            n_diagnostic_rows += len(diagnostics)
            for diagnostic in diagnostics:
                if diagnostic.species_name:
                    species_names.add(str(diagnostic.species_name))
                if diagnostic.evidence_taxid:
                    evidence_taxids.add(str(diagnostic.evidence_taxid))
    return {
        "k_values": sorted(int(k_value) for k_value in panel_index),
        "n_panel_kmer_keys": sum(len(kmer_map) for kmer_map in panel_index.values()),
        "n_diagnostic_rows": n_diagnostic_rows,
        "n_species_names": len(species_names),
        "n_evidence_taxids": len(evidence_taxids),
    }


def write_panel_index_cache(
    *,
    panel_index: dict[int, dict[str, list[DiagnosticKmer]]],
    panel_path: str | Path,
    cache_path: str | Path,
) -> None:
    """Write a pickled panel index cache.

    Parameters
    ----------
    panel_index : dict[int, dict[str, list[DiagnosticKmer]]]
        Loaded KmerSutra panel index.
    panel_path : str or pathlib.Path
        Source panel path used to build the index.
    cache_path : str or pathlib.Path
        Output cache path.
    """
    if not panel_index:
        raise ValueError("panel_index is empty; refusing to write an empty cache")
    cache_file = Path(cache_path)
    cache_file.parent.mkdir(parents=True, exist_ok=True)
    payload: dict[str, Any] = {
        "format_version": CACHE_FORMAT_VERSION,
        "kmersutra_version": __version__,
        "panel_path": str(Path(panel_path)),
        "panel_signature": _panel_signature(panel_path=panel_path),
        "panel_metadata": _panel_index_metadata(panel_index=panel_index),
        "panel_index": _plain_panel_index(panel_index=panel_index),
    }
    with cache_file.open("wb") as handle:
        pickle.dump(payload, handle, protocol=pickle.HIGHEST_PROTOCOL)


def load_panel_index_cache(
    *,
    cache_path: str | Path,
    panel_path: str | Path | None = None,
    require_current: bool = True,
) -> dict[int, dict[str, list[DiagnosticKmer]]]:
    """Load a pickled panel index cache.

    Parameters
    ----------
    cache_path : str or pathlib.Path
        Input cache path.
    panel_path : str or pathlib.Path or None, optional
        Source panel path used for stale-cache checks.
    require_current : bool, optional
        If True, require the cache signature to match the current panel file.

    Returns
    -------
    dict[int, dict[str, list[DiagnosticKmer]]]
        Cached panel index.

    Raises
    ------
    ValueError
        If the cache format is unsupported, stale or malformed.
    """
    with Path(cache_path).open("rb") as handle:
        payload = pickle.load(handle)

    if payload.get("format_version") not in {1, CACHE_FORMAT_VERSION}:
        raise ValueError("Unsupported panel cache format version")

    if payload.get("format_version") != CACHE_FORMAT_VERSION:
        raise ValueError(
            "Panel cache was written by an older cache format; rebuild the cache"
        )

    if payload.get("kmersutra_version") != __version__:
        raise ValueError("Panel cache KmerSutra version does not match this package")

    if panel_path is not None and require_current:
        expected = _panel_signature(panel_path=panel_path)
        observed = payload.get("panel_signature")
        if observed != expected:
            raise ValueError("Panel cache is stale for the supplied panel file")

    panel_index = payload.get("panel_index")
    if not isinstance(panel_index, dict):
        raise ValueError("Panel cache does not contain a valid panel index")
    if not panel_index:
        raise ValueError("Panel cache contains an empty panel index")
    return panel_index


def load_panel_with_cache(
    *,
    panel_path: str | Path,
    cache_path: str | Path | None = None,
    use_cache: bool = False,
    write_cache: bool = False,
) -> tuple[dict[int, dict[str, list[DiagnosticKmer]]], str]:
    """Load a panel, optionally using or writing an index cache.

    Parameters
    ----------
    panel_path : str or pathlib.Path
        Source panel TSV or TSV.GZ file.
    cache_path : str or pathlib.Path or None, optional
        Cache file path. If omitted, a default path next to the panel is used.
    use_cache : bool, optional
        If True, load the cache when it exists and is current.
    write_cache : bool, optional
        If True, write a cache after loading the panel from TSV.

    Returns
    -------
    tuple[dict[int, dict[str, list[DiagnosticKmer]]], str]
        Loaded panel index and a string describing whether it came from
        ``cache`` or ``tsv``.
    """
    resolved_cache = Path(cache_path) if cache_path else get_default_panel_cache_path(
        panel_path=panel_path
    )

    if use_cache and resolved_cache.exists():
        return (
            load_panel_index_cache(
                cache_path=resolved_cache,
                panel_path=panel_path,
                require_current=True,
            ),
            "cache",
        )

    panel_index = load_panel(panel_path=panel_path)
    if write_cache or use_cache:
        write_panel_index_cache(
            panel_index=panel_index,
            panel_path=panel_path,
            cache_path=resolved_cache,
        )
    return panel_index, "tsv"
