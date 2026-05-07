"""Panel index caching utilities for faster KmerSutra screening."""

from __future__ import annotations

import pickle
from pathlib import Path
from typing import Any

from kmersutra.build_panel import DiagnosticKmer, load_panel

CACHE_FORMAT_VERSION = 1


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


def _panel_signature(*, panel_path: str | Path) -> dict[str, int]:
    """Return a compact signature for stale-cache detection.

    Parameters
    ----------
    panel_path : str or pathlib.Path
        Source panel file path.

    Returns
    -------
    dict[str, int]
        File size and modification-time signature.
    """
    stat = Path(panel_path).stat()
    return {"size": stat.st_size, "mtime_ns": stat.st_mtime_ns}


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
    cache_file = Path(cache_path)
    cache_file.parent.mkdir(parents=True, exist_ok=True)
    payload: dict[str, Any] = {
        "format_version": CACHE_FORMAT_VERSION,
        "panel_path": str(Path(panel_path)),
        "panel_signature": _panel_signature(panel_path=panel_path),
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
        If the cache format is unsupported or stale.
    """
    with Path(cache_path).open("rb") as handle:
        payload = pickle.load(handle)

    if payload.get("format_version") != CACHE_FORMAT_VERSION:
        raise ValueError("Unsupported panel cache format version")

    if panel_path is not None and require_current:
        expected = _panel_signature(panel_path=panel_path)
        observed = payload.get("panel_signature")
        if observed != expected:
            raise ValueError("Panel cache is stale for the supplied panel file")

    panel_index = payload.get("panel_index")
    if not isinstance(panel_index, dict):
        raise ValueError("Panel cache does not contain a valid panel index")
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
