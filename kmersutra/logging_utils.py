"""Logging helpers for KmerSutra."""

from __future__ import annotations

import logging
from pathlib import Path


def configure_logging(
    *,
    log_file: str | Path | None = None,
    verbose: bool = False,
) -> logging.Logger:
    """Configure package logging.

    Parameters
    ----------
    log_file : str | pathlib.Path | None, optional
        Optional path to a log file. Parent directories are created if needed.
    verbose : bool, optional
        If true, emit informational messages to stderr. Otherwise only warnings
        and errors are shown on stderr.

    Returns
    -------
    logging.Logger
        Configured KmerSutra logger.
    """
    logger = logging.getLogger("kmersutra")
    logger.setLevel(logging.DEBUG)
    logger.handlers.clear()

    formatter = logging.Formatter(
        "%(asctime)s %(levelname)s [%(name)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.INFO if verbose else logging.WARNING)
    stream_handler.setFormatter(formatter)
    logger.addHandler(stream_handler)

    if log_file is not None:
        log_path = Path(log_file)
        log_path.parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(log_path, mode="w")
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    return logger
