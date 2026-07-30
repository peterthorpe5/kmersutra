"""Create reproducible package file inventories."""

from __future__ import annotations

import argparse
import hashlib
import logging
from collections.abc import Iterable
from pathlib import Path

from kmersutra.table_io import write_records_table

LOGGER = logging.getLogger("kmersutra")

EXCLUDED_DIRECTORY_NAMES = {
    ".git",
    ".mypy_cache",
    ".pytest_cache",
    ".ruff_cache",
    "__pycache__",
    "build",
    "dist",
    "htmlcov",
    "test_results",
}
EXCLUDED_FILE_SUFFIXES = {".pyc", ".pyo"}


def sha256_file(*, path: Path, chunk_size: int = 1024 * 1024) -> str:
    """Return a SHA-256 digest for one file.

    Args:
        path: File to hash.
        chunk_size: Number of bytes read per block.

    Returns:
        Hexadecimal SHA-256 digest.
    """
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def iter_package_files(*, repo_root: Path) -> Iterable[Path]:
    """Yield package files while excluding generated result directories.

    Args:
        repo_root: Repository directory.

    Yields:
        Files in deterministic relative-path order.
    """
    candidates = []
    for path in repo_root.rglob("*"):
        relative = path.relative_to(repo_root)
        if any(part in EXCLUDED_DIRECTORY_NAMES for part in relative.parts[:-1]):
            continue
        if path.name in EXCLUDED_DIRECTORY_NAMES:
            continue
        if path.suffix.lower() in EXCLUDED_FILE_SUFFIXES:
            continue
        if path.is_file():
            candidates.append(path)
    yield from sorted(
        candidates,
        key=lambda candidate: candidate.relative_to(repo_root).as_posix(),
    )


def create_package_inventory(
    *,
    repo_root: str | Path,
    output_path: str | Path,
    logger: logging.Logger | None = None,
) -> int:
    """Write a tab-separated package inventory.

    Args:
        repo_root: Repository directory.
        output_path: Destination TSV path.
        logger: Optional logger.

    Returns:
        Number of inventoried files.

    Raises:
        NotADirectoryError: If ``repo_root`` is not a directory.
    """
    root = Path(repo_root).expanduser().resolve()
    destination = Path(output_path).expanduser().resolve()
    if not root.is_dir():
        raise NotADirectoryError(f"Repository root is not a directory: {root}")

    records = []
    for path in iter_package_files(repo_root=root):
        if path.resolve() == destination:
            continue
        records.append(
            {
                "relative_path": path.relative_to(root).as_posix(),
                "size_bytes": path.stat().st_size,
                "sha256": sha256_file(path=path),
            }
        )
    write_records_table(
        records=records,
        output_path=destination,
        fieldnames=["relative_path", "size_bytes", "sha256"],
        logger=logger or LOGGER,
    )
    (logger or LOGGER).info(
        "Wrote package inventory with %d files: %s",
        len(records),
        destination,
    )
    return len(records)


def build_parser() -> argparse.ArgumentParser:
    """Build the package-inventory argument parser."""
    parser = argparse.ArgumentParser(
        description="Write a deterministic KmerSutra package inventory."
    )
    parser.add_argument("--repo-root", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--verbose", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    """Run the package-inventory command.

    Args:
        argv: Optional explicit command-line arguments.

    Returns:
        Process exit status.
    """
    args = build_parser().parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s [%(name)s] %(message)s",
    )
    create_package_inventory(
        repo_root=args.repo_root,
        output_path=args.output,
        logger=LOGGER,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
