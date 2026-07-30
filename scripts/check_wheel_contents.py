"""Validate that a built KmerSutra wheel contains only intended packages."""

from __future__ import annotations

import argparse
import logging
import zipfile
from pathlib import Path

LOGGER = logging.getLogger("kmersutra.wheel_check")


def build_parser() -> argparse.ArgumentParser:
    """Build the command-line parser."""
    parser = argparse.ArgumentParser(
        description="Check KmerSutra wheel contents for development artefacts."
    )
    parser.add_argument("--dist-dir", required=True)
    parser.add_argument("--verbose", action="store_true")
    return parser


def find_wheel(*, dist_dir: str | Path) -> Path:
    """Resolve exactly one built wheel.

    Args:
        dist_dir: Distribution directory.

    Returns:
        Wheel path.

    Raises:
        ValueError: If the directory does not contain exactly one wheel.
    """
    directory = Path(dist_dir).expanduser().resolve()
    wheels = sorted(directory.glob("kmersutra-*.whl"))
    if len(wheels) != 1:
        raise ValueError(
            f"Expected exactly one KmerSutra wheel in {directory}; "
            f"found {len(wheels)}"
        )
    return wheels[0]


def check_wheel(*, wheel_path: str | Path) -> list[str]:
    """Return prohibited paths found in a wheel.

    Args:
        wheel_path: Wheel archive.

    Returns:
        Prohibited archive members.

    Raises:
        ValueError: If required KmerSutra package files are absent.
    """
    path = Path(wheel_path).expanduser().resolve()
    with zipfile.ZipFile(path) as archive:
        members = archive.namelist()
    if "kmersutra/__init__.py" not in members:
        raise ValueError(f"Wheel does not contain kmersutra/__init__.py: {path}")
    prohibited_prefixes = (
        "tests/",
        "kmersutra_ai_validation/",
        "__pycache__/",
    )
    prohibited_suffixes = (".pyc", ".pyo", ".tmp", ".DS_Store")
    return sorted(
        member
        for member in members
        if member.startswith(prohibited_prefixes)
        or member.endswith(prohibited_suffixes)
    )


def main(argv: list[str] | None = None) -> int:
    """Run the wheel-content check.

    Args:
        argv: Optional explicit argument list.

    Returns:
        Process exit status.
    """
    args = build_parser().parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s [%(name)s] %(message)s",
    )
    wheel = find_wheel(dist_dir=args.dist_dir)
    prohibited = check_wheel(wheel_path=wheel)
    if prohibited:
        raise ValueError(
            "Wheel contains prohibited development artefacts: "
            + ", ".join(prohibited)
        )
    LOGGER.info("Wheel content check passed: %s", wheel)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
