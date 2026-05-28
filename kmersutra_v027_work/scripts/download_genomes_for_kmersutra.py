#!/usr/bin/env python3
"""Compatibility wrapper for kmersutra-download-genomes."""

from __future__ import annotations

from kmersutra.cli.download_genomes_for_kmersutra import main


if __name__ == "__main__":
    raise SystemExit(main())
