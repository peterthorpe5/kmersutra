#!/usr/bin/env python3
"""Download NCBI genomes for KmerSutra.

This wrapper is retained for users who prefer direct script paths.
The recommended installed command is provided by the package entry point.
"""

from __future__ import annotations

from kmersutra.cli.download_ncbi_taxon_genomes import main


if __name__ == "__main__":
    main()
