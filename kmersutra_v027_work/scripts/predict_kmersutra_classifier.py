#!/usr/bin/env python3
"""Predict with a KmerSutra-ML classifier.

This wrapper is retained for users who prefer direct script paths.
The recommended installed command is provided by the package entry point.
"""

from __future__ import annotations

from kmersutra.cli.predict_classifier import main


if __name__ == "__main__":
    main()
