"""Sphinx configuration for the KmerSutra documentation."""

from __future__ import annotations

import sys
from pathlib import Path

REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPOSITORY_ROOT))

from kmersutra import __version__  # noqa: E402

project = "KmerSutra"
author = "Peter Thorpe"
copyright = "2026, Peter Thorpe"
version = __version__
release = __version__

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.intersphinx",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
source_suffix = ".rst"
master_doc = "index"

html_theme = "sphinx_rtd_theme"
html_title = f"KmerSutra {release}"
html_show_sourcelink = True
html_context = {
    "display_github": True,
    "github_user": "peterthorpe5",
    "github_repo": "kmersutra",
    "github_version": "main",
    "conf_py_path": "/docs/source/",
}

napoleon_google_docstring = True
napoleon_numpy_docstring = True
autodoc_typehints = "description"
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
}
