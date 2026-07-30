"""Regression tests for the published KmerSutra command surface."""

from __future__ import annotations

import importlib
import unittest
from pathlib import Path


LEGACY_COMMANDS = {
    "kmersutra-benchmark-postprocess",
    "kmersutra-build-call-training",
    "kmersutra-build-panel",
    "kmersutra-download-genomes",
    "kmersutra-download-ncbi",
    "kmersutra-download-taxonomy",
    "kmersutra-evaluate-call-predictions",
    "kmersutra-extract-features",
    "kmersutra-label-zymo-calls",
    "kmersutra-merge-modules",
    "kmersutra-merge-panels",
    "kmersutra-predict",
    "kmersutra-screen",
    "kmersutra-summarise-comparable-benchmark",
    "kmersutra-summarise-lca",
    "kmersutra-summarise-run",
    "kmersutra-summarise-spikeins",
    "kmersutra-threshold-sweep",
    "kmersutra-train-call-calibrator",
    "kmersutra-train-classifier",
    "kmersutra-validate-call-calibrator",
    "kmersutra-validate-panel",
}

NEW_COMMANDS = {
    "kmersutra-audit-reference-panel",
    "kmersutra-label-mock-calls",
    "kmersutra-run-benchmark",
    "kmersutra-subset-reads",
    "kmersutra-summarise-mock-benchmark",
}


def read_project_scripts(path: Path) -> dict[str, str]:
    """Read the simple ``[project.scripts]`` mapping without TOML dependencies.

    Args:
        path: ``pyproject.toml`` path.

    Returns:
        Script-name to entry-point mapping.
    """
    scripts: dict[str, str] = {}
    in_scripts = False
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if line == "[project.scripts]":
            in_scripts = True
            continue
        if in_scripts and line.startswith("["):
            break
        if not in_scripts or not line or line.startswith("#"):
            continue
        name, separator, value = line.partition("=")
        if separator:
            scripts[name.strip()] = value.strip().strip('"')
    return scripts


class TestCliCompatibilityV051(unittest.TestCase):
    """Protect every v0.50.1 command from accidental removal."""

    def test_all_legacy_and_new_commands_are_registered(self) -> None:
        """The v0.51.0 script table should be a strict compatibility superset."""
        root = Path(__file__).resolve().parents[1]
        scripts = read_project_scripts(root / "pyproject.toml")
        self.assertFalse(LEGACY_COMMANDS.difference(scripts))
        self.assertFalse(NEW_COMMANDS.difference(scripts))

    def test_all_registered_entry_points_import(self) -> None:
        """Every entry-point module and callable should import."""
        root = Path(__file__).resolve().parents[1]
        scripts = read_project_scripts(root / "pyproject.toml")
        for name, entry_point in scripts.items():
            with self.subTest(command=name):
                module_name, callable_name = entry_point.split(":", 1)
                module = importlib.import_module(module_name)
                self.assertTrue(callable(getattr(module, callable_name)))

    def test_legacy_ai_table_aliases_remain_in_source(self) -> None:
        """Published table option aliases should remain available."""
        root = Path(__file__).resolve().parents[1]
        combined = "\n".join(
            (root / relative).read_text(encoding="utf-8")
            for relative in (
                "kmersutra/cli/build_call_training_table.py",
                "kmersutra/cli/train_call_calibrator.py",
            )
        )
        for alias in (
            "--calls_tsv",
            "--out_tsv",
            "--training_tsv",
            "--out_summary_tsv",
            "--out_evaluation_tsv",
        ):
            with self.subTest(alias=alias):
                self.assertIn(alias, combined)


if __name__ == "__main__":
    unittest.main()
