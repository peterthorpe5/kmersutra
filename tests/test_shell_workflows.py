"""Tests for repository shell workflow entry points."""

from __future__ import annotations

import subprocess
import unittest
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]


class TestShellWorkflows(unittest.TestCase):
    """Check syntax and key portable command-line interfaces."""

    def test_all_shell_workflows_pass_bash_syntax_check(self) -> None:
        """Every retained shell entry point should parse with Bash."""
        shell_paths = sorted(
            {
                *REPOSITORY_ROOT.rglob("*.sh"),
                *REPOSITORY_ROOT.rglob("*.slurm"),
            }
        )
        self.assertGreater(len(shell_paths), 0)
        for path in shell_paths:
            with self.subTest(path=path.relative_to(REPOSITORY_ROOT)):
                completed = subprocess.run(
                    ["bash", "-n", str(path)],
                    check=False,
                    capture_output=True,
                    text=True,
                )
                self.assertEqual(completed.returncode, 0, completed.stderr)

    def test_zymo_submit_help_documents_named_paths(self) -> None:
        """The restored Zymo submitter should expose named path options."""
        path = (
            REPOSITORY_ROOT
            / "benchmarks"
            / "zymo_d6300_ont"
            / "submit_zymo_d6300_screen.sh"
        )
        completed = subprocess.run(
            ["bash", str(path), "--help"],
            check=False,
            capture_output=True,
            text=True,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertIn("--reads FILE", completed.stdout)
        self.assertIn("--panel FILE", completed.stdout)
        self.assertIn("--output-root DIR", completed.stdout)

    def test_atcc_submit_help_requires_explicit_model_and_database(self) -> None:
        """The ATCC submitter should expose portable path arguments."""
        path = (
            REPOSITORY_ROOT
            / "benchmarks"
            / "atcc_msa1003_hifi"
            / "submit_atcc_msa1003_benchmark.sh"
        )
        completed = subprocess.run(
            ["bash", str(path), "--help"],
            check=False,
            capture_output=True,
            text=True,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertIn("--database-root DIR", completed.stdout)
        self.assertIn("--ai-model FILE", completed.stdout)


if __name__ == "__main__":
    unittest.main()
