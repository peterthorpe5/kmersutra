"""CLI tests for the SQLite-backed target-evidence build path."""

from __future__ import annotations

import subprocess
import sys
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory


class TestBuildCliTargetEvidence(unittest.TestCase):
    """Tests for target-evidence-only CLI output files."""

    def test_build_cli_writes_target_evidence_outputs_and_ram_log(self) -> None:
        """CLI should produce panel, SQLite, profile and RAM outputs."""
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            target = root / "target.fna"
            neighbour = root / "neighbour.fna"
            target.write_text(">target\nAAAAACCCCCGGGGG\n", encoding="utf-8")
            neighbour.write_text(">neighbour\nAAAAATTTTT\n", encoding="utf-8")
            config = root / "genomes.tsv"
            config.write_text(
                "genome_fasta\tspecies_name\tstrain_name\ttaxid\trole\tclade\n"
                f"{target}\tAlpha target\tstrain1\t11\ttarget_species\tAlphaGenus\n"
                f"{neighbour}\tAlpha neighbour\tstrain2\t12\tnear_neighbour\tAlphaGenus\n",
                encoding="utf-8",
            )
            out_dir = root / "out"
            command = [
                sys.executable,
                "-m",
                "kmersutra.cli.build_clade_kmer_panel",
                "--genome_config",
                str(config),
                "--out_dir",
                str(out_dir),
                "--k_values",
                "5",
                "--target_evidence_only",
                "--profile",
                "--ram_log_interval_seconds",
                "0.01",
                "--verbose",
            ]
            completed = subprocess.run(
                command,
                check=False,
                capture_output=True,
                text=True,
            )
            if completed.returncode != 0:
                raise AssertionError(completed.stderr)
            self.assertTrue((out_dir / "species_kmer_panel.tsv.gz").is_file())
            self.assertTrue((out_dir / "target_evidence_candidates.sqlite").is_file())
            self.assertTrue((out_dir / "target_evidence_build_summary.tsv").is_file())
            self.assertTrue((out_dir / "build_profile_timing.tsv").is_file())
            self.assertTrue((out_dir / "ram_usage.tsv").is_file())


if __name__ == "__main__":
    unittest.main()
