"""Tests for compact scalable panel building."""

from __future__ import annotations

import sys
import unittest
from unittest.mock import patch
from pathlib import Path
from tempfile import TemporaryDirectory

from kmersutra.build_panel import (
    collect_compact_kmer_groups,
    build_panel,
    identify_clade_core_kmers_from_groups,
    identify_species_unique_kmers_from_groups,
    merge_compact_kmer_groups,
)
from kmersutra.config import GenomeConfig
from kmersutra.taxonomy import TaxonomyDatabase


class TestCompactBuild(unittest.TestCase):
    """Tests for compact k-mer grouping and panel building."""

    def test_collect_compact_groups_merges_repeated_kmers(self) -> None:
        """Compact collection should store one group per distinct k-mer key."""
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            fasta = root / "alpha.fna"
            fasta.write_text(">a\nAAAAAAA\n", encoding="utf-8")
            configs = [
                GenomeConfig(
                    genome_fasta=fasta,
                    species_name="Alpha",
                    role="target_species",
                    clade="Demo",
                )
            ]
            groups, summary = collect_compact_kmer_groups(
                genome_configs=configs,
                k_values=[3],
            )
        self.assertEqual(len(groups), 1)
        self.assertEqual(summary[0]["total_observations"], 5)
        self.assertEqual(summary[0]["distinct_kmer_keys"], 1)

    def test_merge_compact_groups_combines_species_metadata(self) -> None:
        """Merging compact groups should preserve source species metadata."""
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            alpha = root / "alpha.fna"
            beta = root / "beta.fna"
            alpha.write_text(">a\nAAAAA\n", encoding="utf-8")
            beta.write_text(">b\nAAAAA\n", encoding="utf-8")
            alpha_groups, _ = collect_compact_kmer_groups(
                genome_configs=[
                    GenomeConfig(
                        genome_fasta=alpha,
                        species_name="Alpha",
                        role="target_species",
                        clade="Demo",
                    )
                ],
                k_values=[5],
            )
            beta_groups, _ = collect_compact_kmer_groups(
                genome_configs=[
                    GenomeConfig(
                        genome_fasta=beta,
                        species_name="Beta",
                        role="target_species",
                        clade="Demo",
                    )
                ],
                k_values=[5],
            )
            merged = merge_compact_kmer_groups(left=alpha_groups, right=beta_groups)
        only_group = next(iter(merged.values()))
        self.assertEqual(only_group.species_names, {"Alpha", "Beta"})

    def test_compact_species_unique_matches_legacy_builder(self) -> None:
        """Compact and legacy builders should retain the same species-level keys."""
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            alpha = root / "alpha.fna"
            beta = root / "beta.fna"
            outgroup = root / "out.fna"
            alpha.write_text(">a\nAAAAACCCCCGGGGG\n", encoding="utf-8")
            beta.write_text(">b\nTTTTTGGGGGCCCCC\n", encoding="utf-8")
            outgroup.write_text(">o\nCCCCCAAAAA\n", encoding="utf-8")
            configs = [
                GenomeConfig(alpha, "Alpha", role="target_species", clade="Demo"),
                GenomeConfig(beta, "Beta", role="target_species", clade="Demo"),
                GenomeConfig(outgroup, "Out", role="outgroup", clade="Other"),
            ]
            legacy, _, _ = build_panel(
                genome_configs=configs,
                k_values=[5],
                compact_build=False,
            )
            compact, _, _ = build_panel(
                genome_configs=configs,
                k_values=[5],
                compact_build=True,
            )
        legacy_keys = {
            (item.panel_type, item.k, item.kmer, item.species_name) for item in legacy
        }
        compact_keys = {
            (item.panel_type, item.k, item.kmer, item.species_name) for item in compact
        }
        self.assertEqual(compact_keys, legacy_keys)

    def test_compact_group_identifies_clade_core(self) -> None:
        """Compact clade-core detection should retain k-mers restricted to a clade."""
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            alpha = root / "alpha.fna"
            beta = root / "beta.fna"
            outgroup = root / "out.fna"
            alpha.write_text(">a\nAAAAACCCCC\n", encoding="utf-8")
            beta.write_text(">b\nAAAAAGGGGG\n", encoding="utf-8")
            outgroup.write_text(">o\nTTTTTGGGGG\n", encoding="utf-8")
            configs = [
                GenomeConfig(alpha, "Alpha", role="target_species", clade="Demo"),
                GenomeConfig(beta, "Beta", role="target_species", clade="Demo"),
                GenomeConfig(outgroup, "Out", role="outgroup", clade="Other"),
            ]
            groups, _ = collect_compact_kmer_groups(
                genome_configs=configs,
                k_values=[5],
            )
            diagnostics = identify_clade_core_kmers_from_groups(
                compact_groups=groups,
                target_clade="Demo",
            )
        self.assertTrue(any(item.panel_type == "clade_core" for item in diagnostics))

    def test_compact_species_unique_function_returns_targets_only(self) -> None:
        """Compact species-unique detection should not retain outgroup-only k-mers."""
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            alpha = root / "alpha.fna"
            outgroup = root / "out.fna"
            alpha.write_text(">a\nAACCGGTTAACC\n", encoding="utf-8")
            outgroup.write_text(">o\nTTTTTTGGGGGG\n", encoding="utf-8")
            configs = [
                GenomeConfig(alpha, "Alpha", role="target_species", clade="Demo"),
                GenomeConfig(outgroup, "Out", role="outgroup", clade="Other"),
            ]
            groups, _ = collect_compact_kmer_groups(
                genome_configs=configs,
                k_values=[5],
            )
            diagnostics = identify_species_unique_kmers_from_groups(
                compact_groups=groups,
                target_species={"Alpha"},
            )
        self.assertTrue(diagnostics)
        self.assertTrue(all(item.species_name == "Alpha" for item in diagnostics))

    def _write_taxdump(self, root: Path) -> None:
        """Write a tiny taxonomy dump for compact taxonomy tests."""
        root.mkdir(parents=True, exist_ok=True)
        (root / "nodes.dmp").write_text(
            "1\t|\t1\t|\tno rank\t|\n"
            "2\t|\t1\t|\tsuperkingdom\t|\n"
            "10\t|\t2\t|\tgenus\t|\n"
            "11\t|\t10\t|\tspecies\t|\n"
            "12\t|\t10\t|\tspecies\t|\n",
            encoding="utf-8",
        )
        (root / "names.dmp").write_text(
            "1\t|\troot\t|\t\t|\tscientific name\t|\n"
            "2\t|\tPathogenia\t|\t\t|\tscientific name\t|\n"
            "10\t|\tAlphaGenus\t|\t\t|\tscientific name\t|\n"
            "11\t|\tAlpha one\t|\t\t|\tscientific name\t|\n"
            "12\t|\tAlpha two\t|\t\t|\tscientific name\t|\n",
            encoding="utf-8",
        )
        (root / "merged.dmp").write_text("", encoding="utf-8")
        (root / "delnodes.dmp").write_text("", encoding="utf-8")

    def test_compact_taxonomic_build_retains_genus_evidence(self) -> None:
        """Compact taxonomy-aware build should classify shared k-mers by genus."""
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            taxdump = root / "taxdump"
            self._write_taxdump(taxdump)
            alpha1 = root / "alpha1.fna"
            alpha2 = root / "alpha2.fna"
            alpha1.write_text(">a1\nAAAAACCCCC\n", encoding="utf-8")
            alpha2.write_text(">a2\nAAAAAGGGGG\n", encoding="utf-8")
            taxonomy = TaxonomyDatabase.from_taxdump(taxonomy_dir=taxdump)
            configs = [
                GenomeConfig(
                    alpha1,
                    "Alpha one",
                    taxid="11",
                    role="target_species",
                    clade="AlphaGenus",
                ),
                GenomeConfig(
                    alpha2,
                    "Alpha two",
                    taxid="12",
                    role="target_species",
                    clade="AlphaGenus",
                ),
            ]
            diagnostics, _, _ = build_panel(
                genome_configs=configs,
                k_values=[5],
                taxonomy_db=taxonomy,
                target_taxid="10",
                compact_build=True,
            )
        self.assertTrue(any(item.evidence_rank == "genus" for item in diagnostics))


class TestCompactBuildCli(unittest.TestCase):
    """Tests for compact build command-line options."""

    def test_build_panel_cli_writes_profile_file(self) -> None:
        """Build CLI should write profile timings when requested."""
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            alpha = root / "alpha.fna"
            beta = root / "beta.fna"
            config = root / "genomes.tsv"
            out_dir = root / "panel"
            alpha.write_text(">a\nAAAAACCCCC\n", encoding="utf-8")
            beta.write_text(">b\nGGGGGTTTTT\n", encoding="utf-8")
            config.write_text(
                "genome_fasta\tspecies_name\trole\tclade\n"
                f"{alpha}\tAlpha\ttarget_species\tDemo\n"
                f"{beta}\tBeta\ttarget_species\tDemo\n",
                encoding="utf-8",
            )
            from kmersutra.cli.build_clade_kmer_panel import main

            argv = [
                "kmersutra-build-panel",
                "--genome_config",
                str(config),
                "--out_dir",
                str(out_dir),
                "--k_values",
                "5",
                "--compact_build",
                "--profile",
            ]
            with patch.object(sys, "argv", argv):
                main()
            self.assertTrue((out_dir / "build_profile_timing.tsv").exists())
            self.assertTrue((out_dir / "species_kmer_panel.tsv.gz").exists())
