"""Tests for panel building."""

import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from kmersutra.build_panel import build_panel
from kmersutra.config import GenomeConfig


class TestBuildPanel(unittest.TestCase):
    """Tests for diagnostic panel construction."""

    def test_build_panel_species_unique(self) -> None:
        """Panel builder should remove k-mers found in outgroups."""
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            alpha = root / "alpha.fna"
            beta = root / "beta.fna"
            outgroup = root / "outgroup.fna"
            alpha.write_text(">a\nAAAAACCCCC\n", encoding="utf-8")
            beta.write_text(">b\nGGGGGCCCCC\n", encoding="utf-8")
            outgroup.write_text(">o\nTTTTTCCCCC\n", encoding="utf-8")
            configs = [
                GenomeConfig(genome_fasta=alpha, species_name="Alpha", role="target_species", clade="Demo"),
                GenomeConfig(genome_fasta=beta, species_name="Beta", role="target_species", clade="Demo"),
                GenomeConfig(genome_fasta=outgroup, species_name="Out", role="outgroup", clade="Out"),
            ]
            diagnostics, summary, collection_summary = build_panel(
                genome_configs=configs,
                k_values=[5],
            )
        species = {item.species_name for item in diagnostics if item.panel_type == "species_unique"}
        self.assertIn("Alpha", species)
        self.assertIn("Beta", species)
        self.assertTrue(all(item.species_name != "Out" for item in diagnostics))
        self.assertTrue(summary)
        self.assertEqual(len(collection_summary), 3)

    def test_build_panel_threads_matches_single_worker(self) -> None:
        """Parallel panel building should match single-worker output."""
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            alpha = root / "alpha.fna"
            beta = root / "beta.fna"
            alpha.write_text(">a\nAAAAACCCCCGGGGG\n", encoding="utf-8")
            beta.write_text(">b\nTTTTTGGGGGCCCCC\n", encoding="utf-8")
            configs = [
                GenomeConfig(genome_fasta=alpha, species_name="Alpha", role="target_species", clade="Demo"),
                GenomeConfig(genome_fasta=beta, species_name="Beta", role="target_species", clade="Demo"),
            ]
            serial, _, _ = build_panel(genome_configs=configs, k_values=[5], threads=1)
            parallel, _, _ = build_panel(genome_configs=configs, k_values=[5], threads=2)
        serial_keys = {(item.k, item.kmer, item.species_name) for item in serial}
        parallel_keys = {(item.k, item.kmer, item.species_name) for item in parallel}
        self.assertEqual(serial_keys, parallel_keys)

    def test_build_panel_thinning_limit(self) -> None:
        """Panel thinning should limit retained k-mers per species and k."""
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            alpha = root / "alpha.fna"
            beta = root / "beta.fna"
            alpha.write_text(">a\nAAAAACCCCCGGGGGTTTTT\n", encoding="utf-8")
            beta.write_text(">b\nCCCCCGGGGGTTTTTAAAAA\n", encoding="utf-8")
            configs = [
                GenomeConfig(genome_fasta=alpha, species_name="Alpha", role="target_species", clade="Demo"),
                GenomeConfig(genome_fasta=beta, species_name="Beta", role="target_species", clade="Demo"),
            ]
            diagnostics, _, _ = build_panel(
                genome_configs=configs,
                k_values=[5],
                max_per_species_per_k=2,
            )
        counts = {}
        for item in diagnostics:
            key = (item.panel_type, item.species_name, item.k)
            counts[key] = counts.get(key, 0) + 1
        self.assertTrue(all(count <= 2 for count in counts.values()))

    def test_build_panel_clade_core_kmers(self) -> None:
        """Panel builder should optionally include clade-core k-mers."""
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            alpha = root / "alpha.fna"
            beta = root / "beta.fna"
            outgroup = root / "outgroup.fna"
            alpha.write_text(">a\nAAAAACCCCC\n", encoding="utf-8")
            beta.write_text(">b\nAAAAAGGGGG\n", encoding="utf-8")
            outgroup.write_text(">o\nTTTTTGGGGG\n", encoding="utf-8")
            configs = [
                GenomeConfig(genome_fasta=alpha, species_name="Alpha", role="target_species", clade="Demo"),
                GenomeConfig(genome_fasta=beta, species_name="Beta", role="target_species", clade="Demo"),
                GenomeConfig(genome_fasta=outgroup, species_name="Out", role="outgroup", clade="Other"),
            ]
            diagnostics, _, _ = build_panel(
                genome_configs=configs,
                k_values=[5],
                target_clade="Demo",
            )
        panel_types = {item.panel_type for item in diagnostics}
        self.assertIn("clade_core", panel_types)

class TestTaxonomicBuildPanel(unittest.TestCase):
    """Tests for taxonomy-aware panel construction."""

    def _write_taxdump(self, root: Path) -> None:
        """Write a tiny taxonomy dump for panel tests."""
        root.mkdir(parents=True, exist_ok=True)
        (root / "nodes.dmp").write_text(
            "1\t|\t1\t|\tno rank\t|\n"
            "2\t|\t1\t|\tsuperkingdom\t|\n"
            "10\t|\t2\t|\tgenus\t|\n"
            "11\t|\t10\t|\tspecies\t|\n"
            "12\t|\t10\t|\tspecies\t|\n"
            "20\t|\t2\t|\tgenus\t|\n"
            "21\t|\t20\t|\tspecies\t|\n",
            encoding="utf-8",
        )
        (root / "names.dmp").write_text(
            "1\t|\troot\t|\t\t|\tscientific name\t|\n"
            "2\t|\tPathogenia\t|\t\t|\tscientific name\t|\n"
            "10\t|\tAlphaGenus\t|\t\t|\tscientific name\t|\n"
            "11\t|\tAlpha one\t|\t\t|\tscientific name\t|\n"
            "12\t|\tAlpha two\t|\t\t|\tscientific name\t|\n"
            "20\t|\tBetaGenus\t|\t\t|\tscientific name\t|\n"
            "21\t|\tBeta one\t|\t\t|\tscientific name\t|\n",
            encoding="utf-8",
        )
        (root / "merged.dmp").write_text("", encoding="utf-8")
        (root / "delnodes.dmp").write_text("", encoding="utf-8")

    def test_taxonomy_aware_panel_builds_genus_evidence(self) -> None:
        """Shared species k-mers should become genus-level evidence."""
        from kmersutra.taxonomy import TaxonomyDatabase

        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            tax_root = root / "taxdump"
            self._write_taxdump(tax_root)
            alpha1 = root / "alpha1.fna"
            alpha2 = root / "alpha2.fna"
            beta = root / "beta.fna"
            alpha1.write_text(">a1\nAAAAACCCCC\n", encoding="utf-8")
            alpha2.write_text(">a2\nAAAAAGGGGG\n", encoding="utf-8")
            beta.write_text(">b\nCATGCCATGC\n", encoding="utf-8")
            taxonomy = TaxonomyDatabase.from_taxdump(taxonomy_dir=tax_root)
            configs = [
                GenomeConfig(genome_fasta=alpha1, species_name="Alpha one", taxid="11", role="target_species", clade="AlphaGenus"),
                GenomeConfig(genome_fasta=alpha2, species_name="Alpha two", taxid="12", role="target_species", clade="AlphaGenus"),
                GenomeConfig(genome_fasta=beta, species_name="Beta one", taxid="21", role="outgroup", clade="BetaGenus"),
            ]
            diagnostics, summary, _ = build_panel(
                genome_configs=configs,
                k_values=[5],
                taxonomy_db=taxonomy,
                target_taxid="10",
            )
        panel_types = {item.panel_type for item in diagnostics}
        self.assertIn("genus_core", panel_types)
        self.assertTrue(any(item.evidence_rank == "genus" for item in diagnostics))
        self.assertTrue(summary)

    def test_taxonomy_aware_panel_filters_outside_target_taxid(self) -> None:
        """Target taxid filtering should remove unrelated outgroup evidence."""
        from kmersutra.taxonomy import TaxonomyDatabase

        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            tax_root = root / "taxdump"
            self._write_taxdump(tax_root)
            alpha = root / "alpha.fna"
            beta = root / "beta.fna"
            alpha.write_text(">a\nAAAAACCCCC\n", encoding="utf-8")
            beta.write_text(">b\nTTTTTGGGGG\n", encoding="utf-8")
            taxonomy = TaxonomyDatabase.from_taxdump(taxonomy_dir=tax_root)
            configs = [
                GenomeConfig(genome_fasta=alpha, species_name="Alpha one", taxid="11", role="target_species", clade="AlphaGenus"),
                GenomeConfig(genome_fasta=beta, species_name="Beta one", taxid="21", role="outgroup", clade="BetaGenus"),
            ]
            diagnostics, _, _ = build_panel(
                genome_configs=configs,
                k_values=[5],
                taxonomy_db=taxonomy,
                target_taxid="10",
            )
        self.assertTrue(all(item.evidence_taxid != "21" for item in diagnostics))
        self.assertTrue(all(item.evidence_taxid in {"10", "11"} for item in diagnostics))
