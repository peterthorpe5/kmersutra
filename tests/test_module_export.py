"""Tests for automatic hierarchical module export."""

from __future__ import annotations

import gzip
import tempfile
import unittest
from pathlib import Path

from kmersutra.hierarchical import load_module_manifest
from kmersutra.io import read_tsv
from kmersutra.module_export import (
    ModuleExportConfig,
    export_hierarchical_modules_from_panel,
    infer_genus_from_record,
    safe_module_id,
)

PANEL_COLUMNS = [
    "kmer",
    "k",
    "panel_type",
    "species_name",
    "clade",
    "source_genomes",
    "source_contigs",
    "example_position",
    "evidence_taxid",
    "evidence_name",
    "evidence_rank",
    "lineage_taxids",
    "source_taxids",
]


def write_panel(path: Path, records: list[dict[str, object]]) -> None:
    """Write a small gzipped KmerSutra panel for testing."""
    with gzip.open(path, "wt") as handle:
        handle.write("\t".join(PANEL_COLUMNS) + "\n")
        for record in records:
            handle.write(
                "\t".join(str(record.get(column, "")) for column in PANEL_COLUMNS)
                + "\n"
            )


class TestAutomaticModuleExport(unittest.TestCase):
    """Tests for automatic module manifest generation."""

    def test_safe_module_id_normalises_labels(self) -> None:
        """Module identifiers should be stable and shell/path safe."""
        observed = safe_module_id(value="Plasmodium cf. malariae", prefix="genus")
        self.assertEqual(observed, "genus_plasmodium_cf._malariae")

    def test_infer_genus_from_record_uses_taxonomic_context(self) -> None:
        """Genus inference should use genus evidence before species labels."""
        self.assertEqual(
            infer_genus_from_record(
                record={
                    "evidence_rank": "genus",
                    "evidence_name": "Plasmodium",
                    "species_name": "",
                    "clade": "Apicomplexa",
                }
            ),
            "Plasmodium",
        )
        self.assertEqual(
            infer_genus_from_record(
                record={
                    "evidence_rank": "species",
                    "evidence_name": "Plasmodium vivax",
                    "species_name": "Plasmodium vivax",
                    "clade": "Plasmodium",
                }
            ),
            "Plasmodium",
        )

    def test_export_hierarchical_modules_writes_loadable_manifest(self) -> None:
        """Export should write gate/detail panels and a loadable manifest."""
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            panel_path = tmp_path / "panel.tsv.gz"
            write_panel(
                panel_path,
                [
                    {
                        "kmer": "A" * 77,
                        "k": 77,
                        "panel_type": "genus_core",
                        "species_name": "",
                        "clade": "Plasmodium",
                        "source_genomes": "g1;g2",
                        "source_contigs": "c1",
                        "example_position": 10,
                        "evidence_taxid": "5820",
                        "evidence_name": "Plasmodium",
                        "evidence_rank": "genus",
                        "lineage_taxids": "1;5820",
                        "source_taxids": "5855;5833",
                    },
                    {
                        "kmer": "C" * 77,
                        "k": 77,
                        "panel_type": "species_unique",
                        "species_name": "Plasmodium vivax",
                        "clade": "Plasmodium",
                        "source_genomes": "g1",
                        "source_contigs": "c1",
                        "example_position": 100,
                        "evidence_taxid": "5855",
                        "evidence_name": "Plasmodium vivax",
                        "evidence_rank": "species",
                        "lineage_taxids": "1;5820;5855",
                        "source_taxids": "5855",
                    },
                    {
                        "kmer": "G" * 77,
                        "k": 77,
                        "panel_type": "species_unique",
                        "species_name": "Babesia bovis",
                        "clade": "Babesia",
                        "source_genomes": "g3",
                        "source_contigs": "c1",
                        "example_position": 20,
                        "evidence_taxid": "5865",
                        "evidence_name": "Babesia bovis",
                        "evidence_rank": "species",
                        "lineage_taxids": "1;5864;5865",
                        "source_taxids": "5865",
                    },
                ],
            )
            result = export_hierarchical_modules_from_panel(
                panel_path=panel_path,
                out_dir=tmp_path / "modules",
                config=ModuleExportConfig(gate_ranks={"genus", "family"}),
            )

            self.assertTrue(result.manifest_path.is_file())
            self.assertTrue(result.summary_path.is_file())
            records = read_tsv(input_path=result.manifest_path)
            module_ids = {record["module_id"] for record in records}
            self.assertIn("clade_plasmodium", module_ids)
            self.assertIn("genus_plasmodium", module_ids)
            self.assertIn("genus_babesia", module_ids)

            loaded = load_module_manifest(manifest_path=result.manifest_path)
            self.assertEqual(len(loaded), len(records))
            plasmodium = [m for m in loaded if m.module_id == "genus_plasmodium"][0]
            self.assertTrue(Path(plasmodium.gate_panel_path).is_file())
            self.assertTrue(Path(plasmodium.module_panel_path).is_file())
            detail_records = read_tsv(input_path=plasmodium.module_panel_path)
            self.assertTrue(
                any(record["species_name"] == "Plasmodium vivax" for record in detail_records)
            )

    def test_export_rejects_bad_gate_thresholds(self) -> None:
        """Invalid module export thresholds should fail before writing output."""
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            panel_path = tmp_path / "panel.tsv.gz"
            write_panel(
                panel_path,
                [
                    {
                        "kmer": "A" * 77,
                        "k": 77,
                        "panel_type": "genus_core",
                        "species_name": "",
                        "clade": "Plasmodium",
                        "source_genomes": "g1",
                        "source_contigs": "c1",
                        "example_position": 10,
                        "evidence_taxid": "5820",
                        "evidence_name": "Plasmodium",
                        "evidence_rank": "genus",
                        "lineage_taxids": "1;5820",
                        "source_taxids": "5855",
                    }
                ],
            )
            with self.assertRaises(ValueError):
                export_hierarchical_modules_from_panel(
                    panel_path=panel_path,
                    out_dir=tmp_path / "modules",
                    config=ModuleExportConfig(
                        gate_ranks={"genus"},
                        min_gate_unique_kmers=0,
                    ),
                )


if __name__ == "__main__":
    unittest.main()
