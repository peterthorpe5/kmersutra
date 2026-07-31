# KmerSutra

KmerSutra is a conservative, taxonomy-aware multi-k framework for
species-resolved metagenomic detection. It combines exact evidence across
multiple k-mer lengths, outgroup-aware marker construction, explicit
species/lineage reporting and optional interpretable evidence calibration.

It was developed in response to a practical problem identified in an ONT
Plasmodium benchmark: Kraken2 and Metabuli detected low-abundance targets
sensitively, but their positive reports also contained numerous non-expected
Plasmodium labels. KmerSutra tests a different operating point, prioritising
clean, inspectable species-level evidence rather than maximal read assignment
or earliest-read detection. It complements general-purpose classifiers rather
than replacing them.

The deterministic evidence and reporting layers remain the scientific core.
The optional AI model calibrates summary evidence; it does not replace the
rule-based calls and it does not perform per-read neural classification.

**Documentation:** [browsable user guide](https://peterthorpe5.github.io/kmersutra/)
· [offline HTML](docs/html/index.html)
· [installation](docs/source/installation.rst)
· [command reference](docs/source/cli_reference.rst)

## Current release

Version 0.51.2 is the ATCC reference-preparation repair release. It preserves
the installed v0.50.1 command surface and adds:

- generic mock-community truth manifests;
- truth-accession exclusion and reference-panel leakage auditing;
- deterministic depth-series generation;
- a restartable benchmark controller;
- a Slurm workflow for ATCC MSA-1003 HiFi, SRR9328980;
- restored, Slurm-native ONT Zymo D6300 screening;
- explicit ATCC `--database-root` and `--ai-model` options;
- early rejection of unresolved configuration placeholders;
- dated test, coverage and package-inventory result directories;
- Read the Docs-style HTML documentation;
- continuous integration, coverage, documentation and wheel-content checks.
- a leakage-controlled ATCC reference builder using held-out target genomes,
  genus-level near neighbours and the established v4 collection as background;
- a tested Python reference gate that prevents screening when any ATCC target
  is missing or any published truth accession leaks into the panel.

See the [changelog](https://github.com/peterthorpe5/kmersutra/blob/main/CHANGELOG.md)
and the [v0.51.2 release notes](docs/releases/V0_51_2_RELEASE_NOTES.md). The
[overlay and cluster handoff](docs/DEPLOYMENT_COMMANDS_v0_51_2.md) uses named
paths and keeps generated validation artefacts outside the repository root.

## Installation

Create the Conda environment:

```bash
conda env create --file environment.yml
conda run --name kmersutra kmersutra-screen --help
```

Install into an existing environment:

```bash
python -m pip install --editable ".[all,dev,docs]"
```

The package supports Python 3.10–3.12. See the
[complete installation guide](docs/source/installation.rst).

## Core workflow

Build a panel:

```bash
kmersutra-build-panel \
    --genome_config genome_config.tsv \
    --out_dir panel_build \
    --k_values 51 77 101 151 \
    --global_candidate_evidence \
    --global_source_index_mode candidate_universe \
    --marker_selection independent_multik_genome_spread \
    --min_cross_k_marker_distance 5000 \
    --threads 8 \
    --profile \
    --verbose
```

Screen a sample using exact matching:

```bash
kmersutra-screen \
    --input sample.fastq.gz \
    --input_format fastq \
    --panel panel_build/species_kmer_panel.tsv.gz \
    --sample_id sample_001 \
    --out_dir results/sample_001 \
    --screen_mode flat \
    --screen_preset exact \
    --call_preset lineage_aware \
    --consolidate_species_calls \
    --same_genus_reportable_min_fraction 0.05 \
    --threads 8 \
    --decompressor auto \
    --use_panel_cache \
    --no_read_level_hits \
    --write_parquet_outputs \
    --profile \
    --verbose
```

All primary tables are TSV, TSV.GZ or Parquet. CSV output is deliberately not
used.

## Independent ATCC benchmark

The benchmark specification is under
[`benchmarks/atcc_msa1003_hifi`](https://github.com/peterthorpe5/kmersutra/tree/main/benchmarks/atcc_msa1003_hifi).
Its primary settings are frozen before examining SRR9328980:

- exact matching;
- k = 51, 77, 101 and 151;
- same-genus reportability fraction = 0.05;
- AI novelty scale = 2.9 without retraining;
- truth reference sequences excluded from the panel;
- flat-panel result primary;
- hierarchical result secondary;
- no primary threshold sweep.

Pass site-specific paths explicitly. The example JSON remains portable and is
not edited to contain private cluster paths:

```bash
DATABASE_ROOT="/absolute/path/to/kmersutra_db"
AI_MODEL="/absolute/path/to/final_internal_calibrator_all_training.json"
OUTPUT_ROOT="/absolute/path/to/kmersutra_atcc_msa1003"

bash benchmarks/atcc_msa1003_hifi/submit_atcc_msa1003_benchmark.sh \
    --config benchmarks/atcc_msa1003_hifi/config.example.json \
    --database-root "${DATABASE_ROOT}" \
    --ai-model "${AI_MODEL}" \
    --output-root "${OUTPUT_ROOT}" \
    --account barton \
    --partition barton \
    --conda-env kmersutra \
    --resume
```

Full instructions are in
[the ATCC guide](docs/source/benchmarks/atcc_msa1003.rst).

## Existing ONT Zymo benchmark

The public ONT Q20 Zymo D6300 benchmark (`ERR5396170`) remains part of the
package and manuscript evidence. Its original qsub wrapper was stored inside
the cluster dataset directory rather than GitHub. Version 0.51.1 provides a
current Slurm replacement:

```bash
bash benchmarks/zymo_d6300_ont/submit_zymo_d6300_screen.sh \
    --reads /absolute/path/to/ERR5396170.fastq.gz \
    --panel /absolute/path/to/species_kmer_panel.tsv.gz \
    --output-root /absolute/path/to/zymo_d6300_screen \
    --min-mixed 0.05 \
    --account barton \
    --partition barton \
    --threads 24 \
    --conda-env kmersutra
```

See the [Zymo benchmark guide](docs/source/benchmarks/zymo_d6300.rst).

## Manifest-driven truth labelling

Label calls from any mock community:

```bash
kmersutra-label-mock-calls \
    --calls_table species_detection_calls.tsv \
    --truth_manifest mock_community_truth.tsv \
    --out_table labelled_calls.tsv.gz \
    --out_category_counts truth_category_counts.tsv \
    --out_coarse_label_counts coarse_label_counts.tsv \
    --verbose
```

The historical `kmersutra-label-zymo-calls` command remains available and
unchanged.

## Reference leakage audit

Before an external benchmark:

```bash
kmersutra-audit-reference-panel \
    --genome_config kmersutra_genome_config.tsv \
    --truth_manifest mock_community_truth.tsv \
    --out_table reference_panel_audit.tsv \
    --fail_on_leakage \
    --verbose
```

The audit checks exact assembly accessions, truth sequence accessions in FASTA
headers and exact FASTA SHA-256 values where the truth manifest provides a
truth FASTA.

## Testing and quality checks

Run the compatibility test suite and keep its logs beneath a dedicated result
directory:

```bash
bash run_tests.sh \
    --results-dir /absolute/project/path/kmersutra_test_results
```

Run release checks, including branch coverage and an HTML coverage browser:

```bash
bash run_quality_checks.sh \
    --results-dir /absolute/project/path/kmersutra_test_results
```

Each run writes its logs, coverage reports and package inventory below one dated
directory. Generated `test_results_*` and `PACKAGE_FILE_INVENTORY_*` files no
longer accumulate in the repository root or home directory.

Historical SGE and AI-validation wrappers are retained for provenance beneath
[`scripts/legacy`](scripts/legacy/README.md). Current Slurm launchers live with
their benchmark specifications under `benchmarks/`.

## Documentation

Build the Read the Docs-style site locally:

```bash
python -m pip install --editable ".[docs]"
bash scripts/build_documentation.sh --output-dir docs/html
```

Then open `docs/html/index.html`. The repository also includes a Read the Docs
configuration and a GitHub Actions documentation build.

The guide covers installation, scientific concepts, panel construction,
screening, outputs, Slurm, reproducibility, troubleshooting, ATCC, Zymo, the
command catalogue and the Python API.

## Citation and licence

Citation metadata are provided in [CITATION.cff](https://github.com/peterthorpe5/kmersutra/blob/main/CITATION.cff). A permanent
software archive identifier should be added after the public release.

An institutional licence decision is still required; see
[LICENCE_SELECTION_REQUIRED.md](https://github.com/peterthorpe5/kmersutra/blob/main/docs/LICENCE_SELECTION_REQUIRED.md). No open-source
licence is granted by this technical rebuild.
