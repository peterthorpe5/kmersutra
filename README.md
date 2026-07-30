# KmerSutra

KmerSutra is a conservative, taxonomy-aware multi-k framework for
species-resolved metagenomic detection. It combines exact evidence across
multiple k-mer lengths, outgroup-aware marker construction, explicit
species/lineage reporting and optional interpretable evidence calibration.

The deterministic evidence and reporting layers remain the scientific core.
The optional AI model calibrates summary evidence; it does not replace the
rule-based calls and it does not perform per-read neural classification.

## Current release

Version 0.51.0 is the benchmark-ready release for locked independent
validation. It preserves all v0.50.1 commands and adds:

- generic mock-community truth manifests;
- truth-accession exclusion and reference-panel leakage auditing;
- deterministic depth-series generation;
- a restartable benchmark controller;
- a Slurm workflow for ATCC MSA-1003 HiFi, SRR9328980;
- continuous integration, coverage and wheel-content checks.

See the [changelog](https://github.com/peterthorpe5/kmersutra/blob/main/CHANGELOG.md)
and the [v0.51.0 release notes](https://github.com/peterthorpe5/kmersutra/blob/main/docs/releases/V0_51_0_RELEASE_NOTES.md).

## Installation

Create the Conda environment:

```bash
conda env create --file environment.yml
conda run --name kmersutra kmersutra-screen --help
```

Install into an existing environment:

```bash
python -m pip install --editable ".[dev,reporting,ncbi,ml,parquet]"
```

The package requires Python 3.10 or newer.

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

Prepare a site-specific configuration, then submit:

```bash
cp benchmarks/atcc_msa1003_hifi/config.example.json \
    benchmarks/atcc_msa1003_hifi/config.cluster.json

bash benchmarks/atcc_msa1003_hifi/submit_atcc_msa1003_benchmark.sh \
    --config benchmarks/atcc_msa1003_hifi/config.cluster.json \
    --output-root /home/pthorpe001/data/benchmarks/kmersutra_atcc_msa1003 \
    --account barton \
    --partition barton \
    --conda-env kmersutra \
    --resume
```

Full instructions are in
[ATCC_MSA1003_HIFI.md](https://github.com/peterthorpe5/kmersutra/blob/main/docs/benchmarking/ATCC_MSA1003_HIFI.md).

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

Run the compatibility test suite:

```bash
bash run_tests.sh
```

Run release checks:

```bash
bash run_quality_checks.sh
```

The release gate includes unit tests, a ratcheted 75% branch-coverage floor,
Ruff static checks, wheel construction and a wheel-content audit. The verified
v0.51.0 rebuild is at 78%; the next engineering target is at least 90%.

## Documentation

- [Installation and cluster deployment](https://github.com/peterthorpe5/kmersutra/blob/main/docs/INSTALLATION.md)
- [v0.51.0 Mac/GitHub/cluster commands](https://github.com/peterthorpe5/kmersutra/blob/main/docs/DEPLOYMENT_COMMANDS_v0_51_0.md)
- [v0.51.0 rebuild handoff](https://github.com/peterthorpe5/kmersutra/blob/main/docs/REBUILD_HANDOFF_v0_51_0.md)
- [Command compatibility](https://github.com/peterthorpe5/kmersutra/blob/main/docs/CLI_COMPATIBILITY.md)
- [Reproducibility model](https://github.com/peterthorpe5/kmersutra/blob/main/docs/REPRODUCIBILITY.md)
- [ATCC benchmark design](https://github.com/peterthorpe5/kmersutra/blob/main/docs/benchmarking/ATCC_MSA1003_HIFI.md)
- [Historical README](https://github.com/peterthorpe5/kmersutra/blob/main/docs/history/README_v0_50_1_legacy.md)

## Citation and licence

Citation metadata are provided in [CITATION.cff](https://github.com/peterthorpe5/kmersutra/blob/main/CITATION.cff). A permanent
software archive identifier should be added after the public v0.51.0 release.

An institutional licence decision is still required; see
[LICENCE_SELECTION_REQUIRED.md](https://github.com/peterthorpe5/kmersutra/blob/main/docs/LICENCE_SELECTION_REQUIRED.md). No open-source
licence is granted by this technical rebuild.
