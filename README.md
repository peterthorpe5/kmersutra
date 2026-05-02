# KmerSutra

KmerSutra is an outgroup-aware k-mer framework for clade- and species-resolved metagenomic detection. The first intended application is species-resolved *Plasmodium* detection in Oxford Nanopore Technologies metagenomic reads and assemblies, but the code is designed for any clade of interest.

The package now includes two linked layers:

1. **KmerSutra rule-based detection**: builds diagnostic k-mer panels and screens reads or assemblies for clade-level and species-level evidence.
2. **KmerSutra-ML**: extracts interpretable feature tables from diagnostic k-mer hits, trains a lightweight open-set classifier, and predicts known-like versus unknown/unresolved sequence evidence.

The ML layer is intentionally dependency-light and interpretable. It is not a black-box deep-learning model. It provides a reproducible baseline for later, more complex models.

## Core idea

KmerSutra builds diagnostic k-mer panels from target genomes and outgroup genomes. It then screens reads or assemblies for evidence that supports:

- clade-level detection
- species-level resolution
- unresolved clade-level signal
- conflicting species evidence
- open-set unknown or unresolved calls

The first implementation uses exact canonical k-mer matching. Optional fuzzy matching by Hamming distance is available for long k-mers, but should be benchmarked carefully before being used for final biological interpretation.

## Installation

For development:

```bash
pip install -e '.[dev]'
```

For all optional extras:

```bash
pip install -e '.[all]'
```

The installed console commands are:

```text
kmersutra-build-panel
kmersutra-screen
kmersutra-download-ncbi
kmersutra-extract-features
kmersutra-train-classifier
kmersutra-predict
```

## Genome configuration

The genome configuration is a tab-separated file.

Required columns:

- `genome_fasta`
- `species_name`
- `role`

Recommended columns:

- `strain_name`
- `taxid`
- `assembly_accession`
- `clade`
- `source`

Supported roles:

- `target_species`
- `target_clade_member`
- `near_neighbour`
- `distant_outgroup`
- `background_pathogen`
- `host`
- `outgroup`
- `target`
- `exclude`

## Rule-based workflow

Build a diagnostic panel:

```bash
kmersutra-build-panel \
  --genome_config examples/example_genome_config.tsv \
  --out_dir panel_out \
  --k_values 71 101 \
  --target_clade Plasmodium \
  --threads 4 \
  --verbose
```

Screen reads or assemblies:

```bash
kmersutra-screen \
  --input sample.fastq.gz \
  --input_format fastq \
  --panel panel_out/species_kmer_panel.tsv.gz \
  --sample_id sample_001 \
  --out_dir screen_out \
  --threads 4 \
  --chunk_size 1000 \
  --max_mismatches 0 \
  --verbose
```

For fuzzy long-k testing:

```bash
kmersutra-screen \
  --input sample.fastq.gz \
  --input_format fastq \
  --panel panel_out/species_kmer_panel.tsv.gz \
  --sample_id sample_001_fuzzy \
  --out_dir screen_out_fuzzy \
  --threads 4 \
  --chunk_size 1000 \
  --max_mismatches 1 \
  --fuzzy_min_k 101 \
  --verbose
```

## KmerSutra-ML workflow

The ML workflow starts from the read-level hit table produced by `kmersutra-screen`.

Extract sequence-level ML features:

```bash
kmersutra-extract-features \
  --hits_tsv screen_out/read_level_species_kmer_hits.tsv.gz \
  --out_tsv screen_out/sequence_ml_features.tsv \
  --verbose
```

The screening command also writes `sequence_ml_features.tsv` automatically.

To train a classifier, the feature table must contain a supervised label column such as:

- `true_species`
- `true_clade`
- `true_label`

Train an open-set prototype classifier:

```bash
kmersutra-train-classifier \
  --features_tsv training_sequence_ml_features.tsv \
  --label_column true_species \
  --out_model_json kmersutra_species_model.json \
  --out_summary_tsv kmersutra_species_model_summary.tsv \
  --distance_quantile 0.95 \
  --verbose
```

Predict labels for new feature records:

```bash
kmersutra-predict \
  --features_tsv new_sequence_ml_features.tsv \
  --model_json kmersutra_species_model.json \
  --out_tsv new_sequence_ml_predictions.tsv \
  --novelty_scale 1.0 \
  --verbose
```

Prediction output includes:

- `prediction`
- `best_label`
- `best_distance`
- `best_threshold`
- `second_label`
- `distance_margin`
- `open_set_status`
- `ml_confidence_score`

The `prediction` is set to `unknown_or_unresolved` when the feature pattern is outside the learned threshold for the closest known class.

## Main outputs

Panel building:

- `species_kmer_panel.tsv.gz`
- `kmer_uniqueness_summary.tsv`
- `kmer_collection_summary.tsv`
- `species_kmer_panel_metadata.json`
- `species_detection_report.html`
- `build_panel.log`

Screening:

- `read_level_species_kmer_hits.tsv.gz`
- `sequence_ml_features.tsv`
- `sample_species_kmer_hits.tsv`
- `sample_species_kmer_evidence.tsv`
- `species_detection_calls.tsv`
- `species_detection_report.html`
- `screen_reads.log`

ML training:

- model JSON
- training summary TSV
- training log

ML prediction:

- prediction TSV
- prediction log

## Rule-based call categories

- `present_high_confidence`
- `present_low_confidence`
- `ambiguous_mixed_signal`
- `not_detected`

Confidence scores are heuristic in this early version. They should be treated as evidence scores until calibrated against spike-in truth data.

## Open-set interpretation

KmerSutra-ML is intended to avoid overclaiming. A read or contig can be classified as known-like when it resembles a trained class, or as `unknown_or_unresolved` when it falls outside the learned novelty threshold.

This supports outputs such as:

- likely target clade, species resolved
- likely target clade, species unresolved
- closest known clade but outside known-species support
- unknown or unresolved pathogen-like sequence

## Important validation note

For serious species-level claims, avoid training and testing on reads from the same genome assemblies. Use genome-level train/validation/test splits wherever possible. This reduces leakage and gives a more honest estimate of performance on novel strains or divergent taxa.

## Testing

Run tests with:

```bash
nose2
```

Run tests with names and docstring summaries:

```bash
nose2 -v
```

The tests use deterministic toy genomes, reads, feature tables, and model inputs so expected outputs are known exactly.

## Parallel execution

KmerSutra supports worker-process parallelism in both the panel-building and screening stages.

Build a panel with multiple workers:

```bash
kmersutra-build-panel \
  --genome_config examples/example_genome_config.tsv \
  --out_dir example_panel_build \
  --k_values 71 101 \
  --target_clade Demo \
  --threads 4 \
  --verbose
```

Screen reads or assemblies with multiple workers:

```bash
kmersutra-screen \
  --input reads.fastq.gz \
  --input_format fastq \
  --panel example_panel_build/species_kmer_panel.tsv.gz \
  --sample_id sample_001 \
  --out_dir sample_001_kmersutra \
  --threads 4 \
  --chunk_size 1000 \
  --max_mismatches 0 \
  --verbose
```

For raw ONT read screening, exact matching should be tested first. Fuzzy matching is currently limited to one or two substitutions and should be restricted to longer k-mers.

## Run-level spike-in summaries

KmerSutra can now create a run-level Excel and HTML summary from a spike-in
summary TSV. This is used by the spike-in shell wrapper after all sample-level
screening jobs have completed.

```bash
kmersutra-summarise-run \
  --summary_tsv spikein_multi_kmersutra_summary.tsv \
  --out_xlsx kmersutra_spikein_summary.xlsx \
  --out_html kmersutra_spikein_summary.html \
  --verbose
```

The Excel workbook includes:

- `Run_Summary`: the raw wide run summary table
- `Species_Long`: one row per replicate, spike level, and species
- `By_Spike`: detection rates and mean evidence by spike level and species
- `Call_Counts`: counts of each call class by species

All sheets use frozen top rows, filter drop-downs, formatted tables, wrapped
text, and sensible column widths.

## Mixed-species calls

By default, KmerSutra now treats multiple species that independently pass the
configured evidence thresholds as `present_in_mixed_sample`. This is intended
for metagenomic and spike-in settings where true mixed-species samples are
possible. To retain the older conservative behaviour, use:

```bash
kmersutra-screen ... --disallow_mixed_species
```

With that option, multiple species passing evidence thresholds are labelled as
`ambiguous_conflicting_signal`.

## Zero-hit samples

For benchmarking, KmerSutra now writes explicit zero-evidence rows for every
species represented in the diagnostic panel. This means zero-spike controls are
reported as `not_detected` rather than producing header-only call tables.
