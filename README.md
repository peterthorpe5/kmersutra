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
kmersutra-download-taxonomy
kmersutra-merge-panels
kmersutra-validate-panel
kmersutra-summarise-run
```


## v0.8 speed and diagnostics options

KmerSutra v0.8 adds practical speed and diagnostic controls for large ONT
metagenomic screens. Exact matching remains the recommended first benchmark.

Useful screening options:

```bash
kmersutra-screen \
  --input sample.fastq.gz \
  --input_format fastq \
  --panel master_kmer_panel.tsv.gz \
  --sample_id sample_001 \
  --out_dir sample_001_kmersutra \
  --threads 12 \
  --chunk_size 5000 \
  --use_panel_cache \
  --profile \
  --max_mismatches 0
```

For very large screens where read-level hit output is not needed:

```bash
kmersutra-screen \
  --input sample.fastq.gz \
  --input_format fastq \
  --panel master_kmer_panel.tsv.gz \
  --sample_id sample_001 \
  --out_dir sample_001_kmersutra \
  --threads 12 \
  --chunk_size 10000 \
  --use_panel_cache \
  --no_read_level_hits \
  --profile
```

These options write `profile_timing.tsv`, which reports wall-clock time for
panel loading, screening, summarisation, and report writing. Panel caching writes
or reuses a pickled index next to the panel by default, avoiding repeated parsing
of large gzip-compressed TSV panels.

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

## Taxonomy-aware evidence levels

KmerSutra can optionally use the NCBI taxonomy dump to assign diagnostic k-mers to the most appropriate evidence level rather than only labelling k-mers as species-specific. This allows a k-mer shared by several species within one genus to be retained as genus-level evidence, while k-mers shared across unrelated groups can be excluded or assigned to broader levels.

Download NCBI taxonomy files:

```bash
kmersutra-download-taxonomy \
  --taxonomy_dir ncbi_taxonomy \
  --verbose
```

The downloader retrieves `taxdmp.zip` from NCBI and extracts the required files:

```text
nodes.dmp
names.dmp
merged.dmp
delnodes.dmp
```

Build a taxonomy-aware panel:

```bash
kmersutra-build-panel \
  --genome_config examples/example_genome_config.tsv \
  --out_dir example_taxonomy_panel \
  --k_values 51 71 101 151 \
  --taxonomy_dir ncbi_taxonomy \
  --download_taxonomy_if_missing \
  --target_taxid 5820 \
  --evidence_ranks species genus family order class phylum superkingdom \
  --threads 12 \
  --verbose
```

The output panel now includes additional taxonomy-aware fields:

```text
evidence_taxid
evidence_name
evidence_rank
lineage_taxids
source_taxids
```

This supports the central KmerSutra question:

```text
What level of taxonomic evidence is supported across the k-mer ladder?
```

For example, a read may be reported as species-level, genus-level, broader clade-level, or unresolved depending on which evidence tiers are supported.

## Merging independently built panels

Large all-pathogen databases can be expensive to build in one pass. KmerSutra
therefore supports a modular workflow where separate panels are built for
host/background genomes, target clades, viruses, bacteria, fungi, or other
pathogen groups, then merged into a globally validated master panel.

Merge panels without taxonomy-aware reassignment:

```bash
kmersutra-merge-panels \
  --panels module_human/species_kmer_panel.tsv.gz module_plasmodium/species_kmer_panel.tsv.gz \
  --out_dir master_kmersutra_panel \
  --verbose
```

Merge panels with NCBI taxonomy-aware evidence reassignment:

```bash
kmersutra-merge-panels \
  --panels module_human/species_kmer_panel.tsv.gz module_plasmodium/species_kmer_panel.tsv.gz module_viral/species_kmer_panel.tsv.gz \
  --out_dir master_kmersutra_panel \
  --taxonomy_dir ncbi_taxonomy \
  --download_taxonomy_if_missing \
  --evidence_ranks species genus family order class phylum superkingdom \
  --verbose
```

The merge step groups identical k-mers across all input panels and assigns the
most specific globally valid evidence level. A k-mer that was species-specific
inside one module may be downgraded to genus-, family-, or broader-level
evidence if it is also found in related taxa. K-mers that cannot be assigned to
a useful evidence rank are written to the removed-conflict table rather than the
master diagnostic panel.

Merge outputs include:

```text
master_kmer_panel.tsv.gz
master_panel_metadata.json
master_validation_summary.tsv
taxonomic_level_summary.tsv
downgraded_kmers.tsv.gz
removed_conflicting_kmers.tsv.gz
```

Validate any panel:

```bash
kmersutra-validate-panel \
  --panel master_kmersutra_panel/master_kmer_panel.tsv.gz \
  --out_dir master_kmersutra_panel/validation \
  --verbose
```

Validation outputs include:

```text
panel_validation_summary.tsv
panel_validation_issues.tsv
taxonomic_level_summary.tsv
```

This module-level design allows KmerSutra to scale to multiple taxonomic spaces
while still checking global specificity before any k-mer is used for detection.

## KmerSutra v0.9.0: scalable panel building

Version 0.9.0 adds a more scalable build path for larger taxonomic databases.
The original builder stored every k-mer occurrence in memory before deciding
whether each k-mer was species-specific, genus-level, clade-level, or
non-diagnostic. That is simple but scales poorly as more pathogen, host, and
outgroup genomes are added.

The default build command now uses compact set-based grouping. This stores one
record per distinct `(k, k-mer)` key and tracks only the metadata needed for
classification, such as source species, source genomes, source clades, and
source taxids. This should reduce memory use substantially for larger database
builds while preserving the same diagnostic output as the legacy builder for
small panels.

### Compact build mode

Compact build mode is enabled by default:

```bash
kmersutra-build-panel \
  --genome_config pathogen_genomes.tsv \
  --out_dir pathogen_panel \
  --k_values 51 71 101 151 \
  --threads 24 \
  --profile \
  --verbose
```

To force the older occurrence-level builder for debugging or regression checks:

```bash
kmersutra-build-panel \
  --genome_config pathogen_genomes.tsv \
  --out_dir pathogen_panel_legacy \
  --k_values 51 71 101 151 \
  --legacy_observation_build \
  --threads 24 \
  --profile \
  --verbose
```

### Build profiling

The `--profile` option writes:

```text
build_profile_timing.tsv
```

This file reports wall-clock time for key build stages, including config
loading, panel construction, panel writing, and summary writing.

### Recommended scalable database strategy

For large pathogen databases, build smaller modules first and then merge and
validate them globally:

```bash
kmersutra-build-panel --genome_config host_primates.tsv --out_dir module_host --k_values 51 71 101 151 --threads 24 --profile
kmersutra-build-panel --genome_config plasmodium.tsv --out_dir module_plasmodium --k_values 51 71 101 151 --threads 24 --profile
kmersutra-build-panel --genome_config viruses.tsv --out_dir module_viruses --k_values 51 71 101 151 --threads 24 --profile

kmersutra-merge-panels \
  --panels module_host/species_kmer_panel.tsv.gz module_plasmodium/species_kmer_panel.tsv.gz module_viruses/species_kmer_panel.tsv.gz \
  --out_dir master_kmersutra_panel \
  --taxonomy_dir ncbi_taxonomy \
  --verbose

kmersutra-validate-panel \
  --panel master_kmersutra_panel/master_kmer_panel.tsv.gz \
  --out_dir master_kmersutra_panel/validation \
  --verbose
```

This keeps module builds manageable while still checking the final master panel
for cross-module conflicts.

## KmerSutra v0.10.0 and v0.11.0 additions

Version 0.10.0 bundled two workflow tools that were previously being used as separate scripts. Version 0.11.0 extends this with broader biological role support and optional assembly quality filters for cleaner database construction.

### Download NCBI genomes into a KmerSutra-ready layout

Use `kmersutra-download-genomes` to download assemblies under one or more NCBI taxids and write organised genome folders plus metadata tables.

```bash
kmersutra-download-genomes \
  --taxid_plan example_taxid_plan.tsv \
  --out_dir ncbi_genomes_for_kmersutra \
  --email your.email@example.ac.uk \
  --source prefer_refseq \
  --formats genomic_fna \
  --decompress \
  --verbose
```

The downloader writes tab-separated metadata and a ready-to-use KmerSutra genome config:

```text
genomes/
ncbi_download_metadata.tsv
kmersutra_genome_config.tsv
query_summary.tsv
run_config.json
logs/download.log
```

A taxid plan should contain columns like:

```text
taxid	role	clade	group_label	max_assemblies	best_per_species	min_total_length	max_total_length	min_scaffold_n50	min_contig_n50
5820	near_neighbour	Plasmodium	Plasmodium		1				
5811	apicomplexan_outgroup	Apicomplexa	Toxoplasma	2					
5807	apicomplexan_outgroup	Apicomplexa	Cryptosporidium	2					
5911	distant_outgroup	Ciliophora	Tetrahymena	2		10000000		500000	100000
```

The resulting `kmersutra_genome_config.tsv` can be passed directly to `kmersutra-build-panel`.

Version 0.11.0 accepts biologically descriptive roles such as:

```text
target_species
target_clade_member
near_neighbour
apicomplexan_outgroup
distant_outgroup
host_or_background
background_pathogen
outgroup
exclude
```

Only target roles are treated as targets during panel construction. Specialist
outgroup and background roles are retained as metadata but treated as non-target
records for uniqueness and conflict checking.

The downloader also supports optional quality filters in the taxid plan or on the
command line:

```text
min_total_length
max_total_length
min_scaffold_n50
min_contig_n50
```

These are useful when a taxid has very small partial assemblies that should not
be used as outgroup representatives.

### Summarise one or more spike-in benchmark run folders

Use `kmersutra-summarise-spikeins` to combine one or more KmerSutra spike-in run folders into TSV, Excel and HTML outputs.

```bash
kmersutra-summarise-spikeins \
  --input_dirs runs \
  --out_dir kmersutra_spikein_summary \
  --run_glob 'spikein_multi_kmersutra*' \
  --expected_replicates 12 \
  --expected_spike_levels '0 1 5 10 25 50 100 250 500 1000 2500 5000' \
  --verbose
```

Outputs include:

```text
combined_run_summary.tsv
species_long_from_wide_summary.tsv
combined_detection_calls.tsv
combined_hit_summary.tsv
authoritative_species_summary.tsv
by_spike_species_summary.tsv
call_counts.tsv
run_qc.tsv
kmersutra_spikein_overall_summary.xlsx
kmersutra_spikein_overall_summary.html
```

## Version 0.12.0: low-memory target-evidence building and RAM logging

Version 0.12.0 adds a new build path for larger near-neighbour and outgroup
panels. The original compact builder is still available and remains useful for
small to moderate panels, but it still holds a global dictionary of every
distinct `(k, k-mer)` key in memory. For long k values and many genomes, most
k-mers are unique, so that dictionary can become too large for normal cluster
jobs.

The new `--target_evidence_only` mode is designed for the immediate KmerSutra
benchmarking question:

> Do the named target species still retain species-level evidence after adding
> near-neighbour and outgroup genomes?

This mode stores candidate k-mers from genomes labelled `target_species` in an
SQLite database and then streams all other genomes against those candidates.
Non-target genomes are therefore used as filters or downgrade evidence, but the
builder does not hold all non-target k-mers in memory and does not create a full
species-level panel for every near-neighbour species.

### Why this is lower memory

The previous compact path keeps a Python dictionary for all observed k-mers from
all genomes. The new target-evidence path keeps only target-candidate k-mers on
disk and updates those records when the same k-mer is seen in near-neighbours or
outgroups. This trades some runtime and disk use for much lower peak RAM.

### What this mode reports

For each target-candidate k-mer, KmerSutra checks whether it is:

- unique to one target species, giving species-level evidence;
- shared with other taxa inside the retained target taxid, giving genus or
  higher taxonomic evidence where supported by the taxonomy database;
- also present outside the retained target taxid, in which case it is not
  retained for the target panel.

This is intentionally conservative. It is not a replacement for a future full
master-panel builder that globally validates every species in a large database.
It is a pragmatic build mode for the current Plasmodium/outgroup benchmark and
for larger target-centred diagnostic panels.

### RAM monitoring

`kmersutra-build-panel` now writes a RAM log by default. The file is:

```bash
ram_usage.tsv
```

inside the output directory, unless `--ram_log_path` is supplied. The log
contains elapsed time, current RSS and peak RSS in bytes and MB. This is useful
on clusters where `/usr/bin/time -v` is not installed or not reliable.

To change the sampling interval:

```bash
--ram_log_interval_seconds 30
```

To disable RAM logging:

```bash
--ram_log_interval_seconds 0
```

### Recommended command for the 33-genome Plasmodium/outgroup panel

For the current Plasmodium benchmark, first build one k value at a time:

```bash
kmersutra-build-panel \
  --genome_config ncbi_genomes_plasmodium_outgroups_v3/kmersutra_genome_config_targets.tsv \
  --out_dir kmersutra_plasmodium_outgroups_v3_target_evidence_k101 \
  --k_values 101 \
  --taxonomy_dir ncbi_taxonomy \
  --download_taxonomy_if_missing \
  --target_taxid 5820 \
  --evidence_ranks species genus family order class phylum superkingdom \
  --target_evidence_only \
  --sqlite_batch_size 50000 \
  --max_per_species_per_k 50000 \
  --profile \
  --verbose
```

Then repeat with `--k_values 71`, or run a combined test only after confirming
that RAM and disk use are acceptable.

`--max_per_species_per_k` is optional. For exploratory builds it is strongly
recommended because a full unthinned long-k panel can contain many millions of
species-level records. Remove or increase the limit only when you are ready for
the corresponding output size.

### New output files

When `--target_evidence_only` is used, the output directory contains the normal
panel files plus:

- `target_evidence_candidates.sqlite`: SQLite candidate database.
- `target_evidence_build_summary.tsv`: target-candidate and overlap counts.
- `ram_usage.tsv`: RAM usage over time.

The standard files are still written:

- `species_kmer_panel.tsv.gz`
- `kmer_uniqueness_summary.tsv`
- `kmer_collection_summary.tsv`
- `species_kmer_panel_metadata.json`
- `build_profile_timing.tsv` when `--profile` is used
- `species_detection_report.html`

### Tests added in v0.12.0

Version 0.12.0 adds tests for:

- SQLite target-candidate database creation.
- Non-target overlap marking.
- species-level and genus-level evidence from the target-evidence path.
- diagnostic stream thinning.
- CLI output generation in `--target_evidence_only` mode.
- RAM helper functions and RAM-log writing.

The full test suite currently contains 139 tests and passed with:

```bash
python -m unittest discover -s tests -v
```

The tests are also compatible with `nose2` when it is installed.

## Version 0.13.0: query-agnostic all-candidate evidence builds

Version 0.13.0 adds a second SQLite-backed low-memory build mode for real-world
unknown-sample screening. Version 0.12.0 introduced `--target_evidence_only`,
which is useful when the benchmark truth is known and the question is whether
named targets retain species-level evidence after near-neighbour and outgroup
filtering. That mode is not sufficient for real samples, because the true target
will not usually be known in advance.

The new mode is:

```bash
kmersutra-build-panel \
  --genome_config kmersutra_genome_config.tsv \
  --out_dir kmersutra_all_candidate_panel \
  --k_values 71 101 \
  --taxonomy_dir ncbi_taxonomy \
  --download_taxonomy_if_missing \
  --evidence_ranks species genus family order class phylum superkingdom \
  --all_candidate_evidence \
  --max_per_species_per_k 50000 \
  --sqlite_batch_size 50000 \
  --ram_log_interval_seconds 30 \
  --profile \
  --verbose
```

This mode is intentionally query-agnostic. It iterates over every eligible
candidate species in the genome configuration, treats that species as the
temporary candidate, streams all other genomes as filters, and retains validated
evidence for the most specific supported taxonomic level. The final panel can
therefore contain markers for many Plasmodium species, selected apicomplexan
outgroups, distant outgroups, and shared genus- or clade-level evidence. This is
the appropriate database design for asking what is present in an unknown sample.

By default, `--all_candidate_evidence` makes all non-host, non-background and
non-excluded genomes reportable candidates. This means roles such as
`target_species`, `near_neighbour`, `outgroup`, `apicomplexan_outgroup`,
`distant_outgroup`, `non_target` and `downloaded` can all contribute reportable
evidence. If a narrower module is required, use `--candidate_roles` to whitelist
specific roles. For example:

```bash
--candidate_roles target_species near_neighbour apicomplexan_outgroup distant_outgroup
```

Do not use `--target_taxid` for a fully query-agnostic broad panel, because it
restricts retained evidence to that subtree. Use `--target_taxid 5820` only when
building a Plasmodium-only module and deliberately excluding reportable outgroup
evidence.

New v0.13.0 outputs include:

- `all_candidate_evidence.sqlite`: retained diagnostic evidence for all selected
  candidate taxa.
- `species_kmer_panel.tsv.gz`: merged panel used by `kmersutra-screen`.
- `target_evidence_build_summary.tsv`: coarse all-candidate build summary. The
  filename is retained for compatibility with v0.12.0 scripts.
- `kmer_collection_summary.tsv`: per-candidate and per-filter collection records.
- `ram_usage.tsv`: RAM usage sampled during the build.
- `build_profile_timing.tsv`: wall-clock timing for major build stages.

The all-candidate build is expected to be slower than `--target_evidence_only`,
because each candidate species is validated against all other genomes. The gain
is that the resulting database is appropriate for unknown-sample screening and
can generate off-target calls during benchmarking. This is the architectural
step needed before comparing KmerSutra fairly with general-purpose classifiers
such as Kraken2 and Metabuli.

Suggested first Plasmodium/outgroup build:

```bash
cd /home/pthorpe001/data/databases/kmersutra_db

kmersutra-build-panel \
  --genome_config ncbi_genomes_plasmodium_outgroups_v3/kmersutra_genome_config_targets_supported_roles.tsv \
  --out_dir kmersutra_builds/kmersutra_plasmodium_outgroups_v3_all_candidate_k71 \
  --k_values 71 \
  --taxonomy_dir ncbi_taxonomy \
  --download_taxonomy_if_missing \
  --evidence_ranks species genus family order class phylum superkingdom \
  --all_candidate_evidence \
  --max_per_species_per_k 50000 \
  --sqlite_batch_size 50000 \
  --ram_log_interval_seconds 30 \
  --profile \
  --verbose
```

Run k=71 and k=101 as separate jobs first. If runtime and disk use are acceptable,
they can later be merged at the final `species_kmer_panel.tsv.gz` level.

Development safeguards in v0.13.0:

- all previous tests are retained;
- new all-candidate unit tests exercise candidate selection, multi-species
  evidence retention, global evidence caps, and CLI output creation;
- the full test suite contains 143 tests and passes with `python -m unittest`;
- all new functions include PEP 8-style docstrings;
- build logging and RAM tracking are retained.

## Version 0.14.0: scalable global all-candidate evidence builds

Version 0.14.0 adds a new global, query-agnostic evidence builder designed to
replace the repeated all-versus-all candidate loop introduced in v0.13.0. The
v0.13 all-candidate mode was biologically correct, because each candidate
species was validated against all other genomes, but it rescanned the same
reference genomes for every candidate species. That is not scalable for larger
panels.

The new mode is:

```bash
--global_candidate_evidence
```

This build path changes the algorithm:

1. Each genome is indexed once for each requested k value.
2. Source metadata for each distinct `(k, kmer)` key is stored in SQLite.
3. The complete set of source taxids for each k-mer is used to assign the most
   specific supported taxonomic evidence level.
4. Evidence is retained for all reportable candidate taxa, including target
   species, near-neighbour species and outgroup species, unless role filters are
   supplied.
5. Optional caps such as `--max_per_species_per_k` are applied after evidence
   assignment.

This is the preferred build mode for unknown-sample panels, where the true
species is not known in advance. It is designed to answer questions such as:

- Which known species has species-level support?
- Is the evidence only genus-level or clade-level?
- Is a non-target or outgroup species supported?
- Does a sample contain conflicting or mixed evidence?

Example broad-panel build:

```bash
kmersutra-build-panel \
  --genome_config kmersutra_genome_config.tsv \
  --out_dir kmersutra_global_candidate_panel_k77_101 \
  --k_values 77 101 \
  --taxonomy_dir ncbi_taxonomy \
  --download_taxonomy_if_missing \
  --evidence_ranks species genus family order class phylum superkingdom \
  --threads 24 \
  --global_candidate_evidence \
  --sqlite_batch_size 50000 \
  --max_per_species_per_k 50000 \
  --ram_log_path ram_usage.tsv \
  --ram_log_interval_seconds 30 \
  --profile \
  --verbose
```

Do not pass `--target_taxid` unless you intentionally want to restrict retained
evidence to one taxonomic subtree. For a fully query-agnostic broad panel,
leave `--target_taxid` unset.

The global candidate build writes the same main screening panel as other build
modes:

```text
species_kmer_panel.tsv.gz
```

It also writes the SQLite global evidence database:

```text
global_candidate_evidence.sqlite
```

For very large production builds, this SQLite database may be treated as a
build intermediate. The final `species_kmer_panel.tsv.gz`, summary files,
profile file and RAM log are the key outputs needed for screening and reporting.

Important distinction between v0.13 and v0.14:

| Mode | Behaviour | Intended use |
| ---- | --------- | ------------ |
| `--all_candidate_evidence` | Validates one candidate species at a time by rescanning all other genomes. | Correct but slow; retained for regression and small panels. |
| `--global_candidate_evidence` | Indexes each genome once, then assigns evidence from the global k-mer source index. | Preferred scalable mode for larger query-agnostic panels. |

Development safeguards in v0.14.0:

- full unit-test suite retained;
- new tests for global candidate indexing, evidence retention, evidence caps and
  CLI output;
- PEP 8-style docstrings on new functions;
- verbose logging for one-pass genome indexing and evidence assignment;
- RAM logging remains available through `--ram_log_path` and
  `--ram_log_interval_seconds`.

## Version 0.15.0: conservative species-call interpretation

Version 0.15.0 adds a stricter interpretation layer for long-read pathogen
screening. The aim is not to maximise ultra-low-abundance read-level
sensitivity. Instead, the new options make species-level calls harder to earn
while keeping weaker or broader taxonomic evidence visible for interpretation.

New screening options include:

- `--call_preset legacy|conservative|strict`
- `--min_best_k`
- `--min_exact_hits`
- `--min_confidence_score`
- `--min_unique_kmer_margin`
- `--min_unique_kmer_ratio`
- `--low_evidence_call present_low_confidence|observed_below_threshold`

The default `legacy` preset preserves previous behaviour for backwards
compatibility. For publication-facing Plasmodium screening with the global
candidate panel, start with:

```bash
kmersutra-screen \
  --input sample.fastq.gz \
  --input_format fastq \
  --sample_id sample1 \
  --panel species_kmer_panel.tsv.gz \
  --out_dir sample1_kmersutra_v015 \
  --threads 4 \
  --chunk_size 10000 \
  --call_preset conservative \
  --no_read_level_hits \
  --profile \
  --verbose
```

The conservative preset requires multi-k support, long-k support, exact-hit
support and a minimum confidence score before a species-level call is reported.
Weak evidence is labelled as `observed_below_threshold` rather than
`present_low_confidence`, so downstream summaries can keep it visible without
counting it as a positive species call.

Version 0.15.0 also writes `sample_taxonomic_kmer_evidence.tsv`, which retains
broader genus, family or phylum evidence from the panel. This is intended to
support calls such as Plasmodium-like signal detected but species unresolved.

All existing tests plus new v0.15 conservative-call and taxonomic-evidence tests
passed with Python `unittest` during packaging.

## v0.16.0 unresolved and possible-novel lineage reporting

Version 0.16.0 adds a sample-level lineage interpretation layer for cases where
species-level evidence is present but does not justify a confident species call.
The intent is to retain weak neighbouring-species evidence without over-reporting
it as a true species detection.

`kmersutra-screen` now writes:

```text
sample_lineage_interpretation.tsv
```

This file separates:

- `species_detected`: one conservative species-level call is supported.
- `mixed_species_detected`: more than one conservative species-level call is
  supported.
- `unresolved_taxonomic_signal`: genus-level or broader evidence is strong, but
  no species-level call is justified.
- `possible_novel_or_unsampled_lineage`: genus-level or broader evidence is
  strong and the sample also has weak evidence spread across multiple related
  neighbouring species. This is intended to flag potentially novel or unsampled
  lineages without forcing a known-species diagnosis.
- `weak_unresolved_neighbour_signal`: weak neighbouring-species evidence is
  present, but the broader taxonomic evidence does not yet pass reporting
  thresholds.
- `no_supported_signal`: no species or taxonomic evidence is supported.

The output also reports the best and second-best species, the absolute and
relative support margin, the best supported taxonomic rank/name, and a heuristic
lineage confidence score. These scores are not calibrated probabilities. They
are intended for conservative ranking and reporting.

Important principle: weak neighbouring-species evidence is evidence, not a final
species diagnosis. v0.16.0 keeps that evidence visible so potentially novel or
unsampled Plasmodium-like signals can be identified without inflating false
species calls.


## Version 0.17.0: lineage-aware mixed-species reporting

KmerSutra v0.17.0 adds a lineage-aware mixed-sample interpretation policy.
The goal is to keep neighbouring-species evidence visible while avoiding the
over-reporting of weak biological neighbours as confident species diagnoses.

The new call preset is:

```bash
--call_preset lineage_aware
```

This preset keeps the conservative v0.15 requirements for multi-k, long-k and
independent-read support, but adds a mixed-species co-dominance check. When
several species pass the basic evidence thresholds, only species with enough
support relative to the strongest species are promoted to reportable mixed
species calls. Weaker passing species are retained as:

```text
neighbour_lineage_evidence
```

This call means that the evidence is not discarded. It can still contribute to
genus-level, unresolved, or possible novel/unsampled lineage reporting, but it
is not promoted to a species-level diagnosis. This is intended for situations
where a real sample contains strong Plasmodium-like evidence plus weak evidence
spread across biologically plausible neighbouring species.

Useful options are:

```bash
--call_preset lineage_aware \
--min_mixed_species_fraction 0.25 \
--low_evidence_call observed_below_threshold
```

The main species-call output now includes additional columns:

```text
reportable_conflicting_unique_kmers
reportable_conflict_ratio
mixed_species_support_fraction
signal_confidence_score
```

The sample-level lineage output remains:

```text
sample_lineage_interpretation.tsv
```

This file should be used for reporting unresolved lineage evidence, including:

```text
species_detected
mixed_species_detected
unresolved_taxonomic_signal
possible_novel_or_unsampled_lineage
weak_unresolved_neighbour_signal
no_supported_signal
```

Recommended test command:

```bash
python -m unittest discover -s tests -v
```

## v0.18.0 genome-spread marker selection

KmerSutra v0.18.0 adds an optional genome-aware marker thinning strategy for
panel builds. This addresses the risk that a simple per-species/per-k cap can
retain a dense block of adjacent, highly overlapping k-mers from one early
scaffold or from one arbitrary lexical region of the k-mer index.

The default remains the legacy behaviour:

```bash
--marker_selection first_seen
```

For larger publication panels, the recommended experimental mode is:

```bash
--marker_selection genome_spread \
--genome_bin_size 10000 \
--max_per_genome_bin 10 \
--max_per_species_per_k 100000
```

In `genome_spread` mode, KmerSutra first assigns taxonomic evidence, then
selects a deterministic subset of reportable markers across source
`genome/contig/bin` groups. This means that a species-level call is supported by
markers distributed across more of the reference genome rather than by many
adjacent sliding-window k-mers from a small region.

For `--global_candidate_evidence` and `--all_candidate_evidence`, the package
avoids applying the per-evidence cap during evidence assignment when
`--marker_selection genome_spread` is used. Instead, the cap is applied during
panel writing, after all retained evidence has been assigned. This avoids
throwing away valid markers before the genome-spread selector can choose a
representative subset.

### Parquet and module-aware databases

The current release still writes the main screening panel as TSV/TSV.GZ for
compatibility with existing KmerSutra screening and summary workflows. For very
large multi-module databases, the recommended next storage layer is an optional
Parquet/Arrow representation of intermediate source-index and retained-evidence
tables. A future module-aware build workflow should support:

1. building local taxonomic modules independently;
2. storing module evidence tables in Parquet or SQLite;
3. merging module evidence tables into a master source index;
4. re-validating every retained marker globally;
5. downgrading or removing markers when cross-module evidence shows they are not
   specific at the originally assigned taxonomic level.

This is important because a marker that is species-specific within a local
Plasmodium panel may become genus-level, clade-level or non-diagnostic after
host, vector, bacterial, fungal, viral or broader eukaryotic modules are merged.

## v0.19.0: genome-spread default and optional Parquet module workflow

KmerSutra v0.19.0 makes `genome_spread` the default marker-selection strategy
for `kmersutra-build-panel`. This avoids the positional bias introduced by a
simple first-seen cap, where many retained k-mers can come from adjacent sliding
windows on the first scaffold encountered. The legacy behaviour remains
available with:

```bash
--marker_selection first_seen
```

For publication-oriented builds, the recommended default is now:

```bash
--marker_selection genome_spread \
--genome_bin_size 10000 \
--max_per_genome_bin 10
```

### Optional Parquet module export

Large KmerSutra databases are expected to be built as taxonomic modules. A
module can now export its global source-index tables as Parquet files when
`pyarrow` is installed:

```bash
kmersutra-build-panel \
  --genome_config module_genome_config.tsv \
  --out_dir plasmodium_module_build \
  --k_values 77 101 151 \
  --taxonomy_dir ncbi_taxonomy \
  --global_candidate_evidence \
  --max_per_species_per_k 100000 \
  --marker_selection genome_spread \
  --write_module_parquet \
  --module_parquet_dir plasmodium_module_build/module_parquet \
  --module_name plasmodium_apicomplexa_v4 \
  --verbose
```

This writes:

```text
global_kmers.parquet
retained_kmers.parquet
build_events.parquet
module_metadata.json
```

The most important table is `global_kmers.parquet`, because it stores source
metadata before the final reportable evidence assignment. This allows modules
built on different genome sets to be merged later and globally revalidated.

### Global revalidation across modules

The new command `kmersutra-merge-modules` merges one or more module source
indexes, assigns evidence levels using the combined source metadata, and writes
a final screenable panel:

```bash
kmersutra-merge-modules \
  --module_dirs plasmodium_module_build/module_parquet host_module/module_parquet \
  --out_dir merged_master_panel \
  --taxonomy_dir ncbi_taxonomy \
  --evidence_ranks species genus family order class phylum superkingdom \
  --max_per_species_per_k 100000 \
  --marker_selection genome_spread \
  --genome_bin_size 10000 \
  --max_per_genome_bin 10 \
  --verbose
```

This workflow is designed for the key KmerSutra use case: a k-mer that appears
species-specific in one local module may need to be downgraded to genus, family,
clade or removed after other modules are merged. The final screenable panel is
therefore based on globally revalidated evidence rather than local-only marker
specificity.

Parquet support is optional. If `pyarrow` is not installed, normal TSV/SQLite
builds still work; only the Parquet module export/import commands require the
optional dependency.
