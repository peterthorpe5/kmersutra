# KmerSutra

KmerSutra is an outgroup-aware k-mer framework for clade- and species-resolved metagenomic detection. The first intended application is species-resolved *Plasmodium* detection in Oxford Nanopore Technologies metagenomic reads and assemblies, but the code is designed for any clade of interest.

## Core idea

KmerSutra builds diagnostic k-mer panels from target genomes and outgroup genomes. It then screens reads or assemblies for evidence that supports:

- clade-level detection
- species-level resolution
- unresolved clade-level signal
- conflicting species evidence

The first implementation uses exact canonical k-mer matching. Optional fuzzy matching by Hamming distance is available for long k-mers, but should be benchmarked carefully before being used for final biological interpretation.

## Typical workflow

```bash
python scripts/build_clade_kmer_panel.py \
  --genome_config examples/example_genome_config.tsv \
  --out_dir panel_out \
  --k_values 31 51 71 101 \
  --target_clade Plasmodium \
  --verbose

python scripts/screen_reads_for_clade_kmers.py \
  --input sample.fastq.gz \
  --input_format fastq \
  --panel panel_out/species_kmer_panel.tsv.gz \
  --sample_id sample_001 \
  --out_dir screen_out \
  --max_mismatches 0 \
  --verbose
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

## Main outputs

Panel building:

- `species_kmer_panel.tsv.gz`
- `kmer_uniqueness_summary.tsv`
- `species_kmer_panel_metadata.json`
- `species_detection_report.html`
- `build_panel.log`

Screening:

- `read_level_species_kmer_hits.tsv.gz`
- `sample_species_kmer_hits.tsv`
- `sample_species_kmer_evidence.tsv`
- `species_detection_calls.tsv`
- `species_detection_report.html`
- `screen_reads.log`

## Call categories

- `present_high_confidence`
- `present_low_confidence`
- `ambiguous_mixed_signal`
- `not_detected`

Confidence scores are heuristic in this early version. They should be treated as evidence scores until calibrated against spike-in truth data.

## Testing

Run tests with:

```bash
nose2
```

The tests use deterministic toy genomes and reads so that expected outputs are known exactly.

## Parallel execution

KmerSutra supports worker-process parallelism in both the panel-building and screening stages.

Build a panel with multiple workers:

```bash
python scripts/build_clade_kmer_panel.py \
  --genome_config examples/example_genome_config.tsv \
  --out_dir example_panel_build \
  --k_values 71 101 \
  --target_clade Demo \
  --threads 4 \
  --verbose
```

Screen reads or assemblies with multiple workers:

```bash
python scripts/screen_reads_for_clade_kmers.py \
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

For raw ONT read screening, exact matching should be tested first. Fuzzy matching is currently limited to one or two substitutions and should be restricted to longer k-mers, for example:

```bash
--max_mismatches 1 --fuzzy_min_k 101
```

## Development tests

KmerSutra tests are written with the standard `unittest` framework and are compatible with `nose2`.

```bash
pip install -e '.[dev]'
nose2 -v
```
