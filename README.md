# KmerSutra

KmerSutra is an outgroup-aware k-mer framework for clade- and species-resolved metagenomic detection. The first intended application is species-resolved *Plasmodium* detection in Oxford Nanopore Technologies metagenomic reads and assemblies, but the code is designed for any clade of interest.


## Requirements

KmerSutra has lightweight Python dependencies and is designed to run in a standard scientific Python environment.

### Core requirements

- Python >= 3.10  
- pandas  
- numpy  
- jinja2  
- biopython  
- openpyxl  
- matplotlib  
- nose2  

### Installation (conda recommended)

```bash
conda create -n kmersutra python=3.10
conda activate kmersutra

pip install pandas numpy jinja2 biopython openpyxl matplotlib nose2
```


### Install KmerSutra

From the repository root:

```bash
pip install -e .
```

### Run tests

```bash
nose2
```


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
