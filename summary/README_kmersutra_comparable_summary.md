# KmerSutra comparable benchmark summary tool

This folder contains a robust Python summary tool for the KmerSutra v0.14 comparable spike-in benchmark. It replaces the earlier development shell script that embedded a large Python block.

The tool is designed for output roots created by the KmerSutra comparable benchmark array, for example:

```bash
/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/runs_kmersutra_v014_global_comparable_20260514_111550
```

## Files

```text
summarise_kmersutra_comparable_benchmark.py
run_kmersutra_comparable_summary.sh
tests/test_summarise_kmersutra_comparable_benchmark.py
README_kmersutra_comparable_summary.md
```

## Main features

The Python summariser:

- reads `kmersutra_v014_comparable_manifest.tsv`;
- scans sample outputs under `samples/<benchmark_family>/<sample_id>/`;
- tolerates partial runs while the SGE array is still running;
- marks missing samples as `missing_calls` rather than silently dropping them;
- collates `species_detection_calls.tsv`;
- collates `sample_species_kmer_evidence.tsv` when present;
- reads `screen_timing.tsv` when present;
- infers expected targets from panel labels and optional panel TSV files;
- treats shuffled samples as negative controls even when their folder name contains a spike count;
- calculates tracked-target sensitivity, specificity, precision, F1 and LOD50/LOD95/LOD100;
- calculates real-world interpretability summaries, including clean sensitivity and off-target burden;
- writes tab-separated output tables;
- writes a formatted Excel workbook;
- writes an HTML report;
- writes a detailed log file.

No comma-separated output files are written.

## Output files

The default output directory is:

```bash
${OUT_ROOT}/summary
```

The main outputs are:

```text
sample_status.tsv
kmersutra_detection_calls_long.tsv
kmersutra_evidence_long.tsv
kmersutra_sample_summary.tsv
progress_by_family.tsv
qc_by_family_spike.tsv
tracked_target_performance.tsv
real_world_by_sample.tsv
real_world_summary.tsv
off_target_summary.tsv
runtime_summary.tsv
kmersutra_comparable_summary.xlsx
kmersutra_comparable_summary.html
kmersutra_comparable_summary.log
```

## Recommended command

Copy the script and wrapper into the project root:

```bash
cd /home/pthorpe001/data/2026_plasmodium_kraken_sensitivity

cp /path/to/summarise_kmersutra_comparable_benchmark.py .
cp /path/to/run_kmersutra_comparable_summary.sh .
chmod +x summarise_kmersutra_comparable_benchmark.py
chmod +x run_kmersutra_comparable_summary.sh
```

Run a partial-safe summary while the array is still running:

```bash
OUT_ROOT="/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/runs_kmersutra_v014_global_comparable_20260514_111550" \
./run_kmersutra_comparable_summary.sh
```

Or call the Python script directly:

```bash
python3 summarise_kmersutra_comparable_benchmark.py \
    --out_root "/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/runs_kmersutra_v014_global_comparable_20260514_111550" \
    --panel1_targets "Plasmodium vivax" \
    --panel2_tsv "/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/PT_nanopore_spike_in_pathogen_detection/configs/pathogen_panel_2.tsv" \
    --panel3_tsv "/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/PT_nanopore_spike_in_pathogen_detection/configs/pathogen_panel_3.tsv" \
    --allow_partial \
    --verbose
```

## Strict final run

Once the SGE array has completed, run in strict mode to ensure no samples are missing:

```bash
OUT_ROOT="/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/runs_kmersutra_v014_global_comparable_20260514_111550" \
ALLOW_PARTIAL=false \
STRICT=true \
./run_kmersutra_comparable_summary.sh
```

## Unit tests

From this folder:

```bash
python -m unittest discover -s tests -v
```

The included tests use synthetic KmerSutra-like outputs and check:

- manifest validation;
- panel target loading;
- shuffled-control negative classification;
- expected versus off-target species counting;
- missing-sample handling for partial summaries;
- target-level sensitivity and LOD calculations;
- real-world/off-target summary generation;
- Excel and HTML report creation;
- strict-mode failure when samples are missing.

## Interpretation notes

The tracked-target table answers:

```text
Was each expected species detected in the samples where it was deliberately spiked?
```

The real-world tables answer the more conservative question:

```text
Was the expected species detected without additional non-expected species calls?
```

This distinction is important. KmerSutra may recover the expected target but still report additional candidate taxa. Those additional calls are captured in `off_target_summary.tsv` and in the clean-sensitivity/off-target-rate metrics.

## v0.15/v0.16 manifest compatibility update

This version auto-detects comparable benchmark manifest files in the output root. It prefers:

1. `kmersutra_comparable_manifest.tsv`
2. `kmersutra_v016_conservative_manifest.tsv`
3. `kmersutra_v015_conservative_manifest.tsv`
4. `kmersutra_v014_comparable_manifest.tsv`

It also normalises version-specific manifest column names. In particular, v0.15-style `spike_reads` is treated as the canonical `spike_n` column used by the summary workflow.

The wrapper can now be called either with `OUT_ROOT`:

```bash
OUT_ROOT="/path/to/runs_kmersutra_v015_conservative_comparable_<STAMP>" \
./run_kmersutra_comparable_summary.sh
```

or with the output root as the first positional argument:

```bash
./run_kmersutra_comparable_summary.sh \
  /path/to/runs_kmersutra_v015_conservative_comparable_<STAMP>
```

If auto-detection is ambiguous, provide an explicit manifest:

```bash
OUT_ROOT="/path/to/run" \
MANIFEST_TSV="/path/to/run/kmersutra_v015_conservative_manifest.tsv" \
./run_kmersutra_comparable_summary.sh
```

This version also supports the v0.16 taxonomic evidence filename `sample_taxonomic_kmer_evidence.tsv` as a fallback when `sample_species_kmer_evidence.tsv` is not present.
