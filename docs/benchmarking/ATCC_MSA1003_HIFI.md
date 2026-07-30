# Locked ATCC MSA-1003 HiFi benchmark

## Scientific purpose

SRR9328980 is a PacBio HiFi shotgun dataset from the 20-organism ATCC MSA-1003
staggered genomic-DNA mock community. It extends validation beyond the existing
Plasmodium and Zymo ONT analyses and tests generalisation across a different
community, sequencing technology and four abundance tiers.

The public dataset catalogue reports 2,419,037 reads, mean read length 8.4 kb,
20.5 Gb total sequence and median quality Q36.

Source records:

- ATCC MSA-1003 product record:
  <https://www.atcc.org/products/msa-1003>
- PacBio public metagenomics dataset catalogue:
  <https://github.com/PacificBiosciences/pb-metagenomics-tools/blob/master/docs/PacBio-Data.md>
- ATCC composition and reference-accession webinar:
  <https://www.atcc.org/-/media/resources/webinar-presentations/2017/introduction-to-atcc-microbiome-standards--an-endtoend-solution-for-the-standardization-of-microbiom.pdf>

## Registered primary design

- Dataset: SRR9328980.
- Exact k-mer matching.
- k-mer ladder: 51, 77, 101 and 151.
- Existing conservative/lineage-aware call framework.
- Same-genus reportability fraction: 0.05.
- Flat-panel full-depth result is primary.
- Exact ATCC truth references are excluded and audited.
- AI novelty scale: 2.9, applied without retraining.
- No threshold sweep may select the primary result.

## Pre-specified secondary analyses

- Deterministic 1%, 5%, 10% and 25% depth subsets.
- Three seeds per depth: 1001, 1002 and 1003.
- Full-depth single-k ablations.
- Hierarchical screening when a module manifest is configured.
- Frozen AI evidence-calibration validation.

Single-k ablations use one positive k value and remove a minimum-best-k
restriction so each individual k can be interrogated. They are secondary
diagnostic analyses, not replacements for the locked multi-k primary result.

## Reference-panel requirement

The exact truth reference sequences in `truth_manifest.tsv` must not appear in
the panel. Same-species alternative assemblies are allowed and expected. Before
screening, the workflow audits:

- exact assembly accession matches;
- truth sequence accessions in every panel FASTA header;
- exact FASTA SHA-256 matches when truth FASTAs/checksums are supplied.

The workflow stops if exact leakage is found.

The reference panel should also contain absent near neighbours and outgroups.
Do not build a 20-strain truth-only panel: that would not provide a meaningful
specificity challenge. The existing expanded query-agnostic panel can be used
only if all 20 expected species are represented by held-out alternatives and it
passes the audit.

## Configure

```bash
REPO="/path/to/kmersutra"
DATABASE_ROOT="/path/to/kmersutra_db"
AI_MODEL="/path/to/final_internal_calibrator_all_training.json"
BENCHMARK_ROOT="/path/to/kmersutra_atcc_msa1003"

cd "${REPO}"
```

The example JSON is kept portable. Do not copy private paths into it and do not
alter the locked scientific settings.

Run the preflight before spending compute:

```bash
conda run --name kmersutra \
    kmersutra-run-benchmark \
    --config benchmarks/atcc_msa1003_hifi/config.example.json \
    --database-root "${DATABASE_ROOT}" \
    --ai-model "${AI_MODEL}" \
    --output-root "${BENCHMARK_ROOT}" \
    --threads 8 \
    --stop-after 00_preflight \
    --resume \
    --verbose
```

## Submit to Slurm

```bash
bash benchmarks/atcc_msa1003_hifi/submit_atcc_msa1003_benchmark.sh \
    --config benchmarks/atcc_msa1003_hifi/config.example.json \
    --database-root "${DATABASE_ROOT}" \
    --ai-model "${AI_MODEL}" \
    --output-root "${BENCHMARK_ROOT}" \
    --account barton \
    --partition barton \
    --threads 8 \
    --memory 128G \
    --time 72:00:00 \
    --conda-env kmersutra \
    --resume
```

If the Slurm time limit is reached, submit the same command again with
`--resume`. Completed stages and completed screen tasks are validated and
skipped.

To rerun one stage deliberately:

```bash
bash benchmarks/atcc_msa1003_hifi/submit_atcc_msa1003_benchmark.sh \
    --config benchmarks/atcc_msa1003_hifi/config.example.json \
    --database-root "${DATABASE_ROOT}" \
    --ai-model "${AI_MODEL}" \
    --output-root "${BENCHMARK_ROOT}" \
    --account barton \
    --partition barton \
    --threads 8 \
    --conda-env kmersutra \
    --resume \
    --start-at 05_screen_single_k \
    --stop-after 05_screen_single_k \
    --force-stage 05_screen_single_k

```

## Primary outputs

```text
stages/00_preflight/reference_panel_leakage_audit.tsv
stages/03_screen_full/SRR9328980/species_detection_calls.tsv
stages/04_screen_depths/depth_screen_tasks.tsv
stages/05_screen_single_k/single_k_screen_tasks.tsv
stages/06_screen_hierarchical/hierarchical_screen_tasks.tsv
stages/07_ai_validation/atcc_ai_metrics.tsv
stages/08_summarise/tables/sample_summary.tsv
stages/08_summarise/tables/abundance_tier_summary.tsv
stages/08_summarise/tables/expected_species_results.tsv.gz
stages/08_summarise/tables/off_target_results.tsv.gz
stages/09_provenance/input_checksums.tsv
```

Raw evidence detection and reportable species detection are kept separate
throughout.
