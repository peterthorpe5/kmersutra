# KmerSutra AI validation scripts, 2026-06-09

This package contains publication-standard wrapper scripts and Python modules for KmerSutra AI/evidence-calibration validation. The shell scripts are deliberately thin: they set up paths, copy inputs to `TMPDIR`, run unit tests, call Python scripts, and copy outputs back. The Python code lives in `scripts/ai_python/` and is unit-tested with `unittest`.

The package replaces the earlier inline-Python shell draft. It is intended for reproducible validation and handover, not for hiding analysis inside a shell heredoc.

## Folder layout

```text
scripts/
  qsub_kmersutra_ai_full_internal_validation_tmpdir.sh
  qsub_kmersutra_ai_external_zymo_validation_tmpdir.sh
  run_kmersutra_ai_validation_tests.sh
  ai_python/
    kmersutra_ai_common.py
    validate_kmersutra_ai_holdouts.py
    apply_internal_model_to_zymo.py
    tests/
      test_kmersutra_ai_common.py
      test_validate_kmersutra_ai_holdouts.py
      test_apply_internal_model_to_zymo.py
HANDOVER_SUMMARY_KmerSutra_AI_validation_20260609.txt
README.md
```

## Install into the project folder

From the benchmark root:

```bash
cd /home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/ONT_ZymoBIOMICS_ENAERR5396170
unzip /path/to/kmersutra_ai_validation_package_20260609.zip
```

This will add or overwrite files under `scripts/` and place the Python code under `scripts/ai_python/`.

## Run unit tests manually

```bash
cd /home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/ONT_ZymoBIOMICS_ENAERR5396170
bash scripts/run_kmersutra_ai_validation_tests.sh
```

The qsub scripts also run the unit tests by default. Set `RUN_TESTS=false` only if deliberately skipping tests.

## Internal Plasmodium validation

This performs the full internal train/test validation using the already extracted internal Plasmodium table.

```bash
cd /home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/ONT_ZymoBIOMICS_ENAERR5396170

qsub \
  -v CALLS_TABLE="${PWD}/tables/final_kmersutra_detection_calls_long.tsv.gz" \
  scripts/qsub_kmersutra_ai_full_internal_validation_tmpdir.sh
```

Main outputs are written under:

```text
ai_validation/runs_ai_full_internal_plasmodium_validation_<timestamp>/outputs/
```

Key outputs:

```text
ai_call_training.tsv.gz
all_training_label_counts.tsv
feature_leakage_audit.tsv
validation_design.tsv
skipped_splits.tsv
holdout_metrics.tsv
holdout_summary_by_split.tsv
holdout_predictions.tsv.gz
final_internal_plasmodium_calibrator_all_training.json
final_internal_model_training_summary.tsv
unit_tests.log
scripts_snapshot/
```

## External Zymo validation

This applies the final internal Plasmodium-trained model to a public Zymo species-detection table without retraining. By default, it tries to find the latest internal model and the `ref_label_minmixed_0p05` Zymo output.

```bash
cd /home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/ONT_ZymoBIOMICS_ENAERR5396170
qsub scripts/qsub_kmersutra_ai_external_zymo_validation_tmpdir.sh
```

If automatic discovery fails, pass paths explicitly:

```bash
qsub \
  -v MODEL_JSON="/path/to/final_internal_plasmodium_calibrator_all_training.json",ZYMO_OUT_DIR="/path/to/ERR5396170" \
  scripts/qsub_kmersutra_ai_external_zymo_validation_tmpdir.sh
```

or:

```bash
qsub \
  -v MODEL_JSON="/path/to/final_internal_plasmodium_calibrator_all_training.json",CALLS_TABLE="/path/to/species_detection_calls.tsv",REFERENCE_LABEL_MAP="/path/to/reference_label_map.tsv" \
  scripts/qsub_kmersutra_ai_external_zymo_validation_tmpdir.sh
```

Key outputs:

```text
external_zymo_feature_table.tsv.gz
external_zymo_predictions.tsv.gz
external_zymo_metrics.tsv
external_zymo_truth_label_counts.tsv
external_zymo_prediction_counts.tsv
external_zymo_expected_target_predictions.tsv
external_zymo_validation_manifest.json
unit_tests.log
scripts_snapshot/
```

## Design notes

The internal validation uses only explicit safe feature columns:

```text
n_hits
n_unique_kmers
n_positive_sequences
n_k_values_positive
best_k
n_exact_hits
n_fuzzy_hits
conflicting_unique_kmers
conflict_ratio
reportable_conflicting_unique_kmers
reportable_conflict_ratio
mixed_species_support_fraction
confidence_score
signal_confidence_score
has_long_k_support
has_multi_k_support
exact_hit_fraction
```

It deliberately avoids benchmark-only or identity/leakage columns such as `spike_n`, `spike_n_per_genome`, `total_spike_n`, `sample_id`, `species_name`, `benchmark_family`, `panel`, `replicate`, truth labels and report labels.

The external Zymo script is a cross-domain stress test. It does not retrain the model. It uses species/reference-label truth at the sample/species-report level, not per-read truth.

## Runtime expectation

These AI validation runs should be very fast compared with KmerSutra screening because they operate on summary tables, not on raw FASTQ reads or the full k-mer panel.
