# KmerSutra AI novelty-scale operating curve script

This small package contains a standalone Python script for summarising external
Zymo AI validation novelty-scale sweeps.

## Files

- `scripts/ai_python/plot_kmersutra_ai_novelty_curve.py`
- `tests/test_plot_kmersutra_ai_novelty_curve.py`

The script scans run folders containing `external_zymo_predictions.tsv.gz` and
uses `external_zymo_predictions.tsv.log` or `run_submission_settings.tsv` to find
the novelty scale for each run.

## Main outputs

- `ai_novelty_sweep_summary.tsv`
- `ai_expected_target_pr_roc_operating_points.tsv`
- `kmersutra_ai_novelty_sweep_operating_curve.png/pdf`
- `kmersutra_ai_novelty_sweep_expected_pr_points.png/pdf`
- `kmersutra_ai_novelty_sweep_expected_roc_points.png/pdf`

## Run tests

```bash
python -m unittest discover -s tests -p "test_plot_kmersutra_ai_novelty_curve.py"
```

## Example command

```bash
python scripts/ai_python/plot_kmersutra_ai_novelty_curve.py \
  --sweep_root /path/to/ai_validation \
  --out_dir /path/to/ai_validation/ai_novelty_sweep_summary_v0501 \
  --exclude_run_regex "20260609_141451" \
  --verbose
```
