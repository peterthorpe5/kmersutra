# Installation and cluster deployment

## Fresh Conda environment

From the repository root:

```bash
conda env create --file environment.yml
conda run --name kmersutra kmersutra-screen --help
```

## Reinstall after a Git pull

```bash
REPO="/path/to/kmersutra"

cd "${REPO}"

git switch main
git pull --ff-only
git status --short

conda run --name kmersutra \
    python -m pip install \
    --no-deps \
    --force-reinstall \
    --editable "${REPO}"

conda run --name kmersutra kmersutra-screen --help
conda run --name kmersutra bash "${REPO}/run_tests.sh"
```

`conda run` is used deliberately so Slurm jobs do not depend on interactive
shell activation.

## Required benchmark tools

Automatic retrieval of SRR9328980 requires `prefetch` and `fasterq-dump` from
SRA Tools. `pigz` is used when available. Supplying a local
`dataset.input_reads` path avoids the SRA download stage.

## Benchmark paths

Keep private cluster paths outside the tracked example configuration. Supply
them as named options:

```bash
kmersutra-run-benchmark \
    --config benchmarks/atcc_msa1003_hifi/config.example.json \
    --database-root /path/to/kmersutra_db \
    --ai-model /path/to/final_internal_calibrator_all_training.json \
    --output-root /path/to/benchmark_outputs \
    --stop-after 00_preflight \
    --resume
```

The environment variables `KMERSUTRA_DB_ROOT` and `KMERSUTRA_AI_MODEL` remain
supported for compatibility. Named options are preferred because the paths are
shown in the recorded command and checked before Slurm submission.
