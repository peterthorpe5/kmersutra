# Installation and cluster deployment

## Fresh Conda environment

From the repository root:

```bash
conda env create --file environment.yml
conda run --name kmersutra kmersutra-screen --help
```

## Reinstall after a Git pull

```bash
REPO="/home/pthorpe001/github_repos/kmersutra"

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

## Database root

The recommended cluster database root is:

```text
/home/pthorpe001/data/databases/kmersutra_db
```

Set it before using the supplied ATCC example configuration:

```bash
export KMERSUTRA_DB_ROOT="/home/pthorpe001/data/databases/kmersutra_db"
export KMERSUTRA_AI_MODEL="/absolute/path/to/final_internal_calibrator_all_training.json"
```
