# v0.51.0 overlay, GitHub and cluster commands

These commands apply the verified rebuild over the existing clone while
preserving the repository's Git history.

## 1. Mac: apply the downloaded overlay and push

```bash
OVERLAY_ZIP="/Users/PThorpe001/Downloads/KmerSutra_v0_51_0_benchmark_ready_repository_20260730.zip"
REPO="/Users/PThorpe001/github_repos/kmersutra"
OVERLAY_DIR="$(mktemp -d /tmp/kmersutra_v0_51_0.XXXXXX)"

cd "${REPO}"

git switch main
git pull --ff-only
git status --short

unzip -q "${OVERLAY_ZIP}" -d "${OVERLAY_DIR}"

OVERLAY_ROOT="${OVERLAY_DIR}/kmersutra_v0_51_0_benchmark_ready_20260730"

rsync -av \
    --exclude ".git/" \
    "${OVERLAY_ROOT}/" \
    "${REPO}/"

bash "${REPO}/scripts/apply_v051_overlay_cleanup.sh" \
    --repo "${REPO}"

git diff --check
git status --short

conda run --name kmersutra \
    python -m pip install \
    --no-deps \
    --force-reinstall \
    --editable "${REPO}"

conda run --name kmersutra \
    bash "${REPO}/run_tests.sh"

git add -A
git status --short

git commit -m "Prepare KmerSutra v0.51.0 locked benchmark release"
git push origin main
```

Review `git status --short` before `git add -A`. The cleanup script stages only
known generated artefacts and root files that have been moved under
`docs/history`.

## 2. Cluster: pull, reinstall and test

Set `REPO` to the actual cluster clone if it differs:

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

conda run --name kmersutra \
    kmersutra-run-benchmark --help

conda run --name kmersutra \
    bash "${REPO}/run_tests.sh"
```

For a fresh environment:

```bash
cd "${REPO}"
conda env create --file environment.yml
```

## 3. Cluster: configure and preflight ATCC

```bash
REPO="/home/pthorpe001/github_repos/kmersutra"
BENCHMARK_ROOT="/home/pthorpe001/data/benchmarks/kmersutra_atcc_msa1003"
CONFIG="${REPO}/benchmarks/atcc_msa1003_hifi/config.cluster.json"

cd "${REPO}"

cp \
    benchmarks/atcc_msa1003_hifi/config.example.json \
    "${CONFIG}"

export KMERSUTRA_DB_ROOT="/home/pthorpe001/data/databases/kmersutra_db"
export KMERSUTRA_AI_MODEL="/absolute/path/to/final_internal_calibrator_all_training.json"

conda run --name kmersutra \
    kmersutra-run-benchmark \
    --config "${CONFIG}" \
    --output-root "${BENCHMARK_ROOT}" \
    --threads 8 \
    --stop-after 00_preflight \
    --resume \
    --verbose
```

Before running the preflight, edit only the site-specific paths in
`config.cluster.json`. The referenced held-out panel and its genome
configuration must already exist. The preflight stops on exact truth-accession
or sequence leakage.

## 4. Cluster: submit the locked Slurm benchmark

```bash
REPO="/home/pthorpe001/github_repos/kmersutra"
BENCHMARK_ROOT="/home/pthorpe001/data/benchmarks/kmersutra_atcc_msa1003"
CONFIG="${REPO}/benchmarks/atcc_msa1003_hifi/config.cluster.json"

cd "${REPO}"

bash benchmarks/atcc_msa1003_hifi/submit_atcc_msa1003_benchmark.sh \
    --config "${CONFIG}" \
    --output-root "${BENCHMARK_ROOT}" \
    --account barton \
    --partition barton \
    --time 72:00:00 \
    --memory 128G \
    --threads 8 \
    --conda-env kmersutra \
    --resume
```

If the job reaches the 72-hour limit, submit the same command again. Validated
completed stages and per-screen task directories are skipped.

To inspect a submission:

```bash
squeue --me
find "${BENCHMARK_ROOT}/submission_logs" -maxdepth 1 -type f -print
find "${BENCHMARK_ROOT}/slurm_runtime" -maxdepth 2 -type f -print
```
