# v0.51.1 overlay, GitHub and cluster commands

These commands apply the v0.51.1 overlay while preserving Git history. Replace
the example paths with local values. No private filesystem path belongs in the
tracked example configuration.

## 1. Mac: apply and push

```bash
OVERLAY_ZIP="/path/to/KmerSutra_v0_51_1_documented_benchmark_ready_repository_20260730.zip"
REPO="/path/to/kmersutra"
OVERLAY_DIR="$(mktemp -d /tmp/kmersutra_v0_51_1.XXXXXX)"

cd "${REPO}"
git switch main
git pull --ff-only
git status --short

unzip -q "${OVERLAY_ZIP}" -d "${OVERLAY_DIR}"
OVERLAY_ROOT="${OVERLAY_DIR}/kmersutra_v0_51_1_documented_benchmark_ready_20260730"

rsync -av --exclude ".git/" "${OVERLAY_ROOT}/" "${REPO}/"

bash "${REPO}/scripts/apply_v0511_overlay_cleanup.sh" --repo "${REPO}"

git diff --check
git status --short

conda run --name kmersutra \
    python -m pip install --no-deps --force-reinstall --editable "${REPO}"

conda run --name kmersutra \
    bash "${REPO}/run_tests.sh" \
    --results-dir "${REPO}/test_results"

git add -A
git status --short
git commit -m "Release KmerSutra v0.51.1 documentation and benchmark repairs"
git push origin main
```

Review `git status --short` before `git add -A`.

## 2. Cluster: pull, reinstall and test

```bash
REPO="/path/to/cluster/kmersutra"
RESULTS_ROOT="/path/to/project/kmersutra_test_results"

cd "${REPO}"
git switch main
git pull --ff-only
git status --short

conda run --name kmersutra \
    python -m pip install --no-deps --force-reinstall --editable "${REPO}"

conda run --name kmersutra kmersutra-screen --help
conda run --name kmersutra kmersutra-run-benchmark --help

conda run --name kmersutra \
    bash "${REPO}/run_quality_checks.sh" \
    --results-dir "${RESULTS_ROOT}"
```

The quality run writes dated logs, coverage HTML and the package inventory
beneath `RESULTS_ROOT`; it creates no root-level test or inventory files.

## 3. ATCC preflight and Slurm submission

```bash
DATABASE_ROOT="/path/to/kmersutra_db"
AI_MODEL="/path/to/final_internal_calibrator_all_training.json"
BENCHMARK_ROOT="/path/to/kmersutra_atcc_msa1003"
CONFIG="${REPO}/benchmarks/atcc_msa1003_hifi/config.example.json"

conda run --name kmersutra \
    kmersutra-run-benchmark \
    --config "${CONFIG}" \
    --database-root "${DATABASE_ROOT}" \
    --ai-model "${AI_MODEL}" \
    --output-root "${BENCHMARK_ROOT}" \
    --threads 8 \
    --stop-after 00_preflight \
    --resume \
    --verbose

bash "${REPO}/benchmarks/atcc_msa1003_hifi/submit_atcc_msa1003_benchmark.sh" \
    --config "${CONFIG}" \
    --database-root "${DATABASE_ROOT}" \
    --ai-model "${AI_MODEL}" \
    --output-root "${BENCHMARK_ROOT}" \
    --account barton \
    --partition barton \
    --time 72:00:00 \
    --memory 128G \
    --threads 8 \
    --conda-env kmersutra \
    --resume
```

The preflight is expected to stop until the held-out ATCC panel and its genome
configuration exist under the selected database root and pass the leakage
audit.
