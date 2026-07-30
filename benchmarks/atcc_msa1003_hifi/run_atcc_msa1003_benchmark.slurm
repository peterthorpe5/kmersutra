#!/usr/bin/env bash
#SBATCH --job-name=KS_ATCC1003
#SBATCH --account=barton
#SBATCH --partition=barton
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G

set -Eeuo pipefail

required_environment=(
    KS_CONFIG
    KS_OUTPUT_ROOT
    KS_DATABASE_ROOT
    KS_AI_MODEL
    KS_THREADS
    KS_CONDA_ENV
)
for variable_name in "${required_environment[@]}"; do
    if [[ -z "${!variable_name:-}" ]]; then
        echo "ERROR: required environment variable is blank: ${variable_name}" >&2
        exit 2
    fi
done

if [[ ! -s "${KS_CONFIG}" ]]; then
    echo "ERROR: configuration is missing or empty: ${KS_CONFIG}" >&2
    exit 1
fi
if [[ ! "${KS_THREADS}" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: KS_THREADS must be a positive integer" >&2
    exit 2
fi

mkdir -p "${KS_OUTPUT_ROOT}/slurm_runtime"
RUNTIME_ROOT="${KS_OUTPUT_ROOT}/slurm_runtime/${SLURM_JOB_ID:-manual}"
mkdir -p "${RUNTIME_ROOT}"
STATUS_TABLE="${RUNTIME_ROOT}/job_exit_status.tsv"

write_status() {
    local status=$?
    {
        printf "setting\tvalue\n"
        printf "exit_status\t%s\n" "${status}"
        printf "finished_at_utc\t%s\n" "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
        printf "hostname\t%s\n" "$(hostname)"
        printf "slurm_job_id\t%s\n" "${SLURM_JOB_ID:-}"
        printf "slurm_cpus_per_task\t%s\n" "${SLURM_CPUS_PER_TASK:-}"
        printf "config\t%s\n" "${KS_CONFIG}"
        printf "output_root\t%s\n" "${KS_OUTPUT_ROOT}"
        printf "database_root\t%s\n" "${KS_DATABASE_ROOT}"
        printf "ai_model\t%s\n" "${KS_AI_MODEL}"
        printf "run_name\t%s\n" "${KS_RUN_NAME:-}"
        printf "conda_env\t%s\n" "${KS_CONDA_ENV}"
    } > "${STATUS_TABLE}"
    if command -v sacct >/dev/null 2>&1 && [[ -n "${SLURM_JOB_ID:-}" ]]; then
        sacct \
            --jobs "${SLURM_JOB_ID}" \
            --format=JobID,JobName,State,Elapsed,AllocCPUS,MaxRSS,MaxVMSize,ExitCode \
            --parsable2 \
            > "${RUNTIME_ROOT}/sacct.tsv" || true
    fi
    exit "${status}"
}
trap write_status EXIT
trap 'echo "ERROR at line ${LINENO}: ${BASH_COMMAND}" >&2' ERR

CONDA_EXE_PATH="$(command -v conda || true)"
if [[ -z "${CONDA_EXE_PATH}" ]]; then
    for candidate in \
        "${HOME}/miniforge3/bin/conda" \
        "${HOME}/mambaforge/bin/conda" \
        "${HOME}/miniconda3/bin/conda" \
        "${HOME}/anaconda3/bin/conda"
    do
        if [[ -x "${candidate}" ]]; then
            CONDA_EXE_PATH="${candidate}"
            break
        fi
    done
fi
if [[ -z "${CONDA_EXE_PATH}" || ! -x "${CONDA_EXE_PATH}" ]]; then
    echo "ERROR: conda executable could not be resolved" >&2
    exit 1
fi

RUN_COMMAND=(
    "${CONDA_EXE_PATH}"
    run
    --no-capture-output
    --name "${KS_CONDA_ENV}"
    kmersutra-run-benchmark
    --config "${KS_CONFIG}"
    --output-root "${KS_OUTPUT_ROOT}"
    --database-root "${KS_DATABASE_ROOT}"
    --ai-model "${KS_AI_MODEL}"
    --threads "${KS_THREADS}"
    --verbose
)

if [[ -n "${KS_RUN_NAME:-}" ]]; then
    RUN_COMMAND+=(--run-name "${KS_RUN_NAME}")
fi
if [[ "${KS_RESUME:-0}" == "1" ]]; then
    RUN_COMMAND+=(--resume)
fi
if [[ -n "${KS_START_AT:-}" ]]; then
    RUN_COMMAND+=(--start-at "${KS_START_AT}")
fi
if [[ -n "${KS_STOP_AFTER:-}" ]]; then
    RUN_COMMAND+=(--stop-after "${KS_STOP_AFTER}")
fi
if [[ -n "${KS_FORCE_STAGES:-}" ]]; then
    IFS=',' read -r -a forced_stages <<< "${KS_FORCE_STAGES}"
    for forced_stage in "${forced_stages[@]}"; do
        RUN_COMMAND+=(--force-stage "${forced_stage}")
    done
fi

{
    printf "setting\tvalue\n"
    printf "started_at_utc\t%s\n" "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    printf "hostname\t%s\n" "$(hostname)"
    printf "slurm_job_id\t%s\n" "${SLURM_JOB_ID:-}"
    printf "slurm_cpus_per_task\t%s\n" "${SLURM_CPUS_PER_TASK:-}"
    printf "conda_executable\t%s\n" "${CONDA_EXE_PATH}"
    printf "conda_env\t%s\n" "${KS_CONDA_ENV}"
    printf "database_root\t%s\n" "${KS_DATABASE_ROOT}"
    printf "ai_model\t%s\n" "${KS_AI_MODEL}"
} > "${RUNTIME_ROOT}/job_start.tsv"

"${CONDA_EXE_PATH}" run --name "${KS_CONDA_ENV}" python --version
"${CONDA_EXE_PATH}" run --name "${KS_CONDA_ENV}" kmersutra-screen --help >/dev/null

if [[ -x /usr/bin/time ]]; then
    /usr/bin/time \
        --verbose \
        --output "${RUNTIME_ROOT}/resource_usage.txt" \
        "${RUN_COMMAND[@]}"
else
    "${RUN_COMMAND[@]}"
fi
