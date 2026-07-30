#!/usr/bin/env bash
#SBATCH --job-name=KS_ZYMO_D6300
#SBATCH --account=barton
#SBATCH --partition=barton
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=24
#SBATCH --mem=96G

set -Eeuo pipefail

required_environment=(
    KS_ZYMO_READS
    KS_ZYMO_PANEL
    KS_ZYMO_OUTPUT_ROOT
    KS_ZYMO_SAMPLE_ID
    KS_ZYMO_MIN_MIXED
    KS_ZYMO_THREADS
    KS_ZYMO_CONDA_ENV
)
for variable_name in "${required_environment[@]}"; do
    if [[ -z "${!variable_name:-}" ]]; then
        echo "ERROR: required environment variable is blank: ${variable_name}" >&2
        exit 2
    fi
done

if [[ ! -s "${KS_ZYMO_READS}" ]]; then
    echo "ERROR: reads are missing or empty: ${KS_ZYMO_READS}" >&2
    exit 1
fi
if [[ ! -s "${KS_ZYMO_PANEL}" ]]; then
    echo "ERROR: panel is missing or empty: ${KS_ZYMO_PANEL}" >&2
    exit 1
fi

mkdir -p \
    "${KS_ZYMO_OUTPUT_ROOT}/results" \
    "${KS_ZYMO_OUTPUT_ROOT}/runtime"
RUNTIME_ROOT="${KS_ZYMO_OUTPUT_ROOT}/runtime/${SLURM_JOB_ID:-manual}"
mkdir -p "${RUNTIME_ROOT}"

CONDA_EXE_PATH="$(command -v conda || true)"
if [[ -z "${CONDA_EXE_PATH}" ]]; then
    echo "ERROR: conda executable is not available" >&2
    exit 1
fi

RESULT_DIR="${KS_ZYMO_OUTPUT_ROOT}/results/${KS_ZYMO_SAMPLE_ID}"
COMMAND=(
    "${CONDA_EXE_PATH}"
    run
    --no-capture-output
    --name "${KS_ZYMO_CONDA_ENV}"
    kmersutra-screen
    --input "${KS_ZYMO_READS}"
    --input_format fastq
    --panel "${KS_ZYMO_PANEL}"
    --sample_id "${KS_ZYMO_SAMPLE_ID}"
    --out_dir "${RESULT_DIR}"
    --screen_mode flat
    --screen_preset exact
    --call_preset raw_ont_sensitive
    --min_mixed_species_fraction "${KS_ZYMO_MIN_MIXED}"
    --same_genus_reportable_min_fraction 0.05
    --consolidate_species_calls
    --threads "${KS_ZYMO_THREADS}"
    --chunk_size 10000
    --decompressor auto
    --use_panel_cache
    --no_read_level_hits
    --write_parquet_outputs
    --profile
    --verbose
)

{
    printf "setting\tvalue\n"
    printf "started_at_utc\t%s\n" "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    printf "slurm_job_id\t%s\n" "${SLURM_JOB_ID:-}"
    printf "reads\t%s\n" "${KS_ZYMO_READS}"
    printf "panel\t%s\n" "${KS_ZYMO_PANEL}"
    printf "sample_id\t%s\n" "${KS_ZYMO_SAMPLE_ID}"
    printf "min_mixed_species_fraction\t%s\n" "${KS_ZYMO_MIN_MIXED}"
    printf "threads\t%s\n" "${KS_ZYMO_THREADS}"
    printf "conda_env\t%s\n" "${KS_ZYMO_CONDA_ENV}"
    printf "command\t"
    printf "%q " "${COMMAND[@]}"
    printf "\n"
} > "${RUNTIME_ROOT}/run_configuration.tsv"

"${CONDA_EXE_PATH}" run --name "${KS_ZYMO_CONDA_ENV}" kmersutra-screen --help \
    >/dev/null

if [[ -x /usr/bin/time ]]; then
    /usr/bin/time \
        --verbose \
        --output "${RUNTIME_ROOT}/resource_usage.txt" \
        "${COMMAND[@]}"
else
    "${COMMAND[@]}"
fi

if command -v sacct >/dev/null 2>&1 && [[ -n "${SLURM_JOB_ID:-}" ]]; then
    sacct \
        --jobs "${SLURM_JOB_ID}" \
        --format=JobID,JobName,State,Elapsed,AllocCPUS,MaxRSS,MaxVMSize,ExitCode \
        --parsable2 \
        > "${RUNTIME_ROOT}/sacct.tsv" || true
fi

printf "setting\tvalue\nstatus\tcomplete\nfinished_at_utc\t%s\n" \
    "$(date -u +%Y-%m-%dT%H:%M:%SZ)" \
    > "${RUNTIME_ROOT}/completion.tsv"
