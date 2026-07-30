#!/usr/bin/env bash

set -Eeuo pipefail

usage() {
    echo "Usage:"
    echo "  bash submit_atcc_msa1003_benchmark.sh \\"
    echo "      --config FILE --output-root DIR [options]"
    echo
    echo "Required:"
    echo "  --config FILE          Locked benchmark JSON configuration."
    echo "  --output-root DIR      Durable benchmark output root."
    echo
    echo "Options:"
    echo "  --run-name NAME        Override config run_name."
    echo "  --account NAME         Slurm account (default: barton)."
    echo "  --partition NAME       Slurm partition (default: barton)."
    echo "  --time HH:MM:SS        Slurm time limit (default: 72:00:00)."
    echo "  --memory SIZE          Slurm memory request (default: 128G)."
    echo "  --threads N            CPUs and KmerSutra workers (default: 8)."
    echo "  --conda-env NAME       Conda environment (default: kmersutra)."
    echo "  --resume               Resume validated completed stages."
    echo "  --start-at STAGE       Start at a named stage."
    echo "  --stop-after STAGE     Stop after a named stage."
    echo "  --force-stage STAGE    Force one stage; may be repeated."
    echo "  --dry-run               Print the sbatch command only."
    echo "  --help                  Show this help."
}

CONFIG=""
OUTPUT_ROOT=""
RUN_NAME=""
ACCOUNT="barton"
PARTITION="barton"
TIME_LIMIT="72:00:00"
MEMORY="128G"
THREADS="8"
CONDA_ENV="kmersutra"
RESUME="0"
START_AT=""
STOP_AFTER=""
FORCE_STAGES=""
DRY_RUN="0"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --config)
            CONFIG="${2:?Missing value for --config}"
            shift 2
            ;;
        --output-root)
            OUTPUT_ROOT="${2:?Missing value for --output-root}"
            shift 2
            ;;
        --run-name)
            RUN_NAME="${2:?Missing value for --run-name}"
            shift 2
            ;;
        --account)
            ACCOUNT="${2:?Missing value for --account}"
            shift 2
            ;;
        --partition)
            PARTITION="${2:?Missing value for --partition}"
            shift 2
            ;;
        --time)
            TIME_LIMIT="${2:?Missing value for --time}"
            shift 2
            ;;
        --memory)
            MEMORY="${2:?Missing value for --memory}"
            shift 2
            ;;
        --threads)
            THREADS="${2:?Missing value for --threads}"
            shift 2
            ;;
        --conda-env)
            CONDA_ENV="${2:?Missing value for --conda-env}"
            shift 2
            ;;
        --resume)
            RESUME="1"
            shift
            ;;
        --start-at)
            START_AT="${2:?Missing value for --start-at}"
            shift 2
            ;;
        --stop-after)
            STOP_AFTER="${2:?Missing value for --stop-after}"
            shift 2
            ;;
        --force-stage)
            FORCE_STAGES="${FORCE_STAGES}${FORCE_STAGES:+,}${2:?Missing value for --force-stage}"
            shift 2
            ;;
        --dry-run)
            DRY_RUN="1"
            shift
            ;;
        --help|-h)
            usage
            exit 0
            ;;
        *)
            echo "ERROR: unknown option: $1" >&2
            usage >&2
            exit 2
            ;;
    esac
done

if [[ -z "${CONFIG}" || -z "${OUTPUT_ROOT}" ]]; then
    echo "ERROR: --config and --output-root are required" >&2
    usage >&2
    exit 2
fi
if [[ ! -s "${CONFIG}" ]]; then
    echo "ERROR: configuration is missing or empty: ${CONFIG}" >&2
    exit 1
fi
if [[ ! "${THREADS}" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: --threads must be a positive integer" >&2
    exit 2
fi
if ! command -v sbatch >/dev/null 2>&1; then
    echo "ERROR: sbatch is not available; run this command on a Slurm login node" >&2
    exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
JOB_SCRIPT="${SCRIPT_DIR}/run_atcc_msa1003_benchmark.slurm"
if [[ ! -s "${JOB_SCRIPT}" ]]; then
    echo "ERROR: Slurm job script is missing: ${JOB_SCRIPT}" >&2
    exit 1
fi

mkdir -p "${OUTPUT_ROOT}/submission_logs"

export KS_CONFIG="$(cd "$(dirname "${CONFIG}")" && pwd)/$(basename "${CONFIG}")"
export KS_OUTPUT_ROOT="$(cd "${OUTPUT_ROOT}" && pwd)"
export KS_RUN_NAME="${RUN_NAME}"
export KS_THREADS="${THREADS}"
export KS_CONDA_ENV="${CONDA_ENV}"
export KS_RESUME="${RESUME}"
export KS_START_AT="${START_AT}"
export KS_STOP_AFTER="${STOP_AFTER}"
export KS_FORCE_STAGES="${FORCE_STAGES}"

SBATCH_COMMAND=(
    sbatch
    --parsable
    --account="${ACCOUNT}"
    --partition="${PARTITION}"
    --time="${TIME_LIMIT}"
    --mem="${MEMORY}"
    --cpus-per-task="${THREADS}"
    --job-name="KS_ATCC1003"
    --output="${OUTPUT_ROOT}/submission_logs/%x.%j.out"
    --error="${OUTPUT_ROOT}/submission_logs/%x.%j.err"
    --export=ALL
    "${JOB_SCRIPT}"
)

echo "Submitting locked ATCC MSA-1003 benchmark"
printf ' %q' "${SBATCH_COMMAND[@]}"
echo

if [[ "${DRY_RUN}" == "1" ]]; then
    exit 0
fi

JOB_ID="$("${SBATCH_COMMAND[@]}")"
echo "Submitted Slurm job: ${JOB_ID}"
echo "Monitor with: squeue -j ${JOB_ID}"
