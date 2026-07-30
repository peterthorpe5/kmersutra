#!/usr/bin/env bash

set -Eeuo pipefail

usage() {
    echo "Usage:"
    echo "  bash submit_zymo_d6300_screen.sh \\"
    echo "      --reads FILE --panel FILE --output-root DIR [options]"
    echo
    echo "Required:"
    echo "  --reads FILE           ERR5396170 FASTQ or FASTQ.GZ."
    echo "  --panel FILE           KmerSutra panel TSV or TSV.GZ."
    echo "  --output-root DIR      Durable screening output directory."
    echo
    echo "Options:"
    echo "  --sample-id NAME       Sample label (default: ERR5396170)."
    echo "  --min-mixed FRACTION   Mixed-species fraction (default: 0.05)."
    echo "  --account NAME         Slurm account (default: barton)."
    echo "  --partition NAME       Slurm partition (default: barton)."
    echo "  --time HH:MM:SS        Time limit (default: 72:00:00)."
    echo "  --memory SIZE          Memory request (default: 96G)."
    echo "  --threads N            CPUs and workers (default: 24)."
    echo "  --conda-env NAME       Conda environment (default: kmersutra)."
    echo "  --dry-run              Print the sbatch command without submitting."
    echo "  --help                 Show this help."
}

READS=""
PANEL=""
OUTPUT_ROOT=""
SAMPLE_ID="ERR5396170"
MIN_MIXED="0.05"
ACCOUNT="barton"
PARTITION="barton"
TIME_LIMIT="72:00:00"
MEMORY="96G"
THREADS="24"
CONDA_ENV="kmersutra"
DRY_RUN="0"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --reads)
            READS="${2:?Missing value for --reads}"
            shift 2
            ;;
        --panel)
            PANEL="${2:?Missing value for --panel}"
            shift 2
            ;;
        --output-root)
            OUTPUT_ROOT="${2:?Missing value for --output-root}"
            shift 2
            ;;
        --sample-id)
            SAMPLE_ID="${2:?Missing value for --sample-id}"
            shift 2
            ;;
        --min-mixed)
            MIN_MIXED="${2:?Missing value for --min-mixed}"
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

if [[ -z "${READS}" || -z "${PANEL}" || -z "${OUTPUT_ROOT}" ]]; then
    echo "ERROR: --reads, --panel and --output-root are required" >&2
    usage >&2
    exit 2
fi
if [[ ! -s "${READS}" ]]; then
    echo "ERROR: reads are missing or empty: ${READS}" >&2
    exit 1
fi
if [[ ! -s "${PANEL}" ]]; then
    echo "ERROR: panel is missing or empty: ${PANEL}" >&2
    exit 1
fi
if [[ ! "${THREADS}" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: --threads must be a positive integer" >&2
    exit 2
fi
if [[ ! "${MIN_MIXED}" =~ ^(0(\.[0-9]+)?|1(\.0+)?)$ ]]; then
    echo "ERROR: --min-mixed must be between 0 and 1" >&2
    exit 2
fi
if ! command -v sbatch >/dev/null 2>&1; then
    echo "ERROR: sbatch is unavailable; submit from a Slurm login node" >&2
    exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
JOB_SCRIPT="${SCRIPT_DIR}/run_zymo_d6300_screen.slurm"
if [[ ! -s "${JOB_SCRIPT}" ]]; then
    echo "ERROR: Slurm job script is missing: ${JOB_SCRIPT}" >&2
    exit 1
fi

mkdir -p "${OUTPUT_ROOT}/submission_logs"
export KS_ZYMO_READS="$(cd "$(dirname "${READS}")" && pwd)/$(basename "${READS}")"
export KS_ZYMO_PANEL="$(cd "$(dirname "${PANEL}")" && pwd)/$(basename "${PANEL}")"
export KS_ZYMO_OUTPUT_ROOT="$(cd "${OUTPUT_ROOT}" && pwd)"
export KS_ZYMO_SAMPLE_ID="${SAMPLE_ID}"
export KS_ZYMO_MIN_MIXED="${MIN_MIXED}"
export KS_ZYMO_THREADS="${THREADS}"
export KS_ZYMO_CONDA_ENV="${CONDA_ENV}"

SBATCH_COMMAND=(
    sbatch
    --parsable
    --account="${ACCOUNT}"
    --partition="${PARTITION}"
    --time="${TIME_LIMIT}"
    --mem="${MEMORY}"
    --cpus-per-task="${THREADS}"
    --job-name="KS_ZYMO_D6300"
    --output="${OUTPUT_ROOT}/submission_logs/%x.%j.out"
    --error="${OUTPUT_ROOT}/submission_logs/%x.%j.err"
    --export=ALL
    "${JOB_SCRIPT}"
)

printf "Submitting Zymo D6300 screen:"
printf " %q" "${SBATCH_COMMAND[@]}"
printf "\n"
if [[ "${DRY_RUN}" == "1" ]]; then
    exit 0
fi

JOB_ID="$("${SBATCH_COMMAND[@]}")"
echo "Submitted Slurm job: ${JOB_ID}"
echo "Monitor with: squeue -j ${JOB_ID}"
