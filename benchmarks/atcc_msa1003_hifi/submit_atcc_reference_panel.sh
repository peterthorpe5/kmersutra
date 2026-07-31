#!/usr/bin/env bash

set -Eeuo pipefail

usage() {
    echo "Usage:"
    echo "  bash submit_atcc_reference_panel.sh \\"
    echo "      --repo DIR --database-root DIR --background-config FILE \\"
    echo "      --truth-manifest FILE --taxonomy-dir DIR --email ADDRESS [options]"
    echo
    echo "Required:"
    echo "  --repo DIR                KmerSutra repository root."
    echo "  --database-root DIR       KmerSutra database root."
    echo "  --background-config FILE  Established v4 genome configuration."
    echo "  --truth-manifest FILE     Locked ATCC truth manifest."
    echo "  --taxonomy-dir DIR        NCBI taxdump directory."
    echo "  --email ADDRESS           Email supplied to NCBI Entrez."
    echo
    echo "Options:"
    echo "  --reference-root DIR      Output panel root (default: DATABASE_ROOT/atcc_msa1003_heldout_v1)."
    echo "  --account NAME            Slurm account (default: barton)."
    echo "  --partition NAME          Slurm partition (default: barton)."
    echo "  --time HH:MM:SS           Slurm time limit (default: 72:00:00)."
    echo "  --memory SIZE             Slurm memory request (default: 128G)."
    echo "  --threads N               CPUs and KmerSutra workers (default: 24)."
    echo "  --conda-env NAME          Conda environment (default: kmersutra)."
    echo "  --target-assemblies N     Alternative references requested per target (default: 5)."
    echo "  --neighbour-assemblies N  Assemblies requested per target genus (default: 25)."
    echo "  --minimum-target-refs N   Required held-out references per target (default: 1)."
    echo "  --max-entrez-records N    Entrez search cap per taxid (default: 500)."
    echo "  --resume                  Reuse validated completed stages."
    echo "  --force-stage STAGE       Force one stage; may be repeated."
    echo "  --dry-run                 Print the sbatch command only."
    echo "  --help                    Show this help."
}

REPO=""
DATABASE_ROOT=""
REFERENCE_ROOT=""
BACKGROUND_CONFIG=""
TRUTH_MANIFEST=""
TAXONOMY_DIR=""
NCBI_EMAIL=""
ACCOUNT="barton"
PARTITION="barton"
TIME_LIMIT="72:00:00"
MEMORY="128G"
THREADS="24"
CONDA_ENV="kmersutra"
TARGET_ASSEMBLIES="5"
NEIGHBOUR_ASSEMBLIES="25"
MINIMUM_TARGET_REFS="1"
MAX_ENTREZ_RECORDS="500"
RESUME="0"
FORCE_STAGES=""
DRY_RUN="0"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --repo)
            REPO="${2:?Missing value for --repo}"
            shift 2
            ;;
        --database-root)
            DATABASE_ROOT="${2:?Missing value for --database-root}"
            shift 2
            ;;
        --reference-root)
            REFERENCE_ROOT="${2:?Missing value for --reference-root}"
            shift 2
            ;;
        --background-config)
            BACKGROUND_CONFIG="${2:?Missing value for --background-config}"
            shift 2
            ;;
        --truth-manifest)
            TRUTH_MANIFEST="${2:?Missing value for --truth-manifest}"
            shift 2
            ;;
        --taxonomy-dir)
            TAXONOMY_DIR="${2:?Missing value for --taxonomy-dir}"
            shift 2
            ;;
        --email)
            NCBI_EMAIL="${2:?Missing value for --email}"
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
        --target-assemblies)
            TARGET_ASSEMBLIES="${2:?Missing value for --target-assemblies}"
            shift 2
            ;;
        --neighbour-assemblies)
            NEIGHBOUR_ASSEMBLIES="${2:?Missing value for --neighbour-assemblies}"
            shift 2
            ;;
        --minimum-target-refs)
            MINIMUM_TARGET_REFS="${2:?Missing value for --minimum-target-refs}"
            shift 2
            ;;
        --max-entrez-records)
            MAX_ENTREZ_RECORDS="${2:?Missing value for --max-entrez-records}"
            shift 2
            ;;
        --resume)
            RESUME="1"
            shift
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

if [[ -z "${REPO}" || -z "${DATABASE_ROOT}" || -z "${BACKGROUND_CONFIG}" \
    || -z "${TRUTH_MANIFEST}" || -z "${TAXONOMY_DIR}" || -z "${NCBI_EMAIL}" ]]; then
    echo "ERROR: all required named options must be supplied" >&2
    usage >&2
    exit 2
fi

for value in "${THREADS}" "${TARGET_ASSEMBLIES}" "${NEIGHBOUR_ASSEMBLIES}" \
    "${MINIMUM_TARGET_REFS}" "${MAX_ENTREZ_RECORDS}"; do
    if [[ ! "${value}" =~ ^[1-9][0-9]*$ ]]; then
        echo "ERROR: numeric options must be positive integers: ${value}" >&2
        exit 2
    fi
done

if [[ ! -d "${REPO}" ]]; then
    echo "ERROR: repository root does not exist: ${REPO}" >&2
    exit 1
fi
if [[ ! -d "${DATABASE_ROOT}" ]]; then
    echo "ERROR: database root does not exist: ${DATABASE_ROOT}" >&2
    exit 1
fi
if [[ ! -s "${BACKGROUND_CONFIG}" ]]; then
    echo "ERROR: background configuration is missing or empty: ${BACKGROUND_CONFIG}" >&2
    exit 1
fi
if [[ ! -s "${TRUTH_MANIFEST}" ]]; then
    echo "ERROR: truth manifest is missing or empty: ${TRUTH_MANIFEST}" >&2
    exit 1
fi
if [[ ! -s "${TAXONOMY_DIR}/nodes.dmp" || ! -s "${TAXONOMY_DIR}/names.dmp" ]]; then
    echo "ERROR: taxonomy directory lacks nodes.dmp or names.dmp: ${TAXONOMY_DIR}" >&2
    exit 1
fi
if ! command -v sbatch >/dev/null 2>&1; then
    echo "ERROR: sbatch is unavailable; submit from a Slurm login node" >&2
    exit 1
fi

REPO="$(cd "${REPO}" && pwd)"
DATABASE_ROOT="$(cd "${DATABASE_ROOT}" && pwd)"
BACKGROUND_CONFIG="$(cd "$(dirname "${BACKGROUND_CONFIG}")" && pwd)/$(basename "${BACKGROUND_CONFIG}")"
TRUTH_MANIFEST="$(cd "$(dirname "${TRUTH_MANIFEST}")" && pwd)/$(basename "${TRUTH_MANIFEST}")"
TAXONOMY_DIR="$(cd "${TAXONOMY_DIR}" && pwd)"
REFERENCE_ROOT="${REFERENCE_ROOT:-${DATABASE_ROOT}/atcc_msa1003_heldout_v1}"
mkdir -p "${REFERENCE_ROOT}/submission_logs"
REFERENCE_ROOT="$(cd "${REFERENCE_ROOT}" && pwd)"

JOB_SCRIPT="${REPO}/benchmarks/atcc_msa1003_hifi/build_atcc_reference_panel.slurm"
if [[ ! -s "${JOB_SCRIPT}" ]]; then
    echo "ERROR: reference-panel Slurm script is missing: ${JOB_SCRIPT}" >&2
    exit 1
fi

export KS_ATCC_REPO="${REPO}"
export KS_ATCC_DATABASE_ROOT="${DATABASE_ROOT}"
export KS_ATCC_REFERENCE_ROOT="${REFERENCE_ROOT}"
export KS_ATCC_BACKGROUND_CONFIG="${BACKGROUND_CONFIG}"
export KS_ATCC_TRUTH_MANIFEST="${TRUTH_MANIFEST}"
export KS_ATCC_TAXONOMY_DIR="${TAXONOMY_DIR}"
export KS_ATCC_NCBI_EMAIL="${NCBI_EMAIL}"
export KS_ATCC_THREADS="${THREADS}"
export KS_ATCC_CONDA_ENV="${CONDA_ENV}"
export KS_ATCC_TARGET_ASSEMBLIES="${TARGET_ASSEMBLIES}"
export KS_ATCC_NEIGHBOUR_ASSEMBLIES="${NEIGHBOUR_ASSEMBLIES}"
export KS_ATCC_MINIMUM_TARGET_REFS="${MINIMUM_TARGET_REFS}"
export KS_ATCC_MAX_ENTREZ_RECORDS="${MAX_ENTREZ_RECORDS}"
export KS_ATCC_RESUME="${RESUME}"
export KS_ATCC_FORCE_STAGES="${FORCE_STAGES}"

SBATCH_COMMAND=(
    sbatch
    --parsable
    --account="${ACCOUNT}"
    --partition="${PARTITION}"
    --time="${TIME_LIMIT}"
    --mem="${MEMORY}"
    --cpus-per-task="${THREADS}"
    --job-name="KS_ATCCref"
    --output="${REFERENCE_ROOT}/submission_logs/%x.%j.out"
    --error="${REFERENCE_ROOT}/submission_logs/%x.%j.err"
    --export=ALL
    "${JOB_SCRIPT}"
)

echo "Submitting leakage-controlled ATCC reference build"
echo "Reference root: ${REFERENCE_ROOT}"
echo "Background configuration: ${BACKGROUND_CONFIG}"
printf ' %q' "${SBATCH_COMMAND[@]}"
echo

if [[ "${DRY_RUN}" == "1" ]]; then
    exit 0
fi

JOB_ID="$("${SBATCH_COMMAND[@]}")"
echo "Submitted Slurm job: ${JOB_ID}"
echo "Monitor with: squeue -j ${JOB_ID}"
