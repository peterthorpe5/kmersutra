#!/bin/bash
#$ -N KSai_ext_zymo
#$ -cwd
#$ -j y
#$ -pe smp 4
#$ -o logs
#$ -e logs

set -Eeuo pipefail
trap 'echo "ERROR at line ${LINENO}: ${BASH_COMMAND}" >&2' ERR

echo "Started at: $(date)"
echo "Host: $(hostname)"
echo "Job ID: ${JOB_ID:-NA}"
echo "NSLOTS: ${NSLOTS:-4}"
echo "TMPDIR: ${TMPDIR:-NA}"

BENCH="${BENCH:-${HOME}/data/2026_plasmodium_kraken_sensitivity/ONT_ZymoBIOMICS_ENAERR5396170}"
MODEL_JSON="${MODEL_JSON:-}"
ZYMO_OUT_DIR="${ZYMO_OUT_DIR:-}"
CALLS_TABLE="${CALLS_TABLE:-}"
REFERENCE_LABEL_MAP="${REFERENCE_LABEL_MAP:-}"
RUN_TESTS="${RUN_TESTS:-true}"

RUN_TAG="ai_external_zymo_validation_$(date +%Y%m%d_%H%M%S)"
FINAL_RUN_DIR="${BENCH}/ai_validation/runs_${RUN_TAG}"
FINAL_OUT_DIR="${FINAL_RUN_DIR}/outputs"

WORK_ROOT="${TMPDIR:-/tmp}/${USER}/kmersutra_${RUN_TAG}"
WORK_INPUT_DIR="${WORK_ROOT}/inputs"
WORK_OUT_DIR="${WORK_ROOT}/outputs"
WORK_SCRIPT_DIR="${WORK_ROOT}/scripts"

PROJECT_SCRIPT_DIR="${BENCH}/scripts"
AI_SCRIPT_DIR="${PROJECT_SCRIPT_DIR}/ai_python"
EXTERNAL_SCRIPT="${AI_SCRIPT_DIR}/apply_internal_model_to_zymo.py"

mkdir -p "${WORK_INPUT_DIR}"
mkdir -p "${WORK_OUT_DIR}"
mkdir -p "${WORK_SCRIPT_DIR}"
mkdir -p "${FINAL_OUT_DIR}"
mkdir -p "${BENCH}/logs"

sync_results() {
    status=$?

    echo "Syncing results back to final output directory..."
    mkdir -p "${FINAL_OUT_DIR}"

    if [ -d "${WORK_OUT_DIR}" ]; then
        rsync -a "${WORK_OUT_DIR}/" "${FINAL_OUT_DIR}/" || true
    fi

    {
        printf "setting\tvalue\n"
        printf "exit_status\t%s\n" "${status}"
        printf "finished_at\t%s\n" "$(date)"
        printf "host\t%s\n" "$(hostname)"
        printf "job_id\t%s\n" "${JOB_ID:-NA}"
        printf "work_root\t%s\n" "${WORK_ROOT}"
        printf "final_out_dir\t%s\n" "${FINAL_OUT_DIR}"
    } > "${FINAL_OUT_DIR}/job_exit_status.tsv" || true

    exit "${status}"
}

trap sync_results EXIT

count_lines() {
    local file_path="$1"
    case "${file_path}" in
        *.gz)
            gzip -cd "${file_path}" | wc -l | tr -d ' '
            ;;
        *)
            wc -l < "${file_path}" | tr -d ' '
            ;;
    esac
}

if [ -z "${MODEL_JSON}" ]; then
    MODEL_JSON="$(find "${BENCH}/ai_validation" -type f \
        -path '*/runs_ai_full_internal_plasmodium_validation_*/outputs/final_internal_plasmodium_calibrator_all_training.json' \
        | sort \
        | tail -n 1)"
fi

if [ -z "${ZYMO_OUT_DIR}" ]; then
    ZYMO_OUT_DIR="$(find "${BENCH}/screening" -maxdepth 3 -type d \
        -name 'ERR5396170' \
        | grep 'ref_label_minmixed_0p05' \
        | sort \
        | tail -n 1 || true)"
fi

if [ -z "${CALLS_TABLE}" ] && [ -n "${ZYMO_OUT_DIR}" ]; then
    CALLS_TABLE="${ZYMO_OUT_DIR}/species_detection_calls.tsv"
fi

if [ -z "${REFERENCE_LABEL_MAP}" ] && [ -n "${ZYMO_OUT_DIR}" ]; then
    REFERENCE_LABEL_MAP="${ZYMO_OUT_DIR}/reference_label_map.tsv"
fi

echo "BENCH: ${BENCH}"
echo "MODEL_JSON: ${MODEL_JSON}"
echo "ZYMO_OUT_DIR: ${ZYMO_OUT_DIR}"
echo "CALLS_TABLE: ${CALLS_TABLE}"
echo "REFERENCE_LABEL_MAP: ${REFERENCE_LABEL_MAP}"
echo "AI_SCRIPT_DIR: ${AI_SCRIPT_DIR}"
echo "EXTERNAL_SCRIPT: ${EXTERNAL_SCRIPT}"
echo "FINAL_OUT_DIR: ${FINAL_OUT_DIR}"
echo "WORK_ROOT: ${WORK_ROOT}"

if [ ! -s "${MODEL_JSON}" ]; then
    echo "ERROR: MODEL_JSON not found or empty: ${MODEL_JSON}" >&2
    exit 1
fi

if [ ! -s "${CALLS_TABLE}" ]; then
    echo "ERROR: CALLS_TABLE not found or empty: ${CALLS_TABLE}" >&2
    exit 1
fi

if [ ! -s "${EXTERNAL_SCRIPT}" ]; then
    echo "ERROR: external validation Python script not found: ${EXTERNAL_SCRIPT}" >&2
    exit 1
fi

command -v python

rsync -a "${AI_SCRIPT_DIR}/" "${WORK_SCRIPT_DIR}/ai_python/" >&2
export PYTHONPATH="${WORK_SCRIPT_DIR}/ai_python:${PYTHONPATH:-}"

if [ "${RUN_TESTS}" = "true" ]; then
    echo "Running AI validation unit tests..."
    python -m unittest discover \
        -s "${WORK_SCRIPT_DIR}/ai_python/tests" \
        -p 'test_*.py' \
        -v \
        > "${WORK_OUT_DIR}/unit_tests.log" 2>&1
    echo "Unit tests passed."
fi

MODEL_JSON_WORK="${WORK_INPUT_DIR}/$(basename "${MODEL_JSON}")"
CALLS_TABLE_WORK="${WORK_INPUT_DIR}/$(basename "${CALLS_TABLE}")"
REFERENCE_LABEL_MAP_WORK=""

rsync -a "${MODEL_JSON}" "${MODEL_JSON_WORK}" >&2
rsync -a "${CALLS_TABLE}" "${CALLS_TABLE_WORK}" >&2

if [ -s "${REFERENCE_LABEL_MAP}" ]; then
    REFERENCE_LABEL_MAP_WORK="${WORK_INPUT_DIR}/$(basename "${REFERENCE_LABEL_MAP}")"
    rsync -a "${REFERENCE_LABEL_MAP}" "${REFERENCE_LABEL_MAP_WORK}" >&2
else
    echo "WARNING: REFERENCE_LABEL_MAP missing or empty; using built-in expected species list only." >&2
fi

python "${WORK_SCRIPT_DIR}/ai_python/apply_internal_model_to_zymo.py" \
    --model_json "${MODEL_JSON_WORK}" \
    --calls_table "${CALLS_TABLE_WORK}" \
    --reference_label_map "${REFERENCE_LABEL_MAP_WORK}" \
    --out_dir "${WORK_OUT_DIR}" \
    --verbose

{
    printf "setting\tvalue\n"
    printf "run_type\texternal_zymo_cross_domain_ai_validation\n"
    printf "training_source\tinternal_plasmodium_full_model\n"
    printf "external_dataset\tERR5396170_ZymoBIOMICS_D6300_ONT_Q20\n"
    printf "bench\t%s\n" "${BENCH}"
    printf "model_json_original\t%s\n" "${MODEL_JSON}"
    printf "model_json_work\t%s\n" "${MODEL_JSON_WORK}"
    printf "zymo_out_dir\t%s\n" "${ZYMO_OUT_DIR}"
    printf "calls_table_original\t%s\n" "${CALLS_TABLE}"
    printf "calls_table_work\t%s\n" "${CALLS_TABLE_WORK}"
    printf "reference_label_map_original\t%s\n" "${REFERENCE_LABEL_MAP:-not_used}"
    printf "reference_label_map_work\t%s\n" "${REFERENCE_LABEL_MAP_WORK:-not_used}"
    printf "run_tests\t%s\n" "${RUN_TESTS}"
    printf "job_id\t%s\n" "${JOB_ID:-NA}"
    printf "nslots\t%s\n" "${NSLOTS:-4}"
    printf "tmpdir\t%s\n" "${TMPDIR:-NA}"
    printf "finished_at\t%s\n" "$(date)"
} > "${WORK_OUT_DIR}/run_submission_settings.tsv"

{
    printf "file\tline_count\n"
    for file in \
        "${WORK_OUT_DIR}/external_zymo_feature_table.tsv.gz" \
        "${WORK_OUT_DIR}/external_zymo_predictions.tsv.gz" \
        "${WORK_OUT_DIR}/external_zymo_metrics.tsv" \
        "${WORK_OUT_DIR}/external_zymo_truth_label_counts.tsv" \
        "${WORK_OUT_DIR}/external_zymo_prediction_counts.tsv" \
        "${WORK_OUT_DIR}/external_zymo_expected_target_predictions.tsv" \
        "${WORK_OUT_DIR}/unit_tests.log"
    do
        if [ -s "${file}" ]; then
            printf "%s\t%s\n" "$(basename "${file}")" "$(count_lines "${file}")"
        fi
    done
} > "${WORK_OUT_DIR}/output_row_counts.tsv"

mkdir -p "${WORK_OUT_DIR}/scripts_snapshot"
rsync -a "${WORK_SCRIPT_DIR}/ai_python/" "${WORK_OUT_DIR}/scripts_snapshot/ai_python/" >&2
cp "$0" "${WORK_OUT_DIR}/scripts_snapshot/$(basename "$0")" || true

echo "Key output files:"
find "${WORK_OUT_DIR}" -maxdepth 2 -type f | sort

echo "Finished at: $(date)"
