#!/bin/bash
#$ -N KSai_full_v501
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
CALLS_TABLE="${CALLS_TABLE:-${BENCH}/tables/final_kmersutra_detection_calls_long.tsv.gz}"
MAX_NOT_DETECTED="${MAX_NOT_DETECTED:-50000}"
DISTANCE_QUANTILE="${DISTANCE_QUANTILE:-0.95}"
FEATURE_PROFILE="${FEATURE_PROFILE:-safe_transformed_bounded}"
SAMPLE_TEST_FRACTION="${SAMPLE_TEST_FRACTION:-0.20}"
RUN_UNIT_TESTS="${RUN_UNIT_TESTS:-0}"
TEST_START_DIR="${TEST_START_DIR:-tests}"

RUN_TAG="ai_full_internal_plasmodium_validation_v0501_${FEATURE_PROFILE}_$(date +%Y%m%d_%H%M%S)"
FINAL_RUN_DIR="${BENCH}/ai_validation/runs_${RUN_TAG}"
FINAL_OUT_DIR="${FINAL_RUN_DIR}/outputs"
WORK_ROOT="${TMPDIR:-/tmp}/${USER}/kmersutra_${RUN_TAG}"
WORK_INPUT_DIR="${WORK_ROOT}/inputs"
WORK_OUT_DIR="${WORK_ROOT}/outputs"

mkdir -p "${WORK_INPUT_DIR}" "${WORK_OUT_DIR}" "${FINAL_OUT_DIR}" "${BENCH}/logs"

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

echo "BENCH: ${BENCH}"
echo "CALLS_TABLE: ${CALLS_TABLE}"
echo "FEATURE_PROFILE: ${FEATURE_PROFILE}"
echo "FINAL_OUT_DIR: ${FINAL_OUT_DIR}"
echo "WORK_ROOT: ${WORK_ROOT}"

if [ ! -s "${CALLS_TABLE}" ]; then
    echo "ERROR: CALLS_TABLE not found or empty: ${CALLS_TABLE}" >&2
    exit 1
fi

command -v kmersutra-build-call-training
command -v kmersutra-validate-call-calibrator

if [ "${RUN_UNIT_TESTS}" = "1" ]; then
    echo "Running unit tests before validation..."
    python -m unittest discover -s "${TEST_START_DIR}" \
        > "${WORK_OUT_DIR}/unit_tests.log" 2>&1
fi

CALLS_TABLE_WORK="${WORK_INPUT_DIR}/$(basename "${CALLS_TABLE}")"
TRAINING_TABLE="${WORK_OUT_DIR}/ai_call_training.tsv.gz"
rsync -a "${CALLS_TABLE}" "${CALLS_TABLE_WORK}" >&2

echo "Building AI-ready call-training table..."
kmersutra-build-call-training \
    --calls_table "${CALLS_TABLE_WORK}" \
    --out_table "${TRAINING_TABLE}" \
    --max_not_detected "${MAX_NOT_DETECTED}" \
    --verbose

echo "Running full internal holdout validation..."
kmersutra-validate-call-calibrator \
    --training_table "${TRAINING_TABLE}" \
    --out_dir "${WORK_OUT_DIR}" \
    --feature_profile "${FEATURE_PROFILE}" \
    --distance_quantile "${DISTANCE_QUANTILE}" \
    --sample_test_fraction "${SAMPLE_TEST_FRACTION}" \
    --verbose

{
    printf "setting\tvalue\n"
    printf "run_type\tfull_internal_plasmodium_ai_holdout_validation_v050\n"
    printf "bench\t%s\n" "${BENCH}"
    printf "calls_table_original\t%s\n" "${CALLS_TABLE}"
    printf "calls_table_work\t%s\n" "${CALLS_TABLE_WORK}"
    printf "training_table\t%s\n" "${TRAINING_TABLE}"
    printf "feature_profile\t%s\n" "${FEATURE_PROFILE}"
    printf "max_not_detected\t%s\n" "${MAX_NOT_DETECTED}"
    printf "distance_quantile\t%s\n" "${DISTANCE_QUANTILE}"
    printf "sample_test_fraction\t%s\n" "${SAMPLE_TEST_FRACTION}"
    printf "final_model\t%s\n" "${WORK_OUT_DIR}/final_internal_calibrator_all_training.json"
    printf "job_id\t%s\n" "${JOB_ID:-NA}"
    printf "nslots\t%s\n" "${NSLOTS:-4}"
    printf "tmpdir\t%s\n" "${TMPDIR:-NA}"
    printf "finished_at\t%s\n" "$(date)"
} > "${WORK_OUT_DIR}/run_submission_settings.tsv"

{
    printf "file\tline_count\n"
    for file in \
        "${TRAINING_TABLE}" \
        "${WORK_OUT_DIR}/all_training_label_counts.tsv" \
        "${WORK_OUT_DIR}/validation_design.tsv" \
        "${WORK_OUT_DIR}/skipped_splits.tsv" \
        "${WORK_OUT_DIR}/holdout_metrics.tsv" \
        "${WORK_OUT_DIR}/holdout_predictions.tsv.gz" \
        "${WORK_OUT_DIR}/holdout_summary_by_split.tsv" \
        "${WORK_OUT_DIR}/final_internal_model_training_summary.tsv"
    do
        if [ -s "${file}" ]; then
            case "${file}" in
                *.gz) count="$(gzip -cd "${file}" | wc -l | tr -d ' ')" ;;
                *) count="$(wc -l < "${file}" | tr -d ' ')" ;;
            esac
            printf "%s\t%s\n" "$(basename "${file}")" "${count}"
        fi
    done
} > "${WORK_OUT_DIR}/output_row_counts.tsv"

find "${WORK_OUT_DIR}" -maxdepth 2 -type f | sort

echo "Finished at: $(date)"
