#!/bin/bash
#$ -N KSai_ext_zymo51
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
NOVELTY_SCALE="${NOVELTY_SCALE:-1.0}"
RUN_UNIT_TESTS="${RUN_UNIT_TESTS:-0}"
TEST_START_DIR="${TEST_START_DIR:-tests}"

RUN_TAG="ai_external_zymo_validation_v0501_safe_transformed_bounded_$(date +%Y%m%d_%H%M%S)"
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

if [ -z "${MODEL_JSON}" ]; then
    MODEL_JSON="$(find "${BENCH}/ai_validation" -type f \
        -path "*/runs_ai_full_internal_plasmodium_validation_v0501_*/outputs/final_internal_calibrator_all_training.json" \
        | sort | tail -n 1)"
fi

if [ -z "${ZYMO_OUT_DIR}" ]; then
    ZYMO_OUT_DIR="$(find "${BENCH}/screening" -maxdepth 2 -type d \
        -name "ERR5396170" \
        | grep "ref_label_minmixed_0p05" \
        | sort | tail -n 1 || true)"
fi

if [ -z "${CALLS_TABLE}" ]; then
    CALLS_TABLE="${ZYMO_OUT_DIR}/species_detection_calls.tsv"
fi

if [ -z "${REFERENCE_LABEL_MAP}" ]; then
    REFERENCE_LABEL_MAP="${ZYMO_OUT_DIR}/reference_label_map.tsv"
fi

echo "BENCH: ${BENCH}"
echo "MODEL_JSON: ${MODEL_JSON}"
echo "ZYMO_OUT_DIR: ${ZYMO_OUT_DIR}"
echo "CALLS_TABLE: ${CALLS_TABLE}"
echo "REFERENCE_LABEL_MAP: ${REFERENCE_LABEL_MAP}"
echo "NOVELTY_SCALE: ${NOVELTY_SCALE}"
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
if [ ! -s "${REFERENCE_LABEL_MAP}" ]; then
    echo "ERROR: REFERENCE_LABEL_MAP not found or empty: ${REFERENCE_LABEL_MAP}" >&2
    exit 1
fi

command -v kmersutra-label-zymo-calls
command -v kmersutra-predict
command -v kmersutra-evaluate-call-predictions

if [ "${RUN_UNIT_TESTS}" = "1" ]; then
    echo "Running unit tests before validation..."
    python -m unittest discover -s "${TEST_START_DIR}" \
        > "${WORK_OUT_DIR}/unit_tests.log" 2>&1
fi

MODEL_JSON_WORK="${WORK_INPUT_DIR}/$(basename "${MODEL_JSON}")"
CALLS_TABLE_WORK="${WORK_INPUT_DIR}/$(basename "${CALLS_TABLE}")"
REFERENCE_LABEL_MAP_WORK="${WORK_INPUT_DIR}/$(basename "${REFERENCE_LABEL_MAP}")"
rsync -a "${MODEL_JSON}" "${MODEL_JSON_WORK}" >&2
rsync -a "${CALLS_TABLE}" "${CALLS_TABLE_WORK}" >&2
rsync -a "${REFERENCE_LABEL_MAP}" "${REFERENCE_LABEL_MAP_WORK}" >&2

FEATURE_TABLE="${WORK_OUT_DIR}/external_zymo_feature_table.tsv.gz"
PREDICTIONS_TABLE="${WORK_OUT_DIR}/external_zymo_predictions.tsv.gz"

echo "Labelling Zymo calls with fine truth categories..."
kmersutra-label-zymo-calls \
    --calls_table "${CALLS_TABLE_WORK}" \
    --reference_label_map "${REFERENCE_LABEL_MAP_WORK}" \
    --out_table "${FEATURE_TABLE}" \
    --out_category_counts "${WORK_OUT_DIR}/external_zymo_truth_category_counts.tsv" \
    --out_coarse_label_counts "${WORK_OUT_DIR}/external_zymo_coarse_label_counts.tsv" \
    --verbose

echo "Applying internal transformed-feature model to Zymo table..."
kmersutra-predict \
    --features_tsv "${FEATURE_TABLE}" \
    --model_json "${MODEL_JSON_WORK}" \
    --out_tsv "${PREDICTIONS_TABLE}" \
    --novelty_scale "${NOVELTY_SCALE}" \
    --verbose

echo "Evaluating Zymo predictions..."
kmersutra-evaluate-call-predictions \
    --predictions_table "${PREDICTIONS_TABLE}" \
    --out_metrics "${WORK_OUT_DIR}/external_zymo_metrics.tsv" \
    --out_prediction_counts "${WORK_OUT_DIR}/external_zymo_prediction_counts.tsv" \
    --out_label_counts "${WORK_OUT_DIR}/external_zymo_ml_label_counts.tsv" \
    --verbose

echo "Evaluating strict Zymo predictions with same-species compatible rows excluded..."
kmersutra-evaluate-call-predictions \
    --predictions_table "${PREDICTIONS_TABLE}" \
    --exclude_label same_species_compatible_reference \
    --out_metrics "${WORK_OUT_DIR}/external_zymo_strict_metrics.tsv" \
    --out_prediction_counts "${WORK_OUT_DIR}/external_zymo_strict_prediction_counts.tsv" \
    --out_label_counts "${WORK_OUT_DIR}/external_zymo_strict_ml_label_counts.tsv" \
    --verbose

{
    printf "setting\tvalue\n"
    printf "run_type\texternal_zymo_cross_domain_ai_validation_v050\n"
    printf "training_source\tinternal_plasmodium_safe_transformed_model\n"
    printf "external_dataset\tERR5396170_ZymoBIOMICS_D6300_ONT_Q20\n"
    printf "bench\t%s\n" "${BENCH}"
    printf "model_json_original\t%s\n" "${MODEL_JSON}"
    printf "model_json_work\t%s\n" "${MODEL_JSON_WORK}"
    printf "zymo_out_dir\t%s\n" "${ZYMO_OUT_DIR}"
    printf "calls_table_original\t%s\n" "${CALLS_TABLE}"
    printf "calls_table_work\t%s\n" "${CALLS_TABLE_WORK}"
    printf "reference_label_map_original\t%s\n" "${REFERENCE_LABEL_MAP}"
    printf "reference_label_map_work\t%s\n" "${REFERENCE_LABEL_MAP_WORK}"
    printf "novelty_scale\t%s\n" "${NOVELTY_SCALE}"
    printf "job_id\t%s\n" "${JOB_ID:-NA}"
    printf "nslots\t%s\n" "${NSLOTS:-4}"
    printf "tmpdir\t%s\n" "${TMPDIR:-NA}"
    printf "finished_at\t%s\n" "$(date)"
} > "${WORK_OUT_DIR}/run_submission_settings.tsv"

{
    printf "file\tline_count\n"
    for file in \
        "${FEATURE_TABLE}" \
        "${PREDICTIONS_TABLE}" \
        "${WORK_OUT_DIR}/external_zymo_metrics.tsv" \
        "${WORK_OUT_DIR}/external_zymo_truth_category_counts.tsv" \
        "${WORK_OUT_DIR}/external_zymo_coarse_label_counts.tsv" \
        "${WORK_OUT_DIR}/external_zymo_prediction_counts.tsv" \
        "${WORK_OUT_DIR}/external_zymo_ml_label_counts.tsv" \
        "${WORK_OUT_DIR}/external_zymo_strict_metrics.tsv" \
        "${WORK_OUT_DIR}/external_zymo_strict_prediction_counts.tsv" \
        "${WORK_OUT_DIR}/external_zymo_strict_ml_label_counts.tsv"
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

find "${WORK_OUT_DIR}" -maxdepth 1 -type f | sort

echo "Finished at: $(date)"
