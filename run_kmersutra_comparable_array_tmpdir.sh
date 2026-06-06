#!/usr/bin/env bash
set -euo pipefail

log() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') INFO  $*"
}

warn() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') WARN  $*" >&2
}

die() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') ERROR $*" >&2
    exit 1
}

require_env() {
    local name="$1"
    if [[ -z "${!name:-}" ]]; then
        die "Required environment variable is unset: ${name}"
    fi
}

require_file() {
    local path="$1"
    local label="$2"
    if [[ ! -s "${path}" ]]; then
        die "${label} missing or empty: ${path}"
    fi
}

require_env MANIFEST
require_env OUT_ROOT
require_env PANEL

TASK_ID="${SGE_TASK_ID:-${SLURM_ARRAY_TASK_ID:-1}}"
THREADS="${THREADS:-1}"
CHUNK_SIZE="${CHUNK_SIZE:-10000}"
MAX_PENDING_CHUNKS="${MAX_PENDING_CHUNKS:-}"
DECOMPRESSOR="${DECOMPRESSOR:-python}"
CALL_PRESET="${CALL_PRESET:-conservative}"
SCREEN_PRESET="${SCREEN_PRESET:-exact}"
MAX_MISMATCHES="${MAX_MISMATCHES:-}"
FUZZY_MIN_K="${FUZZY_MIN_K:-}"
LOW_EVIDENCE_CALL="${LOW_EVIDENCE_CALL:-observed_below_threshold}"
MIN_UNIQUE_KMER_MARGIN="${MIN_UNIQUE_KMER_MARGIN:-0}"
MIN_UNIQUE_KMER_RATIO="${MIN_UNIQUE_KMER_RATIO:-0.0}"
MIN_UNIQUE_KMERS="${MIN_UNIQUE_KMERS:-}"
MIN_POSITIVE_SEQUENCES="${MIN_POSITIVE_SEQUENCES:-}"
MIN_K_VALUES_POSITIVE="${MIN_K_VALUES_POSITIVE:-}"
MAX_CONFLICT_RATIO="${MAX_CONFLICT_RATIO:-}"
MIN_BEST_K="${MIN_BEST_K:-}"
MIN_EXACT_HITS="${MIN_EXACT_HITS:-}"
MIN_TOTAL_HITS="${MIN_TOTAL_HITS:-}"
MIN_CONFIDENCE_SCORE="${MIN_CONFIDENCE_SCORE:-}"
CONSOLIDATE_SPECIES_CALLS="${CONSOLIDATE_SPECIES_CALLS:-false}"
BACKGROUND_CANDIDATE_TAXA="${BACKGROUND_CANDIDATE_TAXA:-}"
BACKGROUND_CANDIDATE_FILE="${BACKGROUND_CANDIDATE_FILE:-}"
DISABLE_SAME_GENUS_NEIGHBOUR_DEMOTION="${DISABLE_SAME_GENUS_NEIGHBOUR_DEMOTION:-false}"
DOMINANT_SPECIES_MIN_MARGIN="${DOMINANT_SPECIES_MIN_MARGIN:-25}"
DOMINANT_SPECIES_MIN_RATIO="${DOMINANT_SPECIES_MIN_RATIO:-2.0}"
WRITE_PARQUET_OUTPUTS="${WRITE_PARQUET_OUTPUTS:-false}"
KEEP_TMP="${KEEP_TMP:-false}"

require_file "${MANIFEST}" "Comparable manifest"
require_file "${PANEL}" "KmerSutra panel"

ROW=$(awk -F '\t' -v task_id="${TASK_ID}" 'NR == task_id + 1 {print; exit}' "${MANIFEST}")
if [[ -z "${ROW}" ]]; then
    die "No manifest row for task ${TASK_ID} in ${MANIFEST}"
fi

IFS=$'\t' read -r SAMPLE_ID INPUT_FASTQ BENCHMARK_FAMILY PANEL_NAME REPLICATE SPIKE_READS SOURCE_RUN_DIR SOURCE_RELATIVE_DIR <<< "${ROW}"
if [[ -z "${SAMPLE_ID}" || -z "${INPUT_FASTQ}" ]]; then
    die "Manifest row ${TASK_ID} has empty sample_id or input_fastq"
fi
require_file "${INPUT_FASTQ}" "Input FASTQ"

SAMPLE_OUT_DIR="${OUT_ROOT}/samples/${BENCHMARK_FAMILY}/${SAMPLE_ID}"
mkdir -p "${SAMPLE_OUT_DIR}"

TMP_BASE="${TMPDIR:-${OUT_ROOT}/tmp}"
TASK_TMP="$(mktemp -d "${TMP_BASE%/}/kmersutra_${SAMPLE_ID}_XXXXXX")"
cleanup() {
    if [[ "${KEEP_TMP}" == "true" ]]; then
        warn "Keeping temporary directory: ${TASK_TMP}"
    else
        rm -rf "${TASK_TMP}"
    fi
}
trap cleanup EXIT

LOCAL_PANEL="${TASK_TMP}/$(basename "${PANEL}")"
LOCAL_FASTQ="${TASK_TMP}/$(basename "${INPUT_FASTQ}")"
cp "${PANEL}" "${LOCAL_PANEL}"
cp "${INPUT_FASTQ}" "${LOCAL_FASTQ}"

ARGS=(
    --input "${LOCAL_FASTQ}"
    --panel "${LOCAL_PANEL}"
    --sample_id "${SAMPLE_ID}"
    --input_format fastq
    --decompressor "${DECOMPRESSOR}"
    --out_dir "${SAMPLE_OUT_DIR}"
    --screen_preset "${SCREEN_PRESET}"
    --threads "${THREADS}"
    --chunk_size "${CHUNK_SIZE}"
    --call_preset "${CALL_PRESET}"
    --low_evidence_call "${LOW_EVIDENCE_CALL}"
    --min_unique_kmer_margin "${MIN_UNIQUE_KMER_MARGIN}"
    --min_unique_kmer_ratio "${MIN_UNIQUE_KMER_RATIO}"
    --dominant_species_min_margin "${DOMINANT_SPECIES_MIN_MARGIN}"
    --dominant_species_min_ratio "${DOMINANT_SPECIES_MIN_RATIO}"
    --profile
    --no_read_level_hits
    --verbose
)

if [[ -n "${MAX_PENDING_CHUNKS}" ]]; then
    ARGS+=(--max_pending_chunks "${MAX_PENDING_CHUNKS}")
fi
if [[ -n "${MAX_MISMATCHES}" ]]; then
    ARGS+=(--max_mismatches "${MAX_MISMATCHES}")
fi
if [[ -n "${FUZZY_MIN_K}" ]]; then
    ARGS+=(--fuzzy_min_k "${FUZZY_MIN_K}")
fi
if [[ -n "${MIN_UNIQUE_KMERS}" ]]; then
    ARGS+=(--min_unique_kmers "${MIN_UNIQUE_KMERS}")
fi
if [[ -n "${MIN_POSITIVE_SEQUENCES}" ]]; then
    ARGS+=(--min_positive_sequences "${MIN_POSITIVE_SEQUENCES}")
fi
if [[ -n "${MIN_K_VALUES_POSITIVE}" ]]; then
    ARGS+=(--min_k_values_positive "${MIN_K_VALUES_POSITIVE}")
fi
if [[ -n "${MAX_CONFLICT_RATIO}" ]]; then
    ARGS+=(--max_conflict_ratio "${MAX_CONFLICT_RATIO}")
fi
if [[ -n "${MIN_BEST_K}" ]]; then
    ARGS+=(--min_best_k "${MIN_BEST_K}")
fi
if [[ -n "${MIN_EXACT_HITS}" ]]; then
    ARGS+=(--min_exact_hits "${MIN_EXACT_HITS}")
fi
if [[ -n "${MIN_TOTAL_HITS}" ]]; then
    ARGS+=(--min_total_hits "${MIN_TOTAL_HITS}")
fi
if [[ -n "${MIN_CONFIDENCE_SCORE}" ]]; then
    ARGS+=(--min_confidence_score "${MIN_CONFIDENCE_SCORE}")
fi
if [[ "${CONSOLIDATE_SPECIES_CALLS}" == "true" ]]; then
    ARGS+=(--consolidate_species_calls)
fi
if [[ -n "${BACKGROUND_CANDIDATE_TAXA}" ]]; then
    # shellcheck disable=SC2206
    BACKGROUND_TAXA_ARRAY=(${BACKGROUND_CANDIDATE_TAXA})
    ARGS+=(--background_candidate_taxa "${BACKGROUND_TAXA_ARRAY[@]}")
fi
if [[ -n "${BACKGROUND_CANDIDATE_FILE}" ]]; then
    ARGS+=(--background_candidate_file "${BACKGROUND_CANDIDATE_FILE}")
fi
if [[ "${DISABLE_SAME_GENUS_NEIGHBOUR_DEMOTION}" == "true" ]]; then
    ARGS+=(--disable_same_genus_neighbour_demotion)
fi
if [[ "${WRITE_PARQUET_OUTPUTS}" == "true" ]]; then
    ARGS+=(--write_parquet_outputs)
fi

COMMAND_TSV="${SAMPLE_OUT_DIR}/kmersutra_screen.command.tsv"
{
    printf 'field\tvalue\n'
    printf 'sample_id\t%s\n' "${SAMPLE_ID}"
    printf 'benchmark_family\t%s\n' "${BENCHMARK_FAMILY}"
    printf 'panel\t%s\n' "${PANEL_NAME}"
    printf 'replicate\t%s\n' "${REPLICATE}"
    printf 'spike_reads\t%s\n' "${SPIKE_READS}"
    printf 'screen_preset\t%s\n' "${SCREEN_PRESET}"
    printf 'call_preset\t%s\n' "${CALL_PRESET}"
    printf 'max_mismatches\t%s\n' "${MAX_MISMATCHES}"
    printf 'fuzzy_min_k\t%s\n' "${FUZZY_MIN_K}"
    printf 'command\t'; printf '%q ' kmersutra-screen "${ARGS[@]}"; printf '\n'
} > "${COMMAND_TSV}"

log "Running sample ${SAMPLE_ID}; screen_preset=${SCREEN_PRESET}; call_preset=${CALL_PRESET}"
START_SECONDS=$(date +%s)
kmersutra-screen "${ARGS[@]}"
END_SECONDS=$(date +%s)

cat > "${SAMPLE_OUT_DIR}/screen_task_timing.tsv" <<EOF
field	value
start_epoch	${START_SECONDS}
end_epoch	${END_SECONDS}
elapsed_seconds	$((END_SECONDS - START_SECONDS))
threads	${THREADS}
chunk_size	${CHUNK_SIZE}
screen_preset	${SCREEN_PRESET}
call_preset	${CALL_PRESET}
EOF

log "Completed sample ${SAMPLE_ID} in $((END_SECONDS - START_SECONDS)) second(s)"
