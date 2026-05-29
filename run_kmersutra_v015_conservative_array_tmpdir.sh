#!/usr/bin/env bash
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -j y
#$ -N KSscreen_v015

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

require_file() {
    local file_path="$1"
    local label="$2"

    if [[ ! -s "${file_path}" ]]; then
        die "${label} does not exist or is empty: ${file_path}"
    fi
}

safe_copy_file() {
    local source_path="$1"
    local destination_path="$2"

    if command -v rsync >/dev/null 2>&1; then
        rsync -a "${source_path}" "${destination_path}"
    else
        cp -f "${source_path}" "${destination_path}"
    fi
}

safe_sync_dir() {
    local source_dir="$1"
    local destination_dir="$2"

    mkdir -p "${destination_dir}"
    if command -v rsync >/dev/null 2>&1; then
        rsync -a "${source_dir}/" "${destination_dir}/"
    else
        cp -a "${source_dir}/." "${destination_dir}/"
    fi
}

MANIFEST="${MANIFEST:?MANIFEST is required}"
OUT_ROOT="${OUT_ROOT:?OUT_ROOT is required}"
PANEL="${PANEL:?PANEL is required}"
THREADS="${THREADS:-4}"
CHUNK_SIZE="${CHUNK_SIZE:-10000}"
MAX_PENDING_CHUNKS="${MAX_PENDING_CHUNKS:-}"
CALL_PRESET="${CALL_PRESET:-conservative}"
MIN_UNIQUE_KMER_MARGIN="${MIN_UNIQUE_KMER_MARGIN:-0}"
MIN_UNIQUE_KMER_RATIO="${MIN_UNIQUE_KMER_RATIO:-0.0}"
LOW_EVIDENCE_CALL="${LOW_EVIDENCE_CALL:-observed_below_threshold}"
MIN_UNIQUE_KMERS="${MIN_UNIQUE_KMERS:-}"
MIN_POSITIVE_SEQUENCES="${MIN_POSITIVE_SEQUENCES:-}"
MIN_K_VALUES_POSITIVE="${MIN_K_VALUES_POSITIVE:-}"
MAX_CONFLICT_RATIO="${MAX_CONFLICT_RATIO:-}"
MIN_BEST_K="${MIN_BEST_K:-}"
MIN_EXACT_HITS="${MIN_EXACT_HITS:-}"
MIN_CONFIDENCE_SCORE="${MIN_CONFIDENCE_SCORE:-}"
WRITE_PANEL_CACHE="${WRITE_PANEL_CACHE:-false}"
USE_PANEL_CACHE="${USE_PANEL_CACHE:-false}"
CONSOLIDATE_SPECIES_CALLS="${CONSOLIDATE_SPECIES_CALLS:-false}"
BACKGROUND_CANDIDATE_TAXA="${BACKGROUND_CANDIDATE_TAXA:-}"
BACKGROUND_CANDIDATE_FILE="${BACKGROUND_CANDIDATE_FILE:-}"
DISABLE_SAME_GENUS_NEIGHBOUR_DEMOTION="${DISABLE_SAME_GENUS_NEIGHBOUR_DEMOTION:-false}"
DOMINANT_SPECIES_MIN_MARGIN="${DOMINANT_SPECIES_MIN_MARGIN:-25}"
DOMINANT_SPECIES_MIN_RATIO="${DOMINANT_SPECIES_MIN_RATIO:-2.0}"
WRITE_PARQUET_OUTPUTS="${WRITE_PARQUET_OUTPUTS:-false}"

require_file "${MANIFEST}" "Manifest"
require_file "${PANEL}" "KmerSutra panel"

if [[ -z "${SGE_TASK_ID:-}" ]]; then
    die "SGE_TASK_ID is not set. This script is intended to run as an SGE array task."
fi

line_number=$((SGE_TASK_ID + 1))
manifest_line="$(awk -v line_number="${line_number}" 'NR == line_number {print; exit}' "${MANIFEST}")"

if [[ -z "${manifest_line}" ]]; then
    die "No manifest row found for SGE_TASK_ID=${SGE_TASK_ID} in ${MANIFEST}"
fi

IFS=$'\t' read -r \
    sample_id \
    input_fastq \
    benchmark_family \
    panel_label \
    replicate \
    spike_reads \
    source_run_dir \
    source_relative_dir <<< "${manifest_line}"

if [[ -z "${sample_id}" || -z "${input_fastq}" ]]; then
    die "Manifest row is malformed for SGE_TASK_ID=${SGE_TASK_ID}: ${manifest_line}"
fi

sample_id="${sample_id/kmersutra_v014/kmersutra_v015}"
require_file "${input_fastq}" "Input FASTQ"

sample_final_dir="${OUT_ROOT}/samples/${benchmark_family}/${sample_id}"
status_dir="${OUT_ROOT}/metrics/task_status"
mkdir -p "${status_dir}" "${OUT_ROOT}/logs"

job_tmp_parent="${TMPDIR:-${OUT_ROOT}/tmp}"
work_root="${job_tmp_parent}/kmersutra_v015_${JOB_ID:-manual}_${SGE_TASK_ID}"
local_fastq_dir="${work_root}/fastq"
local_panel_dir="${work_root}/panel"
sample_work_dir="${work_root}/output"
mkdir -p "${local_fastq_dir}" "${local_panel_dir}" "${sample_work_dir}"

status_file="${status_dir}/task_${SGE_TASK_ID}.tsv"
local_panel="${local_panel_dir}/$(basename "${PANEL}")"
local_fastq="${local_fastq_dir}/${sample_id}.fastq.gz"

cleanup() {
    local exit_code=$?
    if [[ -d "${sample_work_dir}" ]]; then
        safe_sync_dir "${sample_work_dir}" "${sample_final_dir}" || true
    fi
    if [[ "${KEEP_TMP:-false}" != "true" ]]; then
        rm -rf "${work_root}" || true
    else
        warn "Keeping temporary directory because KEEP_TMP=true: ${work_root}"
    fi
    exit "${exit_code}"
}
trap cleanup EXIT

log "Starting KmerSutra v0.15 conservative screen task ${SGE_TASK_ID}"
log "Sample: ${sample_id}"
log "Family: ${benchmark_family}; panel=${panel_label}; replicate=${replicate}; spike_reads=${spike_reads}"
log "Input FASTQ: ${input_fastq}"
log "Panel: ${PANEL}"
log "Output: ${sample_final_dir}"
log "TMPDIR work root: ${work_root}"
log "Threads: ${THREADS}; chunk_size=${CHUNK_SIZE}; call_preset=${CALL_PRESET}"
log "Consolidate species calls: ${CONSOLIDATE_SPECIES_CALLS}; write parquet outputs: ${WRITE_PARQUET_OUTPUTS}"

safe_copy_file "${PANEL}" "${local_panel}"
safe_copy_file "${input_fastq}" "${local_fastq}"

command=(
    kmersutra-screen
    --input "${local_fastq}"
    --panel "${local_panel}"
    --sample_id "${sample_id}"
    --input_format fastq
    --out_dir "${sample_work_dir}"
    --threads "${THREADS}"
    --chunk_size "${CHUNK_SIZE}"
    --call_preset "${CALL_PRESET}"
    --min_unique_kmer_margin "${MIN_UNIQUE_KMER_MARGIN}"
    --min_unique_kmer_ratio "${MIN_UNIQUE_KMER_RATIO}"
    --low_evidence_call "${LOW_EVIDENCE_CALL}"
    --no_read_level_hits
    --profile
    --verbose
)

if [[ -n "${MAX_PENDING_CHUNKS}" ]]; then
    command+=(--max_pending_chunks "${MAX_PENDING_CHUNKS}")
fi
if [[ -n "${MIN_UNIQUE_KMERS}" ]]; then
    command+=(--min_unique_kmers "${MIN_UNIQUE_KMERS}")
fi
if [[ -n "${MIN_POSITIVE_SEQUENCES}" ]]; then
    command+=(--min_positive_sequences "${MIN_POSITIVE_SEQUENCES}")
fi
if [[ -n "${MIN_K_VALUES_POSITIVE}" ]]; then
    command+=(--min_k_values_positive "${MIN_K_VALUES_POSITIVE}")
fi
if [[ -n "${MAX_CONFLICT_RATIO}" ]]; then
    command+=(--max_conflict_ratio "${MAX_CONFLICT_RATIO}")
fi
if [[ -n "${MIN_BEST_K}" ]]; then
    command+=(--min_best_k "${MIN_BEST_K}")
fi
if [[ -n "${MIN_EXACT_HITS}" ]]; then
    command+=(--min_exact_hits "${MIN_EXACT_HITS}")
fi
if [[ -n "${MIN_CONFIDENCE_SCORE}" ]]; then
    command+=(--min_confidence_score "${MIN_CONFIDENCE_SCORE}")
fi
if [[ "${WRITE_PANEL_CACHE}" == "true" ]]; then
    command+=(--write_panel_cache --panel_cache "${local_panel}.cache.pkl")
elif [[ "${USE_PANEL_CACHE}" == "true" ]]; then
    command+=(--use_panel_cache --panel_cache "${local_panel}.cache.pkl")
fi
if [[ "${CONSOLIDATE_SPECIES_CALLS}" == "true" ]]; then
    command+=(--consolidate_species_calls)
fi
if [[ -n "${BACKGROUND_CANDIDATE_TAXA}" ]]; then
    read -r -a BACKGROUND_TAXA_ARRAY <<< "${BACKGROUND_CANDIDATE_TAXA}"
    command+=(--background_candidate_taxa "${BACKGROUND_TAXA_ARRAY[@]}")
fi
if [[ -n "${BACKGROUND_CANDIDATE_FILE}" ]]; then
    command+=(--background_candidate_file "${BACKGROUND_CANDIDATE_FILE}")
fi
if [[ "${DISABLE_SAME_GENUS_NEIGHBOUR_DEMOTION}" == "true" ]]; then
    command+=(--disable_same_genus_neighbour_demotion)
fi
command+=(--dominant_species_min_margin "${DOMINANT_SPECIES_MIN_MARGIN}")
command+=(--dominant_species_min_ratio "${DOMINANT_SPECIES_MIN_RATIO}")
if [[ "${WRITE_PARQUET_OUTPUTS}" == "true" ]]; then
    command+=(--write_parquet_outputs)
fi

printf '%s\t%s\n' "field" "value" > "${sample_work_dir}/kmersutra_v015_screen_metadata.tsv"
printf '%s\t%s\n' "sample_id" "${sample_id}" >> "${sample_work_dir}/kmersutra_v015_screen_metadata.tsv"
printf '%s\t%s\n' "benchmark_family" "${benchmark_family}" >> "${sample_work_dir}/kmersutra_v015_screen_metadata.tsv"
printf '%s\t%s\n' "panel" "${panel_label}" >> "${sample_work_dir}/kmersutra_v015_screen_metadata.tsv"
printf '%s\t%s\n' "replicate" "${replicate}" >> "${sample_work_dir}/kmersutra_v015_screen_metadata.tsv"
printf '%s\t%s\n' "spike_reads" "${spike_reads}" >> "${sample_work_dir}/kmersutra_v015_screen_metadata.tsv"
printf '%s\t%s\n' "input_fastq" "${input_fastq}" >> "${sample_work_dir}/kmersutra_v015_screen_metadata.tsv"
printf '%s\t%s\n' "source_run_dir" "${source_run_dir}" >> "${sample_work_dir}/kmersutra_v015_screen_metadata.tsv"
printf '%s\t%s\n' "source_relative_dir" "${source_relative_dir}" >> "${sample_work_dir}/kmersutra_v015_screen_metadata.tsv"
printf '%s\t%s\n' "call_preset" "${CALL_PRESET}" >> "${sample_work_dir}/kmersutra_v015_screen_metadata.tsv"
printf '%s\t%s\n' "low_evidence_call" "${LOW_EVIDENCE_CALL}" >> "${sample_work_dir}/kmersutra_v015_screen_metadata.tsv"
printf '%s\t%s\n' "min_unique_kmer_margin" "${MIN_UNIQUE_KMER_MARGIN}" >> "${sample_work_dir}/kmersutra_v015_screen_metadata.tsv"
printf '%s\t%s\n' "min_unique_kmer_ratio" "${MIN_UNIQUE_KMER_RATIO}" >> "${sample_work_dir}/kmersutra_v015_screen_metadata.tsv"

start_time="$(date '+%Y-%m-%d %H:%M:%S')"
start_epoch="$(date +%s)"
exit_status=0

log "Command: ${command[*]}"
"${command[@]}" || exit_status=$?

end_time="$(date '+%Y-%m-%d %H:%M:%S')"
end_epoch="$(date +%s)"
runtime_seconds=$((end_epoch - start_epoch))

printf 'label\tstart_time\tend_time\truntime_seconds\texit_status\n' \
    > "${sample_work_dir}/screen_task_timing.tsv"
printf 'kmersutra_screen_v015_conservative\t%s\t%s\t%s\t%s\n' \
    "${start_time}" "${end_time}" "${runtime_seconds}" "${exit_status}" \
    >> "${sample_work_dir}/screen_task_timing.tsv"

printf 'sample_id\tbenchmark_family\tpanel\treplicate\tspike_reads\tinput_fastq\tsample_out\truntime_seconds\texit_status\n' \
    > "${status_file}"
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "${sample_id}" \
    "${benchmark_family}" \
    "${panel_label}" \
    "${replicate}" \
    "${spike_reads}" \
    "${input_fastq}" \
    "${sample_final_dir}" \
    "${runtime_seconds}" \
    "${exit_status}" \
    >> "${status_file}"

safe_sync_dir "${sample_work_dir}" "${sample_final_dir}"

if [[ "${exit_status}" -ne 0 ]]; then
    die "KmerSutra screening failed for ${sample_id}; exit_status=${exit_status}"
fi

log "Finished KmerSutra v0.15 conservative screen for ${sample_id}; runtime=${runtime_seconds}s"
