#!/usr/bin/env bash
#$ -cwd
#$ -j y
#$ -pe smp 24
#$ -jc long
#$ -mods l_hard mfree 100G
#$ -adds l_hard h_vmem 100G
#$ -N KSbuild_global_v028

set -euo pipefail

log_info() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') INFO  $*" >&2
}

log_warn() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') WARN  $*" >&2
}

log_error() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') ERROR $*" >&2
}

fail() {
    log_error "$*"
    exit 1
}

sync_back_outputs() {
    if [ -d "${WORK_OUT_DIR:-}" ]; then
        log_info "Syncing useful outputs back to: ${FINAL_OUT_DIR}"
        mkdir -p "${FINAL_OUT_DIR}"
        if [ "${KEEP_SQLITE:-false}" != "true" ]; then
            find "${WORK_OUT_DIR}" -maxdepth 1 -type f -name '*.sqlite*' -delete || true
        fi
        rsync -av "${WORK_OUT_DIR}/" "${FINAL_OUT_DIR}/" >&2 || true
    fi
}

run_command() {
    local label="$1"
    shift
    local start_time
    local start_epoch
    local end_time
    local end_epoch
    local runtime_seconds
    local exit_status

    log_info "Starting command: ${label}"
    log_info "Command: $*"
    start_time="$(date '+%Y-%m-%d %H:%M:%S')"
    start_epoch="$(date +%s)"
    set +e
    "$@"
    exit_status="$?"
    set -e
    end_time="$(date '+%Y-%m-%d %H:%M:%S')"
    end_epoch="$(date +%s)"
    runtime_seconds=$((end_epoch - start_epoch))
    printf '%s\t%s\t%s\t%s\t%s\n' \
        "${label}" "${start_time}" "${end_time}" "${runtime_seconds}" "${exit_status}" \
        >> "${COMMAND_TIMING_TSV}"
    log_info "Finished command: ${label}; runtime=${runtime_seconds}s; exit_status=${exit_status}"
    if [ "${exit_status}" -ne 0 ]; then
        fail "Command failed: ${label}"
    fi
}

DB_ROOT="${DB_ROOT:-/home/pthorpe001/data/databases/kmersutra_db}"
SOURCE_CONFIG="${SOURCE_CONFIG:-${DB_ROOT}/ncbi_genomes_plasmodium_outgroups_v3/kmersutra_genome_config.tsv}"
RUN_STAMP="${RUN_STAMP:-$(date +%Y%m%d_%H%M%S)}"
BUILD_ROOT="${BUILD_ROOT:-${DB_ROOT}/kmersutra_builds}"
K_VALUES="${K_VALUES:-77 101}"
K_LABEL="$(printf '%s' "${K_VALUES}" | tr ' ' '_')"
FINAL_OUT_DIR="${OUT_DIR:-${BUILD_ROOT}/kmersutra_plasmodium_outgroups_v3_global_candidate_k${K_LABEL}_${RUN_STAMP}}"
TMP_PARENT="${TMPDIR:-/tmp}"
WORK_OUT_DIR="${TMP_PARENT}/$(basename "${FINAL_OUT_DIR}")"
TAXONOMY_DIR="${TAXONOMY_DIR:-${DB_ROOT}/ncbi_taxonomy}"
KMERSUTRA_THREADS="${KMERSUTRA_THREADS:-${NSLOTS:-24}}"
MAX_KMERSUTRA_THREADS="${MAX_KMERSUTRA_THREADS:-24}"
SQLITE_BATCH_SIZE="${SQLITE_BATCH_SIZE:-50000}"
MAX_PER_SPECIES_PER_K="${MAX_PER_SPECIES_PER_K:-100000}"
GLOBAL_SOURCE_INDEX_MODE="${GLOBAL_SOURCE_INDEX_MODE:-candidate_universe}"
GLOBAL_INDEX_PROGRESS_INTERVAL="${GLOBAL_INDEX_PROGRESS_INTERVAL:-1000000}"
MARKER_SELECTION="${MARKER_SELECTION:-genome_spread}"
GENOME_BIN_SIZE="${GENOME_BIN_SIZE:-10000}"
MAX_PER_GENOME_BIN="${MAX_PER_GENOME_BIN:-10}"
WRITE_MODULE_MANIFEST="${WRITE_MODULE_MANIFEST:-true}"
MODULE_MANIFEST_DIR="${MODULE_MANIFEST_DIR:-}"
MODULE_MAX_GATE_RECORDS_PER_K="${MODULE_MAX_GATE_RECORDS_PER_K:-0}"
MODULE_MIN_GATE_UNIQUE_KMERS="${MODULE_MIN_GATE_UNIQUE_KMERS:-1}"
MODULE_MIN_GATE_POSITIVE_SEQUENCES="${MODULE_MIN_GATE_POSITIVE_SEQUENCES:-1}"
MODULE_MIN_GATE_K_VALUES="${MODULE_MIN_GATE_K_VALUES:-1}"
MODULE_MIN_GATE_BEST_K="${MODULE_MIN_GATE_BEST_K:-0}"

WRITE_MODULE_PARQUET="${WRITE_MODULE_PARQUET:-false}"
MODULE_PARQUET_DIR="${MODULE_PARQUET_DIR:-${WORK_OUT_DIR}/module_parquet}"
MODULE_NAME="${MODULE_NAME:-$(basename "${FINAL_OUT_DIR}")}"
PANEL_STORAGE_FORMAT="${PANEL_STORAGE_FORMAT:-auto}"
RAM_LOG_INTERVAL_SECONDS="${RAM_LOG_INTERVAL_SECONDS:-30}"
KEEP_SQLITE="${KEEP_SQLITE:-false}"
EVIDENCE_RANKS="${EVIDENCE_RANKS:-species genus family order class phylum superkingdom}"

if [ "${KMERSUTRA_THREADS}" -gt "${MAX_KMERSUTRA_THREADS}" ]; then
    log_warn "Capping KMERSUTRA_THREADS from ${KMERSUTRA_THREADS} to ${MAX_KMERSUTRA_THREADS}"
    KMERSUTRA_THREADS="${MAX_KMERSUTRA_THREADS}"
fi

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export MALLOC_ARENA_MAX=2

INPUTS_DIR="${WORK_OUT_DIR}/inputs"
METRICS_DIR="${WORK_OUT_DIR}/metrics"
TARGET_CONFIG="${INPUTS_DIR}/kmersutra_genome_config_global_candidate.tsv"
RUN_METADATA_TSV="${WORK_OUT_DIR}/run_metadata.tsv"
COMMAND_TSV="${WORK_OUT_DIR}/kmersutra_build.command.tsv"
COMMAND_TIMING_TSV="${METRICS_DIR}/command_timing.tsv"
POST_BUILD_CHECKS_TSV="${METRICS_DIR}/post_build_checks.tsv"
RAM_LOG_TSV="${WORK_OUT_DIR}/ram_usage.tsv"

trap sync_back_outputs EXIT

rm -rf "${WORK_OUT_DIR}"
mkdir -p "${INPUTS_DIR}" "${METRICS_DIR}" "${FINAL_OUT_DIR}"

log_info "Starting KmerSutra global candidate panel build"
log_info "Host: $(hostname)"
log_info "DB root: ${DB_ROOT}"
log_info "Source config: ${SOURCE_CONFIG}"
log_info "Work output directory: ${WORK_OUT_DIR}"
log_info "Final output directory: ${FINAL_OUT_DIR}"
log_info "TMPDIR: ${TMPDIR:-unset}"
log_info "K values: ${K_VALUES}"
log_info "Threads: ${KMERSUTRA_THREADS}"
log_info "Keep SQLite: ${KEEP_SQLITE}"
log_info "Global source-index mode: ${GLOBAL_SOURCE_INDEX_MODE}"
log_info "Candidate-universe mode samples genome-spread candidates before conflict annotation"
log_info "Global index progress interval: ${GLOBAL_INDEX_PROGRESS_INTERVAL}"
log_info "Marker selection: ${MARKER_SELECTION}"
log_info "Genome bin size: ${GENOME_BIN_SIZE}"
log_info "Max per genome bin: ${MAX_PER_GENOME_BIN}"
log_info "Write module Parquet: ${WRITE_MODULE_PARQUET}"
log_info "Write module manifest: ${WRITE_MODULE_MANIFEST}"
log_info "Module manifest dir: ${MODULE_MANIFEST_DIR:-${WORK_OUT_DIR}/hierarchical_modules}"
log_info "Panel storage format: ${PANEL_STORAGE_FORMAT}"

df -h "${TMP_PARENT}" >&2 || true

command -v kmersutra-build-panel >/dev/null 2>&1 || fail "kmersutra-build-panel not found on PATH"
[ -s "${SOURCE_CONFIG}" ] || fail "Source config missing or empty: ${SOURCE_CONFIG}"
[ -d "${TAXONOMY_DIR}" ] || fail "Taxonomy directory missing: ${TAXONOMY_DIR}"

{
    printf 'field\tvalue\n'
    printf 'run_stamp\t%s\n' "${RUN_STAMP}"
    printf 'db_root\t%s\n' "${DB_ROOT}"
    printf 'source_config\t%s\n' "${SOURCE_CONFIG}"
    printf 'target_config\t%s\n' "${TARGET_CONFIG}"
    printf 'work_out_dir\t%s\n' "${WORK_OUT_DIR}"
    printf 'final_out_dir\t%s\n' "${FINAL_OUT_DIR}"
    printf 'taxonomy_dir\t%s\n' "${TAXONOMY_DIR}"
    printf 'k_values\t%s\n' "${K_VALUES}"
    printf 'kmersutra_threads\t%s\n' "${KMERSUTRA_THREADS}"
    printf 'sqlite_batch_size\t%s\n' "${SQLITE_BATCH_SIZE}"
    printf 'max_per_species_per_k\t%s\n' "${MAX_PER_SPECIES_PER_K}"
    printf 'global_source_index_mode\t%s\n' "${GLOBAL_SOURCE_INDEX_MODE}"
    printf 'global_index_progress_interval\t%s\n' "${GLOBAL_INDEX_PROGRESS_INTERVAL}"
    printf 'ram_log_interval_seconds\t%s\n' "${RAM_LOG_INTERVAL_SECONDS}"
    printf 'keep_sqlite\t%s\n' "${KEEP_SQLITE}"
    printf 'job_id\t%s\n' "${JOB_ID:-NA}"
    printf 'nslots\t%s\n' "${NSLOTS:-NA}"
    printf 'start_time\t%s\n' "$(date '+%Y-%m-%d %H:%M:%S')"
} > "${RUN_METADATA_TSV}"

printf 'label\tstart_time\tend_time\truntime_seconds\texit_status\n' > "${COMMAND_TIMING_TSV}"

log_info "Creating global-candidate genome config"
awk -F '\t' 'BEGIN {OFS="\t"}
NR == 1 {print; next}
$2 == "Plasmodium falciparum" || $2 == "Plasmodium vivax" || $2 == "Plasmodium knowlesi" {$5 = "target_species"}
{print}' "${SOURCE_CONFIG}" > "${TARGET_CONFIG}"

[ -s "${TARGET_CONFIG}" ] || fail "Target config was not created: ${TARGET_CONFIG}"

log_info "Checking FASTA files exist and are non-empty"
MISSING_FASTA_COUNT="$(
    awk -F '\t' 'NR > 1 {print $1}' "${TARGET_CONFIG}" \
        | while read -r fasta_path; do
            if [ ! -s "${fasta_path}" ]; then
                printf 'MISSING_OR_EMPTY\t%s\n' "${fasta_path}"
            fi
        done \
        | tee "${METRICS_DIR}/missing_or_empty_fastas.tsv" \
        | wc -l
)"

if [ "${MISSING_FASTA_COUNT}" -ne 0 ]; then
    fail "Found ${MISSING_FASTA_COUNT} missing or empty FASTA files"
fi

read -r -a K_VALUE_ARRAY <<< "${K_VALUES}"
read -r -a EVIDENCE_RANK_ARRAY <<< "${EVIDENCE_RANKS}"

BUILD_COMMAND=(
    kmersutra-build-panel
    --genome_config "${TARGET_CONFIG}"
    --out_dir "${WORK_OUT_DIR}"
    --k_values "${K_VALUE_ARRAY[@]}"
    --taxonomy_dir "${TAXONOMY_DIR}"
    --download_taxonomy_if_missing
    --evidence_ranks "${EVIDENCE_RANK_ARRAY[@]}"
    --threads "${KMERSUTRA_THREADS}"
    --global_candidate_evidence
    --sqlite_batch_size "${SQLITE_BATCH_SIZE}"
    --max_per_species_per_k "${MAX_PER_SPECIES_PER_K}"
    --global_source_index_mode "${GLOBAL_SOURCE_INDEX_MODE}"
    --global_index_progress_interval "${GLOBAL_INDEX_PROGRESS_INTERVAL}"
    --marker_selection "${MARKER_SELECTION}"
    --genome_bin_size "${GENOME_BIN_SIZE}"
    --max_per_genome_bin "${MAX_PER_GENOME_BIN}"
    --module_min_gate_unique_kmers "${MODULE_MIN_GATE_UNIQUE_KMERS}"
    --module_min_gate_positive_sequences "${MODULE_MIN_GATE_POSITIVE_SEQUENCES}"
    --module_min_gate_k_values "${MODULE_MIN_GATE_K_VALUES}"
    --module_min_gate_best_k "${MODULE_MIN_GATE_BEST_K}"
    --module_max_gate_records_per_k "${MODULE_MAX_GATE_RECORDS_PER_K}"
    --panel_storage_format "${PANEL_STORAGE_FORMAT}"
    --ram_log_path "${RAM_LOG_TSV}"
    --ram_log_interval_seconds "${RAM_LOG_INTERVAL_SECONDS}"
    --profile
    --verbose
)

if [ "${WRITE_MODULE_MANIFEST}" != "true" ]; then
    BUILD_COMMAND+=(--no_write_module_manifest)
elif [ -n "${MODULE_MANIFEST_DIR}" ]; then
    BUILD_COMMAND+=(--module_manifest_dir "${MODULE_MANIFEST_DIR}")
fi

if [ "${WRITE_MODULE_PARQUET}" = "true" ]; then
    BUILD_COMMAND+=(
        --write_module_parquet
        --module_parquet_dir "${MODULE_PARQUET_DIR}"
        --module_name "${MODULE_NAME}"
    )
fi

printf '%q ' "${BUILD_COMMAND[@]}" > "${COMMAND_TSV}"
printf '\n' >> "${COMMAND_TSV}"

run_command "kmersutra_build_panel_global_candidate" "${BUILD_COMMAND[@]}"

log_info "Running post-build checks"
PANEL_TSV_GZ="${WORK_OUT_DIR}/species_kmer_panel.tsv.gz"
PROFILE_TSV="${WORK_OUT_DIR}/build_profile_timing.tsv"
SUMMARY_TSV="${WORK_OUT_DIR}/target_evidence_build_summary.tsv"
COLLECTION_TSV="${WORK_OUT_DIR}/kmer_collection_summary.tsv"

{
    printf 'check\tstatus\tvalue\n'
    for path in "${PANEL_TSV_GZ}" "${PROFILE_TSV}" "${SUMMARY_TSV}" "${COLLECTION_TSV}" "${RAM_LOG_TSV}"; do
        if [ -s "${path}" ]; then
            printf '%s\texists\t%s\n' "$(basename "${path}")" "${path}"
        else
            printf '%s\tmissing_or_empty\t%s\n' "$(basename "${path}")" "${path}"
        fi
    done
} > "${POST_BUILD_CHECKS_TSV}"

if [ ! -s "${PANEL_TSV_GZ}" ]; then
    fail "Expected panel missing or empty: ${PANEL_TSV_GZ}"
fi

log_info "Summarising evidence ranks"
{
    printf 'n_records\tevidence_rank\tevidence_name\tspecies_name\tk\n'
    zcat "${PANEL_TSV_GZ}" \
        | awk -F '\t' 'NR == 1 {
            for (i = 1; i <= NF; i++) {
                if ($i == "evidence_rank") rank = i
                if ($i == "evidence_name") name = i
                if ($i == "species_name") species = i
                if ($i == "k") kval = i
            }
            next
        }
        {print $rank "\t" $name "\t" $species "\t" $kval}' \
        | sort \
        | uniq -c \
        | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5}'
} > "${METRICS_DIR}/panel_evidence_rank_summary.tsv"

{
    printf 'end_time\t%s\n' "$(date '+%Y-%m-%d %H:%M:%S')"
    printf 'panel_tsv_gz\t%s\n' "${PANEL_TSV_GZ}"
    printf 'profile_tsv\t%s\n' "${PROFILE_TSV}"
    printf 'ram_tsv\t%s\n' "${RAM_LOG_TSV}"
    printf 'summary_tsv\t%s\n' "${SUMMARY_TSV}"
} >> "${RUN_METADATA_TSV}"

log_info "Done"
log_info "Final output directory: ${FINAL_OUT_DIR}"
