#!/usr/bin/env bash
#SBATCH --job-name=KS_ATCCref
#SBATCH --account=barton
#SBATCH --partition=barton
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=24
#SBATCH --mem=128G

set -Eeuo pipefail

required_environment=(
    KS_ATCC_REPO
    KS_ATCC_DATABASE_ROOT
    KS_ATCC_REFERENCE_ROOT
    KS_ATCC_BACKGROUND_CONFIG
    KS_ATCC_TRUTH_MANIFEST
    KS_ATCC_TAXONOMY_DIR
    KS_ATCC_NCBI_EMAIL
    KS_ATCC_THREADS
    KS_ATCC_CONDA_ENV
    KS_ATCC_TARGET_ASSEMBLIES
    KS_ATCC_NEIGHBOUR_ASSEMBLIES
    KS_ATCC_MINIMUM_TARGET_REFS
    KS_ATCC_MAX_ENTREZ_RECORDS
)

for variable_name in "${required_environment[@]}"; do
    if [[ -z "${!variable_name:-}" ]]; then
        echo "ERROR: required environment variable is blank: ${variable_name}" >&2
        exit 2
    fi
done

if [[ ! -s "${KS_ATCC_BACKGROUND_CONFIG}" ]]; then
    echo "ERROR: background configuration is missing: ${KS_ATCC_BACKGROUND_CONFIG}" >&2
    exit 1
fi
if [[ ! -s "${KS_ATCC_TRUTH_MANIFEST}" ]]; then
    echo "ERROR: truth manifest is missing: ${KS_ATCC_TRUTH_MANIFEST}" >&2
    exit 1
fi

CONDA_EXE_PATH="$(command -v conda || true)"
if [[ -z "${CONDA_EXE_PATH}" ]]; then
    for candidate in \
        "${HOME}/miniforge3/bin/conda" \
        "${HOME}/mambaforge/bin/conda" \
        "${HOME}/miniconda3/bin/conda" \
        "${HOME}/anaconda3/bin/conda"
    do
        if [[ -x "${candidate}" ]]; then
            CONDA_EXE_PATH="${candidate}"
            break
        fi
    done
fi
if [[ -z "${CONDA_EXE_PATH}" || ! -x "${CONDA_EXE_PATH}" ]]; then
    echo "ERROR: conda executable could not be resolved" >&2
    exit 1
fi

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export MALLOC_ARENA_MAX=2

REFERENCE_ROOT="${KS_ATCC_REFERENCE_ROOT}"
PLAN_DIR="${REFERENCE_ROOT}/reference_plan"
DOWNLOAD_DIR="${REFERENCE_ROOT}/downloads"
AUDIT_DIR="${REFERENCE_ROOT}/reference_audit"
PANEL_BUILD_DIR="${REFERENCE_ROOT}/panel_build"
STATE_DIR="${REFERENCE_ROOT}/workflow_control/stage_tokens"
RUNTIME_DIR="${REFERENCE_ROOT}/slurm_runtime/${SLURM_JOB_ID:-manual}"
TAXID_PLAN="${PLAN_DIR}/atcc_taxid_plan.tsv"
DOWNLOADED_CONFIG="${DOWNLOAD_DIR}/kmersutra_genome_config.tsv"
FINAL_CONFIG="${REFERENCE_ROOT}/kmersutra_genome_config.tsv"
PREBUILD_COVERAGE="${AUDIT_DIR}/prebuild_target_coverage.tsv"
AUDIT_TABLE="${AUDIT_DIR}/reference_panel_leakage_audit.tsv"
FINAL_COVERAGE="${AUDIT_DIR}/atcc_species_coverage.tsv"
GATE_SUMMARY="${AUDIT_DIR}/atcc_reference_gate_summary.tsv"
PANEL="${PANEL_BUILD_DIR}/species_kmer_panel.tsv.gz"
MODULE_MANIFEST="${PANEL_BUILD_DIR}/hierarchical_modules/kmersutra_module_manifest.tsv"
STATUS_TABLE="${RUNTIME_DIR}/job_exit_status.tsv"

mkdir -p "${PLAN_DIR}" "${DOWNLOAD_DIR}" "${AUDIT_DIR}" "${STATE_DIR}" "${RUNTIME_DIR}"

write_status() {
    local status=$?
    {
        printf 'setting\tvalue\n'
        printf 'exit_status\t%s\n' "${status}"
        printf 'finished_at_utc\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
        printf 'hostname\t%s\n' "$(hostname)"
        printf 'slurm_job_id\t%s\n' "${SLURM_JOB_ID:-}"
        printf 'reference_root\t%s\n' "${REFERENCE_ROOT}"
        printf 'background_config\t%s\n' "${KS_ATCC_BACKGROUND_CONFIG}"
        printf 'truth_manifest\t%s\n' "${KS_ATCC_TRUTH_MANIFEST}"
    } > "${STATUS_TABLE}"
    if command -v sacct >/dev/null 2>&1 && [[ -n "${SLURM_JOB_ID:-}" ]]; then
        sacct \
            --jobs "${SLURM_JOB_ID}" \
            --format=JobID,JobName,State,Elapsed,AllocCPUS,MaxRSS,MaxVMSize,ExitCode \
            --parsable2 \
            > "${RUNTIME_DIR}/sacct.tsv" || true
    fi
    exit "${status}"
}
trap write_status EXIT
trap 'echo "ERROR at line ${LINENO}: ${BASH_COMMAND}" >&2' ERR

RUNNER=(
    "${CONDA_EXE_PATH}"
    run
    --no-capture-output
    --name "${KS_ATCC_CONDA_ENV}"
)

SIGNATURE="$({
    sha256sum "${KS_ATCC_TRUTH_MANIFEST}" "${KS_ATCC_BACKGROUND_CONFIG}"
    sha256sum \
        "${KS_ATCC_TAXONOMY_DIR}/nodes.dmp" \
        "${KS_ATCC_TAXONOMY_DIR}/names.dmp" \
        "${KS_ATCC_TAXONOMY_DIR}/merged.dmp" \
        "${KS_ATCC_TAXONOMY_DIR}/delnodes.dmp"
    printf '%s\n' \
        "${KS_ATCC_TARGET_ASSEMBLIES}" \
        "${KS_ATCC_NEIGHBOUR_ASSEMBLIES}" \
        "${KS_ATCC_MINIMUM_TARGET_REFS}" \
        "${KS_ATCC_MAX_ENTREZ_RECORDS}" \
        "0.51.2"
} | sha256sum | cut -d ' ' -f 1)"

stage_forced() {
    local stage="$1"
    local forced_stage
    IFS=',' read -r -a forced_stages <<< "${KS_ATCC_FORCE_STAGES:-}"
    for forced_stage in "${forced_stages[@]}"; do
        if [[ "${forced_stage}" == "${stage}" ]]; then
            return 0
        fi
    done
    return 1
}

token_matches() {
    local token="$1"
    local key
    local value
    [[ -s "${token}" ]] || return 1
    while IFS=$'\t' read -r key value; do
        if [[ "${key}" == "signature" && "${value}" == "${SIGNATURE}" ]]; then
            return 0
        fi
    done < "${token}"
    return 1
}

stage_complete() {
    local stage="$1"
    shift
    local token="${STATE_DIR}/${stage}.token"
    local expected_path
    [[ "${KS_ATCC_RESUME:-0}" == "1" ]] || return 1
    stage_forced "${stage}" && return 1
    token_matches "${token}" || return 1
    for expected_path in "$@"; do
        [[ -s "${expected_path}" ]] || return 1
    done
    return 0
}

write_token() {
    local stage="$1"
    local token="${STATE_DIR}/${stage}.token"
    local temporary_token="${token}.tmp.${SLURM_JOB_ID:-manual}"
    {
        printf 'stage\t%s\n' "${stage}"
        printf 'signature\t%s\n' "${SIGNATURE}"
        printf 'completed_at_utc\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
        printf 'slurm_job_id\t%s\n' "${SLURM_JOB_ID:-}"
    } > "${temporary_token}"
    mv "${temporary_token}" "${token}"
}

run_stage() {
    local stage="$1"
    shift
    local function_name="$1"
    shift
    if stage_complete "${stage}" "$@"; then
        echo "SKIP ${stage}: validated completion token and outputs"
        return 0
    fi
    echo "START ${stage}"
    "${function_name}"
    local expected_path
    for expected_path in "$@"; do
        if [[ ! -s "${expected_path}" ]]; then
            echo "ERROR: ${stage} did not produce ${expected_path}" >&2
            return 1
        fi
    done
    write_token "${stage}"
    echo "COMPLETE ${stage}"
}

stage_plan() {
    "${RUNNER[@]}" kmersutra-prepare-atcc-reference \
        --verbose \
        plan \
        --truth-manifest "${KS_ATCC_TRUTH_MANIFEST}" \
        --taxonomy-dir "${KS_ATCC_TAXONOMY_DIR}" \
        --output "${TAXID_PLAN}" \
        --target-assemblies "${KS_ATCC_TARGET_ASSEMBLIES}" \
        --near-neighbour-assemblies-per-genus "${KS_ATCC_NEIGHBOUR_ASSEMBLIES}"
}

stage_download() {
    "${RUNNER[@]}" kmersutra-download-genomes \
        --taxid_plan "${TAXID_PLAN}" \
        --exclude_accessions_table "${KS_ATCC_TRUTH_MANIFEST}" \
        --out_dir "${DOWNLOAD_DIR}" \
        --email "${KS_ATCC_NCBI_EMAIL}" \
        --source prefer_refseq \
        --formats genomic_fna \
        --decompress \
        --delete_archives \
        --max_entrez_records_per_taxid "${KS_ATCC_MAX_ENTREZ_RECORDS}" \
        --batch_size 100 \
        --retries 5 \
        --sleep_seconds 0.34 \
        --verbose
}

stage_finalise() {
    "${RUNNER[@]}" kmersutra-prepare-atcc-reference \
        --verbose \
        finalise \
        --downloaded-config "${DOWNLOADED_CONFIG}" \
        --background-config "${KS_ATCC_BACKGROUND_CONFIG}" \
        --truth-manifest "${KS_ATCC_TRUTH_MANIFEST}" \
        --taxonomy-dir "${KS_ATCC_TAXONOMY_DIR}" \
        --output-config "${FINAL_CONFIG}" \
        --coverage-table "${PREBUILD_COVERAGE}" \
        --minimum-target-references "${KS_ATCC_MINIMUM_TARGET_REFS}"
}

stage_audit() {
    "${RUNNER[@]}" kmersutra-audit-reference-panel \
        --genome_config "${FINAL_CONFIG}" \
        --truth_manifest "${KS_ATCC_TRUTH_MANIFEST}" \
        --out_table "${AUDIT_TABLE}" \
        --fail_on_leakage \
        --verbose

    "${RUNNER[@]}" kmersutra-prepare-atcc-reference \
        --verbose \
        gate \
        --audit-table "${AUDIT_TABLE}" \
        --truth-manifest "${KS_ATCC_TRUTH_MANIFEST}" \
        --coverage-table "${FINAL_COVERAGE}" \
        --summary-table "${GATE_SUMMARY}" \
        --minimum-target-references "${KS_ATCC_MINIMUM_TARGET_REFS}"
}

stage_build() {
    if [[ -d "${PANEL_BUILD_DIR}" ]]; then
        local partial_dir="${REFERENCE_ROOT}/partial_builds/panel_build.$(date -u +%Y%m%dT%H%M%SZ)"
        mkdir -p "$(dirname "${partial_dir}")"
        mv "${PANEL_BUILD_DIR}" "${partial_dir}"
    fi
    mkdir -p "${PANEL_BUILD_DIR}"

    "${RUNNER[@]}" kmersutra-build-panel \
        --genome_config "${FINAL_CONFIG}" \
        --out_dir "${PANEL_BUILD_DIR}" \
        --k_values 51 77 101 151 \
        --taxonomy_dir "${KS_ATCC_TAXONOMY_DIR}" \
        --evidence_ranks species genus family order class phylum superkingdom \
        --threads "${KS_ATCC_THREADS}" \
        --global_candidate_evidence \
        --marker_profile raw_ont_lod_balanced \
        --candidate_k_order input \
        --global_source_index_mode candidate_universe \
        --global_index_progress_interval 5000000 \
        --sqlite_batch_size 50000 \
        --max_per_species_per_k 100000 \
        --marker_selection independent_multik_genome_spread \
        --min_cross_k_marker_distance 150 \
        --genome_bin_size 10000 \
        --max_per_genome_bin 10 \
        --max_per_genome_bin_by_k 51:4,77:6,101:4,151:1 \
        --assembly_aware_binning \
        --panel_storage_format auto \
        --module_min_gate_unique_kmers 1 \
        --module_min_gate_positive_sequences 1 \
        --module_min_gate_k_values 1 \
        --module_min_gate_best_k 0 \
        --module_max_gate_records_per_k 0 \
        --ram_log_path "${PANEL_BUILD_DIR}/ram_usage.tsv" \
        --ram_log_interval_seconds 30 \
        --profile \
        --verbose
}

stage_publish() {
    ln -sfn "panel_build/species_kmer_panel.tsv.gz" \
        "${REFERENCE_ROOT}/species_kmer_panel.tsv.gz"
    ln -sfn "panel_build/hierarchical_modules/kmersutra_module_manifest.tsv" \
        "${REFERENCE_ROOT}/kmersutra_module_manifest.tsv"
    sha256sum \
        "${FINAL_CONFIG}" \
        "${AUDIT_TABLE}" \
        "${GATE_SUMMARY}" \
        "${PANEL}" \
        "${MODULE_MANIFEST}" \
        > "${REFERENCE_ROOT}/SHA256SUMS.txt"
}

{
    printf 'setting\tvalue\n'
    printf 'started_at_utc\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    printf 'hostname\t%s\n' "$(hostname)"
    printf 'slurm_job_id\t%s\n' "${SLURM_JOB_ID:-}"
    printf 'signature\t%s\n' "${SIGNATURE}"
    printf 'reference_root\t%s\n' "${REFERENCE_ROOT}"
    printf 'threads\t%s\n' "${KS_ATCC_THREADS}"
} > "${RUNTIME_DIR}/job_start.tsv"

"${RUNNER[@]}" python --version
"${RUNNER[@]}" kmersutra-prepare-atcc-reference --help >/dev/null

run_stage "00_plan" stage_plan "${TAXID_PLAN}"
run_stage "01_download" stage_download \
    "${DOWNLOADED_CONFIG}" \
    "${DOWNLOAD_DIR}/ncbi_assembly_candidate_audit.tsv"
run_stage "02_finalise" stage_finalise "${FINAL_CONFIG}" "${PREBUILD_COVERAGE}"
run_stage "03_audit" stage_audit "${AUDIT_TABLE}" "${FINAL_COVERAGE}" "${GATE_SUMMARY}"
run_stage "04_build" stage_build "${PANEL}" "${MODULE_MANIFEST}"
run_stage "05_publish" stage_publish \
    "${REFERENCE_ROOT}/species_kmer_panel.tsv.gz" \
    "${REFERENCE_ROOT}/kmersutra_module_manifest.tsv" \
    "${REFERENCE_ROOT}/SHA256SUMS.txt"

echo "ATCC reference panel completed: ${REFERENCE_ROOT}"
