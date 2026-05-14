#!/usr/bin/env bash
#$ -j y
#$ -cwd
#$ -V
#$ -pe smp 2
#$ -jc long
#$ -N KSsummary_v014

set -euo pipefail

log_info() {
    printf '%s INFO  %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" >&2
}

log_error() {
    printf '%s ERROR %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" >&2
}

fail() {
    log_error "$*"
    exit 1
}

require_file() {
    [[ -s "$1" ]] || fail "Required file missing or empty: $1"
}

require_dir() {
    [[ -d "$1" ]] || fail "Required directory missing: $1"
}

PROJECT_DIR="${PROJECT_DIR:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity}"
REPO_DIR="${REPO_DIR:-${PROJECT_DIR}/PT_nanopore_spike_in_pathogen_detection}"
SUMMARY_SCRIPT="${SUMMARY_SCRIPT:-${PROJECT_DIR}/summarise_kmersutra_comparable_benchmark.py}"
OUT_ROOT="${OUT_ROOT:?OUT_ROOT is required}"
OUT_DIR="${OUT_DIR:-${OUT_ROOT}/summary}"
PANEL2_TSV="${PANEL2_TSV:-${REPO_DIR}/configs/pathogen_panel_2.tsv}"
PANEL3_TSV="${PANEL3_TSV:-${REPO_DIR}/configs/pathogen_panel_3.tsv}"
PANEL1_TARGET="${PANEL1_TARGET:-Plasmodium vivax}"
ALLOW_PARTIAL="${ALLOW_PARTIAL:-true}"
STRICT="${STRICT:-false}"
VERBOSE="${VERBOSE:-true}"

require_dir "${PROJECT_DIR}"
require_dir "${OUT_ROOT}"
require_file "${SUMMARY_SCRIPT}"
require_file "${OUT_ROOT}/kmersutra_v014_comparable_manifest.tsv"

if [[ -f "${PANEL2_TSV}" ]]; then
    PANEL2_ARGS=(--panel2_tsv "${PANEL2_TSV}")
else
    PANEL2_ARGS=()
    log_info "Panel 2 TSV not found; proceeding without it: ${PANEL2_TSV}"
fi

if [[ -f "${PANEL3_TSV}" ]]; then
    PANEL3_ARGS=(--panel3_tsv "${PANEL3_TSV}")
else
    PANEL3_ARGS=()
    log_info "Panel 3 TSV not found; proceeding without it: ${PANEL3_TSV}"
fi

PARTIAL_ARGS=()
if [[ "${ALLOW_PARTIAL}" == "true" ]]; then
    PARTIAL_ARGS+=(--allow_partial)
fi
if [[ "${STRICT}" == "true" ]]; then
    PARTIAL_ARGS+=(--strict)
fi
if [[ "${VERBOSE}" == "true" ]]; then
    PARTIAL_ARGS+=(--verbose)
fi

log_info "Project directory: ${PROJECT_DIR}"
log_info "Output root: ${OUT_ROOT}"
log_info "Summary output: ${OUT_DIR}"
log_info "Summary script: ${SUMMARY_SCRIPT}"
log_info "Panel 1 target: ${PANEL1_TARGET}"
log_info "Panel 2 TSV: ${PANEL2_TSV}"
log_info "Panel 3 TSV: ${PANEL3_TSV}"

python3 "${SUMMARY_SCRIPT}" \
    --out_root "${OUT_ROOT}" \
    --out_dir "${OUT_DIR}" \
    --panel1_targets "${PANEL1_TARGET}" \
    "${PANEL2_ARGS[@]}" \
    "${PANEL3_ARGS[@]}" \
    "${PARTIAL_ARGS[@]}"

log_info "KmerSutra comparable summary complete"
log_info "Summary directory: ${OUT_DIR}"
