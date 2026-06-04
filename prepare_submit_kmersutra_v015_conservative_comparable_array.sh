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

require_file() {
    local file_path="$1"
    local label="$2"

    if [[ ! -s "${file_path}" ]]; then
        die "${label} does not exist or is empty: ${file_path}"
    fi
}

SCRIPT_DIR="$(
    cd "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1
    pwd
)"

PROJECT_DIR="${PROJECT_DIR:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity}"
DB_ROOT="${DB_ROOT:-/home/pthorpe001/data/databases/kmersutra_db}"
SOURCE_RUN_ROOT="${SOURCE_RUN_ROOT:-${PROJECT_DIR}/runs_kmersutra_v014_global_comparable_20260514_111550}"
PANEL="${PANEL:-${DB_ROOT}/kmersutra_builds/kmersutra_plasmodium_outgroups_v3_global_candidate_k77_101_20260512_113441/species_kmer_panel.tsv.gz}"
ARRAY_SCRIPT="${ARRAY_SCRIPT:-${SCRIPT_DIR}/run_kmersutra_v015_conservative_array_tmpdir.sh}"
RUN_STAMP="${RUN_STAMP:-$(date '+%Y%m%d_%H%M%S')}"
OUT_ROOT="${OUT_ROOT:-${PROJECT_DIR}/runs_kmersutra_v015_conservative_comparable_${RUN_STAMP}}"
SOURCE_MANIFEST="${SOURCE_MANIFEST:-}"
TC_LIMIT="${TC_LIMIT:-100}"
THREADS="${THREADS:-4}"
CHUNK_SIZE="${CHUNK_SIZE:-10000}"
MAX_PENDING_CHUNKS="${MAX_PENDING_CHUNKS:-}"
CALL_PRESET="${CALL_PRESET:-conservative}"
LOW_EVIDENCE_CALL="${LOW_EVIDENCE_CALL:-observed_below_threshold}"
MIN_UNIQUE_KMER_MARGIN="${MIN_UNIQUE_KMER_MARGIN:-0}"
MIN_UNIQUE_KMER_RATIO="${MIN_UNIQUE_KMER_RATIO:-0.0}"
H_RT="${H_RT:-}"
H_VMEM="${H_VMEM:-120G}"
M_MEM_FREE="${M_MEM_FREE:-120G}"
JOB_NAME="${JOB_NAME:-KSscreen_v015}"
CONSOLIDATE_SPECIES_CALLS="${CONSOLIDATE_SPECIES_CALLS:-false}"
BACKGROUND_CANDIDATE_TAXA="${BACKGROUND_CANDIDATE_TAXA:-}"
BACKGROUND_CANDIDATE_FILE="${BACKGROUND_CANDIDATE_FILE:-}"
DISABLE_SAME_GENUS_NEIGHBOUR_DEMOTION="${DISABLE_SAME_GENUS_NEIGHBOUR_DEMOTION:-false}"
DOMINANT_SPECIES_MIN_MARGIN="${DOMINANT_SPECIES_MIN_MARGIN:-25}"
DOMINANT_SPECIES_MIN_RATIO="${DOMINANT_SPECIES_MIN_RATIO:-2.0}"
WRITE_PARQUET_OUTPUTS="${WRITE_PARQUET_OUTPUTS:-false}"

# Optional threshold overrides for kmersutra-screen. Leave empty to use the preset.
MIN_UNIQUE_KMERS="${MIN_UNIQUE_KMERS:-}"
MIN_POSITIVE_SEQUENCES="${MIN_POSITIVE_SEQUENCES:-}"
MIN_K_VALUES_POSITIVE="${MIN_K_VALUES_POSITIVE:-}"
MAX_CONFLICT_RATIO="${MAX_CONFLICT_RATIO:-}"
MIN_BEST_K="${MIN_BEST_K:-}"
MIN_EXACT_HITS="${MIN_EXACT_HITS:-}"
MIN_CONFIDENCE_SCORE="${MIN_CONFIDENCE_SCORE:-}"

require_file "${PANEL}" "KmerSutra panel"
require_file "${ARRAY_SCRIPT}" "Array worker script"

if [[ -z "${SOURCE_MANIFEST}" ]]; then
    for candidate in \
        "${SOURCE_RUN_ROOT}/kmersutra_v014_comparable_manifest.tsv" \
        "${SOURCE_RUN_ROOT}/kmersutra_comparable_manifest.tsv" \
        "${SOURCE_RUN_ROOT}/summary/sample_status.tsv" \
        "${SOURCE_RUN_ROOT}/sample_status.tsv" \
        "${SOURCE_RUN_ROOT}/screen_status_summary.tsv" \
        "${SOURCE_RUN_ROOT}/summary/screen_status_summary.tsv"
    do
        if [[ -s "${candidate}" ]]; then
            SOURCE_MANIFEST="${candidate}"
            break
        fi
    done
fi

require_file "${SOURCE_MANIFEST}" "Source manifest"

mkdir -p "${OUT_ROOT}/inputs" \
         "${OUT_ROOT}/logs" \
         "${OUT_ROOT}/metrics/task_status" \
         "${OUT_ROOT}/samples"

MANIFEST="${OUT_ROOT}/kmersutra_v015_conservative_manifest.tsv"
COMPAT_MANIFEST="${OUT_ROOT}/kmersutra_comparable_manifest.tsv"
MISSING_FASTQS="${OUT_ROOT}/metrics/missing_or_empty_fastqs.tsv"

log "Project dir: ${PROJECT_DIR}"
log "Source run root: ${SOURCE_RUN_ROOT}"
log "Source manifest: ${SOURCE_MANIFEST}"
log "Panel: ${PANEL}"
log "Output root: ${OUT_ROOT}"
log "Array script: ${ARRAY_SCRIPT}"
log "Throttle: ${TC_LIMIT} concurrent array task(s)"
log "Threads per task: ${THREADS}"
log "Call preset: ${CALL_PRESET}; low evidence call: ${LOW_EVIDENCE_CALL}"
log "Consolidate species calls: ${CONSOLIDATE_SPECIES_CALLS}; background candidate taxa: ${BACKGROUND_CANDIDATE_TAXA:-none}"
log "h_rt request: ${H_RT:-queue_default}"

awk -F '\t' '
BEGIN {
    OFS = "\t"
}
NR == 1 {
    for (i = 1; i <= NF; i++) {
        header[$i] = i
    }
    required["sample_id"] = 1
    required["input_fastq"] = 1
    required["benchmark_family"] = 1
    required["panel"] = 1
    required["replicate"] = 1
    missing = ""
    for (name in required) {
        if (!(name in header)) {
            missing = missing " " name
        }
    }
    if (!("spike_reads" in header) && !("spike_n" in header)) {
        missing = missing " spike_reads_or_spike_n"
    }
    if (missing != "") {
        print "Missing required manifest columns:" missing > "/dev/stderr"
        exit 2
    }
    print "sample_id", "input_fastq", "benchmark_family", "panel", \
          "replicate", "spike_reads", "source_run_dir", "source_relative_dir"
    next
}
{
    sample_id = $(header["sample_id"])
    input_fastq = $(header["input_fastq"])
    benchmark_family = $(header["benchmark_family"])
    panel = $(header["panel"])
    replicate = $(header["replicate"])
    if ("spike_reads" in header) {
        spike_reads = $(header["spike_reads"])
    } else {
        spike_reads = $(header["spike_n"])
    }
    if ("source_run_dir" in header) {
        source_run_dir = $(header["source_run_dir"])
    } else {
        source_run_dir = "NA"
    }
    if ("source_relative_dir" in header) {
        source_relative_dir = $(header["source_relative_dir"])
    } else {
        source_relative_dir = "NA"
    }

    gsub(/^kmersutra_v014/, "kmersutra_v015", sample_id)
    gsub(/kmersutra_v014/, "kmersutra_v015", sample_id)
    gsub(/^\/gpfs\/uod-scale-01\/cluster\/gjb_lab\/pthorpe001\//, "/home/pthorpe001/data/", input_fastq)
    gsub(/^\/gpfs\/uod-scale-01\/cluster\/gjb_lab\/pthorpe001\//, "/home/pthorpe001/data/", source_run_dir)

    if (sample_id == "" || input_fastq == "") {
        next
    }
    print sample_id, input_fastq, benchmark_family, panel, replicate, \
          spike_reads, source_run_dir, source_relative_dir
}' "${SOURCE_MANIFEST}" > "${MANIFEST}"

cp "${MANIFEST}" "${COMPAT_MANIFEST}"
cp "${SOURCE_MANIFEST}" "${OUT_ROOT}/inputs/$(basename "${SOURCE_MANIFEST}")"

: > "${MISSING_FASTQS}"
awk -F '\t' 'NR > 1 {print $2}' "${MANIFEST}" | while IFS= read -r fastq_path
 do
    if [[ ! -s "${fastq_path}" ]]; then
        printf '%s\n' "${fastq_path}" >> "${MISSING_FASTQS}"
    fi
done

if [[ -s "${MISSING_FASTQS}" ]]; then
    warn "Some FASTQs are missing or empty. See: ${MISSING_FASTQS}"
    head "${MISSING_FASTQS}" >&2 || true
    die "Refusing to submit until missing FASTQs are resolved."
fi

N_TASKS=$(( $(wc -l < "${MANIFEST}") - 1 ))
if [[ "${N_TASKS}" -le 0 ]]; then
    die "Manifest contains no samples: ${MANIFEST}"
fi

cat > "${OUT_ROOT}/run_metadata.tsv" <<EOF
field\tvalue
run_stamp\t${RUN_STAMP}
project_dir\t${PROJECT_DIR}
source_run_root\t${SOURCE_RUN_ROOT}
source_manifest\t${SOURCE_MANIFEST}
manifest\t${MANIFEST}
panel\t${PANEL}
out_root\t${OUT_ROOT}
array_script\t${ARRAY_SCRIPT}
call_preset\t${CALL_PRESET}
low_evidence_call\t${LOW_EVIDENCE_CALL}
min_unique_kmer_margin\t${MIN_UNIQUE_KMER_MARGIN}
min_unique_kmer_ratio\t${MIN_UNIQUE_KMER_RATIO}
threads_per_task\t${THREADS}
chunk_size\t${CHUNK_SIZE}
tc_limit\t${TC_LIMIT}
n_tasks\t${N_TASKS}
EOF

log "Prepared manifest: ${MANIFEST}"
log "Number of samples: ${N_TASKS}"
log "Submitting SGE array with -tc ${TC_LIMIT}"

export MANIFEST OUT_ROOT PANEL THREADS CHUNK_SIZE MAX_PENDING_CHUNKS
export CALL_PRESET LOW_EVIDENCE_CALL MIN_UNIQUE_KMER_MARGIN MIN_UNIQUE_KMER_RATIO
export MIN_UNIQUE_KMERS MIN_POSITIVE_SEQUENCES MIN_K_VALUES_POSITIVE
export MAX_CONFLICT_RATIO MIN_BEST_K MIN_EXACT_HITS MIN_CONFIDENCE_SCORE
export CONSOLIDATE_SPECIES_CALLS BACKGROUND_CANDIDATE_TAXA BACKGROUND_CANDIDATE_FILE
export DISABLE_SAME_GENUS_NEIGHBOUR_DEMOTION DOMINANT_SPECIES_MIN_MARGIN
export DOMINANT_SPECIES_MIN_RATIO WRITE_PARQUET_OUTPUTS

QSUB_ARGS=(
    -V
    -cwd
    -j y
    -N "${JOB_NAME}"
    -pe smp "${THREADS}"
    -l "h_vmem=${H_VMEM}"
    -l "m_mem_free=${M_MEM_FREE}"
    -t "1-${N_TASKS}"
    -tc "${TC_LIMIT}"
    -o "${OUT_ROOT}/logs"
)
if [[ -n "${H_RT}" ]]; then
    QSUB_ARGS+=(-l "h_rt=${H_RT}")
fi

qsub "${QSUB_ARGS[@]}" "${ARRAY_SCRIPT}"

log "Submitted KmerSutra v0.15 conservative comparable array"
log "Output root: ${OUT_ROOT}"
log "After completion, rerun your KmerSutra comparable summary with OUT_ROOT=${OUT_ROOT}"
