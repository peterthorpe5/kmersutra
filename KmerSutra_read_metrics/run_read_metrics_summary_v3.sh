#!/usr/bin/env bash
set -eo pipefail

PROJECT_DIR="${PROJECT_DIR:-/home/${USER}/data/2026_plasmodium_kraken_sensitivity}"
RUN_ROOT="${RUN_ROOT:-${PROJECT_DIR}/runs_kmersutra_v042_locus_balanced_sensitive_k1_u10_p2_e10_20260605_053947}"
MANIFEST="${MANIFEST:-${RUN_ROOT}/kmersutra_comparable_manifest.tsv}"
OUT_DIR="${OUT_DIR:-${RUN_ROOT}/read_metrics}"
MAX_READS_PER_FASTQ="${MAX_READS_PER_FASTQ:-100000}"

mkdir -p "${OUT_DIR}"

python "${PROJECT_DIR}/scripts/summarise_fastq_read_metrics.py" \
  --manifest "${MANIFEST}" \
  --out_dir "${OUT_DIR}" \
  --mode simulated_sources \
  --max_reads_per_fastq "${MAX_READS_PER_FASTQ}" \
  --group_columns benchmark_family panel fastq_file_type \
  --verbose

echo
echo "Wrote read metrics to:"
echo "  ${OUT_DIR}"
echo
echo "Created files:"
find "${OUT_DIR}" -maxdepth 1 -type f -printf "  %f\n" | sort

echo
echo "First grouped simulated-source rows:"
if [ -s "${OUT_DIR}/simulated_source_group_fastq_read_metrics.tsv" ]; then
    column -t -s $'\t' "${OUT_DIR}/simulated_source_group_fastq_read_metrics.tsv" | head -n 30
else
    echo "No simulated source grouped summary found."
fi
