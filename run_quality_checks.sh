#!/usr/bin/env bash

set -Eeuo pipefail
trap 'echo "ERROR at line ${LINENO}: ${BASH_COMMAND}" >&2' ERR

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${REPO_ROOT}"

usage() {
    echo "Usage:"
    echo "  bash run_quality_checks.sh [--results-dir DIR] [--run-label LABEL]"
}

RESULTS_DIR="${TEST_RESULTS_ROOT:-${REPO_ROOT}/test_results}"
RUN_LABEL="$(date -u +%Y%m%dT%H%M%SZ)"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --results-dir)
            RESULTS_DIR="${2:?Missing value for --results-dir}"
            shift 2
            ;;
        --run-label)
            RUN_LABEL="${2:?Missing value for --run-label}"
            shift 2
            ;;
        --help|-h)
            usage
            exit 0
            ;;
        *)
            echo "ERROR: unknown option: $1" >&2
            usage >&2
            exit 2
            ;;
    esac
done

RUN_DIR="${RESULTS_DIR%/}/${RUN_LABEL}"
COVERAGE_FILE_PATH="${RUN_DIR}/.coverage"
mkdir -p \
    "${RUN_DIR}/coverage_html" \
    "${RUN_DIR}/package_inventories"

python -m ruff check \
    kmersutra \
    tests \
    scripts/check_wheel_contents.py \
    2>&1 | tee "${RUN_DIR}/ruff.log"

COVERAGE_FILE="${COVERAGE_FILE_PATH}" python -m coverage erase
COVERAGE_FILE="${COVERAGE_FILE_PATH}" python -m coverage run \
    -m unittest discover -s tests -p 'test_*.py' \
    2>&1 | tee "${RUN_DIR}/main_unittest.log"
COVERAGE_FILE="${COVERAGE_FILE_PATH}" python -m coverage report \
    2>&1 | tee "${RUN_DIR}/coverage_report.txt"
COVERAGE_FILE="${COVERAGE_FILE_PATH}" python -m coverage xml \
    -o "${RUN_DIR}/coverage.xml"
COVERAGE_FILE="${COVERAGE_FILE_PATH}" python -m coverage html \
    -d "${RUN_DIR}/coverage_html"

bash kmersutra_ai_validation/scripts/run_kmersutra_ai_validation_tests.sh \
    2>&1 | tee "${RUN_DIR}/ai_validation_unittest.log"

python -m build --outdir "${RUN_DIR}/dist" \
    2>&1 | tee "${RUN_DIR}/package_build.log"
python scripts/check_wheel_contents.py --dist-dir "${RUN_DIR}/dist" \
    2>&1 | tee "${RUN_DIR}/wheel_contents.log"
python -m kmersutra.package_inventory \
    --repo-root "${REPO_ROOT}" \
    --output \
    "${RUN_DIR}/package_inventories/package_file_inventory.tsv"
git diff --check 2>&1 | tee "${RUN_DIR}/git_diff_check.log"

{
    printf "setting\tvalue\n"
    printf "completed_at_utc\t%s\n" "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    printf "results_directory\t%s\n" "${RUN_DIR}"
    printf "coverage_report\t%s\n" "${RUN_DIR}/coverage_report.txt"
    printf "coverage_html\t%s\n" "${RUN_DIR}/coverage_html/index.html"
    printf "package_inventory\t%s\n" \
        "${RUN_DIR}/package_inventories/package_file_inventory.tsv"
} > "${RUN_DIR}/quality_check_summary.tsv"

echo "Quality-check results: ${RUN_DIR}"
