#!/usr/bin/env bash

set -Eeuo pipefail
trap 'echo "ERROR at line ${LINENO}: ${BASH_COMMAND}" >&2' ERR

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${REPO_ROOT}"

usage() {
    echo "Usage:"
    echo "  bash run_tests.sh [--results-dir DIR] [--run-label LABEL]"
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
mkdir -p "${RUN_DIR}"

set +e
python -m unittest discover -s tests -p 'test_*.py' -v \
    2>&1 | tee "${RUN_DIR}/main_unittest.log"
MAIN_STATUS="${PIPESTATUS[0]}"

bash kmersutra_ai_validation/scripts/run_kmersutra_ai_validation_tests.sh \
    2>&1 | tee "${RUN_DIR}/ai_validation_unittest.log"
AI_STATUS="${PIPESTATUS[0]}"
set -e

{
    printf "suite\texit_status\tlog\n"
    printf "main_unittest\t%s\t%s\n" \
        "${MAIN_STATUS}" "${RUN_DIR}/main_unittest.log"
    printf "ai_validation_unittest\t%s\t%s\n" \
        "${AI_STATUS}" "${RUN_DIR}/ai_validation_unittest.log"
} > "${RUN_DIR}/test_status.tsv"

echo "Test results: ${RUN_DIR}"
if [[ "${MAIN_STATUS}" -ne 0 || "${AI_STATUS}" -ne 0 ]]; then
    exit 1
fi
