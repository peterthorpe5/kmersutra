#!/usr/bin/env bash

set -Eeuo pipefail
trap 'echo "ERROR at line ${LINENO}: ${BASH_COMMAND}" >&2' ERR

usage() {
    echo "Usage:"
    echo "  bash scripts/build_documentation.sh [--output-dir DIR]"
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
OUTPUT_DIR="${REPO_ROOT}/docs/html"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --output-dir)
            OUTPUT_DIR="${2:?Missing value for --output-dir}"
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

mkdir -p "${OUTPUT_DIR}"
python -m sphinx \
    -W \
    --keep-going \
    -b html \
    "${REPO_ROOT}/docs/source" \
    "${OUTPUT_DIR}"
touch "${OUTPUT_DIR}/.nojekyll"
echo "Documentation: ${OUTPUT_DIR}/index.html"
