#!/usr/bin/env bash

set -Eeuo pipefail
trap 'echo "ERROR at line ${LINENO}: ${BASH_COMMAND}" >&2' ERR

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${REPO_ROOT}"

python -m unittest discover -s tests -p 'test_*.py' -v
bash kmersutra_ai_validation/scripts/run_kmersutra_ai_validation_tests.sh
