#!/usr/bin/env bash

set -Eeuo pipefail
trap 'echo "ERROR at line ${LINENO}: ${BASH_COMMAND}" >&2' ERR

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${REPO_ROOT}"

python -m ruff check kmersutra tests scripts/check_wheel_contents.py
python -m coverage erase
python -m coverage run -m unittest discover -s tests -p 'test_*.py'
python -m coverage report
bash kmersutra_ai_validation/scripts/run_kmersutra_ai_validation_tests.sh

rm -rf build dist ./*.egg-info
python -m build
python scripts/check_wheel_contents.py --dist-dir dist
git diff --check
