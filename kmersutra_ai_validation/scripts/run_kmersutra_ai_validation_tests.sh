#!/bin/bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AI_DIR="${SCRIPT_DIR}/ai_python"

export PYTHONPATH="${AI_DIR}:${PYTHONPATH:-}"
python -m unittest discover \
    -s "${AI_DIR}/tests" \
    -p 'test_*.py' \
    -v
