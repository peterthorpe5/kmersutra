#!/usr/bin/env bash
set -euo pipefail

printf 'KmerSutra v0.32.1 repair validation
'
printf 'Working directory: %s
' "$(pwd)"

if [[ ! -f summarise_kmersutra_comparable_benchmark.py ]]; then
    printf 'ERROR: summarise_kmersutra_comparable_benchmark.py is missing from the repository root.
' >&2
    exit 1
fi

if ! grep -q 'independent_multik_genome_spread' tests/test_cli.py; then
    printf 'ERROR: tests/test_cli.py does not contain the v0.32 default marker-selection expectation.
' >&2
    exit 1
fi

python -m unittest discover -s tests -v
if command -v nose2 >/dev/null 2>&1; then
    nose2 -v
else
    printf 'WARNING: nose2 is not installed in this environment; unittest completed.
' >&2
fi
