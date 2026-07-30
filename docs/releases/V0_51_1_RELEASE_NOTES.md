# KmerSutra v0.51.1 release notes

Date: 30 July 2026

## Purpose

Version 0.51.1 repairs the operational gaps found during the first ATCC
preflight attempt and restores the public Zymo benchmark workflow to a visible,
supported location. It also supplies the publication-facing documentation and
test-result organisation needed for a maintainable release.

## Benchmark corrections

- `kmersutra-run-benchmark` accepts `--database-root` and `--ai-model`.
- The ATCC Slurm submitter requires both paths and validates them before
  submission.
- Configuration loading fails explicitly when `${KMERSUTRA_DB_ROOT}`,
  `${KMERSUTRA_AI_MODEL}` or another placeholder remains unresolved.
- The portable example JSON is used directly; private cluster paths are not
  committed.
- A Slurm-native `ERR5396170` Zymo D6300 screen is supplied under
  `benchmarks/zymo_d6300_ont`.

The historical Zymo qsub wrapper was stored beneath the dataset's local
`${BENCH}/scripts` directory and was never committed to GitHub. The new wrapper
replaces that operational role without pretending to reproduce an unversioned
file byte for byte.

## Documentation

The release contains Sphinx source configured with the Read the Docs theme,
an offline HTML build, a Read the Docs configuration and an automated
documentation build. The guide covers the complete user journey from
installation to benchmark interpretation.

## Test and audit outputs

`run_tests.sh` and `run_quality_checks.sh` accept `--results-dir` and
`--run-label`. They write dated logs, branch coverage, HTML coverage, package
build records, wheel-content checks and package inventories beneath one result
directory.

New package inventories are tab-separated and live under
`package_inventories/`; they are no longer placed in the repository root.

The main suite now runs 465 tests and whole-package branch coverage is 80%, up
from 78%. The enforced coverage floor is 80%; 90% remains the next development
target.

## Historical scripts

Former root-level Grid Engine wrappers are retained under
`scripts/legacy/sge/`. Historical v0.50-v0.50.1 AI-validation wrappers are
under `scripts/legacy/ai_validation/`. Current Slurm launchers live with the
ATCC and Zymo benchmark specifications under `benchmarks/`.

## Compatibility

All v0.50.1 installed console commands remain available. Existing environment
variables continue to work for benchmark configuration, while explicit named
path options provide a safer route for cluster submission.
