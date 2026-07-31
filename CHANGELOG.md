# Changelog

All notable user-facing changes are recorded here. Historical release notes are
retained under [`docs/history/release_notes`](docs/history/release_notes).

## 0.51.2 - 2026-07-31

### Fixed

- Replaced the unusable ATCC placeholder panel path with a complete,
  leakage-controlled reference-preparation workflow.
- Removed the shell-specific AWK reference gate that failed on implementations
  where `index` is reserved; the gate is now implemented and tested in Python.
- Prevented the established Plasmodium/outgroup v4 panel from being treated as
  an ATCC target panel. It remains background evidence in the extended build.

### Added

- `kmersutra-prepare-atcc-reference` commands for taxid planning, genome-config
  finalisation and reference-gate evaluation.
- Restartable Slurm stages for target and genus-neighbour acquisition, truth
  accession exclusion, reference audit, v0.46-profile panel construction and
  atomic completion tokens.
- Unit tests for target/genus planning, role and accession deduplication,
  held-out target coverage, truth-sequence leakage and gate blocking.

### Changed

- Package version increased from 0.51.1 to 0.51.2.
- The ATCC example configuration now enables the hierarchical module manifest
  produced by the reference build.

## 0.51.1 - 2026-07-30

### Fixed

- ATCC benchmark submission now requires explicit `--database-root` and
  `--ai-model` paths, while retaining environment-variable compatibility.
- Unresolved `${VARIABLE}` placeholders fail before stage execution with the
  missing variable names and corrective guidance.
- Removed the invalid `nose2` table that caused third-party `nose2` to fail
  while reading `pyproject.toml`; KmerSutra continues to use `unittest`.
- Restored the public ONT Zymo D6300 screening workflow as a current
  Slurm-native wrapper.

### Added

- Read the Docs-style Sphinx source, offline HTML build and automated
  documentation validation.
- A browsable guide covering installation, concepts, panel construction,
  screening, outputs, Slurm, reproducibility, troubleshooting, ATCC and Zymo.
- Named-output test and quality scripts that collect all generated artefacts
  beneath a dated result directory.
- Deterministic `kmersutra-package-inventory` TSV generation.
- Shell syntax and benchmark-interface regression tests.
- Expanded compressed-I/O, deterministic depth-subset, Parquet and taxonomy
  tests, raising whole-package branch coverage from 78% to 80%.

### Changed

- Package version increased from 0.51.0 to 0.51.1.
- Generated test logs, coverage outputs and package inventories no longer
  belong in the repository root or the user's home directory.
- Historical shell workflows are documented separately from current benchmark
  entry points.
- The enforced branch-coverage floor is ratcheted from 75% to 80%.

## 0.51.0 - 2026-07-30

### Added

- Generic, manifest-driven mock-community truth labelling.
- Reference-panel leakage audit for truth assembly accessions, sequence
  accessions in FASTA headers and exact FASTA checksums.
- Downloader support for explicit accession exclusions.
- Deterministic, one-pass FASTQ/FASTA depth-series generation.
- Restartable benchmark controller with checksummed stage manifests, atomic
  output publication, safe resume and controlled stage reruns.
- Locked Slurm workflow for ATCC MSA-1003 PacBio HiFi dataset SRR9328980.
- Publication-facing README, citation metadata, continuous integration,
  coverage configuration and packaging checks.
- Compatibility regression tests for every v0.50.1 console command.

### Changed

- Package version increased from 0.50.1 to 0.51.0.
- Python package discovery now excludes the separate historical AI validation
  scripts and all tests from built wheels.
- Historical inventories, release notes and test logs moved under
  `docs/history`.

### Compatibility

- Every v0.50.1 console command is retained unchanged.
- Historical option aliases such as `--calls_tsv`, `--out_tsv`,
  `--training_tsv`, `--out_summary_tsv` and `--out_evaluation_tsv` remain
  supported.
- The Zymo-specific truth command remains available. New benchmarks should use
  `kmersutra-label-mock-calls`.

## 0.50.1 - 2026-06-09

- Canonical bounded transformed-feature AI validation release.
- See
  [`docs/history/release_notes/V0_50_1_RELEASE_NOTES.txt`](docs/history/release_notes/V0_50_1_RELEASE_NOTES.txt).
