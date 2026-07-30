# Changelog

All notable user-facing changes are recorded here. Historical release notes are
retained under [`docs/history/release_notes`](docs/history/release_notes).

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
