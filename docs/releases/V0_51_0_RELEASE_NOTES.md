# KmerSutra v0.51.0 release notes

Release date: 30 July 2026

## Purpose

Version 0.51.0 prepares KmerSutra for a locked independent benchmark while
preserving the full v0.50.1 command surface.

## Scientific additions

- Generic mock-community truth manifests replace benchmark-specific truth code
  for new validations.
- Exact truth assembly accessions, truth sequence accessions in FASTA headers
  and exact FASTA checksums can be audited before screening.
- NCBI downloads can exclude truth accessions at candidate-assembly and
  post-download FASTA-header stages.
- Deterministic depth subsets are produced in one input pass.
- The ATCC MSA-1003 HiFi benchmark is registered before results are examined.

The primary ATCC settings are exact matching, k = 51/77/101/151, same-genus
reportability fraction 0.05 and AI novelty scale 2.9. The existing AI model is
applied without retraining. Primary threshold tuning is prohibited.

## Workflow additions

`kmersutra-run-benchmark` provides:

- named command-line options;
- staged preflight, acquisition, subsampling, screening, ablation, hierarchical,
  AI, summary and provenance steps;
- configuration digests;
- durable success/failed state manifests;
- output validation before a stage is skipped;
- `--resume`, `--start-at`, `--stop-after` and repeatable `--force-stage`;
- per-screen atomic publication;
- full input checksums and runtime metadata.

The Slurm launcher is:

```text
benchmarks/atcc_msa1003_hifi/submit_atcc_msa1003_benchmark.sh
```

## Release engineering

- Root development artefacts were removed or archived.
- The README is now publication-facing.
- GitHub Actions tests Python 3.10, 3.11 and 3.12.
- Coverage, Ruff, wheel construction and wheel-content checks are configured.
- Full branch coverage increased to 78%, with a 75% ratcheted CI floor and a
  documented 90% target.
- `CITATION.cff` is supplied.
- The built wheel excludes historical AI validation scripts and tests.

## Compatibility

All v0.50.1 console commands remain installed. The Zymo-specific command
`kmersutra-label-zymo-calls` remains available for reproduction of the existing
public validation. New mock communities should use
`kmersutra-label-mock-calls`.

## Manual release decisions

An institutional software licence and permanent archive DOI must still be
selected before the first public v0.51.0 release.
