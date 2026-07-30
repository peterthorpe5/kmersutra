# KmerSutra v0.51.0 rebuild handoff

## Outcome

This repository is a backward-compatible rebuild of v0.50.1 for a locked
independent benchmark. The deterministic rule-based caller remains the
scientific core. No existing console command was removed or renamed.

## Verified release state

- 416 main tests pass; one environment-dependent `pigz` test is skipped.
- 13 historical AI-validation tests pass.
- All 27 installed console commands return successful `--help`.
- Ruff static analysis passes.
- Branch coverage is 78%, above the ratcheted 75% release floor.
- Resource warnings are treated as errors in an additional successful test run.
- The complete benchmark controller passes a synthetic end-to-end and resume
  test.
- Source and wheel builds succeed.
- The wheel contains only `kmersutra` and `kmersutra.cli`; it excludes tests,
  cache files and the historical `kmersutra_ai_validation` package.

## ATCC workflow

The locked benchmark is:

```text
ATCC MSA-1003 PacBio HiFi
SRR9328980
k = 51, 77, 101, 151
exact matching
same-genus reportability fraction = 0.05
AI novelty scale = 2.9 without retraining
```

The workflow is restartable across preflight, acquisition, depth subsets,
flat screening, single-k ablations, hierarchical screening, frozen AI
validation, summary and provenance stages. Stage success requires the current
configuration digest and non-empty declared outputs.

## Required site-specific inputs

The repository deliberately does not invent or bundle:

- the held-out ATCC reference panel;
- its `kmersutra_genome_config.tsv`;
- the canonical frozen AI model path;
- institutional licence terms.

Set those paths in `config.cluster.json`. The preflight leakage audit must pass
before the full Slurm job is submitted.

## Handoff commands

Use
[`DEPLOYMENT_COMMANDS_v0_51_0.md`](DEPLOYMENT_COMMANDS_v0_51_0.md) for the exact
Mac overlay, GitHub push, cluster pull/reinstall and Slurm submission commands.
