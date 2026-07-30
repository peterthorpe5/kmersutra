# Reproducibility model

## Stage state

A benchmark stage is skipped only when its JSON control manifest:

- records `status=success`;
- has the current canonical configuration digest;
- declares at least one output;
- points to outputs that still exist;
- points to non-empty files where a file is expected.

File existence alone is not treated as evidence of completion.

## Partial and failed work

The controller writes `running`, `failed` or `success` stage state. Individual
screening tasks publish through a temporary directory and are renamed into
place only after `species_detection_calls.tsv` has been validated.

Partial task outputs are preserved outside successful task paths for
diagnostics, while subsequent resume runs can safely rerun them.

## Frozen scientific settings

The registered ATCC benchmark rejects changes to:

- SRR9328980;
- exact matching;
- k = 51, 77, 101 and 151;
- maximum mismatches = 0;
- same-genus reportability fraction = 0.05;
- AI novelty scale = 2.9;
- prohibition of primary threshold tuning.

An exploratory analysis must use a distinct benchmark identifier and be
reported separately from the locked primary analysis.

## Provenance

The final stage records:

- SHA-256 digests and sizes for reads, truth, panel, genome configuration and
  frozen AI model;
- KmerSutra and Python versions;
- executable and platform information;
- hostname and Slurm job identifiers;
- per-stage success state.

Screening itself writes timing and resource tables. The Slurm wrapper also
records `/usr/bin/time -v` and `sacct` output when those tools are available.
