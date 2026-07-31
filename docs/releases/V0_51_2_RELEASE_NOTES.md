# KmerSutra v0.51.2 release notes

Release date: 31 July 2026

## Purpose

Version 0.51.2 repairs the missing reference-preparation stage for the locked
ATCC MSA-1003 PacBio HiFi benchmark. The earlier v4/v0.46 panel is a valid
Plasmodium/outgroup build, but it represents none of the 20 ATCC targets. The
scientific gate correctly blocked screening rather than producing an
all-negative, uninterpretable result.

## Reference design

The new workflow:

- requests alternative NCBI assemblies for all 20 expected ATCC species;
- samples additional species from every distinct target genus as absent near
  neighbours;
- excludes the published truth sequence and assembly accessions during
  acquisition and repeats the audit before panel construction;
- removes duplicate target rows returned through genus queries;
- keeps the established v4 genome collection as background evidence;
- builds k = 51, 77, 101 and 151 markers with the locked v0.46
  `raw_ont_lod_balanced` profile;
- publishes flat and hierarchical panels only after the reference gate passes.

## Execution

The reference preparation is split into restartable stages with configuration
signatures and non-empty output validation. It can be resubmitted with
`--resume` after the cluster time limit. Failed partial panel builds are moved
under `partial_builds` rather than treated as complete.

The complete command sequence is in
[`DEPLOYMENT_COMMANDS_v0_51_2.md`](../DEPLOYMENT_COMMANDS_v0_51_2.md).

## Validation

- 469 repository tests pass; four optional dependency tests are skipped.
- New leakage-control, coverage, deduplication and taxonomic-plan tests pass.
- All retained shell and Slurm files pass `bash -n` syntax validation.
