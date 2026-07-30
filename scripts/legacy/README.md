# Historical cluster wrappers

These scripts are retained for provenance and exact reconstruction of older
KmerSutra analyses. They are not the recommended launchers for new work.

- `sge/` contains the former root-level Grid Engine (`qsub`) and validation
  wrappers from the v0.15-v0.49 development series.
- `ai_validation/` contains the v0.50-v0.50.1 Grid Engine AI-validation
  wrappers.

Current restartable benchmark launchers live with their specifications:

- ATCC MSA-1003 HiFi:
  `benchmarks/atcc_msa1003_hifi/submit_atcc_msa1003_benchmark.sh`
- ZymoBIOMICS D6300 ONT:
  `benchmarks/zymo_d6300_ont/submit_zymo_d6300_screen.sh`

The relocation does not change installed KmerSutra commands. If an older
notebook or job record cites a former root-level script, use the same filename
under `scripts/legacy/sge/`. The former `scripts/qsub_*_v050*.sh` wrappers use
the same filename under `scripts/legacy/ai_validation/`.

The legacy scripts may contain historical cluster paths and scheduler
directives. Review them before reuse. New analysis should use the maintained
Slurm wrappers and named command-line options.
