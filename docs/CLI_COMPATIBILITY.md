# Command compatibility

Version 0.51.0 retains every console command distributed in v0.50.1:

```text
kmersutra-benchmark-postprocess
kmersutra-build-call-training
kmersutra-build-panel
kmersutra-download-genomes
kmersutra-download-ncbi
kmersutra-download-taxonomy
kmersutra-evaluate-call-predictions
kmersutra-extract-features
kmersutra-label-zymo-calls
kmersutra-merge-modules
kmersutra-merge-panels
kmersutra-predict
kmersutra-screen
kmersutra-summarise-comparable-benchmark
kmersutra-summarise-lca
kmersutra-summarise-run
kmersutra-summarise-spikeins
kmersutra-threshold-sweep
kmersutra-train-call-calibrator
kmersutra-train-classifier
kmersutra-validate-call-calibrator
kmersutra-validate-panel
```

The v0.32 generic table names remain preferred, while historical aliases remain
valid:

```text
--calls_tsv
--out_tsv
--training_tsv
--out_summary_tsv
--out_evaluation_tsv
```

Version 0.51.0 adds:

```text
kmersutra-audit-reference-panel
kmersutra-label-mock-calls
kmersutra-run-benchmark
kmersutra-subset-reads
kmersutra-summarise-mock-benchmark
```

Compatibility is enforced by `tests/test_cli_compatibility_v051.py`.
