# ONT ZymoBIOMICS D6300 benchmark

This directory restores a visible, current Slurm entry point for the public
ONT Q20 ZymoBIOMICS D6300 benchmark (`ERR5396170`).

The historical cluster workflow used:

```text
qsub_screen_zymo_same_species_reference_labels_v046_exact_minmixed_tmpdir.sh
```

That file lived beneath the dataset-specific `${BENCH}/scripts/` directory and
was never committed to the KmerSutra Git repository. It was therefore absent
from the v0.51.0 repository overlay even though the KmerSutra commands and
Zymo truth/reporting code remained present.

The replacement is Slurm-native and takes named paths:

```bash
bash benchmarks/zymo_d6300_ont/submit_zymo_d6300_screen.sh \
    --reads /absolute/path/to/ERR5396170.fastq.gz \
    --panel /absolute/path/to/species_kmer_panel.tsv.gz \
    --output-root /absolute/path/to/zymo_d6300_screen \
    --min-mixed 0.05 \
    --account barton \
    --partition barton \
    --threads 24 \
    --conda-env kmersutra
```

The `0.05` setting reproduces the manuscript-facing report-layer operating
point that recovered all ten expected organisms with no reportable
non-expected species in the same-species reference-label challenge.

The wrapper records the complete command, input paths, Slurm allocation,
observed `/usr/bin/time` metrics and `sacct` metrics when available. It writes
TSV, TSV.GZ and Parquet outputs through `kmersutra-screen`.

This benchmark is historical evidence and may be rerun for reproducibility. It
must not be used to retune the locked ATCC MSA-1003 primary analysis.
