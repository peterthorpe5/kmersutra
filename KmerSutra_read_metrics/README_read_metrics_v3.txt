# KmerSutra read metrics helper v3

This helper summarises read length and quality metrics for the KmerSutra comparable benchmark.

The previous helper wrote outputs successfully but the wrapper looked for the wrong group-summary filename. This version writes and prints explicit filenames.

Default mode is `simulated_sources`, because the manuscript question is about the simulated pathogen read quality rather than only the final mixed sample FASTQs. It looks in each manifest `source_run_dir` for files such as:

- `simulated_pathogen_unaligned_reads.fastq.gz`
- `simulated_pathogen_aligned_reads.fastq.gz`
- `sim_pool.fastq.gz`
- `train_reads.fastq.gz`

Outputs:

- `simulated_source_fastq_read_metrics.tsv.gz`
- `simulated_source_group_fastq_read_metrics.tsv`
- `all_fastq_read_metrics.tsv.gz`
- `all_group_fastq_read_metrics.tsv`

All outputs are tab-separated.
