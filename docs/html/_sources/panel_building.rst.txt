Building reference panels
=========================

Genome configuration
--------------------

Panel construction is driven by a tab-separated genome configuration. Each row
identifies a FASTA, species label, taxid, role and any clade or module metadata.
Use absolute FASTA paths on the cluster and retain assembly accessions for
provenance and leakage auditing.

Typical build
-------------

.. code-block:: bash

   kmersutra-build-panel \
       --genome_config genome_config.tsv \
       --out_dir panel_build \
       --k_values 51 77 101 151 \
       --global_candidate_evidence \
       --global_source_index_mode candidate_universe \
       --marker_selection independent_multik_genome_spread \
       --min_cross_k_marker_distance 5000 \
       --threads 8 \
       --profile \
       --verbose

The global candidate universe is important when modules or independently built
panels will later be merged. A marker that appears unique inside one small
module may be shared elsewhere in the full database.

Reference leakage
-----------------

For an external benchmark, exclude the truth assemblies before building the
panel and run:

.. code-block:: bash

   kmersutra-audit-reference-panel \
       --genome_config heldout_genome_config.tsv \
       --truth_manifest mock_truth.tsv \
       --out_table reference_panel_audit.tsv \
       --fail_on_leakage \
       --verbose

The audit checks configured assembly accessions, truth sequence accessions in
FASTA headers and exact FASTA SHA-256 values when truth FASTAs are supplied.

Validation
----------

Use ``kmersutra-validate-panel`` after construction and before screening. Keep
the build configuration, logs, marker summaries, timings and checksums with the
panel.
