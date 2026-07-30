Output guide
============

Primary screening tables
------------------------

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - File
     - Purpose
   * - ``species_detection_calls.tsv``
     - Final per-species report-layer calls and thresholds.
   * - ``sample_species_kmer_evidence.tsv``
     - Aggregated raw species marker evidence by k value.
   * - ``sample_species_kmer_hits.tsv``
     - Marker-hit detail used to support aggregated evidence.
   * - ``sample_taxonomic_kmer_evidence.tsv``
     - Evidence retained at the most specific supported taxonomic rank.
   * - ``sample_lineage_interpretation.tsv``
     - Unresolved, near-neighbour and broader-lineage interpretation.
   * - ``module_activation.tsv``
     - Gate decisions for hierarchical screening, when populated.
   * - ``profile_timing.tsv``
     - Observed timings for major steps.

TSV is the portable default. Large tables may also be emitted as TSV.GZ or
Parquet. KmerSutra does not use comma-separated result files.

Interpretation order
--------------------

Start with ``species_detection_calls.tsv`` for the reported conclusion, then
trace any surprising call through species evidence, taxonomic evidence and the
lineage table. Do not infer absence solely from a species not being promoted to
``reportable``.
