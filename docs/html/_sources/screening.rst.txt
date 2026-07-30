Screening reads
===============

Exact and rescue modes
----------------------

``--screen_preset exact`` is the primary validated mode. It searches exact
marker matches and was used for the locked ATCC primary benchmark.

``--screen_preset raw_ont_sensitive`` adds a candidate-restricted rescue path
for noisy ONT data. Use it only where that operating point has been
pre-specified or is clearly labelled exploratory.

Flat and hierarchical modes
---------------------------

Flat mode screens one complete panel. Hierarchical mode first evaluates broad
gate panels and opens detailed modules when configured evidence thresholds are
met. If no module manifest is supplied, hierarchical mode falls back to the
single-panel compatibility behaviour.

Report-layer controls
---------------------

``--call_preset``
   Selects a coherent threshold set. ``lineage_aware`` is a useful conservative
   default for exact screening.

``--min_mixed_species_fraction``
   Requires each reportable organism in a mixed community to retain a fraction
   of the strongest species support.

``--same_genus_reportable_min_fraction``
   Controls whether a weaker same-genus species remains reportable or is
   retained as neighbour-lineage evidence.

``--consolidate_species_calls``
   Applies same-genus dominance and background-candidate interpretation.

Performance
-----------

For large FASTQ files, use ``--decompressor auto``, panel caching,
``--no_read_level_hits`` and a practical chunk size. Record observed wall time
and peak memory; do not report requested Slurm memory as tool memory.
