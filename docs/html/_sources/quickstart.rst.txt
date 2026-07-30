Quick start
===========

Install
-------

Clone the repository and create the supplied Conda environment:

.. code-block:: bash

   git clone https://github.com/peterthorpe5/kmersutra.git
   cd kmersutra
   conda env create --file environment.yml
   conda run --name kmersutra kmersutra-screen --help

The environment includes the screening, reporting, NCBI, Parquet, development
and documentation dependencies used by the release checks.

Screen an existing panel
------------------------

.. code-block:: bash

   conda run --name kmersutra \
       kmersutra-screen \
       --input sample.fastq.gz \
       --input_format fastq \
       --panel panel/species_kmer_panel.tsv.gz \
       --sample_id sample_001 \
       --out_dir results/sample_001 \
       --screen_mode flat \
       --screen_preset exact \
       --call_preset lineage_aware \
       --consolidate_species_calls \
       --same_genus_reportable_min_fraction 0.05 \
       --threads 8 \
       --decompressor auto \
       --use_panel_cache \
       --no_read_level_hits \
       --write_parquet_outputs \
       --profile \
       --verbose

The principal outputs are:

``species_detection_calls.tsv``
   Report-layer species decisions. Start here when asking what KmerSutra would
   report.

``sample_species_kmer_evidence.tsv``
   Raw species-level marker evidence before report-layer interpretation.

``sample_taxonomic_kmer_evidence.tsv``
   Evidence retained at genus, family or broader ranks when the sequence does
   not justify a species-specific assignment.

``sample_lineage_interpretation.tsv``
   Rank-aware interpretation of unresolved or near-neighbour evidence.

``profile_timing.tsv``
   Observed wall-clock timings for major screening phases.

Next steps
----------

Read :doc:`concepts` before interpreting weak or same-genus evidence. Use
:doc:`panel_building` if you need to build a new reference panel, and
:doc:`hpc` for Slurm execution.
