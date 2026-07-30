ONT ZymoBIOMICS D6300
=====================

The public ONT Q20 Zymo D6300 dataset is ``ERR5396170``. The historical
same-species reference-label wrapper lived inside the dataset workspace rather
than the Git repository:

``qsub_screen_zymo_same_species_reference_labels_v046_exact_minmixed_tmpdir.sh``

KmerSutra v0.51.1 supplies a Slurm replacement:

.. code-block:: bash

   bash benchmarks/zymo_d6300_ont/submit_zymo_d6300_screen.sh \
       --reads /absolute/path/to/ERR5396170.fastq.gz \
       --panel /absolute/path/to/species_kmer_panel.tsv.gz \
       --output-root /absolute/path/to/zymo_d6300_screen \
       --min-mixed 0.05 \
       --account barton \
       --partition barton \
       --threads 24 \
       --conda-env kmersutra

The ``0.05`` setting is the manuscript-facing operating point from the
same-species reference-label challenge. It recovered all ten expected
organisms with no reportable non-expected species. Other threshold points were
exploratory and must be labelled as such.

Use ``kmersutra-label-zymo-calls`` and the existing benchmark-summary commands
to separate expected organisms, same-species alternatives and non-expected
species in the resulting report.
