KmerSutra documentation
=======================

KmerSutra is a conservative, taxonomy-aware, multi-k-mer evidence framework
for species-resolved interpretation of long-read metagenomic data. It was
developed in response to a practical problem: sensitive general-purpose
classifiers can detect low-abundance targets while still producing numerous
plausible-looking off-target species labels. KmerSutra occupies a different
operating point, prioritising auditable evidence and clean reporting.

The deterministic evidence and reporting layers are the scientific core. The
optional calibration model is downstream of those rules and does not replace
them.

Start here
----------

* :doc:`quickstart` — install KmerSutra and run a small screen.
* :doc:`concepts` — understand marker ranks, raw evidence and reportable calls.
* :doc:`panel_building` — construct an outgroup-aware multi-k panel.
* :doc:`screening` — screen reads and interpret the output tables.
* :doc:`hpc` — run reproducibly on Slurm.
* :doc:`benchmarks/index` — reproduce the Zymo and locked ATCC benchmarks.

.. note::

   KmerSutra is not intended to maximise read assignment or claim the earliest
   possible detection. Use it when conservative, inspectable species-level
   evidence is more important than permissive classification.

.. toctree::
   :maxdepth: 2
   :caption: User guide

   quickstart
   concepts
   installation
   panel_building
   screening
   outputs
   hpc
   reproducibility
   troubleshooting

.. toctree::
   :maxdepth: 2
   :caption: Benchmarks

   benchmarks/index

.. toctree::
   :maxdepth: 2
   :caption: Reference

   cli_reference
   api
   development
