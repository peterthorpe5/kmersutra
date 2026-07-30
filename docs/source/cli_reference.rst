Command reference
=================

Core commands
-------------

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Command
     - Purpose
   * - ``kmersutra-build-panel``
     - Build an outgroup-aware, multi-k marker panel.
   * - ``kmersutra-screen``
     - Screen FASTQ or FASTA reads and produce rank-aware calls.
   * - ``kmersutra-run-benchmark``
     - Run a checksummed, restartable benchmark workflow.
   * - ``kmersutra-audit-reference-panel``
     - Detect truth-assembly or sequence leakage.
   * - ``kmersutra-subset-reads``
     - Produce deterministic depth subsets.
   * - ``kmersutra-summarise-mock-benchmark``
     - Summarise truth-labelled mock-community results.
   * - ``kmersutra-package-inventory``
     - Write a deterministic TSV file inventory outside the repository root.

Panel and module commands
-------------------------

``kmersutra-download-ncbi``, ``kmersutra-download-genomes``,
``kmersutra-validate-panel``, ``kmersutra-merge-panels`` and
``kmersutra-merge-modules`` support reference acquisition, validation and
modular scaling.

Reporting and calibration commands
----------------------------------

The package retains the historical Zymo, comparable-benchmark, threshold-sweep,
feature-extraction and optional calibration commands. Run ``COMMAND --help``
for the complete current option list. The compatibility inventory is tested in
``tests/test_cli_compatibility_v051.py``.
