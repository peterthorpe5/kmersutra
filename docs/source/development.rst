Development and testing
=======================

Code changes require unit tests, defensive validation, logging and Google-style
docstrings. Python code follows a 100-character line length. Result tables are
tab-separated; compressed TSV and Parquet are used for large outputs.

Run the complete quality workflow:

.. code-block:: bash

   bash run_quality_checks.sh \
       --results-dir /absolute/project/path/kmersutra_test_results \
       --run-label local_validation

The command writes test logs, text/XML/HTML coverage, build logs, wheel audit
and package inventory under one results directory. Version 0.51.1 enforces 80%
whole-package branch coverage; 90% is the next development target.

Build the documentation:

.. code-block:: bash

   bash scripts/build_documentation.sh \
       --output-dir docs/html

Open ``docs/html/index.html`` locally or publish the same build through the
included GitHub Pages workflow or Read the Docs configuration.
