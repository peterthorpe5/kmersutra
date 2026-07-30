Installation
============

Supported environment
---------------------

KmerSutra supports Python 3.10 through 3.12. Python 3.11 is used in the supplied
Conda environment and the benchmark instructions.

Conda installation
------------------

.. code-block:: bash

   git clone https://github.com/peterthorpe5/kmersutra.git
   cd kmersutra
   conda env create --file environment.yml

Verify the installed commands:

.. code-block:: bash

   conda run --name kmersutra kmersutra-build-panel --help
   conda run --name kmersutra kmersutra-screen --help
   conda run --name kmersutra kmersutra-run-benchmark --help

Update an existing environment after pulling a release:

.. code-block:: bash

   conda run --name kmersutra \
       python -m pip install \
       --no-deps \
       --force-reinstall \
       --editable .

Pip installation
----------------

Install the core package:

.. code-block:: bash

   python -m pip install .

Install all optional features and development tools:

.. code-block:: bash

   python -m pip install ".[all,dev,docs]"

Run the verified quality workflow:

.. code-block:: bash

   bash run_quality_checks.sh \
       --results-dir /absolute/project/path/kmersutra_test_results

All generated logs, coverage reports and package inventories are written
beneath the chosen results directory.
