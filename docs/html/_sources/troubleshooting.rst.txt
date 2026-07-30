Troubleshooting
===============

Unresolved configuration variable
---------------------------------

If preflight reports ``KMERSUTRA_AI_MODEL`` or ``KMERSUTRA_DB_ROOT`` as
unresolved, pass explicit named paths:

.. code-block:: bash

   kmersutra-run-benchmark \
       --config benchmarks/atcc_msa1003_hifi/config.example.json \
       --database-root /absolute/path/to/kmersutra_db \
       --ai-model /absolute/path/to/final_internal_calibrator_all_training.json \
       --output-root /absolute/path/to/benchmark_outputs \
       --stop-after 00_preflight \
       --resume

Reference leakage failure
-------------------------

Do not disable the audit. Inspect
``reference_panel_leakage_audit.tsv``, remove truth assemblies or sequences
from the panel source configuration, rebuild and rerun preflight.

Only 13 tests appear
--------------------

The 13 tests are the separate historical AI-validation suite. The main suite
runs first. Use:

.. code-block:: bash

   bash run_quality_checks.sh \
       --results-dir /absolute/project/path/kmersutra_test_results

The main test count and branch coverage are recorded in the dated result
directory.

``nose2`` configuration error
-----------------------------

KmerSutra uses Python ``unittest`` and does not require ``nose2``. Run the
supplied scripts rather than invoking ``nose2``.

Slurm command cannot find Conda
-------------------------------

Run submission from a login shell where ``conda`` is on ``PATH``. The supplied
batch wrappers then use the resolved executable through ``conda run``.
