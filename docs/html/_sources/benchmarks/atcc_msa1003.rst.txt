ATCC MSA-1003 PacBio HiFi
=========================

Dataset
-------

The locked external benchmark uses ``SRR9328980`` and a 20-organism ATCC
MSA-1003 community with five organisms at each of four abundance tiers:
18%, 1.8%, 0.18% and 0.02%.

Frozen primary settings
-----------------------

* exact matching;
* k = 51, 77, 101 and 151;
* same-genus reportability fraction = 0.05;
* AI novelty scale = 2.9 without retraining;
* truth accessions excluded from the reference panel;
* flat-panel result primary;
* hierarchical result secondary;
* no primary threshold sweep.

Preflight
---------

Use the example configuration directly and supply site-specific paths as named
options:

.. code-block:: bash

   REPO="/absolute/path/to/kmersutra"
   DATABASE_ROOT="/absolute/path/to/kmersutra_db"
   AI_MODEL="/absolute/path/to/final_internal_calibrator_all_training.json"
   OUTPUT_ROOT="/absolute/path/to/kmersutra_atcc_msa1003"

   cd "${REPO}"

   conda run --name kmersutra \
       kmersutra-run-benchmark \
       --config benchmarks/atcc_msa1003_hifi/config.example.json \
       --database-root "${DATABASE_ROOT}" \
       --ai-model "${AI_MODEL}" \
       --output-root "${OUTPUT_ROOT}" \
       --threads 8 \
       --stop-after 00_preflight \
       --resume \
       --verbose

Submission
----------

.. code-block:: bash

   bash benchmarks/atcc_msa1003_hifi/submit_atcc_msa1003_benchmark.sh \
       --config benchmarks/atcc_msa1003_hifi/config.example.json \
       --database-root "${DATABASE_ROOT}" \
       --ai-model "${AI_MODEL}" \
       --output-root "${OUTPUT_ROOT}" \
       --account barton \
       --partition barton \
       --time 72:00:00 \
       --memory 128G \
       --threads 8 \
       --conda-env kmersutra \
       --resume

The preflight will still stop until the held-out panel and genome configuration
exist beneath ``${DATABASE_ROOT}/atcc_msa1003_heldout_v1`` and pass the leakage
audit.
