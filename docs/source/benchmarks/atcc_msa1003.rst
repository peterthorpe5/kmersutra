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

The established v4/v0.46 Plasmodium/outgroup panel contains none of the 20
ATCC target species and must not be used as the ATCC target panel. Version
0.51.2 first builds an extended, leakage-controlled panel. It requests held-out
assemblies for every expected species, samples additional species from each
target genus, excludes the published truth accessions, and retains the v4
genome collection as background evidence.

Build the reference panel before submitting the benchmark:

.. code-block:: bash

   REPO="/absolute/path/to/kmersutra"
   DATABASE_ROOT="/absolute/path/to/kmersutra_db"
   V4_CONFIG="${DATABASE_ROOT}/kmersutra_builds/selected_v4_build/inputs/kmersutra_genome_config_global_candidate.tsv"
   TRUTH_MANIFEST="${REPO}/benchmarks/atcc_msa1003_hifi/truth_manifest.tsv"
   TAXONOMY_DIR="${DATABASE_ROOT}/ncbi_taxonomy"

   bash "${REPO}/benchmarks/atcc_msa1003_hifi/submit_atcc_reference_panel.sh" \
       --repo "${REPO}" \
       --database-root "${DATABASE_ROOT}" \
       --background-config "${V4_CONFIG}" \
       --truth-manifest "${TRUTH_MANIFEST}" \
       --taxonomy-dir "${TAXONOMY_DIR}" \
       --email "name@example.org" \
       --account barton \
       --partition barton \
       --time 72:00:00 \
       --memory 128G \
       --threads 24 \
       --resume

The completed job writes
``${DATABASE_ROOT}/atcc_msa1003_heldout_v1/reference_audit/atcc_reference_gate_summary.tsv``.
Do not submit the benchmark unless ``gate_status`` is ``PASS`` and all 20
species are represented with zero leakage and zero incomplete FASTAs.

Benchmark preflight
-------------------

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

The full controller repeats the leakage audit during stage ``00_preflight``.
The older site-specific ``config.v4_v046.cluster.json`` pointed directly to the
wrong v4 panel and must not be reused.
