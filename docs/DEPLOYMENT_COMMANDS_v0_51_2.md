# KmerSutra v0.51.2 cluster handoff

## Install the pulled release

```bash
REPO="/gpfs/uod-scale-01/cluster/gjb_lab/pthorpe001/2026_plasmodium_kraken_sensitivity/kmersutra"

cd "${REPO}"
git switch main
git pull --ff-only

conda run --name kmersutra \
    python -m pip install \
    --no-deps \
    --force-reinstall \
    --editable "${REPO}"

conda run --name kmersutra \
    kmersutra-prepare-atcc-reference --help >/dev/null
```

## Submit the ATCC reference build

Set `NCBI_EMAIL` to the email address that should be supplied to NCBI Entrez.

```bash
PROJECT_ROOT="/gpfs/uod-scale-01/cluster/gjb_lab/pthorpe001/2026_plasmodium_kraken_sensitivity"
REPO="${PROJECT_ROOT}/kmersutra"
DATABASE_ROOT="/home/pthorpe001/data/databases/kmersutra_db"
REFERENCE_ROOT="${DATABASE_ROOT}/atcc_msa1003_heldout_v1"
V4_BUILD="${DATABASE_ROOT}/kmersutra_builds/kmersutra_plasmodium_outgroups_v4_global_candidate_k51_77_101_151_rawontlodbalanced_v046_20260605_150721"
V4_CONFIG="${V4_BUILD}/inputs/kmersutra_genome_config_global_candidate.tsv"
TRUTH_MANIFEST="${REPO}/benchmarks/atcc_msa1003_hifi/truth_manifest.tsv"
TAXONOMY_DIR="${DATABASE_ROOT}/ncbi_taxonomy"
NCBI_EMAIL="YOUR_NCBI_EMAIL_ADDRESS"

cd "${REPO}"

bash benchmarks/atcc_msa1003_hifi/submit_atcc_reference_panel.sh \
    --repo "${REPO}" \
    --database-root "${DATABASE_ROOT}" \
    --reference-root "${REFERENCE_ROOT}" \
    --background-config "${V4_CONFIG}" \
    --truth-manifest "${TRUTH_MANIFEST}" \
    --taxonomy-dir "${TAXONOMY_DIR}" \
    --email "${NCBI_EMAIL}" \
    --account barton \
    --partition barton \
    --time 72:00:00 \
    --memory 128G \
    --threads 24 \
    --target-assemblies 5 \
    --neighbour-assemblies 25 \
    --minimum-target-refs 1 \
    --max-entrez-records 500 \
    --conda-env kmersutra \
    --resume
```

If the job reaches 72 hours, submit the identical command again. Stage tokens
are reused only when the configuration signature and declared outputs match.

## Verify the gate

```bash
GATE_SUMMARY="${REFERENCE_ROOT}/reference_audit/atcc_reference_gate_summary.tsv"
COVERAGE_TABLE="${REFERENCE_ROOT}/reference_audit/atcc_species_coverage.tsv"

test -s "${GATE_SUMMARY}" && cat "${GATE_SUMMARY}"
test -s "${COVERAGE_TABLE}" && cat "${COVERAGE_TABLE}"
```

The required result is 20 represented species, no missing species, no leakage,
no incomplete FASTAs and `gate_status` equal to `PASS`.

## Submit the complete locked benchmark

Run this only after the reference gate passes:

```bash
BENCHMARK_ROOT="${PROJECT_ROOT}/benchmarks/kmersutra_atcc_msa1003"
AI_MODEL="${PROJECT_ROOT}/ONT_ZymoBIOMICS_ENAERR5396170/ai_validation/runs_ai_full_internal_plasmodium_validation_v0501_safe_transformed_bounded_20260609_141127/outputs/final_internal_calibrator_all_training.json"
RUN_NAME="atcc_msa1003_hifi_srr9328980_heldout_v0512_20260731"

cd "${REPO}"

bash benchmarks/atcc_msa1003_hifi/submit_atcc_msa1003_benchmark.sh \
    --config benchmarks/atcc_msa1003_hifi/config.example.json \
    --database-root "${DATABASE_ROOT}" \
    --ai-model "${AI_MODEL}" \
    --output-root "${BENCHMARK_ROOT}" \
    --run-name "${RUN_NAME}" \
    --account barton \
    --partition barton \
    --threads 24 \
    --memory 128G \
    --time 72:00:00 \
    --conda-env kmersutra \
    --resume
```

Do not reuse `config.v4_v046.cluster.json`; it points to the panel that contains
zero ATCC target species.
