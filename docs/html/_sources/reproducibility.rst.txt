Reproducibility and provenance
==============================

A defensible KmerSutra run retains:

* the resolved configuration;
* input paths, sizes and SHA-256 digests;
* reference leakage audit;
* complete commands;
* software and Python versions;
* per-stage success manifests;
* declared output checksums;
* observed timing and memory;
* final run-status table.

Benchmark settings
------------------

Primary benchmark thresholds must be registered before examining the result.
If the data are used to select a threshold, that analysis is exploratory and
cannot also be presented as independent validation.

Testing artefacts
-----------------

``run_quality_checks.sh`` writes a dated result bundle beneath
``test_results/`` or a named project-storage path. Package inventories are
written under ``package_inventories/`` inside that run. Generated audits are
not placed in the repository root.
