Slurm execution
===============

KmerSutra benchmark wrappers submit from a login node and invoke the installed
commands through ``conda run``. They do not assume that ``conda activate`` works
inside a non-interactive batch shell.

ATCC benchmark
--------------

Use :doc:`benchmarks/atcc_msa1003` for the locked, restartable workflow.

Zymo benchmark
--------------

Use :doc:`benchmarks/zymo_d6300` for the restored Slurm screening wrapper.

Resuming
--------

The benchmark controller validates configuration digests, declared outputs and
checksums before skipping a stage. A file merely existing is not evidence that
the stage completed.

If a 72-hour job reaches its limit, resubmit the same ATCC command with
``--resume``. Use ``--start-at``, ``--stop-after`` and repeatable
``--force-stage`` only when a controlled partial run is needed.

Resource reporting
------------------

The wrappers retain ``/usr/bin/time --verbose`` and ``sacct`` output. Report
observed peak RSS/MaxRSS and wall time. The Slurm ``--mem`` request is an
allocation, not measured KmerSutra memory use.
