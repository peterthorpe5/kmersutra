Benchmark workflows
===================

KmerSutra includes two distinct public mock-community workflows:

* :doc:`atcc_msa1003` is the locked independent PacBio HiFi validation. Its
  thresholds must not be retuned after examining the result.
* :doc:`zymo_d6300` restores the historical ONT Zymo D6300 screening workflow
  as a Slurm-native, named-option entry point.

The internal Plasmodium benchmark remains part of the scientific evidence
base, but it is not rerun by either public benchmark wrapper.

.. toctree::
   :maxdepth: 1

   atcc_msa1003
   zymo_d6300
