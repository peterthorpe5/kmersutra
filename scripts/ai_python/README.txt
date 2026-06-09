KmerSutra v0.50.0 AI validation Python code now lives inside the installed
kmersutra package rather than being embedded in qsub shell scripts.

The qsub wrappers in ../ call the following packaged command-line interfaces:

- kmersutra-build-call-training
- kmersutra-validate-call-calibrator
- kmersutra-label-zymo-calls
- kmersutra-predict
- kmersutra-evaluate-call-predictions

Keeping the Python implementation inside the package allows normal unit tests,
versioning, import checks and reproducible installation. This directory is kept
as a placeholder for project-specific Python helpers if they are needed later.
