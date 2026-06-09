#!/bin/bash
set -Eeuo pipefail
python -m unittest \
    tests.test_ai_calibration \
    tests.test_ai_transformed_features \
    tests.test_zymo_truth \
    tests.test_call_validation
