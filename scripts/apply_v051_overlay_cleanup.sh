#!/usr/bin/env bash

set -Eeuo pipefail
trap 'echo "ERROR at line ${LINENO}: ${BASH_COMMAND}" >&2' ERR

usage() {
    echo "Usage:"
    echo "  bash scripts/apply_v051_overlay_cleanup.sh --repo REPOSITORY"
}

REPO=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --repo)
            REPO="${2:?Missing value for --repo}"
            shift 2
            ;;
        --help|-h)
            usage
            exit 0
            ;;
        *)
            echo "ERROR: unknown option: $1" >&2
            usage >&2
            exit 2
            ;;
    esac
done

if [[ -z "${REPO}" ]]; then
    echo "ERROR: --repo is required" >&2
    usage >&2
    exit 2
fi
if [[ ! -d "${REPO}/.git" ]]; then
    echo "ERROR: not a Git repository: ${REPO}" >&2
    exit 1
fi

cd "${REPO}"

historical_root_files=(
    PACKAGE_FILE_INVENTORY_v0_39_0.txt
    PACKAGE_FILE_INVENTORY_v0_43_0.txt
    PACKAGE_FILE_INVENTORY_v0_45_0.txt
    PACKAGE_FILE_INVENTORY_v0_47_0.txt
    PACKAGE_FILE_INVENTORY_v0_48_0.txt
    PACKAGE_FILE_INVENTORY_v0_49_0.txt
    PACKAGE_FILE_INVENTORY_v0_50_0.txt
    test_results_v0_17_0.txt
    test_results_v0_18_0.txt
    test_results_v0_19_0.txt
    test_results_v0_20_0.txt
    test_results_v0_21_0.txt
    test_results_v0_22_0.txt
    test_results_v0_23_0.txt
    test_results_v0_25_0.txt
    test_results_v0_27_0.txt
    test_results_v0_28_0.txt
    test_results_v0_28_1.txt
    test_results_v0_29_0.txt
    test_results_v0_30_0.txt
    test_results_v0_31_0.txt
    test_results_v0_32_0.txt
    test_results_v0_32_0_nose2.txt
    test_results_v0_33_0.txt
    test_results_v0_34_0.txt
    test_results_v0_35_0.txt
    test_results_v0_36_0.txt
    test_results_v0_37_0.txt
    test_results_v0_38_0.txt
    test_results_v0_39_0.txt
    test_results_v0_39_0_full.txt
    test_results_v0_40_0.txt
    test_results_v0_41_0.txt
    test_results_v0_42_0.txt
    test_results_v0_43_0.txt
    test_results_v0_44_0.txt
    test_results_v0_45_0.txt
    test_results_v0_46_0.txt
    test_results_v0_47_0.txt
    test_results_v0_49_0.txt
    test_results_v0_50_0.txt
    test_results_v0_50_1.txt
    V0_26_RELEASE_NOTES.txt
    V0_27_RELEASE_NOTES.txt
    V0_29_RELEASE_NOTES.txt
    V0_30_RELEASE_NOTES.txt
    V0_31_RELEASE_NOTES.txt
    V0_32_1_NOSE2_REPAIR_NOTES.txt
    V0_32_RELEASE_NOTES.txt
    V0_33_RELEASE_NOTES.txt
    V0_34_RELEASE_NOTES.txt
    V0_35_RELEASE_NOTES.txt
    V0_36_RELEASE_NOTES.txt
    V0_37_RELEASE_NOTES.txt
    V0_38_RELEASE_NOTES.txt
    V0_39_RELEASE_NOTES.txt
    V0_40_RELEASE_NOTES.txt
    V0_41_RELEASE_NOTES.txt
    V0_42_RELEASE_NOTES.txt
    V0_43_RELEASE_NOTES.txt
    V0_44_RELEASE_NOTES.txt
    V0_45_RELEASE_NOTES.txt
    V0_46_RELEASE_NOTES.txt
    V0_47_RELEASE_NOTES.txt
    V0_48_RELEASE_NOTES.txt
    V0_49_RELEASE_NOTES.txt
    V0_50_1_RELEASE_NOTES.txt
    V0_50_RELEASE_NOTES.txt
    summary/test_results.txt
)

git rm -r --ignore-unmatch -- \
    "${historical_root_files[@]}" \
    .DS_Store \
    __pycache__ \
    kmersutra/__pycache__ \
    kmersutra/cli/__pycache__ \
    kmersutra_ai_validation/.DS_Store \
    summary.zip \
    tests/__pycache__ \
    tests/test_call_consolidation.py.tmp

echo "Repository cleanup staged successfully."
