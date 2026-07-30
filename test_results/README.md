# Generated test and package-audit results

`run_tests.sh` and `run_quality_checks.sh` write each run beneath a dated
subdirectory here by default. Generated logs, coverage reports, HTML coverage
and package inventories are intentionally ignored by Git.

To place results on project storage instead, use the named option:

```bash
bash run_quality_checks.sh \
    --results-dir /absolute/project/path/kmersutra_test_results
```

Each quality run contains:

- `main_unittest.log`;
- `ai_validation_unittest.log`;
- `coverage_report.txt`, `coverage.xml` and `coverage_html/index.html`;
- package-build and wheel-content logs;
- `package_inventories/package_file_inventory.tsv`;
- `quality_check_summary.tsv`.
