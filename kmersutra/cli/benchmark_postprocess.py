"""Command-line interface for KmerSutra benchmark post-processing."""

from __future__ import annotations

import argparse
from pathlib import Path

from kmersutra.benchmark_postprocess import run_benchmark_postprocess
from kmersutra.logging_utils import configure_logging


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Run KmerSutra post-screening benchmark interpretation steps: LCA "
            "summary, AI-ready training table construction, optional AI "
            "calibrator training, and combined Excel/HTML reporting."
        )
    )
    parser.add_argument(
        "--summary_dir",
        default=None,
        help=(
            "Comparable-summary directory containing kmersutra_detection_calls_long.*. "
            "Used when --calls_table is omitted."
        ),
    )
    parser.add_argument(
        "--calls_table",
        default=None,
        help="Explicit long detection/call table (.tsv, .tsv.gz or .parquet).",
    )
    parser.add_argument(
        "--out_dir",
        required=True,
        help="Output directory for v0.35 benchmark post-processing files.",
    )
    parser.add_argument(
        "--taxonomy_dir",
        required=True,
        help="Directory containing NCBI taxonomy dump files.",
    )
    parser.add_argument(
        "--taxon_map_table",
        required=True,
        help="Taxon name-to-taxid mapping table, usually kmersutra_genome_config.tsv.",
    )
    parser.add_argument("--taxid_column", default=None)
    parser.add_argument("--taxon_name_column", default=None)
    parser.add_argument("--taxon_map_name_column", default=None)
    parser.add_argument("--taxon_map_taxid_column", default=None)
    parser.add_argument("--min_unique_kmers", type=int, default=1)
    parser.add_argument("--min_positive_sequences", type=int, default=1)
    parser.add_argument("--min_best_k", type=int, default=0)
    parser.add_argument("--max_not_detected", type=int, default=50000)
    parser.add_argument(
        "--train_calibrator",
        action="store_true",
        help="Also train/evaluate the AI call calibrator from the generated training table.",
    )
    parser.add_argument("--test_fraction", type=float, default=0.2)
    parser.add_argument("--distance_quantile", type=float, default=0.95)
    parser.add_argument("--group_columns", nargs="+", default=["sample_id"])
    parser.add_argument(
        "--no_excel",
        action="store_true",
        help="Do not write the formatted Excel workbook.",
    )
    parser.add_argument(
        "--no_html",
        action="store_true",
        help="Do not write the HTML report.",
    )
    parser.add_argument("--verbose", action="store_true")
    return parser.parse_args()


def main() -> None:
    """Run benchmark post-processing."""
    args = parse_args()
    out_dir = Path(args.out_dir).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    logger = configure_logging(
        log_file=out_dir / "kmersutra_benchmark_postprocess.log",
        verbose=args.verbose,
    )
    logger.info("Starting KmerSutra benchmark post-processing")
    outputs = run_benchmark_postprocess(
        summary_dir=args.summary_dir,
        calls_table=args.calls_table,
        out_dir=out_dir,
        taxonomy_dir=args.taxonomy_dir,
        taxon_map_table=args.taxon_map_table,
        taxid_column=args.taxid_column,
        taxon_name_column=args.taxon_name_column,
        taxon_map_name_column=args.taxon_map_name_column,
        taxon_map_taxid_column=args.taxon_map_taxid_column,
        min_unique_kmers=args.min_unique_kmers,
        min_positive_sequences=args.min_positive_sequences,
        min_best_k=args.min_best_k,
        max_not_detected=args.max_not_detected,
        train_calibrator=args.train_calibrator,
        test_fraction=args.test_fraction,
        distance_quantile=args.distance_quantile,
        group_columns=args.group_columns,
        write_excel=not args.no_excel,
        write_html=not args.no_html,
        logger=logger,
    )
    for name, path in sorted(outputs.items()):
        logger.info("Output %s: %s", name, path)
    logger.info("Finished KmerSutra benchmark post-processing")


if __name__ == "__main__":
    main()
