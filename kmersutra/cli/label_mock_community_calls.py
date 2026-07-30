"""Label KmerSutra calls using a generic mock-community truth manifest."""

from __future__ import annotations

import argparse
from pathlib import Path

from kmersutra.logging_utils import configure_logging
from kmersutra.mock_community import write_mock_ai_feature_table


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments.

    Args:
        argv: Optional explicit argument list.

    Returns:
        Parsed arguments.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Label KmerSutra species-call rows using a mock-community truth manifest."
        )
    )
    parser.add_argument("--calls_table", required=True)
    parser.add_argument("--truth_manifest", required=True)
    parser.add_argument("--out_table", required=True)
    parser.add_argument("--out_category_counts", required=True)
    parser.add_argument("--out_coarse_label_counts", required=True)
    parser.add_argument("--verbose", action="store_true")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    """Run manifest-driven truth labelling.

    Args:
        argv: Optional explicit argument list.

    Returns:
        Process exit status.
    """
    args = parse_args(argv)
    output_path = Path(args.out_table).expanduser().resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    logger = configure_logging(
        log_file=output_path.with_suffix(".log"),
        verbose=args.verbose,
    )
    records = write_mock_ai_feature_table(
        calls_table=args.calls_table,
        truth_manifest=args.truth_manifest,
        output_table=output_path,
        category_counts_table=args.out_category_counts,
        coarse_label_counts_table=args.out_coarse_label_counts,
        logger=logger,
    )
    logger.info("Wrote %d labelled mock-community records", len(records))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
