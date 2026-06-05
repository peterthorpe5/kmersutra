#!/usr/bin/env python3
"""Summarise FASTQ read length and quality metrics for KmerSutra benchmarks.

This script supports two related use cases.

1. Mixed benchmark sample FASTQs listed in the comparable benchmark manifest.
2. Source simulated pathogen FASTQs stored in the source run directories under
   the project ``runs/`` folder.

Outputs are tab-separated. Per-file outputs are written as ``.tsv.gz`` and
group summaries are written as ``.tsv``. No comma-separated outputs are written.

Examples
--------
Summarise only source simulated pathogen reads::

    python summarise_fastq_read_metrics_v3.py \
        --manifest RUN_ROOT/kmersutra_comparable_manifest.tsv \
        --out_dir RUN_ROOT/read_metrics \
        --mode simulated_sources \
        --max_reads_per_fastq 100000 \
        --verbose

Summarise mixed sample FASTQs from the manifest::

    python summarise_fastq_read_metrics_v3.py \
        --manifest RUN_ROOT/kmersutra_comparable_manifest.tsv \
        --out_dir RUN_ROOT/read_metrics_mixed \
        --mode mixed_samples \
        --max_reads_per_fastq 100000 \
        --verbose
"""

from __future__ import annotations

import argparse
import gzip
import logging
import math
import statistics
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Tuple


LOGGER = logging.getLogger("summarise_fastq_read_metrics")


COLUMN_ALIASES = {
    "spike_n": ["spike_n", "spike_reads", "spike_n_per_genome"],
    "input_fastq": ["input_fastq", "source_path", "fastq_path", "sample_fastq"],
    "source_run_dir": ["source_run_dir", "source_run_root", "run_dir"],
    "benchmark_family": ["benchmark_family", "family"],
    "panel": ["panel", "panel_name"],
    "sample_id": ["sample_id", "sample"],
    "replicate": ["replicate", "rep"],
}


SOURCE_FASTQ_CANDIDATES = [
    "simulated_pathogen_unaligned_reads.fastq.gz",
    "simulated_pathogen_aligned_reads.fastq.gz",
    "simulated_pathogen_aligned_error_profile",
    "sim_pool.fastq.gz",
    "sim_pool.fastq",
    "train_reads.fastq.gz",
]


def configure_logging(verbose: bool) -> None:
    """Configure console logging.

    Parameters
    ----------
    verbose:
        If true, use INFO logging; otherwise use WARNING logging.
    """
    level = logging.INFO if verbose else logging.WARNING
    logging.basicConfig(
        level=level,
        format="%(asctime)s %(levelname)s [%(name)s] %(message)s",
    )


def open_text(path: Path):
    """Open plain or gzip-compressed text.

    Parameters
    ----------
    path:
        Input file path.

    Returns
    -------
    Text file handle
        Readable text handle.
    """
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("rt", encoding="utf-8", errors="replace")


def open_output_text(path: Path):
    """Open plain or gzip-compressed output text.

    Parameters
    ----------
    path:
        Output file path.

    Returns
    -------
    Text file handle
        Writable text handle.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    if str(path).endswith(".gz"):
        return gzip.open(path, "wt", encoding="utf-8")
    return path.open("wt", encoding="utf-8")


def read_tsv(path: Path) -> List[Dict[str, str]]:
    """Read a tab-separated file into dictionaries.

    Parameters
    ----------
    path:
        TSV file path.

    Returns
    -------
    list of dict
        Rows keyed by header names.

    Raises
    ------
    ValueError
        If the file is empty or has no header.
    """
    with open_text(path=path) as handle:
        header_line = handle.readline()
        if not header_line:
            raise ValueError(f"Input TSV is empty: {path}")
        header = header_line.rstrip("\n").split("\t")
        rows = []
        for line in handle:
            values = line.rstrip("\n").split("\t")
            row = {
                key: values[index] if index < len(values) else ""
                for index, key in enumerate(header)
            }
            rows.append(row)
    return rows


def normalise_manifest_columns(rows: List[Dict[str, str]]) -> List[Dict[str, str]]:
    """Normalise common comparable-manifest column aliases.

    Parameters
    ----------
    rows:
        Manifest rows.

    Returns
    -------
    list of dict
        Rows with canonical aliases added where possible.
    """
    if not rows:
        return rows

    available = set(rows[0])
    alias_map = {}
    for canonical, aliases in COLUMN_ALIASES.items():
        if canonical in available:
            continue
        for alias in aliases:
            if alias in available:
                alias_map[canonical] = alias
                LOGGER.info(
                    "Using manifest column alias %s for requested column %s",
                    alias,
                    canonical,
                )
                break

    if not alias_map:
        return rows

    normalised = []
    for row in rows:
        copied = dict(row)
        for canonical, alias in alias_map.items():
            copied[canonical] = row.get(alias, "")
        normalised.append(copied)
    return normalised


def parse_fastq(path: Path) -> Iterator[Tuple[str, str, str]]:
    """Yield FASTQ records as title, sequence and quality.

    Parameters
    ----------
    path:
        FASTQ or FASTQ.GZ path.

    Yields
    ------
    tuple
        Title, sequence and quality string.
    """
    with open_text(path=path) as handle:
        while True:
            title = handle.readline()
            if not title:
                break
            sequence = handle.readline()
            plus = handle.readline()
            quality = handle.readline()
            if not quality:
                LOGGER.warning("Truncated FASTQ record encountered in %s", path)
                break
            if not title.startswith("@"):
                LOGGER.warning("Unexpected FASTQ title in %s: %s", path, title[:80])
            yield title.rstrip("\n"), sequence.rstrip("\n"), quality.rstrip("\n")


def calculate_n50(lengths: Sequence[int]) -> float:
    """Calculate read N50 from read lengths.

    Parameters
    ----------
    lengths:
        Read lengths.

    Returns
    -------
    float
        N50 read length or NaN when no lengths are available.
    """
    if not lengths:
        return math.nan
    total = sum(lengths)
    threshold = total / 2
    running = 0
    for length in sorted(lengths, reverse=True):
        running += length
        if running >= threshold:
            return float(length)
    return float(sorted(lengths)[-1])


def safe_median(values: Sequence[float]) -> float:
    """Return median or NaN for an empty sequence.

    Parameters
    ----------
    values:
        Numeric values.

    Returns
    -------
    float
        Median value, or NaN.
    """
    if not values:
        return math.nan
    return float(statistics.median(values))


def summarise_fastq(path: Path, max_reads: int) -> Dict[str, object]:
    """Summarise read length and base quality metrics for a FASTQ file.

    Parameters
    ----------
    path:
        FASTQ or FASTQ.GZ path.
    max_reads:
        Maximum reads to inspect. Use 0 to read all records.

    Returns
    -------
    dict
        Summary metrics for one FASTQ file.
    """
    lengths: List[int] = []
    read_mean_qualities: List[float] = []
    total_bases = 0
    total_quality = 0
    q10_bases = 0
    q20_bases = 0
    q30_bases = 0

    for index, (_, sequence, quality) in enumerate(parse_fastq(path=path), start=1):
        if max_reads and index > max_reads:
            break
        length = len(sequence)
        lengths.append(length)
        total_bases += length
        if quality:
            qualities = [ord(character) - 33 for character in quality]
            quality_sum = sum(qualities)
            total_quality += quality_sum
            read_mean_qualities.append(quality_sum / len(qualities))
            q10_bases += sum(1 for value in qualities if value >= 10)
            q20_bases += sum(1 for value in qualities if value >= 20)
            q30_bases += sum(1 for value in qualities if value >= 30)

    n_reads = len(lengths)
    mean_read_length = total_bases / n_reads if n_reads else math.nan
    mean_q = total_quality / total_bases if total_bases else math.nan

    return {
        "fastq_path": str(path),
        "fastq_exists": path.exists(),
        "n_reads_inspected": n_reads,
        "max_reads_per_fastq": max_reads,
        "total_bases_inspected": total_bases,
        "mean_read_length": mean_read_length,
        "median_read_length": safe_median(lengths),
        "read_n50": calculate_n50(lengths),
        "min_read_length": min(lengths) if lengths else math.nan,
        "max_read_length": max(lengths) if lengths else math.nan,
        "mean_base_quality": mean_q,
        "median_read_mean_quality": safe_median(read_mean_qualities),
        "q10_base_fraction": q10_bases / total_bases if total_bases else math.nan,
        "q20_base_fraction": q20_bases / total_bases if total_bases else math.nan,
        "q30_base_fraction": q30_bases / total_bases if total_bases else math.nan,
    }


def write_tsv(rows: List[Dict[str, object]], path: Path) -> None:
    """Write dictionaries as a tab-separated table.

    Parameters
    ----------
    rows:
        Rows to write.
    path:
        Output path.
    """
    if rows:
        fieldnames = list(rows[0].keys())
    else:
        fieldnames = ["status"]
        rows = [{"status": "no_rows"}]

    with open_output_text(path=path) as handle:
        handle.write("\t".join(fieldnames) + "\n")
        for row in rows:
            handle.write(
                "\t".join(format_value(row.get(field, "")) for field in fieldnames)
                + "\n"
            )


def format_value(value: object) -> str:
    """Format a value for TSV output.

    Parameters
    ----------
    value:
        Value to format.

    Returns
    -------
    str
        TSV-safe value.
    """
    if value is None:
        return ""
    if isinstance(value, float):
        if math.isnan(value):
            return "NaN"
        return f"{value:.6g}"
    return str(value).replace("\t", " ").replace("\n", " ")


def infer_file_type(path: Path) -> str:
    """Infer source FASTQ type from filename.

    Parameters
    ----------
    path:
        FASTQ path.

    Returns
    -------
    str
        File type label.
    """
    name = path.name
    if "unaligned" in name:
        return "simulated_pathogen_unaligned"
    if "aligned" in name:
        return "simulated_pathogen_aligned"
    if name.startswith("sim_pool"):
        return "sim_pool"
    if name.startswith("train_reads"):
        return "train_reads"
    if name == "mixed.fastq.gz":
        return "mixed_sample"
    return "fastq"


def build_mixed_sample_targets(rows: List[Dict[str, str]]) -> List[Dict[str, str]]:
    """Build FASTQ targets for mixed benchmark samples.

    Parameters
    ----------
    rows:
        Normalised manifest rows.

    Returns
    -------
    list of dict
        Target metadata rows.
    """
    targets = []
    for row in rows:
        fastq = row.get("input_fastq", "")
        if not fastq:
            continue
        target = dict(row)
        target["fastq_path"] = fastq
        target["fastq_file_type"] = "mixed_sample"
        targets.append(target)
    return targets


def build_simulated_source_targets(rows: List[Dict[str, str]]) -> List[Dict[str, str]]:
    """Build FASTQ targets for source simulated pathogen read files.

    Parameters
    ----------
    rows:
        Normalised manifest rows.

    Returns
    -------
    list of dict
        Target metadata rows, deduplicated by source directory and file type.
    """
    seen = set()
    targets = []

    for row in rows:
        source_run_dir = row.get("source_run_dir", "")
        if not source_run_dir:
            continue
        run_dir = Path(source_run_dir)
        for filename in SOURCE_FASTQ_CANDIDATES:
            path = run_dir / filename
            if not path.exists() or path.is_dir():
                continue
            key = (str(path), row.get("benchmark_family", ""), row.get("panel", ""))
            if key in seen:
                continue
            seen.add(key)
            target = dict(row)
            target["fastq_path"] = str(path)
            target["fastq_file_type"] = infer_file_type(path=path)
            target["source_run_name"] = run_dir.name
            targets.append(target)

    return targets


def summarise_targets(targets: List[Dict[str, str]], max_reads: int) -> List[Dict[str, object]]:
    """Summarise all target FASTQ files.

    Parameters
    ----------
    targets:
        Target metadata rows.
    max_reads:
        Maximum reads per FASTQ. Use 0 for all reads.

    Returns
    -------
    list of dict
        Per-target summary rows.
    """
    summary_rows = []
    for index, target in enumerate(targets, start=1):
        path = Path(target["fastq_path"])
        LOGGER.info("Summarising %s/%s: %s", index, len(targets), path)

        base = {
            "sample_id": target.get("sample_id", ""),
            "benchmark_family": target.get("benchmark_family", ""),
            "panel": target.get("panel", ""),
            "replicate": target.get("replicate", ""),
            "spike_n": target.get("spike_n", ""),
            "source_run_name": target.get("source_run_name", ""),
            "fastq_file_type": target.get("fastq_file_type", infer_file_type(path)),
        }

        if not path.exists():
            row = {
                **base,
                "fastq_path": str(path),
                "fastq_exists": False,
                "n_reads_inspected": 0,
                "max_reads_per_fastq": max_reads,
                "total_bases_inspected": 0,
                "mean_read_length": math.nan,
                "median_read_length": math.nan,
                "read_n50": math.nan,
                "min_read_length": math.nan,
                "max_read_length": math.nan,
                "mean_base_quality": math.nan,
                "median_read_mean_quality": math.nan,
                "q10_base_fraction": math.nan,
                "q20_base_fraction": math.nan,
                "q30_base_fraction": math.nan,
            }
        else:
            row = {**base, **summarise_fastq(path=path, max_reads=max_reads)}
        summary_rows.append(row)
    return summary_rows


def numeric_mean(values: List[float]) -> float:
    """Calculate a mean excluding NaN values.

    Parameters
    ----------
    values:
        Numeric values.

    Returns
    -------
    float
        Mean or NaN.
    """
    clean = [value for value in values if not math.isnan(value)]
    if not clean:
        return math.nan
    return sum(clean) / len(clean)


def group_summary(rows: List[Dict[str, object]], group_columns: Sequence[str]) -> List[Dict[str, object]]:
    """Summarise per-FASTQ metrics by grouping columns.

    Parameters
    ----------
    rows:
        Per-FASTQ summary rows.
    group_columns:
        Columns used for grouping.

    Returns
    -------
    list of dict
        Grouped summary rows.
    """
    grouped: Dict[Tuple[str, ...], List[Dict[str, object]]] = defaultdict(list)
    for row in rows:
        key = tuple(str(row.get(column, "")) for column in group_columns)
        grouped[key].append(row)

    output = []
    numeric_columns = [
        "n_reads_inspected",
        "total_bases_inspected",
        "mean_read_length",
        "median_read_length",
        "read_n50",
        "mean_base_quality",
        "median_read_mean_quality",
        "q10_base_fraction",
        "q20_base_fraction",
        "q30_base_fraction",
    ]

    for key, group_rows in sorted(grouped.items()):
        row = {column: key[index] for index, column in enumerate(group_columns)}
        row["n_fastq_files"] = len(group_rows)
        row["n_existing_fastq_files"] = sum(
            1 for group_row in group_rows if str(group_row.get("fastq_exists")) == "True"
        )
        for column in numeric_columns:
            values = []
            for group_row in group_rows:
                value = group_row.get(column, math.nan)
                try:
                    values.append(float(value))
                except (TypeError, ValueError):
                    values.append(math.nan)
            row[f"mean_{column}"] = numeric_mean(values)
            row[f"median_{column}"] = safe_median(
                [value for value in values if not math.isnan(value)]
            )
        output.append(row)
    return output


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed arguments.
    """
    parser = argparse.ArgumentParser(
        description="Summarise KmerSutra benchmark FASTQ read metrics."
    )
    parser.add_argument("--manifest", required=True, help="Comparable manifest TSV.")
    parser.add_argument("--out_dir", required=True, help="Output directory.")
    parser.add_argument(
        "--mode",
        choices=["mixed_samples", "simulated_sources", "both"],
        default="simulated_sources",
        help="Which FASTQ files to summarise.",
    )
    parser.add_argument(
        "--max_reads_per_fastq",
        type=int,
        default=100000,
        help="Maximum reads to inspect per FASTQ; use 0 for all reads.",
    )
    parser.add_argument(
        "--group_columns",
        nargs="+",
        default=["benchmark_family", "panel", "spike_n", "fastq_file_type"],
        help="Columns for grouped summaries.",
    )
    parser.add_argument("--verbose", action="store_true", help="Enable INFO logging.")
    return parser.parse_args()


def main() -> None:
    """Run FASTQ read metric summarisation."""
    args = parse_args()
    configure_logging(verbose=args.verbose)

    manifest_path = Path(args.manifest)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    rows = normalise_manifest_columns(read_tsv(path=manifest_path))
    if not rows:
        raise ValueError(f"Manifest contains no rows: {manifest_path}")

    all_sample_rows: List[Dict[str, object]] = []
    all_group_rows: List[Dict[str, object]] = []

    if args.mode in {"mixed_samples", "both"}:
        mixed_targets = build_mixed_sample_targets(rows=rows)
        LOGGER.info("Found %s mixed sample FASTQ targets", len(mixed_targets))
        mixed_rows = summarise_targets(
            targets=mixed_targets,
            max_reads=args.max_reads_per_fastq,
        )
        mixed_group_rows = group_summary(
            rows=mixed_rows,
            group_columns=args.group_columns,
        )
        write_tsv(mixed_rows, out_dir / "mixed_sample_fastq_read_metrics.tsv.gz")
        write_tsv(mixed_group_rows, out_dir / "mixed_sample_group_fastq_read_metrics.tsv")
        all_sample_rows.extend(mixed_rows)
        all_group_rows.extend(mixed_group_rows)

    if args.mode in {"simulated_sources", "both"}:
        source_targets = build_simulated_source_targets(rows=rows)
        LOGGER.info("Found %s simulated source FASTQ targets", len(source_targets))
        source_rows = summarise_targets(
            targets=source_targets,
            max_reads=args.max_reads_per_fastq,
        )
        source_group_rows = group_summary(
            rows=source_rows,
            group_columns=args.group_columns,
        )
        write_tsv(source_rows, out_dir / "simulated_source_fastq_read_metrics.tsv.gz")
        write_tsv(
            source_group_rows,
            out_dir / "simulated_source_group_fastq_read_metrics.tsv",
        )
        all_sample_rows.extend(source_rows)
        all_group_rows.extend(source_group_rows)

    write_tsv(all_sample_rows, out_dir / "all_fastq_read_metrics.tsv.gz")
    write_tsv(all_group_rows, out_dir / "all_group_fastq_read_metrics.tsv")

    LOGGER.info("Wrote FASTQ metrics to %s", out_dir)


if __name__ == "__main__":
    main()
