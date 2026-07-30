"""Restartable, checksummed KmerSutra benchmark workflow."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import logging
import os
import platform
import re
import shutil
import subprocess
import sys
import time
import uuid
from collections.abc import Callable, Iterable, Mapping
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from kmersutra import __version__
from kmersutra.depth_subsets import create_depth_subsets
from kmersutra.mock_benchmark_summary import summarise_mock_benchmark
from kmersutra.mock_community import load_truth_manifest, write_mock_ai_feature_table
from kmersutra.reference_audit import write_reference_panel_audit
from kmersutra.table_io import read_records_table, write_records_table

STAGES = [
    "00_preflight",
    "01_acquire_reads",
    "02_depth_subsets",
    "03_screen_full",
    "04_screen_depths",
    "05_screen_single_k",
    "06_screen_hierarchical",
    "07_ai_validation",
    "08_summarise",
    "09_provenance",
]

LOCKED_ATCC_BENCHMARK_ID = "atcc_msa1003_hifi_srr9328980_locked_v1"
LOCKED_ATCC_K_VALUES = [51, 77, 101, 151]
LOCKED_ATCC_MINMIXED = 0.05
LOCKED_ATCC_NOVELTY_SCALE = 2.9
UNRESOLVED_VARIABLE_PATTERN = re.compile(
    r"\$\{(?P<name>[A-Za-z_][A-Za-z0-9_]*)\}"
)


@dataclass(frozen=True)
class StageResult:
    """Result returned by one completed benchmark stage.

    Attributes:
        outputs: Files required to validate stage completion.
        detail: Optional compact stage detail.
    """

    outputs: tuple[Path, ...]
    detail: str = ""


@dataclass(frozen=True)
class WorkflowOptions:
    """Runtime options that are not part of the frozen scientific config."""

    output_root: Path
    run_name: str
    threads: int
    resume: bool
    start_at: str | None
    stop_after: str | None
    force_stages: frozenset[str]
    dry_run: bool


def utc_now() -> str:
    """Return a second-resolution UTC timestamp."""
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def sha256_bytes(value: bytes) -> str:
    """Return a SHA-256 digest for bytes."""
    return hashlib.sha256(value).hexdigest()


def sha256_path(path: Path, *, chunk_size: int = 1024 * 1024) -> str:
    """Return a SHA-256 digest for one file.

    Args:
        path: File to hash.
        chunk_size: Read size in bytes.

    Returns:
        Hexadecimal digest.
    """
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_write_json(*, path: Path, value: Any) -> None:
    """Write JSON through a same-directory temporary file.

    Args:
        path: Output JSON path.
        value: JSON-serialisable value.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.parent / f".{path.name}.{uuid.uuid4().hex}.tmp"
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")
    temporary.replace(path)


def _expand_string(
    *,
    value: str,
    config_path: Path,
    environment_overrides: Mapping[str, str] | None = None,
) -> str:
    """Expand environment and repository placeholders in config text.

    Args:
        value: Configuration string.
        config_path: Source configuration path.
        environment_overrides: Optional explicit values that take precedence
            over process environment variables.

    Returns:
        Expanded configuration string.
    """
    repo_root = Path(__file__).resolve().parents[1]
    replaced = (
        value.replace("${CONFIG_DIR}", str(config_path.parent))
        .replace("${REPO_ROOT}", str(repo_root))
    )
    for name, replacement in (environment_overrides or {}).items():
        replaced = replaced.replace(f"${{{name}}}", replacement)
    return os.path.expandvars(os.path.expanduser(replaced))


def _expand_config_value(
    *,
    value: Any,
    config_path: Path,
    environment_overrides: Mapping[str, str] | None = None,
) -> Any:
    """Recursively expand string values in a JSON configuration."""
    if isinstance(value, str):
        return _expand_string(
            value=value,
            config_path=config_path,
            environment_overrides=environment_overrides,
        )
    if isinstance(value, list):
        return [
            _expand_config_value(
                value=item,
                config_path=config_path,
                environment_overrides=environment_overrides,
            )
            for item in value
        ]
    if isinstance(value, dict):
        return {
            key: _expand_config_value(
                value=item,
                config_path=config_path,
                environment_overrides=environment_overrides,
            )
            for key, item in value.items()
        }
    return value


def _find_unresolved_variables(*, value: Any) -> set[str]:
    """Return unresolved ``${VARIABLE}`` names from nested configuration."""
    if isinstance(value, str):
        return {
            match.group("name")
            for match in UNRESOLVED_VARIABLE_PATTERN.finditer(value)
        }
    if isinstance(value, list):
        return set().union(
            *(_find_unresolved_variables(value=item) for item in value),
            set(),
        )
    if isinstance(value, dict):
        return set().union(
            *(_find_unresolved_variables(value=item) for item in value.values()),
            set(),
        )
    return set()


def load_benchmark_config(
    *,
    config_path: str | Path,
    environment_overrides: Mapping[str, str] | None = None,
) -> tuple[dict[str, Any], str]:
    """Load, expand and digest a benchmark JSON configuration.

    Args:
        config_path: JSON configuration path.
        environment_overrides: Optional explicit placeholder values.

    Returns:
        Expanded configuration and canonical SHA-256 digest.

    Raises:
        FileNotFoundError: If the configuration is missing.
        ValueError: If the root is not a JSON object.
    """
    path = Path(config_path).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(f"Benchmark configuration not found: {path}")
    with path.open("r", encoding="utf-8") as handle:
        raw = json.load(handle)
    if not isinstance(raw, dict):
        raise ValueError("Benchmark configuration root must be a JSON object")
    config = _expand_config_value(
        value=raw,
        config_path=path,
        environment_overrides=environment_overrides,
    )
    unresolved = sorted(_find_unresolved_variables(value=config))
    if unresolved:
        exports = "\n".join(
            f"  export {name}=/absolute/path" for name in unresolved
        )
        raise ValueError(
            "Unresolved benchmark configuration variable(s): "
            f"{', '.join(unresolved)}.\n"
            "Set them in the environment or pass the corresponding named "
            "path option. For example:\n"
            f"{exports}"
        )
    canonical = json.dumps(
        config,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return config, sha256_bytes(canonical)


def require_mapping(
    *,
    config: dict[str, Any],
    key: str,
) -> dict[str, Any]:
    """Return a required mapping from configuration.

    Args:
        config: Configuration mapping.
        key: Required key.

    Returns:
        Nested mapping.

    Raises:
        ValueError: If the key is absent or not an object.
    """
    value = config.get(key)
    if not isinstance(value, dict):
        raise ValueError(f"Benchmark configuration requires object {key!r}")
    return value


def require_path(
    *,
    config: dict[str, Any],
    key: str,
    must_exist: bool = True,
) -> Path:
    """Resolve a required path from configuration.

    Args:
        config: Configuration mapping.
        key: Required key.
        must_exist: Require an existing non-empty file.

    Returns:
        Resolved path.

    Raises:
        ValueError: If the setting is blank.
        FileNotFoundError: If the required file is missing or empty.
    """
    text = str(config.get(key, "") or "").strip()
    if not text:
        raise ValueError(f"Benchmark configuration requires path {key!r}")
    path = Path(text).expanduser().resolve()
    if must_exist and (not path.is_file() or path.stat().st_size <= 0):
        raise FileNotFoundError(f"Configured file is missing or empty: {path}")
    return path


def validate_locked_atcc_config(*, config: dict[str, Any]) -> None:
    """Validate the pre-specified ATCC primary analysis settings.

    Args:
        config: Expanded benchmark configuration.

    Raises:
        ValueError: If a locked setting differs from the registered design.
    """
    benchmark_id = str(config.get("benchmark_id", ""))
    if benchmark_id != LOCKED_ATCC_BENCHMARK_ID:
        return
    dataset = require_mapping(config=config, key="dataset")
    screen = require_mapping(config=config, key="screen")
    ai = require_mapping(config=config, key="ai")
    if str(dataset.get("sra_accession", "")) != "SRR9328980":
        raise ValueError("Locked ATCC benchmark must use SRR9328980")
    if str(screen.get("screen_preset", "")) != "exact":
        raise ValueError("Locked ATCC benchmark must use exact screening")
    observed_k = [int(value) for value in screen.get("k_values", [])]
    if observed_k != LOCKED_ATCC_K_VALUES:
        raise ValueError(
            f"Locked ATCC k values must be {LOCKED_ATCC_K_VALUES}, "
            f"not {observed_k}"
        )
    observed_minmixed = float(
        screen.get("same_genus_reportable_min_fraction", -1.0)
    )
    if abs(observed_minmixed - LOCKED_ATCC_MINMIXED) > 1e-12:
        raise ValueError(
            "Locked ATCC same-genus reportability fraction must be 0.05"
        )
    if int(screen.get("max_mismatches", -1)) != 0:
        raise ValueError("Locked ATCC benchmark must use max_mismatches=0")
    if abs(float(ai.get("novelty_scale", -1.0)) - LOCKED_ATCC_NOVELTY_SCALE) > 1e-12:
        raise ValueError("Locked ATCC AI novelty scale must be 2.9")
    if bool(config.get("allow_primary_threshold_tuning", True)):
        raise ValueError(
            "Locked ATCC benchmark prohibits primary threshold tuning; "
            "set allow_primary_threshold_tuning=false"
        )


def validate_stage_name(value: str | None, *, option: str) -> None:
    """Validate an optional stage name."""
    if value is not None and value not in STAGES:
        raise ValueError(
            f"{option} must be one of {', '.join(STAGES)}; received {value!r}"
        )


def selected_stages(
    *,
    start_at: str | None,
    stop_after: str | None,
) -> list[str]:
    """Resolve an inclusive stage range.

    Args:
        start_at: Optional first stage.
        stop_after: Optional final stage.

    Returns:
        Ordered selected stages.

    Raises:
        ValueError: If the range is reversed or a stage is invalid.
    """
    validate_stage_name(start_at, option="--start-at")
    validate_stage_name(stop_after, option="--stop-after")
    start_index = STAGES.index(start_at) if start_at else 0
    stop_index = STAGES.index(stop_after) if stop_after else len(STAGES) - 1
    if start_index > stop_index:
        raise ValueError("--start-at occurs after --stop-after")
    return STAGES[start_index : stop_index + 1]


def command_display(command: Iterable[str]) -> str:
    """Return a shell-readable command without executing shell syntax."""
    import shlex

    return " ".join(shlex.quote(str(part)) for part in command)


def run_command(
    *,
    command: list[str],
    log_path: Path,
    logger: logging.Logger,
    cwd: Path | None = None,
    environment: dict[str, str] | None = None,
) -> None:
    """Run a command with combined file logging.

    Args:
        command: Command and arguments.
        log_path: Combined stdout/stderr log.
        logger: Workflow logger.
        cwd: Optional working directory.
        environment: Optional complete environment.

    Raises:
        RuntimeError: If the command exits unsuccessfully.
    """
    log_path.parent.mkdir(parents=True, exist_ok=True)
    logger.info("Running: %s", command_display(command))
    with log_path.open("w", encoding="utf-8") as handle:
        handle.write(f"command\t{command_display(command)}\n")
        handle.flush()
        completed = subprocess.run(
            command,
            cwd=cwd,
            env=environment,
            stdout=handle,
            stderr=subprocess.STDOUT,
            check=False,
            text=True,
        )
    if completed.returncode != 0:
        raise RuntimeError(
            f"Command failed with exit status {completed.returncode}: "
            f"{command_display(command)}. See {log_path}"
        )


def _open_panel_text(path: Path, mode: str):
    """Open a plain or gzip-compressed panel."""
    if path.suffix.lower() == ".gz":
        return gzip.open(path, mode, encoding="utf-8", newline="")
    return path.open(mode, encoding="utf-8", newline="")


def filter_panel_by_k(
    *,
    panel_path: str | Path,
    output_path: str | Path,
    k_value: int,
) -> int:
    """Stream one k value from a TSV/TSV.GZ panel.

    Args:
        panel_path: Input panel.
        output_path: Filtered output panel.
        k_value: Required k value.

    Returns:
        Number of retained rows.

    Raises:
        ValueError: If the panel has no recognised k column or retains no rows.
    """
    source = Path(panel_path).expanduser().resolve()
    destination = Path(output_path).expanduser().resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.parent / (
        f".{destination.name}.{uuid.uuid4().hex}.tmp{destination.suffix}"
    )
    retained = 0
    try:
        with _open_panel_text(source, "rt") as input_handle:
            reader = csv.DictReader(input_handle, delimiter="\t")
            fieldnames = list(reader.fieldnames or [])
            k_column = next(
                (column for column in ("k", "k_value") if column in fieldnames),
                None,
            )
            if k_column is None:
                raise ValueError(
                    f"Panel has no k or k_value column: {source}"
                )
            with _open_panel_text(temporary, "wt") as output_handle:
                writer = csv.DictWriter(
                    output_handle,
                    delimiter="\t",
                    fieldnames=fieldnames,
                    lineterminator="\n",
                )
                writer.writeheader()
                for row in reader:
                    try:
                        observed = int(str(row.get(k_column, "")))
                    except ValueError:
                        continue
                    if observed != int(k_value):
                        continue
                    writer.writerow(row)
                    retained += 1
        if retained == 0:
            raise ValueError(f"Panel contains no rows for k={k_value}: {source}")
        temporary.replace(destination)
    except Exception:
        temporary.unlink(missing_ok=True)
        raise
    return retained


class BenchmarkWorkflow:
    """Run restartable benchmark stages."""

    def __init__(
        self,
        *,
        config: dict[str, Any],
        config_digest: str,
        options: WorkflowOptions,
        logger: logging.Logger,
    ) -> None:
        """Initialise workflow state."""
        self.config = config
        self.config_digest = config_digest
        self.options = options
        self.logger = logger
        self.run_root = options.output_root / options.run_name
        self.stage_root = self.run_root / "stages"
        self.manifest_root = self.run_root / "workflow_control" / "stage_manifests"
        self.cache_root = self.run_root / "cache"
        self.current_stage = ""
        self.stage_methods: dict[str, Callable[[Path], StageResult]] = {
            "00_preflight": self.stage_preflight,
            "01_acquire_reads": self.stage_acquire_reads,
            "02_depth_subsets": self.stage_depth_subsets,
            "03_screen_full": self.stage_screen_full,
            "04_screen_depths": self.stage_screen_depths,
            "05_screen_single_k": self.stage_screen_single_k,
            "06_screen_hierarchical": self.stage_screen_hierarchical,
            "07_ai_validation": self.stage_ai_validation,
            "08_summarise": self.stage_summarise,
            "09_provenance": self.stage_provenance,
        }

    def prepare_run_root(self) -> None:
        """Create or validate the benchmark run root."""
        if self.run_root.exists() and not self.options.resume:
            raise FileExistsError(
                f"Run directory already exists; use --resume: {self.run_root}"
            )
        self.stage_root.mkdir(parents=True, exist_ok=True)
        self.manifest_root.mkdir(parents=True, exist_ok=True)
        self.cache_root.mkdir(parents=True, exist_ok=True)

    def manifest_path(self, stage: str) -> Path:
        """Return the control manifest path for a stage."""
        return self.manifest_root / f"{stage}.json"

    def stage_dir(self, stage: str) -> Path:
        """Return the durable directory for a stage."""
        return self.stage_root / stage

    def stage_is_complete(self, stage: str) -> bool:
        """Return whether a stage manifest and all declared outputs are valid."""
        manifest_path = self.manifest_path(stage)
        if not manifest_path.is_file():
            return False
        try:
            with manifest_path.open("r", encoding="utf-8") as handle:
                manifest = json.load(handle)
        except (OSError, json.JSONDecodeError):
            return False
        if manifest.get("status") != "success":
            return False
        if manifest.get("configuration_digest") != self.config_digest:
            return False
        outputs = manifest.get("outputs", [])
        if not isinstance(outputs, list) or not outputs:
            return False
        for output in outputs:
            path = Path(str(output))
            if not path.exists():
                return False
            if path.is_file() and path.stat().st_size <= 0:
                return False
        return True

    def write_stage_state(
        self,
        *,
        stage: str,
        status: str,
        started_at: str,
        outputs: Iterable[Path] = (),
        detail: str = "",
        error: str = "",
    ) -> None:
        """Write a stage control manifest."""
        atomic_write_json(
            path=self.manifest_path(stage),
            value={
                "stage": stage,
                "status": status,
                "configuration_digest": self.config_digest,
                "started_at_utc": started_at,
                "updated_at_utc": utc_now(),
                "outputs": [str(path) for path in outputs],
                "detail": detail,
                "error": error,
                "kmersutra_version": __version__,
            },
        )

    def run(self) -> None:
        """Run the selected stage range."""
        validate_locked_atcc_config(config=self.config)
        selected = selected_stages(
            start_at=self.options.start_at,
            stop_after=self.options.stop_after,
        )
        unknown_forced = self.options.force_stages.difference(STAGES)
        if unknown_forced:
            raise ValueError(
                "Unknown forced stage(s): " + ", ".join(sorted(unknown_forced))
            )
        self.prepare_run_root()
        for stage in selected:
            if (
                stage not in self.options.force_stages
                and self.options.resume
                and self.stage_is_complete(stage)
            ):
                self.logger.info("Skipping completed stage %s", stage)
                continue
            if self.options.dry_run:
                self.logger.info("Dry run: would execute stage %s", stage)
                continue
            started_at = utc_now()
            stage_directory = self.stage_dir(stage)
            stage_directory.mkdir(parents=True, exist_ok=True)
            self.write_stage_state(
                stage=stage,
                status="running",
                started_at=started_at,
            )
            self.logger.info("Starting stage %s", stage)
            start_time = time.perf_counter()
            try:
                self.current_stage = stage
                result = self.stage_methods[stage](stage_directory)
                elapsed = time.perf_counter() - start_time
                for output in result.outputs:
                    if not output.exists():
                        raise FileNotFoundError(
                            f"Stage {stage} did not create expected output: {output}"
                        )
                    if output.is_file() and output.stat().st_size <= 0:
                        raise ValueError(
                            f"Stage {stage} created an empty output: {output}"
                        )
                detail = (
                    f"elapsed_seconds={elapsed:.6f}"
                    + (f";{result.detail}" if result.detail else "")
                )
                self.write_stage_state(
                    stage=stage,
                    status="success",
                    started_at=started_at,
                    outputs=result.outputs,
                    detail=detail,
                )
                self.logger.info(
                    "Completed stage %s in %.3f seconds",
                    stage,
                    elapsed,
                )
            except Exception as exc:
                self.write_stage_state(
                    stage=stage,
                    status="failed",
                    started_at=started_at,
                    error=f"{type(exc).__name__}: {exc}",
                )
                self.logger.exception("Stage %s failed", stage)
                raise

    def _reference_config(self) -> dict[str, Any]:
        """Return the reference section."""
        return require_mapping(config=self.config, key="reference")

    def _dataset_config(self) -> dict[str, Any]:
        """Return the dataset section."""
        return require_mapping(config=self.config, key="dataset")

    def _screen_config(self) -> dict[str, Any]:
        """Return the screen section."""
        return require_mapping(config=self.config, key="screen")

    def _depth_config(self) -> dict[str, Any]:
        """Return the depth section."""
        return require_mapping(config=self.config, key="depth")

    def _ai_config(self) -> dict[str, Any]:
        """Return the AI section."""
        return require_mapping(config=self.config, key="ai")

    def truth_manifest_path(self) -> Path:
        """Return the truth manifest path."""
        return require_path(
            config=self._reference_config(),
            key="truth_manifest",
        )

    def flat_panel_path(self) -> Path:
        """Return the primary flat panel path."""
        return require_path(
            config=self._reference_config(),
            key="flat_panel",
        )

    def panel_genome_config_path(self) -> Path:
        """Return the panel genome configuration path."""
        return require_path(
            config=self._reference_config(),
            key="panel_genome_config",
        )

    def resolve_reads_path(self) -> Path:
        """Return local reads from config or the acquisition stage."""
        configured = str(self._dataset_config().get("input_reads", "") or "").strip()
        if configured:
            path = Path(configured).expanduser().resolve()
        else:
            path_file = self.stage_dir("01_acquire_reads") / "reads_path.txt"
            if not path_file.is_file():
                raise FileNotFoundError(
                    f"Reads acquisition path record is missing: {path_file}"
                )
            path = Path(path_file.read_text(encoding="utf-8").strip())
        if not path.is_file() or path.stat().st_size <= 0:
            raise FileNotFoundError(f"Benchmark reads are missing or empty: {path}")
        return path

    def stage_preflight(self, stage_dir: Path) -> StageResult:
        """Validate frozen configuration and reference independence."""
        truth_manifest = self.truth_manifest_path()
        load_truth_manifest(
            manifest_path=truth_manifest,
            logger=self.logger,
        )
        flat_panel = self.flat_panel_path()
        genome_config = self.panel_genome_config_path()
        audit_path = stage_dir / "reference_panel_leakage_audit.tsv"
        write_reference_panel_audit(
            genome_config=genome_config,
            truth_manifest=truth_manifest,
            output_table=audit_path,
            fail_on_leakage=True,
            logger=self.logger,
        )
        resolved_config = stage_dir / "resolved_config.json"
        atomic_write_json(path=resolved_config, value=self.config)
        input_rows = []
        for label, path in (
            ("truth_manifest", truth_manifest),
            ("flat_panel", flat_panel),
            ("panel_genome_config", genome_config),
        ):
            input_rows.append(
                {
                    "input_name": label,
                    "path": str(path),
                    "size_bytes": path.stat().st_size,
                    "sha256": sha256_path(path),
                }
            )
        ai = self._ai_config()
        if bool(ai.get("enabled", False)):
            model_path = require_path(config=ai, key="model_json")
            input_rows.append(
                {
                    "input_name": "ai_model_json",
                    "path": str(model_path),
                    "size_bytes": model_path.stat().st_size,
                    "sha256": sha256_path(model_path),
                }
            )
        input_manifest = stage_dir / "frozen_input_manifest.tsv"
        write_records_table(
            records=input_rows,
            output_path=input_manifest,
            fieldnames=["input_name", "path", "size_bytes", "sha256"],
            logger=self.logger,
        )
        return StageResult(
            outputs=(audit_path, resolved_config, input_manifest),
            detail=f"configuration_digest={self.config_digest}",
        )

    def stage_acquire_reads(self, stage_dir: Path) -> StageResult:
        """Resolve a local input or acquire an SRA FASTQ."""
        dataset = self._dataset_config()
        configured = str(dataset.get("input_reads", "") or "").strip()
        path_record = stage_dir / "reads_path.txt"
        if configured:
            reads = Path(configured).expanduser().resolve()
            if not reads.is_file() or reads.stat().st_size <= 0:
                raise FileNotFoundError(
                    f"Configured input_reads is missing or empty: {reads}"
                )
            path_record.write_text(str(reads) + "\n", encoding="utf-8")
            return StageResult(
                outputs=(path_record,),
                detail="input_source=configured_local_file",
            )

        accession = str(dataset.get("sra_accession", "") or "").strip()
        if not accession:
            raise ValueError(
                "dataset.input_reads or dataset.sra_accession is required"
            )
        for executable in ("prefetch", "fasterq-dump"):
            if shutil.which(executable) is None:
                raise RuntimeError(
                    f"{executable} is required to acquire {accession}; "
                    "install sra-tools or supply dataset.input_reads"
                )
        sra_root = stage_dir / "sra"
        fastq_root = stage_dir / "fastq"
        fasterq_tmp = stage_dir / "fasterq_tmp"
        sra_root.mkdir(parents=True, exist_ok=True)
        fastq_root.mkdir(parents=True, exist_ok=True)
        fasterq_tmp.mkdir(parents=True, exist_ok=True)
        run_command(
            command=[
                "prefetch",
                "--output-directory",
                str(sra_root),
                accession,
            ],
            log_path=stage_dir / "prefetch.log",
            logger=self.logger,
        )
        sra_candidates = sorted(sra_root.rglob(f"{accession}.sra"))
        sra_input = str(sra_candidates[0]) if len(sra_candidates) == 1 else accession
        run_command(
            command=[
                "fasterq-dump",
                sra_input,
                "--outdir",
                str(fastq_root),
                "--threads",
                str(self.options.threads),
                "--temp",
                str(fasterq_tmp),
            ],
            log_path=stage_dir / "fasterq_dump.log",
            logger=self.logger,
        )
        fastq_candidates = sorted(fastq_root.glob(f"{accession}*.fastq"))
        if len(fastq_candidates) != 1:
            raise ValueError(
                f"Expected one FASTQ for {accession}; found {len(fastq_candidates)}"
            )
        source_fastq = fastq_candidates[0]
        compressed = fastq_root / f"{accession}.fastq.gz"
        if shutil.which("pigz"):
            run_command(
                command=[
                    "pigz",
                    "--processes",
                    str(self.options.threads),
                    "--keep",
                    str(source_fastq),
                ],
                log_path=stage_dir / "pigz.log",
                logger=self.logger,
            )
            generated = Path(str(source_fastq) + ".gz")
            if generated != compressed:
                generated.replace(compressed)
        else:
            with source_fastq.open("rb") as source, gzip.open(
                compressed,
                "wb",
            ) as target:
                shutil.copyfileobj(source, target)
        source_fastq.unlink()
        path_record.write_text(str(compressed) + "\n", encoding="utf-8")
        return StageResult(
            outputs=(path_record, compressed),
            detail=f"sra_accession={accession}",
        )

    def stage_depth_subsets(self, stage_dir: Path) -> StageResult:
        """Create deterministic depth subsets."""
        depth = self._depth_config()
        manifest = stage_dir / "depth_subset_manifest.tsv"
        create_depth_subsets(
            input_path=self.resolve_reads_path(),
            input_format=str(self._dataset_config().get("input_format", "fastq")),
            output_dir=stage_dir / "reads",
            fractions=[float(value) for value in depth.get("fractions", [])],
            seeds=[int(value) for value in depth.get("seeds", [])],
            sample_prefix=str(self.config.get("sample_id", "benchmark")),
            compress=True,
            decompressor=str(
                self._screen_config().get("decompressor", "auto")
            ),
            manifest_path=manifest,
            logger=self.logger,
        )
        return StageResult(outputs=(manifest,))

    def _screen_command(
        self,
        *,
        reads: Path,
        panel: Path,
        sample_id: str,
        output_dir: Path,
        hierarchical: bool = False,
        module_manifest: Path | None = None,
        single_k: bool = False,
    ) -> list[str]:
        """Build one exact KmerSutra screen command."""
        screen = self._screen_config()
        command = [
            sys.executable,
            "-m",
            "kmersutra.cli.screen_reads_for_clade_kmers",
            "--input",
            str(reads),
            "--input_format",
            str(self._dataset_config().get("input_format", "fastq")),
            "--sample_id",
            sample_id,
            "--out_dir",
            str(output_dir),
            "--screen_mode",
            "hierarchical" if hierarchical else "flat",
            "--screen_preset",
            str(screen.get("screen_preset", "exact")),
            "--max_mismatches",
            str(int(screen.get("max_mismatches", 0))),
            "--call_preset",
            str(screen.get("call_preset", "lineage_aware")),
            "--same_genus_reportable_min_fraction",
            str(
                float(
                    screen.get(
                        "same_genus_reportable_min_fraction",
                        0.05,
                    )
                )
            ),
            "--threads",
            str(self.options.threads),
            "--chunk_size",
            str(int(screen.get("chunk_size", 10_000))),
            "--decompressor",
            str(screen.get("decompressor", "auto")),
            "--panel_cache",
            str(self.cache_root / f"{panel.name}.pickle"),
            "--use_panel_cache",
            "--write_panel_cache",
            "--consolidate_species_calls",
            "--no_read_level_hits",
            "--profile",
            "--verbose",
        ]
        if bool(screen.get("write_parquet_outputs", True)):
            command.append("--write_parquet_outputs")
        if hierarchical:
            if module_manifest is None:
                raise ValueError("Hierarchical screening requires a module manifest")
            command.extend(["--module_manifest", str(module_manifest)])
        else:
            command.extend(["--panel", str(panel)])
        if single_k:
            command.extend(
                [
                    "--min_k_values_positive",
                    str(int(screen.get("single_k_min_k_values_positive", 1))),
                    "--min_best_k",
                    str(int(screen.get("single_k_min_best_k", 0))),
                ]
            )
        return command

    def _run_screen_task(
        self,
        *,
        task_root: Path,
        reads: Path,
        panel: Path,
        sample_id: str,
        hierarchical: bool = False,
        module_manifest: Path | None = None,
        single_k: bool = False,
    ) -> Path:
        """Run one screen through an atomic task directory."""
        calls = task_root / "species_detection_calls.tsv"
        force_task = self.current_stage in self.options.force_stages
        if calls.is_file() and calls.stat().st_size > 0 and not force_task:
            self.logger.info("Skipping completed screen task %s", sample_id)
            return calls
        temporary = task_root.parent / f".{task_root.name}.{uuid.uuid4().hex}.tmp"
        if temporary.exists():
            shutil.rmtree(temporary)
        temporary.mkdir(parents=True, exist_ok=False)
        try:
            command = self._screen_command(
                reads=reads,
                panel=panel,
                sample_id=sample_id,
                output_dir=temporary,
                hierarchical=hierarchical,
                module_manifest=module_manifest,
                single_k=single_k,
            )
            run_command(
                command=command,
                log_path=temporary / "command.log",
                logger=self.logger,
            )
            temporary_calls = temporary / "species_detection_calls.tsv"
            if (
                not temporary_calls.is_file()
                or temporary_calls.stat().st_size <= 0
            ):
                raise ValueError(
                    f"Screen task did not create species calls: {sample_id}"
                )
            if task_root.exists():
                shutil.rmtree(task_root)
            temporary.replace(task_root)
        except Exception:
            shutil.rmtree(temporary, ignore_errors=True)
            raise
        return calls

    def stage_screen_full(self, stage_dir: Path) -> StageResult:
        """Run the primary full-depth flat screen."""
        sample_id = str(self.config.get("sample_id", "benchmark"))
        task_root = stage_dir / sample_id
        calls = self._run_screen_task(
            task_root=task_root,
            reads=self.resolve_reads_path(),
            panel=self.flat_panel_path(),
            sample_id=sample_id,
        )
        return StageResult(outputs=(calls,))

    def stage_screen_depths(self, stage_dir: Path) -> StageResult:
        """Screen every deterministic depth subset."""
        subset_manifest = self.stage_dir("02_depth_subsets") / "depth_subset_manifest.tsv"
        rows = read_records_table(
            input_path=subset_manifest,
            required_columns=["sample_id", "reads_path", "fraction", "seed"],
            logger=self.logger,
        )
        task_rows: list[dict[str, object]] = []
        calls_paths: list[Path] = []
        for row in rows:
            sample_id = str(row["sample_id"])
            calls = self._run_screen_task(
                task_root=stage_dir / sample_id,
                reads=Path(row["reads_path"]),
                panel=self.flat_panel_path(),
                sample_id=sample_id,
            )
            calls_paths.append(calls)
            task_rows.append(
                {
                    "sample_id": sample_id,
                    "analysis_type": "depth_series",
                    "fraction": row["fraction"],
                    "seed": row["seed"],
                    "k_value": "multi_k",
                    "calls_table": str(calls),
                }
            )
        task_manifest = stage_dir / "depth_screen_tasks.tsv"
        write_records_table(
            records=task_rows,
            output_path=task_manifest,
            fieldnames=[
                "sample_id",
                "analysis_type",
                "fraction",
                "seed",
                "k_value",
                "calls_table",
            ],
            logger=self.logger,
        )
        return StageResult(outputs=(task_manifest, *calls_paths))

    def stage_screen_single_k(self, stage_dir: Path) -> StageResult:
        """Run full-depth single-k ablations."""
        k_values = [int(value) for value in self._screen_config()["k_values"]]
        task_rows: list[dict[str, object]] = []
        calls_paths: list[Path] = []
        for k_value in k_values:
            filtered_panel = stage_dir / "panels" / f"panel_k{k_value}.tsv.gz"
            if (
                self.current_stage in self.options.force_stages
                or not filtered_panel.is_file()
                or filtered_panel.stat().st_size <= 0
            ):
                filter_panel_by_k(
                    panel_path=self.flat_panel_path(),
                    output_path=filtered_panel,
                    k_value=k_value,
                )
            sample_id = (
                f"{self.config.get('sample_id', 'benchmark')}_single_k_{k_value}"
            )
            calls = self._run_screen_task(
                task_root=stage_dir / sample_id,
                reads=self.resolve_reads_path(),
                panel=filtered_panel,
                sample_id=sample_id,
                single_k=True,
            )
            calls_paths.append(calls)
            task_rows.append(
                {
                    "sample_id": sample_id,
                    "analysis_type": "single_k_ablation",
                    "fraction": "1.00000000",
                    "seed": "",
                    "k_value": k_value,
                    "calls_table": str(calls),
                }
            )
        task_manifest = stage_dir / "single_k_screen_tasks.tsv"
        write_records_table(
            records=task_rows,
            output_path=task_manifest,
            fieldnames=[
                "sample_id",
                "analysis_type",
                "fraction",
                "seed",
                "k_value",
                "calls_table",
            ],
            logger=self.logger,
        )
        return StageResult(outputs=(task_manifest, *calls_paths))

    def stage_screen_hierarchical(self, stage_dir: Path) -> StageResult:
        """Run the pre-specified secondary hierarchical screen."""
        reference = self._reference_config()
        manifest_text = str(reference.get("module_manifest", "") or "").strip()
        skipped = stage_dir / "stage_skipped.tsv"
        if not manifest_text:
            write_records_table(
                records=[
                    {
                        "stage": "06_screen_hierarchical",
                        "status": "skipped",
                        "reason": "reference.module_manifest_not_configured",
                    }
                ],
                output_path=skipped,
                fieldnames=["stage", "status", "reason"],
                logger=self.logger,
            )
            return StageResult(outputs=(skipped,), detail="skipped=true")
        module_manifest = Path(manifest_text).expanduser().resolve()
        if not module_manifest.is_file():
            raise FileNotFoundError(
                f"Configured module manifest is missing: {module_manifest}"
            )
        sample_id = f"{self.config.get('sample_id', 'benchmark')}_hierarchical"
        calls = self._run_screen_task(
            task_root=stage_dir / sample_id,
            reads=self.resolve_reads_path(),
            panel=self.flat_panel_path(),
            sample_id=sample_id,
            hierarchical=True,
            module_manifest=module_manifest,
        )
        task_manifest = stage_dir / "hierarchical_screen_tasks.tsv"
        write_records_table(
            records=[
                {
                    "sample_id": sample_id,
                    "analysis_type": "hierarchical_secondary",
                    "fraction": "1.00000000",
                    "seed": "",
                    "k_value": "multi_k",
                    "calls_table": str(calls),
                }
            ],
            output_path=task_manifest,
            fieldnames=[
                "sample_id",
                "analysis_type",
                "fraction",
                "seed",
                "k_value",
                "calls_table",
            ],
            logger=self.logger,
        )
        return StageResult(outputs=(task_manifest, calls))

    def stage_ai_validation(self, stage_dir: Path) -> StageResult:
        """Apply the frozen AI model without retraining."""
        ai = self._ai_config()
        skipped = stage_dir / "stage_skipped.tsv"
        if not bool(ai.get("enabled", False)):
            write_records_table(
                records=[
                    {
                        "stage": "07_ai_validation",
                        "status": "skipped",
                        "reason": "ai.enabled=false",
                    }
                ],
                output_path=skipped,
                fieldnames=["stage", "status", "reason"],
                logger=self.logger,
            )
            return StageResult(outputs=(skipped,), detail="skipped=true")
        model = require_path(config=ai, key="model_json")
        calls = (
            self.stage_dir("03_screen_full")
            / str(self.config.get("sample_id", "benchmark"))
            / "species_detection_calls.tsv"
        )
        features = stage_dir / "atcc_mock_features.tsv.gz"
        category_counts = stage_dir / "atcc_truth_category_counts.tsv"
        label_counts = stage_dir / "atcc_coarse_label_counts.tsv"
        write_mock_ai_feature_table(
            calls_table=calls,
            truth_manifest=self.truth_manifest_path(),
            output_table=features,
            category_counts_table=category_counts,
            coarse_label_counts_table=label_counts,
            logger=self.logger,
        )
        predictions = stage_dir / "atcc_ai_predictions.tsv.gz"
        run_command(
            command=[
                sys.executable,
                "-m",
                "kmersutra.cli.predict_classifier",
                "--features_tsv",
                str(features),
                "--model_json",
                str(model),
                "--out_tsv",
                str(predictions),
                "--novelty_scale",
                str(float(ai.get("novelty_scale", 2.9))),
                "--verbose",
            ],
            log_path=stage_dir / "predict.log",
            logger=self.logger,
        )
        metrics = stage_dir / "atcc_ai_metrics.tsv"
        prediction_counts = stage_dir / "atcc_ai_prediction_counts.tsv"
        truth_counts = stage_dir / "atcc_ai_truth_counts.tsv"
        run_command(
            command=[
                sys.executable,
                "-m",
                "kmersutra.cli.evaluate_call_predictions",
                "--predictions_table",
                str(predictions),
                "--out_metrics",
                str(metrics),
                "--out_prediction_counts",
                str(prediction_counts),
                "--out_label_counts",
                str(truth_counts),
                "--verbose",
            ],
            log_path=stage_dir / "evaluate.log",
            logger=self.logger,
        )
        return StageResult(
            outputs=(
                features,
                category_counts,
                label_counts,
                predictions,
                metrics,
                prediction_counts,
                truth_counts,
            )
        )

    def _collect_task_rows(self) -> list[dict[str, str]]:
        """Collect completed screen task manifests."""
        rows = [
            {
                "sample_id": str(self.config.get("sample_id", "benchmark")),
                "analysis_type": "full_flat_primary",
                "fraction": "1.00000000",
                "seed": "",
                "k_value": "multi_k",
                "calls_table": str(
                    self.stage_dir("03_screen_full")
                    / str(self.config.get("sample_id", "benchmark"))
                    / "species_detection_calls.tsv"
                ),
            }
        ]
        for manifest in (
            self.stage_dir("04_screen_depths") / "depth_screen_tasks.tsv",
            self.stage_dir("05_screen_single_k") / "single_k_screen_tasks.tsv",
            self.stage_dir("06_screen_hierarchical")
            / "hierarchical_screen_tasks.tsv",
        ):
            if manifest.is_file() and manifest.stat().st_size > 0:
                rows.extend(read_records_table(input_path=manifest))
        return rows

    def stage_summarise(self, stage_dir: Path) -> StageResult:
        """Summarise primary, depth, ablation and hierarchical tasks."""
        task_manifest = stage_dir / "benchmark_tasks.tsv"
        rows = self._collect_task_rows()
        write_records_table(
            records=rows,
            output_path=task_manifest,
            fieldnames=[
                "sample_id",
                "analysis_type",
                "fraction",
                "seed",
                "k_value",
                "calls_table",
            ],
            logger=self.logger,
        )
        paths = summarise_mock_benchmark(
            task_manifest=task_manifest,
            truth_manifest=self.truth_manifest_path(),
            output_dir=stage_dir / "tables",
            logger=self.logger,
        )
        return StageResult(outputs=(task_manifest, *paths.values()))

    def stage_provenance(self, stage_dir: Path) -> StageResult:
        """Write final checksums, environment and run-state summary."""
        input_paths = [
            self.truth_manifest_path(),
            self.flat_panel_path(),
            self.panel_genome_config_path(),
            self.resolve_reads_path(),
        ]
        ai = self._ai_config()
        if bool(ai.get("enabled", False)):
            input_paths.append(require_path(config=ai, key="model_json"))
        checksum_rows = []
        for path in input_paths:
            self.logger.info("Calculating final SHA-256 for %s", path)
            checksum_rows.append(
                {
                    "path": str(path),
                    "size_bytes": path.stat().st_size,
                    "sha256": sha256_path(path),
                }
            )
        checksums = stage_dir / "input_checksums.tsv"
        write_records_table(
            records=checksum_rows,
            output_path=checksums,
            fieldnames=["path", "size_bytes", "sha256"],
            logger=self.logger,
        )
        environment_path = stage_dir / "runtime_environment.json"
        atomic_write_json(
            path=environment_path,
            value={
                "generated_at_utc": utc_now(),
                "kmersutra_version": __version__,
                "python_version": sys.version,
                "python_executable": sys.executable,
                "platform": platform.platform(),
                "hostname": platform.node(),
                "slurm_job_id": os.environ.get("SLURM_JOB_ID", ""),
                "slurm_array_task_id": os.environ.get(
                    "SLURM_ARRAY_TASK_ID",
                    "",
                ),
                "threads": self.options.threads,
                "configuration_digest": self.config_digest,
            },
        )
        run_status = stage_dir / "run_status.tsv"
        status_rows = []
        for stage in STAGES:
            manifest = self.manifest_path(stage)
            status = "not_run"
            if stage == self.current_stage:
                # This file is written before the controller can atomically mark
                # the provenance stage successful. Reaching this branch means
                # every preceding operation in the stage has succeeded.
                status = "success"
            elif manifest.is_file():
                with manifest.open("r", encoding="utf-8") as handle:
                    status = str(json.load(handle).get("status", "unknown"))
            status_rows.append({"stage": stage, "status": status})
        write_records_table(
            records=status_rows,
            output_path=run_status,
            fieldnames=["stage", "status"],
            logger=self.logger,
        )
        return StageResult(outputs=(checksums, environment_path, run_status))


def build_parser() -> argparse.ArgumentParser:
    """Build the benchmark-controller argument parser."""
    parser = argparse.ArgumentParser(
        description="Run a restartable KmerSutra benchmark workflow."
    )
    parser.add_argument("--config", required=True)
    parser.add_argument("--output-root", required=True)
    parser.add_argument(
        "--database-root",
        default=None,
        help=(
            "Override ${KMERSUTRA_DB_ROOT} in the configuration with this "
            "absolute database root."
        ),
    )
    parser.add_argument(
        "--ai-model",
        default=None,
        help=(
            "Override ${KMERSUTRA_AI_MODEL} in the configuration with this "
            "frozen calibrator JSON."
        ),
    )
    parser.add_argument("--run-name", default=None)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--start-at", choices=STAGES, default=None)
    parser.add_argument("--stop-after", choices=STAGES, default=None)
    parser.add_argument(
        "--force-stage",
        action="append",
        choices=STAGES,
        default=[],
    )
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--verbose", action="store_true")
    return parser


def run_benchmark(argv: list[str] | None = None) -> int:
    """Parse arguments and run the benchmark controller.

    Args:
        argv: Optional explicit argument list.

    Returns:
        Process exit status.
    """
    args = build_parser().parse_args(argv)
    if args.threads <= 0:
        raise ValueError("--threads must be positive")
    environment_overrides = {}
    if args.database_root:
        environment_overrides["KMERSUTRA_DB_ROOT"] = str(
            Path(args.database_root).expanduser().resolve()
        )
    if args.ai_model:
        environment_overrides["KMERSUTRA_AI_MODEL"] = str(
            Path(args.ai_model).expanduser().resolve()
        )
    config, digest = load_benchmark_config(
        config_path=args.config,
        environment_overrides=environment_overrides,
    )
    validate_locked_atcc_config(config=config)
    output_root = Path(args.output_root).expanduser().resolve()
    run_name = (
        str(args.run_name).strip()
        if args.run_name
        else str(config.get("run_name", "")).strip()
    )
    if not run_name:
        raise ValueError("--run-name or config run_name is required")
    options = WorkflowOptions(
        output_root=output_root,
        run_name=run_name,
        threads=args.threads,
        resume=args.resume,
        start_at=args.start_at,
        stop_after=args.stop_after,
        force_stages=frozenset(args.force_stage),
        dry_run=args.dry_run,
    )
    run_root = output_root / run_name
    from kmersutra.logging_utils import configure_logging

    logger = configure_logging(
        log_file=run_root / "logs" / "benchmark_controller.log",
        verbose=args.verbose,
    )
    logger.info("KmerSutra benchmark controller version %s", __version__)
    logger.info("Configuration digest: %s", digest)
    workflow = BenchmarkWorkflow(
        config=config,
        config_digest=digest,
        options=options,
        logger=logger,
    )
    workflow.run()
    return 0
