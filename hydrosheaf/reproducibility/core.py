"""Small, dependency-light primitives for reproducible output runs.

This module does not make an arbitrary historical analysis reproducible by
itself.  It fixes process-level sources of variation and records enough
information to distinguish an exact replay from a newly generated rerun.
"""

from __future__ import annotations

import hashlib
import json
import os
import random
from pathlib import Path
from typing import Any, Iterable, Mapping


DEFAULT_SEED = 1729
DEFAULT_EPOCH = "2026-09-15T00:00:00+00:00"


def apply_deterministic_environment(
    *, seed: int = DEFAULT_SEED, epoch: str = DEFAULT_EPOCH
) -> dict[str, str]:
    """Set process environment variables used by reproducible generators.

    ``PYTHONHASHSEED`` must be present before the interpreter starts to affect
    hash randomisation.  The runner therefore also supplies it in the child
    process environment; this function remains useful for NumPy, plotting and
    numeric-library controls after startup.
    """

    env = {
        "PYTHONHASHSEED": str(int(seed)),
        "SOURCE_DATE_EPOCH": str(_epoch_seconds(epoch)),
        "TZ": "UTC",
        "LC_ALL": "C",
        "LANG": "C",
        "MPLBACKEND": "Agg",
        "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
        "OMP_NUM_THREADS": "1",
        "NUMEXPR_NUM_THREADS": "1",
        "HYDROSHEAF_REPRODUCIBLE": "1",
        "HYDROSHEAF_REPRODUCIBLE_SEED": str(int(seed)),
        "HYDROSHEAF_REPRODUCIBLE_EPOCH": epoch,
    }
    os.environ.update(env)
    return env


def seed_all(seed: int = DEFAULT_SEED) -> None:
    """Seed Python and legacy NumPy global RNGs for older scripts."""

    random.seed(int(seed))
    try:
        import numpy as np

        np.random.seed(int(seed))
    except ImportError:  # pragma: no cover - NumPy is a project dependency
        pass


def reproducible_now_iso() -> str:
    """Return a stable UTC timestamp when reproducible mode is enabled.

    Ordinary interactive runs retain wall-clock timestamps.  A reproducible
    run opts in through ``HYDROSHEAF_REPRODUCIBLE=1`` and receives the frozen
    epoch recorded in its run manifest.
    """

    if os.environ.get("HYDROSHEAF_REPRODUCIBLE") == "1":
        return os.environ.get("HYDROSHEAF_REPRODUCIBLE_EPOCH", DEFAULT_EPOCH)
    from datetime import datetime, timezone

    return datetime.now(timezone.utc).isoformat()


def file_sha256(path: Path) -> str:
    """Return the SHA-256 digest of *path* as lowercase hexadecimal."""

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def manifest_for_tree(
    root: Path,
    *,
    exclude_names: Iterable[str] = (),
    exclude_prefixes: Iterable[str] = (),
    include_suffixes: Iterable[str] | None = None,
) -> dict[str, Any]:
    """Build a stable relative-path/size/hash manifest for a file tree."""

    root = root.resolve()
    excluded_names = set(exclude_names)
    included_suffixes = (
        {suffix.lower() for suffix in include_suffixes}
        if include_suffixes is not None
        else None
    )
    excluded_prefixes = tuple(
        prefix.replace("\\", "/").rstrip("/") + "/"
        for prefix in exclude_prefixes
    )
    files: list[dict[str, Any]] = []
    for path in sorted(root.rglob("*")):
        if not path.is_file():
            continue
        if included_suffixes is not None and path.suffix.lower() not in included_suffixes:
            continue
        rel = path.relative_to(root).as_posix()
        if path.name in excluded_names or any(
            rel.startswith(prefix) for prefix in excluded_prefixes
        ):
            continue
        files.append(
            {
                "path": rel,
                "bytes": path.stat().st_size,
                "sha256": file_sha256(path),
            }
        )
    return {
        "schema": "hydrosheaf-output-tree-manifest-v1",
        "files": files,
        "file_count": len(files),
        "total_bytes": sum(int(item["bytes"]) for item in files),
    }


def compare_trees(reference: Mapping[str, Any], candidate: Mapping[str, Any]) -> dict[str, Any]:
    """Compare two manifests without discarding size or hash evidence."""

    def index(value: Mapping[str, Any]) -> dict[str, Mapping[str, Any]]:
        return {str(item["path"]): item for item in value.get("files", [])}

    ref = index(reference)
    got = index(candidate)
    missing = sorted(set(ref) - set(got))
    extra = sorted(set(got) - set(ref))
    changed = []
    for path in sorted(set(ref) & set(got)):
        if ref[path].get("sha256") != got[path].get("sha256"):
            changed.append(
                {
                    "path": path,
                    "reference": ref[path],
                    "candidate": got[path],
                }
            )
    return {
        "exact": not missing and not extra and not changed,
        "missing": missing,
        "extra": extra,
        "changed": changed,
        "reference_file_count": len(ref),
        "candidate_file_count": len(got),
    }


def write_json(path: Path, value: Mapping[str, Any]) -> None:
    """Write stable, UTF-8 JSON with a final newline."""

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True, ensure_ascii=False) + "\n",
        encoding="utf-8",
        newline="\n",
    )


def _epoch_seconds(value: str) -> int:
    from datetime import datetime, timezone

    parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
    if parsed.tzinfo is None:
        parsed = parsed.replace(tzinfo=timezone.utc)
    return int(parsed.timestamp())
