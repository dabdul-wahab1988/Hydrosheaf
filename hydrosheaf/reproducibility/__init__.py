"""Reproducible-output helpers for HydroSheaf research runs.

The public API is deliberately small.  Output-producing scripts should use
the helpers here for run metadata, hashing, and deterministic process setup;
the orchestration policy lives in :mod:`scripts.reproduce_outputs`.
"""

from .core import (
    DEFAULT_EPOCH,
    DEFAULT_SEED,
    apply_deterministic_environment,
    compare_trees,
    file_sha256,
    manifest_for_tree,
    reproducible_now_iso,
    seed_all,
    write_json,
)

__all__ = [
    "DEFAULT_EPOCH",
    "DEFAULT_SEED",
    "apply_deterministic_environment",
    "compare_trees",
    "file_sha256",
    "manifest_for_tree",
    "reproducible_now_iso",
    "seed_all",
    "write_json",
]
