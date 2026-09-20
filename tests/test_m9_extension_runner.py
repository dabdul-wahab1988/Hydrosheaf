"""Tests for the HydroSheaf M9.1–M9.3 master benchmark runner and artifact generation."""

from __future__ import annotations

import json
from pathlib import Path
import pytest


def test_m9_extension_runner_artifacts_exist():
    out_dir = Path(".codex_work/m9_extension_verification_v1")
    required_files = [
        "run_manifest.json",
        "claim_readiness.json",
        "generator_provenance.json",
        "inference_records.json",
        "truth_scoring_only.json",
        "loss_decomposition.json",
        "dynamic_kernel_summaries.json",
        "per_tracer_results.json",
        "held_out_metrics.json",
        "abstention_diagnostics.json",
        "artifact_hashes.json",
    ]

    for fname in required_files:
        fpath = out_dir / fname
        assert fpath.exists(), f"Missing required runner output: {fname}"
        assert fpath.stat().st_size > 0, f"Empty output artifact: {fname}"

    # Verify claim readiness boundaries
    with open(out_dir / "claim_readiness.json", "r", encoding="utf-8") as f:
        readiness = json.load(f)

    assert readiness["controlled_synthetic_ready"] is True
    assert readiness["field_validation_claimed"] is False
    assert readiness["universal_superiority_claimed"] is False
    assert "Not field validation" in readiness["claim_boundary"]

    # Verify manifest
    with open(out_dir / "run_manifest.json", "r", encoding="utf-8") as f:
        manifest = json.load(f)

    assert manifest["status"] == "COMPLETED"
    assert len(manifest["artifact_hashes"]) >= 9
