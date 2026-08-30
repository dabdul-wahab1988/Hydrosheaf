"""Tests for benchmark bootstrap uncertainty estimation.

Covers:
- _compute_metrics_from_ids correctness
- bootstrap_benchmark_metrics: CI ordering, delta CIs, probability bounds,
  seed reproducibility, small-set warnings, edge-exclusion rules
- Integration: bootstrap enabled/disabled in run_assumption_benchmark,
  CSV CI columns, MD Uncertainty section, manuscript_claim_allowed unchanged
"""

import json
import os
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from hydrosheaf.calibration.adapters import TopologyCalibrationObservation
from hydrosheaf.calibration.benchmark import run_assumption_benchmark
from hydrosheaf.calibration.benchmark_bootstrap import (
    _compute_metrics_from_ids,
    bootstrap_benchmark_metrics,
)
from hydrosheaf.calibration.config import CalibrationConfig


# ── helpers ──────────────────────────────────────────────────────────

def _write_csv(filepath, rows, columns):
    df = pd.DataFrame(rows, columns=columns)
    df.to_csv(filepath, index=False)


def _make_temp_dir():
    return tempfile.TemporaryDirectory()


def _build_test_files(temp_dir):
    """Same 4-edge topology as test_assumption_benchmark.py."""
    base = Path(temp_dir)

    edges_file = base / "candidate_edges.csv"
    _write_csv(
        edges_file,
        [
            {"edge_id": "AtoB", "u": "A", "v": "B", "edge_confidence": 0.7},
            {"edge_id": "AtoC", "u": "A", "v": "C", "edge_confidence": 0.6},
            {"edge_id": "BtoC", "u": "B", "v": "C", "edge_confidence": 0.3},
            {"edge_id": "BtoD", "u": "B", "v": "D", "edge_confidence": 0.4},
        ],
        ["edge_id", "u", "v", "edge_confidence"],
    )

    cal_labels_file = base / "calibration_labels.csv"
    _write_csv(
        cal_labels_file,
        [
            {"edge_id": "AtoB", "observed_present": 1, "weight": 1.0},
            {"edge_id": "AtoC", "observed_present": 1, "weight": 1.0},
            {"edge_id": "BtoC", "observed_present": 0, "weight": 1.0},
        ],
        ["edge_id", "observed_present", "weight"],
    )

    val_labels_file = base / "validation_labels.csv"
    _write_csv(
        val_labels_file,
        [
            {"edge_id": "AtoB", "observed_present": 1, "weight": 1.0},
            {"edge_id": "BtoC", "observed_present": 0, "weight": 1.0},
            {"edge_id": "BtoD", "observed_present": 1, "weight": 1.0},
        ],
        ["edge_id", "observed_present", "weight"],
    )

    return {
        "edges": str(edges_file),
        "cal_labels": str(cal_labels_file),
        "val_labels": str(val_labels_file),
    }


def _make_config(files, output_dir, extra_settings=None):
    settings = {
        "candidates_file": files["edges"],
        "observations_file": files["cal_labels"],
        "validation_observations_file": files["val_labels"],
        "max_iterations": 5,
        "output_dir": output_dir,
        "engine": "internal",
    }
    if extra_settings:
        settings.update(extra_settings)
    return CalibrationConfig(
        problem_type="topology",
        max_nfev=5,
        output_dir=output_dir,
        engine="internal",
        adapter_settings=settings,
    )


def _make_observations(pairs):
    """Build a list of TopologyCalibrationObservation from (edge_id, present) pairs."""
    return [TopologyCalibrationObservation(eid, float(p)) for eid, p in pairs]


# ── _compute_metrics_from_ids unit tests ─────────────────────────────

def test_compute_metrics_perfect():
    """All selected edges match observed_present=1, all unselected match 0."""
    selected = {"A", "B"}
    metrics = _compute_metrics_from_ids(
        selected,
        edge_ids=["A", "B", "C", "D"],
        observed_present=[1.0, 1.0, 0.0, 0.0],
    )
    assert metrics["precision"] == 1.0
    assert metrics["recall"] == 1.0
    assert metrics["f1"] == 1.0
    assert metrics["accuracy"] == 1.0


def test_compute_metrics_mixed():
    """Known TP/FP/TN/FN values."""
    selected = {"A", "C"}  # C selected but absent -> FP
    metrics = _compute_metrics_from_ids(
        selected,
        edge_ids=["A", "B", "C", "D"],
        observed_present=[1.0, 0.0, 0.0, 1.0],
    )
    # TP: A (present, selected)
    # FP: C (absent, selected)
    # TN: B (absent, not selected)
    # FN: D (present, not selected)
    assert metrics["precision"] == 0.5  # 1 TP / (1+1)
    assert metrics["recall"] == 0.5     # 1 TP / (1+1)
    assert metrics["f1"] == 0.5
    assert metrics["accuracy"] == 0.5  # 2/4 correct


def test_compute_metrics_no_selected():
    """Nothing selected, all presents -> 0 recall, NaN/0 precision."""
    selected = set()
    metrics = _compute_metrics_from_ids(
        selected,
        edge_ids=["A", "B"],
        observed_present=[1.0, 1.0],
    )
    assert metrics["precision"] == 0.0
    assert metrics["recall"] == 0.0
    assert metrics["f1"] == 0.0


def test_compute_metrics_excludes_unlabeled():
    """Selected edge IDs not present in observation list must NOT affect metrics."""
    selected = {"A", "EXTRA"}
    metrics = _compute_metrics_from_ids(
        selected,
        edge_ids=["A", "B"],
        observed_present=[1.0, 0.0],
    )
    # A: present=1, selected -> TP
    # B: present=0, not selected -> TN
    # EXTRA: selected but not in obs -> excluded
    assert metrics["precision"] == 1.0
    assert metrics["recall"] == 1.0
    assert metrics["f1"] == 1.0
    assert metrics["accuracy"] == 1.0


# ── bootstrap_benchmark_metrics unit tests ───────────────────────────

def test_bootstrap_variant_cis_ordering():
    """For each variant, CI low <= CI high for all four metrics."""
    val_obs = _make_observations([("A", 1.0), ("B", 0.0), ("C", 1.0), ("D", 0.0)])
    selected_ids = {
        "baseline": {"A", "B"},
        "assumption_calibrated": {"A", "C"},
    }

    result = bootstrap_benchmark_metrics(
        selected_ids_by_variant=selected_ids,
        val_observations=val_obs,
        n_boot=500,
        seed=42,
    )

    for vname, ci in result["variant_cis"].items():
        for metric in ("precision", "recall", "f1", "accuracy"):
            lo, hi = getattr(ci, metric)
            assert lo <= hi, f"{vname}.{metric}: {lo} > {hi}"
            assert 0.0 <= lo <= 1.0
            assert 0.0 <= hi <= 1.0


def test_bootstrap_delta_cis_ordering():
    """Delta CI bounds must be ordered low <= high."""
    val_obs = _make_observations([("A", 1.0), ("B", 0.0), ("C", 1.0), ("D", 0.0)])
    selected_ids = {
        "baseline": {"A", "B"},
        "assumption_calibrated": {"A", "C"},
    }

    result = bootstrap_benchmark_metrics(
        selected_ids_by_variant=selected_ids,
        val_observations=val_obs,
        n_boot=500,
        seed=42,
    )

    dc = result["delta_cis"]
    for metric in ("delta_precision", "delta_recall", "delta_f1", "delta_accuracy"):
        lo, hi = getattr(dc, metric)
        assert lo <= hi, f"{metric}: {lo} > {hi}"
        assert -1.0 <= lo <= 1.0
        assert -1.0 <= hi <= 1.0


def test_bootstrap_probability_bounds():
    """probability_delta_*_gt_0 must be in [0, 1]."""
    val_obs = _make_observations([("A", 1.0), ("B", 0.0), ("C", 1.0)])
    selected_ids = {
        "baseline": {"A"},
        "assumption_calibrated": {"A", "B"},
    }

    result = bootstrap_benchmark_metrics(
        selected_ids_by_variant=selected_ids,
        val_observations=val_obs,
        n_boot=200,
        seed=7,
    )

    dc = result["delta_cis"]
    assert 0.0 <= dc.probability_delta_f1_gt_0 <= 1.0
    assert 0.0 <= dc.probability_delta_precision_gt_0 <= 1.0


def test_bootstrap_seed_reproducibility():
    """Same seed must produce identical CI bounds."""
    val_obs = _make_observations([("A", 1.0), ("B", 0.0), ("C", 1.0), ("D", 0.0)])
    selected_ids = {
        "baseline": {"A", "C"},
        "assumption_calibrated": {"A", "B", "C"},
    }

    r1 = bootstrap_benchmark_metrics(selected_ids, val_obs, n_boot=300, seed=99)
    r2 = bootstrap_benchmark_metrics(selected_ids, val_obs, n_boot=300, seed=99)

    for vname in r1["variant_cis"]:
        ci1 = r1["variant_cis"][vname]
        ci2 = r2["variant_cis"][vname]
        for metric in ("precision", "recall", "f1", "accuracy"):
            assert getattr(ci1, metric) == getattr(ci2, metric), (
                f"{vname}.{metric} differs"
            )

    dc1 = r1["delta_cis"]
    dc2 = r2["delta_cis"]
    for metric in ("delta_precision", "delta_recall", "delta_f1", "delta_accuracy"):
        assert getattr(dc1, metric) == getattr(dc2, metric), f"{metric} differs"
    assert dc1.probability_delta_f1_gt_0 == dc2.probability_delta_f1_gt_0


def test_bootstrap_small_validation_warning():
    """Fewer than 20 observations must produce a bootstrap_warning."""
    val_obs = _make_observations([("A", 1.0), ("B", 0.0)])
    selected_ids = {"baseline": {"A"}, "assumption_calibrated": {"A", "B"}}

    result = bootstrap_benchmark_metrics(selected_ids, val_obs, n_boot=100, seed=1)
    assert "bootstrap_warning" in result
    assert "unreliable" in result["bootstrap_warning"].lower()


def test_bootstrap_large_validation_no_warning():
    """20+ observations must NOT produce a warning."""
    val_obs = _make_observations(
        [(f"E{i}", 1.0 if i % 2 == 0 else 0.0) for i in range(20)]
    )
    selected_ids = {
        "baseline": {f"E{i}" for i in range(0, 20, 2)},
        "assumption_calibrated": {f"E{i}" for i in range(0, 20, 2)},
    }

    result = bootstrap_benchmark_metrics(
        selected_ids, val_obs, n_boot=100, seed=1,
    )
    assert "bootstrap_warning" not in result


def test_bootstrap_excludes_unlabeled_edges():
    """Selected edges not in validation labels must not affect resampled metrics."""
    # 5 validation edges: A-E, all present
    val_obs = _make_observations([(c, 1.0) for c in "ABCDE"])
    # baseline selects A,B,C but also EXTRA which is not in val_obs
    selected_ids = {
        "baseline": {"A", "B", "C", "EXTRA"},
        "assumption_calibrated": {"A", "B", "C", "EXTRA"},
    }

    result = bootstrap_benchmark_metrics(
        selected_ids, val_obs, n_boot=500, seed=42,
    )

    # baseline: 3 TP (A,B,C), 2 FN (D,E), 0 FP → precision=1.0, recall=0.6
    # This should hold in every bootstrap sample since EXTRA is never in resampled obs
    for vname in ("baseline", "assumption_calibrated"):
        ci = result["variant_cis"][vname]
        # Precision must always be 1.0 (no FPs possible because EXTRA is excluded)
        assert ci.precision == (1.0, 1.0), f"{vname} precision CI: {ci.precision}"


# ── integration tests (full benchmark pipeline) ──────────────────────

def test_benchmark_bootstrap_enabled_adds_uncertainty():
    """When benchmark_bootstrap: true, report must contain 'uncertainty'."""
    with _make_temp_dir() as tmp:
        files = _build_test_files(tmp)
        out_dir = str(Path(tmp) / "out")
        config = _make_config(files, out_dir, extra_settings={
            "benchmark_bootstrap": True,
            "benchmark_bootstrap_n": 200,
            "benchmark_bootstrap_seed": 42,
        })

        report = run_assumption_benchmark(config)

        assert "uncertainty" in report
        u = report["uncertainty"]
        assert u["n_boot"] == 200
        assert u["seed"] == 42
        assert "method" in u
        assert "variant_cis" in u
        assert "improvement_summary_ci" in u

        # All three variants must have CIs
        for vname in ("baseline", "null_model_defaults", "assumption_calibrated"):
            assert vname in u["variant_cis"]
            for m in ("precision", "recall", "f1", "accuracy"):
                ci = u["variant_cis"][vname][m]
                assert len(ci) == 2
                assert ci[0] <= ci[1]

        # Delta CIs must have expected keys
        dci = u["improvement_summary_ci"]
        for key in ("delta_precision", "delta_recall", "delta_f1", "delta_accuracy",
                     "probability_delta_f1_gt_0", "probability_delta_precision_gt_0"):
            assert key in dci


def test_benchmark_bootstrap_disabled_no_uncertainty():
    """Backward compatible: no 'uncertainty' key when bootstrap is disabled."""
    with _make_temp_dir() as tmp:
        files = _build_test_files(tmp)
        out_dir = str(Path(tmp) / "out")
        config = _make_config(files, out_dir)

        report = run_assumption_benchmark(config)
        assert "uncertainty" not in report


def test_benchmark_bootstrap_csv_ci_columns():
    """When bootstrap enabled, CSV must have CI columns with float values."""
    with _make_temp_dir() as tmp:
        files = _build_test_files(tmp)
        out_dir = str(Path(tmp) / "out")
        config = _make_config(files, out_dir, extra_settings={
            "benchmark_bootstrap": True,
            "benchmark_bootstrap_n": 200,
            "benchmark_bootstrap_seed": 99,
        })

        run_assumption_benchmark(config)

        csv_path = Path(out_dir) / "assumption_benchmark_summary.csv"
        assert csv_path.exists()

        with open(csv_path) as f:
            lines = f.readlines()

        header = lines[0].strip().split(",")
        for col in ("precision_low95", "precision_high95", "recall_low95",
                     "recall_high95", "f1_low95", "f1_high95",
                     "accuracy_low95", "accuracy_high95"):
            assert col in header, f"Missing CSV column: {col}"

        # Variant rows must have non-empty CI values
        for line in lines[1:4]:  # 3 variant rows
            parts = line.strip().split(",")
            # CI columns start after accuracy and before tp
            ci_idx = header.index("precision_low95")
            for i in range(ci_idx, ci_idx + 8):
                val = parts[i].strip()
                assert val != "", f"Empty CI value at row col {header[i]}"
                float(val)  # must be parseable

        # Improvement row must also have CI values
        imp_line = lines[4].strip().split(",")
        ci_idx = header.index("precision_low95")
        for i in range(ci_idx, ci_idx + 8):
            val = imp_line[i].strip()
            assert val != "", f"Empty CI value at improvement row col {header[i]}"
            float(val)


def test_benchmark_bootstrap_md_uncertainty_section():
    """When bootstrap enabled, MD report must contain 'Uncertainty' section."""
    with _make_temp_dir() as tmp:
        files = _build_test_files(tmp)
        out_dir = str(Path(tmp) / "out")
        config = _make_config(files, out_dir, extra_settings={
            "benchmark_bootstrap": True,
            "benchmark_bootstrap_n": 200,
            "benchmark_bootstrap_seed": 7,
        })

        run_assumption_benchmark(config)

        md_path = Path(out_dir) / "assumption_benchmark_report.md"
        md_text = md_path.read_text()

        assert "## Uncertainty" in md_text
        assert "Percentile bootstrap" in md_text
        assert "95% Bootstrap CIs" in md_text
        assert "Paired Delta" in md_text
        assert "P(delta > 0)" in md_text
        assert "Bootstrap interval for delta F1" in md_text


def test_benchmark_bootstrap_disabled_md_no_uncertainty():
    """When bootstrap disabled, MD report must NOT contain Uncertainty section."""
    with _make_temp_dir() as tmp:
        files = _build_test_files(tmp)
        out_dir = str(Path(tmp) / "out")
        config = _make_config(files, out_dir)

        run_assumption_benchmark(config)

        md_path = Path(out_dir) / "assumption_benchmark_report.md"
        md_text = md_path.read_text()
        assert "## Uncertainty" not in md_text


def test_manuscript_claim_allowed_unchanged_by_bootstrap():
    """Bootstrap must not alter manuscript_claim_allowed logic."""
    with _make_temp_dir() as tmp:
        files = _build_test_files(tmp)
        out_dir = str(Path(tmp) / "out")

        # Without calibration_result
        config = _make_config(files, out_dir, extra_settings={
            "benchmark_bootstrap": True,
            "benchmark_bootstrap_n": 100,
            "benchmark_bootstrap_seed": 1,
        })
        report = run_assumption_benchmark(config)
        assert report["independent_validation"] is False
        assert report["manuscript_claim_allowed"] is False  # no calibration_result

        # With calibration_result
        config2 = _make_config(files, out_dir + "_2", extra_settings={
            "benchmark_bootstrap": True,
            "benchmark_bootstrap_n": 100,
            "benchmark_bootstrap_seed": 1,
        })
        cal_result = {
            "success": True,
            "phi": 0.01,
            "optimal_parameters": {"edge_logit__AtoB": 1.0},
        }
        report2 = run_assumption_benchmark(config2, calibration_result=cal_result)
        assert report2["independent_validation"] is False
        assert report2["manuscript_claim_allowed"] is False
        assert "uncertainty" in report2  # bootstrap still present


def test_benchmark_bootstrap_metrics_only_validation_edges():
    """Bootstrap CI must reflect the fact that unlabeled edges cannot
    create false positives. The point-estimate precision for baseline
    is 1.0 (AtoB selected=TP, AtoC selected but excluded from val),
    and the precision CI upper bound must be 1.0 since FPs are
    structurally impossible when selected edges not in validation
    labels are excluded."""
    with _make_temp_dir() as tmp:
        files = _build_test_files(tmp)
        out_dir = str(Path(tmp) / "out")
        config = _make_config(files, out_dir, extra_settings={
            "benchmark_bootstrap": True,
            "benchmark_bootstrap_n": 500,
            "benchmark_bootstrap_seed": 42,
        })

        report = run_assumption_benchmark(config)

        # Baseline selects AtoB + AtoC (edge_confidence 0.7, 0.6)
        # AtoC is NOT in validation labels → excluded, never creates FP
        baseline_pt = report["variants"]["baseline"]
        assert baseline_pt["metrics"]["precision"] == 1.0
        assert baseline_pt["confusion_matrix"]["false_positives"] == 0

        # Precision CI upper bound must be 1.0: FPs structurally impossible
        baseline_ci = report["uncertainty"]["variant_cis"]["baseline"]
        assert baseline_ci["precision"][1] == pytest.approx(1.0)
