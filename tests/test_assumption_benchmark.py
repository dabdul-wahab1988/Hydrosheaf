"""Integration tests for the assumption calibration benchmark workflow.

Tests:
- Benchmark evaluates only validation-labeled edges (unlabeled selected edges excluded)
- All three variants (baseline, null_model_defaults, assumption_calibrated) in report
- Improvement deltas computed correctly
- manuscript_claim_allowed matches independent_validation
- Benchmark reuses calibration_result without re-running calibration
- CSV, JSON, and MD report files are written
"""

import json
import os
import tempfile
from pathlib import Path

import pandas as pd
import pytest

from hydrosheaf.calibration.benchmark import (
    _evaluate_variant,
    _write_benchmark_csv,
    run_assumption_benchmark,
)
from hydrosheaf.calibration.config import CalibrationConfig
from hydrosheaf.calibration.validation_workflow import (
    _load_topology_observations,
    _resolve_path,
)
from hydrosheaf.config import Config as HConfig
from hydrosheaf.graph.types import Edge


# ── helpers ──────────────────────────────────────────────────────────

def _write_csv(filepath, rows, columns):
    """Write a list of dicts to CSV."""
    df = pd.DataFrame(rows, columns=columns)
    df.to_csv(filepath, index=False)


def _make_temp_dir():
    return tempfile.TemporaryDirectory()


def _build_test_files(temp_dir):
    """Build candidate edges, calibration labels, and validation labels
    CSV files for a 4-edge topology."""
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
    """Build a CalibrationConfig for topology benchmark."""
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


# ── _evaluate_variant unit tests ─────────────────────────────────────

def test_evaluate_variant_metrics_only_validation_edges():
    """Metrics must only be computed on validation-labeled edges.
    Selected edges absent from val_observations must not contribute
    to FP or TN."""
    from hydrosheaf.calibration.adapters import TopologyCalibrationObservation

    candidates = [
        Edge(edge_id="A_B", u="A", v="B", attrs={"edge_confidence": 0.7}),
        Edge(edge_id="X_Y", u="X", v="Y", attrs={"edge_confidence": 0.6}),
        Edge(edge_id="C_D", u="C", v="D", attrs={"edge_confidence": 0.3}),
    ]
    val_obs = [
        TopologyCalibrationObservation("A_B", 1.0),
        TopologyCalibrationObservation("C_D", 0.0),
    ]
    cfg = HConfig()

    result = _evaluate_variant(None, candidates, cfg, val_obs, "test")

    # A_B selected (edge_confidence 0.7 >= 0.5) → TP (present=1, selected)
    # C_D NOT selected (edge_confidence 0.3 < 0.5) → TN (present=0, not selected)
    # X_Y selected but NOT in val_obs → excluded
    cm = result["confusion_matrix"]
    assert cm["true_positives"] == 1
    assert cm["false_positives"] == 0
    assert cm["true_negatives"] == 1
    assert cm["false_negatives"] == 0
    assert result["n_evaluated_validation_labels"] == 2
    assert result["n_selected_edges"] == 2  # A_B and X_Y
    assert result["metrics"]["precision"] == 1.0
    assert result["metrics"]["recall"] == 1.0
    assert result["metrics"]["f1"] == 1.0


def test_evaluate_variant_no_samples_falls_back_to_threshold():
    """When samples is None, _evaluate_variant thresholds on edge_confidence >= 0.5."""
    from hydrosheaf.calibration.adapters import TopologyCalibrationObservation

    candidates = [
        Edge(edge_id="E1", u="A", v="B", attrs={"edge_confidence": 0.9}),
        Edge(edge_id="E2", u="A", v="C", attrs={"edge_confidence": 0.2}),
    ]
    val_obs = [
        TopologyCalibrationObservation("E1", 1.0),
        TopologyCalibrationObservation("E2", 0.0),
    ]
    cfg = HConfig()

    result = _evaluate_variant(None, candidates, cfg, val_obs, "test")
    cm = result["confusion_matrix"]
    assert cm["true_positives"] == 1   # E1 selected, present=1
    assert cm["true_negatives"] == 1   # E2 not selected, present=0
    assert cm["false_positives"] == 0
    assert cm["false_negatives"] == 0


# ── benchmark integration tests ──────────────────────────────────────

def test_benchmark_evaluates_only_validation_edges():
    """Benchmark must compute metrics on validation-labeled edges only.
    A->C is selected in baseline (edge_confidence=0.6) but absent from
    validation labels → must not create a FP."""
    with _make_temp_dir() as tmp:
        files = _build_test_files(tmp)
        out_dir = str(Path(tmp) / "out")
        config = _make_config(files, out_dir)

        report = run_assumption_benchmark(config)

        # Validation labels: AtoB=1, BtoC=0, BtoD=1
        # Baseline: AtoB + AtoC selected (0.7, 0.6 >= 0.5)
        #   BtoC (0.3) and BtoD (0.4) not selected
        #   On val labels: AtoB → TP, BtoC → TN, BtoD → FN
        #   AtoC excluded → FP=0
        baseline = report["variants"]["baseline"]
        cm = baseline["confusion_matrix"]
        assert cm["false_positives"] == 0
        assert cm["true_positives"] == 1
        assert cm["true_negatives"] == 1
        assert cm["false_negatives"] == 1

        # All variants should evaluate exactly 3 validation labels
        for vname in ("baseline", "null_model_defaults", "assumption_calibrated"):
            assert report["variants"][vname]["n_evaluated_validation_labels"] == 3


def test_benchmark_all_three_variants_present():
    """Report must include baseline, null_model_defaults, and assumption_calibrated."""
    with _make_temp_dir() as tmp:
        files = _build_test_files(tmp)
        out_dir = str(Path(tmp) / "out")
        config = _make_config(files, out_dir)

        report = run_assumption_benchmark(config)

        assert set(report["variants"].keys()) == {
            "baseline", "null_model_defaults", "assumption_calibrated",
        }
        for vname in ("baseline", "null_model_defaults", "assumption_calibrated"):
            v = report["variants"][vname]
            assert "metrics" in v
            assert "confusion_matrix" in v
            assert "n_selected_edges" in v


def test_benchmark_improvement_deltas():
    """delta values = assumption_calibrated metric - baseline metric."""
    with _make_temp_dir() as tmp:
        files = _build_test_files(tmp)
        out_dir = str(Path(tmp) / "out")
        config = _make_config(files, out_dir)

        report = run_assumption_benchmark(config)

        imp = report["improvement_summary"]
        baseline_m = report["variants"]["baseline"]["metrics"]
        calibrated_m = report["variants"]["assumption_calibrated"]["metrics"]

        for key in ("precision", "recall", "f1", "accuracy"):
            expected = round(calibrated_m[key] - baseline_m[key], 6)
            assert imp[f"delta_{key}"] == pytest.approx(expected, abs=1e-6)


def test_benchmark_manuscript_claim_allowed():
    """manuscript_claim_allowed requires independent_validation AND actual
    calibration parameters. Without calibration_result, it must be False."""
    with _make_temp_dir() as tmp:
        files = _build_test_files(tmp)
        out_dir = str(Path(tmp) / "out")
        config = _make_config(files, out_dir)

        # No calibration_result → manuscript_claim_allowed is False
        report = run_assumption_benchmark(config)
        # The fixture reuses edge IDs across calibration and validation.  A
        # different filename is not an independent validation unit.
        assert report["independent_validation"] is False
        assert "overlapping edge IDs" in report["independence_reason"]
        assert report["manuscript_claim_allowed"] is False

        # With calibration_result containing optimal_parameters → True
        cal_result = {
            "success": True,
            "phi": 0.01,
            "optimal_parameters": {"edge_logit__AtoB": 1.0},
        }
        report2 = run_assumption_benchmark(config, calibration_result=cal_result)
        assert report2["independent_validation"] is False
        assert report2["manuscript_claim_allowed"] is False


def test_benchmark_same_file_raises():
    """Passing the same file for calibration and validation must raise ValueError."""
    with _make_temp_dir() as tmp:
        base = Path(tmp)
        edges_file = base / "edges.csv"
        labels_file = base / "labels.csv"
        _write_csv(
            edges_file,
            [{"edge_id": "A_B", "u": "A", "v": "B", "edge_confidence": 0.5}],
            ["edge_id", "u", "v", "edge_confidence"],
        )
        _write_csv(
            labels_file,
            [{"edge_id": "A_B", "observed_present": 1, "weight": 1.0}],
            ["edge_id", "observed_present", "weight"],
        )
        out_dir = str(base / "out")
        config = CalibrationConfig(
            problem_type="topology",
            max_nfev=5,
            output_dir=out_dir,
            engine="internal",
            adapter_settings={
                "candidates_file": str(edges_file),
                "observations_file": str(labels_file),
                "validation_observations_file": str(labels_file),
            },
        )
        with pytest.raises(ValueError, match="must be different"):
            run_assumption_benchmark(config)


def test_benchmark_missing_validation_file_raises():
    """Benchmark should raise when validation_observations_file is missing."""
    with _make_temp_dir() as tmp:
        base = Path(tmp)
        edges_file = base / "edges.csv"
        labels_file = base / "cal_labels.csv"
        _write_csv(
            edges_file,
            [{"edge_id": "A_B", "u": "A", "v": "B", "edge_confidence": 0.5}],
            ["edge_id", "u", "v", "edge_confidence"],
        )
        _write_csv(
            labels_file,
            [{"edge_id": "A_B", "observed_present": 1, "weight": 1.0}],
            ["edge_id", "observed_present", "weight"],
        )
        out_dir = str(base / "out")
        config = CalibrationConfig(
            problem_type="topology",
            max_nfev=5,
            output_dir=out_dir,
            engine="internal",
            adapter_settings={
                "candidates_file": str(edges_file),
                "observations_file": str(labels_file),
            },
        )
        with pytest.raises(ValueError, match="validation_observations_file"):
            run_assumption_benchmark(config)


def test_benchmark_reuses_calibration_result():
    """When calibration_result is provided, the calibrated variant must
    use its optimal parameters without running calibration."""
    with _make_temp_dir() as tmp:
        files = _build_test_files(tmp)
        out_dir = str(Path(tmp) / "out")
        config = _make_config(files, out_dir)

        # Build a calibration_result that boosts BtoD confidence
        # so it goes from 0.4 (unselected) to > 0.5 (selected)
        calibration_result = {
            "success": True,
            "phi": 0.01,
            "optimal_parameters": {
                "edge_logit__AtoB": 2.0,   # high confidence → selected
                "edge_logit__AtoC": 1.0,   # moderate → selected
                "edge_logit__BtoC": -1.0,  # low → not selected
                "edge_logit__BtoD": 1.0,   # boosted above threshold → selected
            },
        }

        report = run_assumption_benchmark(
            config, calibration_result=calibration_result,
        )

        calibrated = report["variants"]["assumption_calibrated"]
        cm = calibrated["confusion_matrix"]

        # With BtoD boosted (selected) and in validation as present:
        # AtoB → TP, BtoC → TN, BtoD → TP
        # AtoC selected but not in validation → excluded
        assert cm["true_positives"] == 2
        assert cm["false_positives"] == 0
        assert cm["true_negatives"] == 1
        assert cm["false_negatives"] == 0

        assert calibrated["metrics"]["precision"] == 1.0
        assert calibrated["metrics"]["recall"] == 1.0
        assert calibrated["metrics"]["f1"] == 1.0


def test_benchmark_writes_output_files():
    """Benchmark must write JSON, CSV, and MD report files."""
    with _make_temp_dir() as tmp:
        files = _build_test_files(tmp)
        out_dir = str(Path(tmp) / "out")
        config = _make_config(files, out_dir)

        run_assumption_benchmark(config)

        json_path = Path(out_dir) / "assumption_benchmark_results.json"
        csv_path = Path(out_dir) / "assumption_benchmark_summary.csv"
        md_path = Path(out_dir) / "assumption_benchmark_report.md"

        assert json_path.exists()
        assert csv_path.exists()
        assert md_path.exists()

        # JSON must be parseable and contain expected fields
        with open(json_path) as f:
            data = json.load(f)
        assert "variants" in data
        assert "improvement_summary" in data
        assert "manuscript_claim_allowed" in data

        # CSV must have header and at least 4 rows (3 variants + improvement)
        with open(csv_path) as f:
            lines = f.readlines()
        assert len(lines) >= 5  # header + 3 variants + improvement row
        assert "variant" in lines[0]

        # MD must reference the variants
        md_text = md_path.read_text()
        assert "baseline" in md_text
        assert "assumption_calibrated" in md_text
        assert "Improvement" in md_text


def test_benchmark_json_structure():
    """Verify the benchmark JSON report contains all top-level fields."""
    with _make_temp_dir() as tmp:
        files = _build_test_files(tmp)
        out_dir = str(Path(tmp) / "out")
        config = _make_config(files, out_dir)

        report = run_assumption_benchmark(config)

        for key in (
            "calibration_label_file",
            "validation_label_file",
            "independent_validation",
            "manuscript_claim_allowed",
            "n_calibration_labels",
            "n_validation_labels",
            "n_overlapping_edge_ids",
            "calibrated_assumption_parameters",
            "fixed_config_only_thresholds",
            "variants",
            "improvement_summary",
        ):
            assert key in report, f"Missing key: {key}"

        # Check improvement summary keys
        for key in ("delta_precision", "delta_recall", "delta_f1", "delta_accuracy"):
            assert key in report["improvement_summary"]


def test_benchmark_variant_selected_edge_ids():
    """Each variant should report its selected_edge_ids."""
    with _make_temp_dir() as tmp:
        files = _build_test_files(tmp)
        out_dir = str(Path(tmp) / "out")
        config = _make_config(files, out_dir)

        report = run_assumption_benchmark(config)

        for vname, vdata in report["variants"].items():
            assert "selected_edge_ids" in vdata
            assert isinstance(vdata["selected_edge_ids"], list)
            assert vdata["n_selected_edges"] == len(vdata["selected_edge_ids"])


def test_benchmark_overlapping_edge_ids_tracked():
    """n_overlapping_edge_ids should count edges in both cal and val files."""
    with _make_temp_dir() as tmp:
        files = _build_test_files(tmp)
        out_dir = str(Path(tmp) / "out")
        config = _make_config(files, out_dir)

        report = run_assumption_benchmark(config)

        # cal has: AtoB, AtoC, BtoC (3 edges)
        # val has: AtoB, BtoC, BtoD (3 edges)
        # overlap: AtoB, BtoC (2 edges)
        assert report["n_calibration_labels"] == 3
        assert report["n_validation_labels"] == 3
        assert report["n_overlapping_edge_ids"] == 2


def test_benchmark_with_config_only_thresholds():
    """Config-only evidence thresholds should appear in
    fixed_config_only_thresholds, not calibrated_assumption_parameters."""
    with _make_temp_dir() as tmp:
        files = _build_test_files(tmp)
        out_dir = str(Path(tmp) / "out")
        extra = {
            "evidence_threshold_probable": 0.75,
            "evidence_threshold_validated": 0.95,
        }
        config = _make_config(files, out_dir, extra_settings=extra)

        report = run_assumption_benchmark(config)

        fixed = report["fixed_config_only_thresholds"]
        assert fixed["evidence_threshold_probable"] == 0.75
        assert fixed["evidence_threshold_validated"] == 0.95

        cal_params = report["calibrated_assumption_parameters"]
        assert "evidence_threshold_probable" not in cal_params
        assert "evidence_threshold_validated" not in cal_params
