"""Integration tests for the assumption calibration validation workflow.

Tests:
- Calibration and validation label files must be distinct
- Overlapping edge IDs are not independent without an explicit grouped contract
- Validation metrics are computed correctly from known edge selections
- Results JSON includes calibrated params, validation metrics, independent_validation
- Config-only thresholds are reported but not listed as calibrated params
- Backward compatibility: workflow raises clearly when no validation file is present
- Invalid assumption_params still fail fast (delegated to TopologyCalibrationAdapter)
"""

import json
import os
import tempfile
from pathlib import Path

import pandas as pd
import pytest

from hydrosheaf.calibration.config import CalibrationConfig
from hydrosheaf.calibration.validation_workflow import (
    _extract_assumption_params_from_optimal,
    _extract_config_only_thresholds,
    _load_topology_observations,
    _resolve_path,
    run_assumption_calibration_validation,
)
from hydrosheaf.config import Config as HConfig


# ── helpers ──────────────────────────────────────────────────────────

def _write_csv(filepath, rows, columns):
    """Write a list of dicts to CSV."""
    df = pd.DataFrame(rows, columns=columns)
    df.to_csv(filepath, index=False)


def _make_temp_dir():
    return tempfile.TemporaryDirectory()


def _build_test_files(temp_dir):
    """Build candidate edges, calibration labels, and validation labels
    CSV files for a 4-edge topology. Returns a dict of file paths."""
    base = Path(temp_dir)

    # 4 candidate edges: A->B (high prior), A->C (medium), B->C (low), B->D (low)
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

    # Calibration labels: only A->B and A->C are present, B->C is absent
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

    # Validation labels: A->B present (agrees), B->C absent (agrees),
    # B->D present (model never saw this edge in cal → tests generalisation)
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
    """Build a CalibrationConfig for topology calibration with validation."""
    settings = {
        "candidates_file": files["edges"],
        "observations_file": files["cal_labels"],
        "validation_observations_file": files["val_labels"],
        "max_iterations": 60,
        "output_dir": output_dir,
        "engine": "internal",
    }
    if extra_settings:
        settings.update(extra_settings)
    return CalibrationConfig(
        problem_type="topology",
        max_nfev=60,
        output_dir=output_dir,
        engine="internal",
        adapter_settings=settings,
    )


# ── helper function tests ────────────────────────────────────────────

def test_resolve_path_normalises_abs():
    """_resolve_path should return an absolute normalised path."""
    cwd = os.getcwd()
    result = _resolve_path("relative/path/file.csv")
    assert os.path.isabs(result)
    assert result == os.path.abspath(os.path.join(cwd, "relative/path/file.csv"))


def test_resolve_path_none_returns_empty():
    assert _resolve_path(None) == ""


def test_load_topology_observations():
    with _make_temp_dir() as tmp:
        p = Path(tmp) / "obs.csv"
        _write_csv(
            p,
            [
                {"edge_id": "X_Y", "observed_present": 1, "weight": 1.0},
                {"edge_id": "Z_W", "observed_present": 0, "weight": 0.5},
            ],
            ["edge_id", "observed_present", "weight"],
        )
        obs = _load_topology_observations(str(p))
        assert len(obs) == 2
        assert obs[0].edge_id == "X_Y"
        assert obs[0].observed_present == 1.0
        assert obs[0].weight == 1.0
        assert obs[1].edge_id == "Z_W"
        assert obs[1].observed_present == 0.0
        assert obs[1].weight == 0.5


def test_load_topology_observations_skips_nan_edge_id():
    with _make_temp_dir() as tmp:
        p = Path(tmp) / "obs.csv"
        _write_csv(
            p,
            [
                {"edge_id": "", "observed_present": 1, "weight": 1.0},
                {"edge_id": "valid", "observed_present": 1, "weight": 1.0},
            ],
            ["edge_id", "observed_present", "weight"],
        )
        obs = _load_topology_observations(str(p))
        assert len(obs) == 1
        assert obs[0].edge_id == "valid"


def test_load_topology_observations_accepts_present_column():
    """Should accept 'present' as an alternative to 'observed_present'."""
    with _make_temp_dir() as tmp:
        p = Path(tmp) / "obs.csv"
        _write_csv(
            p,
            [
                {"edge_id": "A_B", "present": 1, "weight": 2.0},
                {"edge_id": "C_D", "present": 0, "weight": 1.0},
            ],
            ["edge_id", "present", "weight"],
        )
        obs = _load_topology_observations(str(p))
        assert len(obs) == 2
        assert obs[0].observed_present == 1.0
        assert obs[1].observed_present == 0.0


def test_load_topology_observations_accepts_value_column():
    """Should accept 'value' as an alternative to 'observed_present'."""
    with _make_temp_dir() as tmp:
        p = Path(tmp) / "obs.csv"
        _write_csv(
            p,
            [
                {"edge_id": "X_Y", "value": 1.0, "weight": 1.0},
                {"edge_id": "Z_W", "value": 0.0, "weight": 1.0},
            ],
            ["edge_id", "value", "weight"],
        )
        obs = _load_topology_observations(str(p))
        assert len(obs) == 2
        assert obs[0].observed_present == 1.0
        assert obs[1].observed_present == 0.0


def test_extract_assumption_params_filters_correctly():
    optimal = {
        "edge_logit__A_B": 2.0,
        "null_model_weight": 0.8,
        "sheaf_weight_isotope": 3.0,
        "other_param": 1.0,
    }
    result = _extract_assumption_params_from_optimal(
        optimal, ["null_model_weight", "sheaf_weight_isotope"],
    )
    assert result == {"null_model_weight": 0.8, "sheaf_weight_isotope": 3.0}
    assert "edge_logit__A_B" not in result
    assert "other_param" not in result


def test_extract_config_only_thresholds():
    cfg = HConfig()
    cfg.evidence_threshold_probable = 0.7
    cfg.evidence_threshold_validated = 0.9
    result = _extract_config_only_thresholds(cfg)
    assert result["evidence_threshold_probable"] == 0.7
    assert result["evidence_threshold_validated"] == 0.9
    # These should NOT include calibratable assumption params
    assert "null_model_weight" not in result
    assert "sheaf_weight_isotope" not in result


# ── validation workflow tests ────────────────────────────────────────

def test_validation_workflow_raises_on_same_file():
    """Passing the same file for calibration and validation labels must
    raise ValueError before any calibration runs."""
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
            run_assumption_calibration_validation(config)


def test_validation_workflow_missing_validation_file_raises():
    """Validation workflow should raise when validation_observations_file
    is not provided."""
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
                # No validation_observations_file
            },
        )
        with pytest.raises(ValueError, match="validation_observations_file"):
            run_assumption_calibration_validation(config)


def test_validation_workflow_metrics():
    """End-to-end test: calibrate on one label set, validate on another,
    verify metrics match expected values.

    Setup:
      Calibration labels:  A->B=1,  A->C=1,  B->C=0  (3 labels)
      Validation labels:   A->B=1,  B->C=0,  B->D=1  (3 labels, B->D novel)

    After calibration:
      - A->B  (prior 0.7, cal=1) → selected
      - A->C  (prior 0.6, cal=1) → selected
      - B->C  (prior 0.3, cal=0) → NOT selected
      - B->D  (prior 0.4, no cal) → NOT selected (p ~ 0.4 < 0.5)

    Metrics are computed only on edges in the validation labels file
    (A->C is excluded because it does not appear there):

      Validation file:    A->B=1       B->C=0       B->D=1
      Selected?           yes (TP)     no (TN)      no (FN)
      TP=1  FP=0  TN=1  FN=1
      precision=1.0  recall=0.5  f1≈0.667  accuracy=2/3≈0.667
    """
    with _make_temp_dir() as tmp:
        files = _build_test_files(tmp)
        out_dir = str(Path(tmp) / "out")
        config = _make_config(files, out_dir)

        report = run_assumption_calibration_validation(config)

        metrics = report["validation_metrics"]
        assert metrics["precision"] == pytest.approx(1.0, abs=0.1)
        assert metrics["recall"] == pytest.approx(0.5, abs=0.1)
        assert metrics["f1"] == pytest.approx(2.0 / 3.0, abs=0.1)
        assert metrics["accuracy"] == pytest.approx(2.0 / 3.0, abs=0.1)

        # FP must be 0 because A->C (selected but not in validation) is excluded
        cm = report["validation_confusion_matrix"]
        assert cm["false_positives"] == 0
        assert cm["true_positives"] == 1
        assert cm["true_negatives"] == 1
        assert cm["false_negatives"] == 1

        # Confusion matrix should have sensible counts
        cm = report["validation_confusion_matrix"]
        assert cm["true_positives"] >= 0
        assert cm["false_positives"] >= 0
        assert cm["false_negatives"] >= 0
        assert cm["true_negatives"] >= 0


def test_validation_report_json_structure():
    """Verify the report JSON contains all required fields with correct types."""
    with _make_temp_dir() as tmp:
        files = _build_test_files(tmp)
        out_dir = str(Path(tmp) / "out")
        config = _make_config(files, out_dir)

        report = run_assumption_calibration_validation(config)

        # Required top-level fields
        assert "calibrated_assumption_parameters" in report
        assert "fixed_config_only_thresholds" in report
        assert "selected_edges_all" in report
        assert "evaluated_validation_edge_ids" in report
        assert "validation_confusion_matrix" in report
        assert "validation_metrics" in report
        assert "calibration_label_file" in report
        assert "validation_label_file" in report
        assert "independent_validation" in report
        assert "independence_reason" in report
        assert "n_calibration_labels" in report
        assert "n_validation_labels" in report
        assert "n_overlapping_edge_ids" in report
        assert "calibration_phi" in report
        assert "calibration_success" in report

        # Separate filenames do not make overlapping labelled edges independent.
        assert report["independent_validation"] is False
        assert "overlapping edge IDs" in report["independence_reason"]

        # Label counts must be correct
        assert report["n_calibration_labels"] == 3
        assert report["n_validation_labels"] == 3

        # Validation metrics must have expected keys
        for key in ("precision", "recall", "f1", "accuracy"):
            assert key in report["validation_metrics"]

        # Confusion matrix fields
        for key in ("true_positives", "false_positives",
                     "true_negatives", "false_negatives"):
            assert key in report["validation_confusion_matrix"]

        # Report file should exist on disk
        report_path = Path(out_dir) / "assumption_validation_results.json"
        assert report_path.exists()
        with open(report_path) as f:
            disk_report = json.load(f)
        assert disk_report["independent_validation"] is False


def test_config_only_thresholds_in_report_not_calibrated():
    """evidence_threshold_probable / _validated should appear in
    fixed_config_only_thresholds, not in calibrated_assumption_parameters."""
    with _make_temp_dir() as tmp:
        files = _build_test_files(tmp)
        out_dir = str(Path(tmp) / "out")
        extra = {
            "evidence_threshold_probable": 0.7,
            "evidence_threshold_validated": 0.9,
        }
        config = _make_config(files, out_dir, extra_settings=extra)

        report = run_assumption_calibration_validation(config)

        fixed = report["fixed_config_only_thresholds"]
        assert fixed["evidence_threshold_probable"] == 0.7
        assert fixed["evidence_threshold_validated"] == 0.9

        calibrated = report["calibrated_assumption_parameters"]
        assert "evidence_threshold_probable" not in calibrated
        assert "evidence_threshold_validated" not in calibrated


def test_validation_workflow_emits_report_file():
    """Verify the report JSON file is written to the output directory."""
    with _make_temp_dir() as tmp:
        files = _build_test_files(tmp)
        out_dir = str(Path(tmp) / "out")
        config = _make_config(files, out_dir)

        run_assumption_calibration_validation(config)

        report_path = Path(out_dir) / "assumption_validation_results.json"
        assert report_path.exists()
        with open(report_path) as f:
            data = json.load(f)
        assert "validation_metrics" in data


def test_invalid_assumption_params_still_fail_fast():
    """TopologyCalibrationAdapter.__init__ raises ValueError for unknown
    assumption_params. This is tested via the adapter directly rather than
    the full workflow."""
    from hydrosheaf.calibration.adapters import (
        TopologyCalibrationAdapter,
        TopologyCalibrationObservation,
    )
    from hydrosheaf.graph.types import Edge

    edge = Edge(edge_id="X_Y", u="X", v="Y", attrs={"edge_confidence": 0.6})
    with pytest.raises(ValueError, match="Unknown assumption_params"):
        TopologyCalibrationAdapter(
            candidate_edges=[edge],
            observations=[TopologyCalibrationObservation("X_Y", 1.0)],
            assumption_params=["nonexistent_param"],
        )


def test_validation_workflow_overlapping_edge_ids_warning():
    """Reversed labels in different files must not pass as independent."""
    with _make_temp_dir() as tmp:
        base = Path(tmp)
        edges_file = base / "edges.csv"
        _write_csv(
            edges_file,
            [
                {"edge_id": "A_B", "u": "A", "v": "B", "edge_confidence": 0.7},
                {"edge_id": "C_D", "u": "C", "v": "D", "edge_confidence": 0.5},
            ],
            ["edge_id", "u", "v", "edge_confidence"],
        )
        cal_file = base / "cal.csv"
        _write_csv(
            cal_file,
            [
                {"edge_id": "A_B", "observed_present": 1, "weight": 1.0},
                {"edge_id": "C_D", "observed_present": 0, "weight": 1.0},
            ],
            ["edge_id", "observed_present", "weight"],
        )
        val_file = base / "val.csv"
        _write_csv(
            val_file,
            [
                {"edge_id": "A_B", "observed_present": 0, "weight": 1.0},
                {"edge_id": "C_D", "observed_present": 1, "weight": 1.0},
            ],
            ["edge_id", "observed_present", "weight"],
        )
        out_dir = str(base / "out")
        config = CalibrationConfig(
            problem_type="topology",
            max_nfev=30,
            output_dir=out_dir,
            engine="internal",
            adapter_settings={
                "candidates_file": str(edges_file),
                "observations_file": str(cal_file),
                "validation_observations_file": str(val_file),
            },
        )

        report = run_assumption_calibration_validation(config)
        # Both files share all edge IDs → overlap = 2
        assert report["n_overlapping_edge_ids"] == 2
        assert report["independent_validation"] is False
        assert "reversed labels" in report["independence_reason"]


def test_validation_workflow_non_overlapping_edge_ids_remains_independent():
    """A genuinely disjoint validation label set remains independent."""
    with _make_temp_dir() as tmp:
        files = _build_test_files(tmp)
        base = Path(tmp)
        disjoint_val_file = base / "validation_labels_disjoint.csv"
        _write_csv(
            disjoint_val_file,
            [{"edge_id": "BtoD", "observed_present": 1, "weight": 1.0}],
            ["edge_id", "observed_present", "weight"],
        )
        files["val_labels"] = str(disjoint_val_file)

        report = run_assumption_calibration_validation(
            _make_config(files, str(base / "out")),
        )

        assert report["n_overlapping_edge_ids"] == 0
        assert report["independent_validation"] is True
        assert "disjoint" in report["independence_reason"]


def test_validation_workflow_explicit_grouped_contract_allows_overlap():
    """A verified disjoint group assignment is the explicit overlap exception."""
    with _make_temp_dir() as tmp:
        files = _build_test_files(tmp)
        base = Path(tmp)
        grouped_cal_file = base / "calibration_grouped.csv"
        grouped_val_file = base / "validation_grouped.csv"
        grouped_columns = [
            "edge_id", "observed_present", "weight", "split_group",
        ]
        _write_csv(
            grouped_cal_file,
            [
                {"edge_id": "AtoB", "observed_present": 1, "weight": 1.0,
                 "split_group": "calibration"},
                {"edge_id": "AtoC", "observed_present": 1, "weight": 1.0,
                 "split_group": "calibration"},
                {"edge_id": "BtoC", "observed_present": 0, "weight": 1.0,
                 "split_group": "calibration"},
            ],
            grouped_columns,
        )
        _write_csv(
            grouped_val_file,
            [
                {"edge_id": "AtoB", "observed_present": 0, "weight": 1.0,
                 "split_group": "validation"},
                {"edge_id": "AtoC", "observed_present": 0, "weight": 1.0,
                 "split_group": "validation"},
                {"edge_id": "BtoC", "observed_present": 1, "weight": 1.0,
                 "split_group": "validation"},
            ],
            grouped_columns,
        )
        files["cal_labels"] = str(grouped_cal_file)
        files["val_labels"] = str(grouped_val_file)
        contract = {
            "kind": "grouped",
            "group_column": "split_group",
            "calibration_groups": ["calibration"],
            "validation_groups": ["validation"],
            "allow_overlapping_edge_ids": True,
        }

        report = run_assumption_calibration_validation(
            _make_config(
                files,
                str(base / "out"),
                extra_settings={
                    "grouped_independence_contract": contract,
                },
            ),
        )

        assert report["n_overlapping_edge_ids"] == 3
        assert report["independent_validation"] is True
        assert "explicit grouped_independence_contract" in report[
            "independence_reason"
        ]


def test_validation_workflow_no_empty_cal_observations():
    """When calibration observations are empty, a clear error is raised."""
    with _make_temp_dir() as tmp:
        base = Path(tmp)
        edges_file = base / "edges.csv"
        _write_csv(
            edges_file,
            [{"edge_id": "A_B", "u": "A", "v": "B", "edge_confidence": 0.5}],
            ["edge_id", "u", "v", "edge_confidence"],
        )
        empty_cal = base / "empty_cal.csv"
        _write_csv(empty_cal, [], ["edge_id", "observed_present", "weight"])
        val_file = base / "val.csv"
        _write_csv(
            val_file,
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
                "observations_file": str(empty_cal),
                "validation_observations_file": str(val_file),
            },
        )
        with pytest.raises(ValueError, match="calibration observations"):
            run_assumption_calibration_validation(config)


def test_validation_workflow_no_empty_val_observations():
    """When validation observations are empty, a clear error is raised."""
    with _make_temp_dir() as tmp:
        base = Path(tmp)
        edges_file = base / "edges.csv"
        _write_csv(
            edges_file,
            [{"edge_id": "A_B", "u": "A", "v": "B", "edge_confidence": 0.5}],
            ["edge_id", "u", "v", "edge_confidence"],
        )
        cal_file = base / "cal.csv"
        _write_csv(
            cal_file,
            [{"edge_id": "A_B", "observed_present": 1, "weight": 1.0}],
            ["edge_id", "observed_present", "weight"],
        )
        empty_val = base / "empty_val.csv"
        _write_csv(empty_val, [], ["edge_id", "observed_present", "weight"])
        out_dir = str(base / "out")
        config = CalibrationConfig(
            problem_type="topology",
            max_nfev=5,
            output_dir=out_dir,
            engine="internal",
            adapter_settings={
                "candidates_file": str(edges_file),
                "observations_file": str(cal_file),
                "validation_observations_file": str(empty_val),
            },
        )
        with pytest.raises(ValueError, match="validation observations"):
            run_assumption_calibration_validation(config)
