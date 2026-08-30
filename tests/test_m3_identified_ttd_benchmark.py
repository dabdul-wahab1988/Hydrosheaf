from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pandas as pd
import pytest

from hydrosheaf.nuclear.joint_lpm import predict_lpm_tracers

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = (
    REPO_ROOT
    / "M3"
    / "m3_age_benchmark"
    / "scripts"
    / "run_m3_identified_ttd_benchmark.py"
)
PROTOCOL = (
    REPO_ROOT / "M3" / "m3_age_benchmark" / "configs" / "identified_ttd_protocol.yaml"
)


def _load_module():
    spec = importlib.util.spec_from_file_location("m3_identified_ttd", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    # Register before exec_module: dataclasses resolves string annotations
    # (the script uses `from __future__ import annotations`) via
    # sys.modules[cls.__module__].__dict__, which is None for an unregistered
    # module. Without this, any dataclass in the script fails to construct.
    sys.modules[spec.name] = module
    try:
        spec.loader.exec_module(module)
    except Exception:
        sys.modules.pop(spec.name, None)
        raise
    return module


def _synthetic_row():
    sample_year = 2020.0
    predictions = predict_lpm_tracers(
        "PFM",
        {"mean_age_years": 20.0},
        sample_year,
        ["3H", "3H/3He", "SF6", "14C"],
        max_age_years=100.0,
    )
    return {
        "site_id": "synthetic-1",
        "StudyUnit": "test",
        "AqGroup": "test-aquifer",
        "sample_year": sample_year,
        "lat": 40.0,
        "lon": -100.0,
        "tritium_TU": predictions["3H"],
        "tritium_sigma_TU": 0.2,
        "he3_trit_TU": predictions["3H/3He"],
        "he3_trit_sigma_TU": 0.2,
        "sf6_pptv": predictions["SF6"],
        "sf6_sigma_pptv": 0.2,
        "c14_pmc": predictions["14C"],
        "c14_sigma_pmc": 1.0,
        "reference_age_years": 999999.0,
        "FracAnthropocene": 0.0,
    }


def test_protocol_age_grid_is_locked_increasing_and_contains_cutoffs():
    module = _load_module()
    protocol = module.load_protocol(PROTOCOL)
    ages = module.build_protocol_age_grid(protocol)

    assert ages[0] == 0.0
    assert ages[-1] == 50000.0
    assert all(right > left for left, right in zip(ages, ages[1:]))
    assert 70.0 in ages


def test_split_removes_tritium_and_its_daughter_equivalent():
    module = _load_module()
    protocol = module.load_protocol(PROTOCOL)
    observations, _ = module.prepare_row_observations(_synthetic_row(), protocol)
    conditioning, held_out = module.split_held_out_observation(observations, "3H")

    assert held_out is not None and held_out.tracer == "3H"
    assert "3H" not in {item.tracer for item in conditioning}
    assert "3H/3He" not in {item.tracer for item in conditioning}


def test_held_out_row_uses_concentration_not_reference_age():
    module = _load_module()
    protocol = module.load_protocol(PROTOCOL)
    first = module.evaluate_held_out_row(_synthetic_row(), "SF6", protocol)
    changed = _synthetic_row()
    changed["reference_age_years"] = 1.0
    changed["FracAnthropocene"] = 1.0
    second = module.evaluate_held_out_row(changed, "SF6", protocol)

    assert first["reference_fields_used"] is False
    assert first["conditioning_tracers"] == second["conditioning_tracers"]
    assert first["prediction_lower"] == pytest.approx(second["prediction_lower"])
    assert first["prediction_upper"] == pytest.approx(second["prediction_upper"])
    assert first["conditioned_on_held_out_observation"] is False


def test_missing_held_out_tracer_is_not_eligible():
    module = _load_module()
    protocol = module.load_protocol(PROTOCOL)
    row = _synthetic_row()
    row["sf6_pptv"] = None
    result = module.evaluate_held_out_row(row, "SF6", protocol)

    assert result["status"] == "NOT_ELIGIBLE"
    assert result["reason"] == "held_out_observation_unavailable"


def test_empty_development_subset_has_a_machine_readable_summary():
    module = _load_module()
    summary = module.summarize_results(pd.DataFrame())

    assert summary["n_rows"] == 0
    assert summary["held_out_coverage"] is None


def test_unsupported_held_out_likelihood_is_preserved_as_abstention(monkeypatch):
    module = _load_module()
    protocol = module.load_protocol(PROTOCOL)
    monkeypatch.setattr(
        module,
        "assess_held_out_tracer",
        lambda *args, **kwargs: {
            "status": "ABSTAIN",
            "reason": "no_declared_linear_interval_semantics",
        },
    )

    result = module.evaluate_held_out_row(_synthetic_row(), "SF6", protocol)

    assert result["held_out_status"] == "ABSTAIN"
    assert result["held_out_compatible"] != result["held_out_compatible"]
    assert result["held_out_reason"] == "no_declared_linear_interval_semantics"
    assert result["conditioned_on_held_out_observation"] is False


def test_unexpected_row_failure_is_audited_instead_of_stopping_the_run(monkeypatch):
    module = _load_module()
    protocol = module.load_protocol(PROTOCOL)
    monkeypatch.setattr(
        module,
        "evaluate_held_out_row",
        lambda *args, **kwargs: (_ for _ in ()).throw(RuntimeError("synthetic")),
    )

    result = module.evaluate_held_out_row_safely(_synthetic_row(), "SF6", protocol)

    assert result["status"] == "ABSTAIN"
    assert result["reason"] == "evaluation_error:RuntimeError:synthetic"
    assert result["reference_fields_used"] is False
