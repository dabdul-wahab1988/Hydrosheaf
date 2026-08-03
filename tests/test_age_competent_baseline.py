from __future__ import annotations

import ast
from pathlib import Path
import math

import pytest

from hydrosheaf.validation.age_competent_baseline import (
    ABSTAIN,
    SELECT,
    AgeBaselineConfig,
    AgeCompetentBaseline,
    enumerate_age_candidates,
)


def _age_channel(age_years: float = 20.0) -> dict[str, object]:
    half_lives = {
        "tritium_TU": (8.0, 12.32),
        "argon39_pmc": (100.0, 269.0),
        "c14_pmc": (100.0, 5730.0),
    }
    features = {
        field: reference * math.exp(-math.log(2.0) * age_years / half_life)
        for field, (reference, half_life) in half_lives.items()
    }
    features.update(
        {
            "tritium_sigma_TU": 0.15,
            "argon39_sigma_pmc": 2.0,
            "c14_sigma_pmc": 0.8,
            "sample_date": 2025.5,
        }
    )
    return {"tracer_age_history": {"nodes": {"site-a": features}}}


def test_candidate_enumeration_is_deterministic_and_bounded() -> None:
    config = AgeBaselineConfig(
        max_age_years=20.0,
        candidate_step_years=5.0,
        quadrature_step_years=1.0,
        gamma_shapes=(2.0, 4.0),
    )
    first = enumerate_age_candidates(config)
    second = enumerate_age_candidates(config)

    assert first.candidate_hash == second.candidate_hash
    assert [candidate.candidate_id for candidate in first.candidates] == [
        candidate.candidate_id for candidate in second.candidates
    ]
    assert first.age_grid_years == (0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0)
    assert {candidate.family for candidate in first.candidates} == {
        "exponential",
        "gamma",
    }
    assert min(candidate.mean_age_years for candidate in first.candidates) == 0.0
    assert max(candidate.mean_age_years for candidate in first.candidates) == 20.0
    assert first.to_audit_record()["truth_blind"] is True


def test_source_has_no_local_inference_or_existing_specialist_import() -> None:
    source = Path(__file__).resolve().parents[1] / "hydrosheaf" / "validation" / "age_competent_baseline.py"
    tree = ast.parse(source.read_text(encoding="utf-8"))
    imported_modules: list[str] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported_modules.extend(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            imported_modules.append(str(node.module or ""))
    assert not any(module.startswith("hydrosheaf.inference") for module in imported_modules)
    assert "hydrosheaf.validation.specialist_baselines" not in imported_modules
    assert "hydrosheaf.validation.specialist_candidate_generation" not in imported_modules


def test_prediction_is_truth_blind_and_ignores_non_tracer_observations() -> None:
    baseline = AgeCompetentBaseline()
    ordinary = _age_channel()
    ordinary["tracer_age_history"]["nodes"]["site-a"].update(
        {"head_meas": 100.0, "x_m": 12.0, "y_m": 48.0, "aquifer_layer": 2}
    )
    clean = baseline.predict(_age_channel())
    with_non_tracer = baseline.predict(ordinary)

    for field in ("mean_age_years", "age_95_low", "age_95_high", "decision"):
        assert with_non_tracer["site-a"][field] == clean["site-a"][field]

    leaked = _age_channel()
    leaked["tracer_age_history"]["nodes"]["site-a"]["true_age_years"] = 20.0
    with pytest.raises(ValueError, match="Truth/reference field"):
        baseline.predict(leaked)

    posterior_leak = _age_channel()
    posterior_leak["tracer_age_history"]["nodes"]["site-a"]["hydrosheaf_posterior"] = {
        "mean_age_years": 20.0
    }
    with pytest.raises(ValueError, match="Truth/reference field"):
        baseline.predict(posterior_leak)


def test_supported_age_output_has_valid_point_and_interval_uncertainty() -> None:
    prediction = AgeCompetentBaseline().predict(_age_channel())["site-a"]

    assert prediction["decision"] == SELECT
    assert prediction["age_95_low"] <= prediction["mean_age_years"]
    assert prediction["mean_age_years"] <= prediction["age_95_high"]
    assert 0.0 <= prediction["age_95_low"]
    assert prediction["age_95_high"] <= 120.0
    assert prediction["uncertainty_years"] >= 0.0
    assert prediction["tracers_used"] == ["argon39_pmc", "c14_pmc", "tritium_TU"]
    assert prediction["interval_kind"] == "model_based_95_percentile_not_calibrated"


def test_unsupported_and_ambiguous_inputs_abstain() -> None:
    baseline = AgeCompetentBaseline()

    missing = baseline.predict(
        {"tracer_age_history": {"nodes": {"missing": {"sample_date": 2025.5}}}}
    )["missing"]
    assert missing["decision"] == ABSTAIN
    assert missing["reason"] == "insufficient_supported_age_tracers"
    assert missing["age_95_low"] == 0.0
    assert missing["age_95_high"] == 120.0

    ambiguous = _age_channel(age_years=20.0)
    ambiguous_features = ambiguous["tracer_age_history"]["nodes"]["site-a"]
    ambiguous_features["argon39_pmc"] = 20.0
    result = baseline.predict(ambiguous)["site-a"]
    assert result["decision"] == ABSTAIN
    assert "tracer_age_disagreement" in str(result["reason"])

