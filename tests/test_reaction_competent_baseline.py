"""Deterministic tests for the independent reaction competence baseline."""

from __future__ import annotations

import ast
import importlib.util
import math
import sys
from pathlib import Path

import pytest

_MODULE_PATH = (
    Path(__file__).parents[1]
    / "hydrosheaf"
    / "validation"
    / "reaction_competent_baseline.py"
)
_SPEC = importlib.util.spec_from_file_location(
    "reaction_competent_baseline_under_test", _MODULE_PATH
)
assert _SPEC is not None and _SPEC.loader is not None
_MODULE = importlib.util.module_from_spec(_SPEC)
sys.modules[_SPEC.name] = _MODULE
_SPEC.loader.exec_module(_MODULE)

ABSTAIN = _MODULE.ABSTAIN
NULL_FAMILY = _MODULE.NULL_FAMILY
SELECT = _MODULE.SELECT
ReactionCompetentBaseline = _MODULE.ReactionCompetentBaseline


def _observations() -> dict[str, object]:
    return {
        "reaction_chemistry": {
            "edges": {
                "carbonate_edge": {
                    "upstream": {"Ca": 1.00, "HCO3": 2.00},
                    "downstream": {"Ca": 1.20, "HCO3": 2.20},
                },
                "unchanged_edge": {
                    "upstream": {"Ca": 1.00, "HCO3": 2.00},
                    "downstream": {"Ca": 1.00, "HCO3": 2.00},
                },
                "noisy_edge": {
                    "upstream": {"Ca": 1.00, "HCO3": 2.00},
                    "downstream": {"Ca": 1.20, "HCO3": 2.10},
                },
                "unsupported_edge": {
                    "upstream": {"DOC": 1.00},
                    "downstream": {"DOC": 1.10},
                },
            }
        }
    }


def test_candidate_generation_is_deterministic_and_includes_null() -> None:
    baseline = ReactionCompetentBaseline()

    first = baseline.generate_candidates(_observations())
    second = baseline.generate_candidates(_observations())

    assert {
        edge_id: [candidate.to_dict() for candidate in candidates]
        for edge_id, candidates in first.items()
    } == {
        edge_id: [candidate.to_dict() for candidate in candidates]
        for edge_id, candidates in second.items()
    }
    families = {candidate.family for candidate in first["carbonate_edge"]}
    assert NULL_FAMILY in families
    assert "carbonate" in families
    carbonate = next(
        candidate
        for candidate in first["carbonate_edge"]
        if candidate.family == "carbonate"
    )
    assert carbonate.extent == pytest.approx(4.0)
    assert carbonate.residual_norm == pytest.approx(0.0)
    assert all(candidate.extent >= 0.0 for candidate in first["carbonate_edge"])
    assert "unsupported_edge" not in first


def test_null_candidate_wins_on_unchanged_chemistry_with_explicit_penalty() -> None:
    baseline = ReactionCompetentBaseline()

    candidates = baseline.generate_candidates(_observations())["unchanged_edge"]
    null = next(
        candidate for candidate in candidates if candidate.family == NULL_FAMILY
    )
    reaction = next(
        candidate for candidate in candidates if candidate.family != NULL_FAMILY
    )
    prediction = baseline.predict(_observations())["unchanged_edge"]

    assert null.complexity_penalty == pytest.approx(0.0)
    assert reaction.complexity_penalty > 0.0
    assert prediction["family"] == NULL_FAMILY
    assert prediction["decision"] == SELECT
    assert prediction["reason"] == "no_reaction_supported_by_delta_tolerance"
    assert prediction["probabilities"][NULL_FAMILY] == pytest.approx(1.0)
    assert prediction["diagnostics"]["null_complexity_penalty"] == pytest.approx(0.0)


def test_residual_and_extent_uncertainty_diagnostics_are_finite_and_explicit() -> None:
    baseline = ReactionCompetentBaseline()
    prediction = baseline.predict(_observations())["noisy_edge"]
    carbonate = next(
        candidate
        for candidate in prediction["candidate_explanations"]
        if candidate["family"] == "carbonate"
    )

    assert math.isfinite(carbonate["extent"])
    assert math.isfinite(carbonate["residual_norm"])
    assert set(carbonate["residuals"]) == {"Ca", "HCO3"}
    assert carbonate["residual_norm"] > 0.0
    assert carbonate["uncertainty"]["extent_std"] is not None
    lower, upper = carbonate["uncertainty"]["extent_ci95"]
    assert 0.0 <= lower <= carbonate["extent"] <= upper
    assert carbonate["uncertainty"]["calibrated"] is False
    assert prediction["uncertainty"]["selection_probability_calibrated"] is False


def test_conflicting_signature_abstains_with_residual_diagnostic() -> None:
    observations = {
        "reaction_chemistry": {
            "edges": {
                "conflict": {
                    "upstream": {"Ca": 1.00, "HCO3": 2.00},
                    "downstream": {"Ca": 1.20, "HCO3": 1.80},
                }
            }
        }
    }

    prediction = ReactionCompetentBaseline().predict(observations)["conflict"]

    assert prediction["decision"] == ABSTAIN
    assert prediction["reason"] in {
        "ambiguous_reaction_explanation",
        "high_residual",
    }
    assert prediction["family"] == NULL_FAMILY
    assert prediction["residual_norm"] is not None
    assert prediction["diagnostics"]["selection_probability_is_uncalibrated"] is True


def test_unsupported_input_abstains_without_inventing_a_candidate() -> None:
    observations = {
        "reaction_chemistry": {
            "edges": {
                "unsupported": {
                    "upstream": {"DOC": 1.0},
                    "downstream": {"DOC": 1.2},
                },
                "negative": {
                    "upstream": {"Ca": 1.0, "HCO3": 2.0},
                    "downstream": {"Ca": -1.0, "HCO3": 2.0},
                },
            }
        }
    }

    predictions = ReactionCompetentBaseline().predict(observations)

    unsupported = predictions["unsupported"]
    assert unsupported["decision"] == ABSTAIN
    assert unsupported["family"] == NULL_FAMILY
    assert unsupported["reason"] == "unsupported_chemistry_fields"
    assert "DOC" in unsupported["diagnostics"]["unsupported_fields"]
    assert unsupported["candidate_explanations"] == []

    negative = predictions["negative"]
    assert negative["decision"] == ABSTAIN
    assert negative["reason"] == "invalid_chemistry_values"
    assert "Ca" in negative["diagnostics"]["invalid_fields"]


def test_truth_reference_fields_are_rejected() -> None:
    observations = _observations()
    observations["reaction_chemistry"]["reference_edges"] = {}  # type: ignore[index]

    with pytest.raises(ValueError, match="Truth/reference field"):
        ReactionCompetentBaseline().predict(observations)


def test_audit_record_declares_independence_and_uncalibrated_uncertainty() -> None:
    record = ReactionCompetentBaseline().to_audit_record()

    assert record["control"]["truth_blind"] is True
    assert record["control"]["uses_hydrosheaf_inference"] is False
    assert record["control"]["uses_hydrosheaf_posterior"] is False
    assert record["control"]["uses_phreeqc_outputs_as_truth"] is False
    assert record["uncertainty"]["calibrated"] is False


def test_module_has_no_hydrosheaf_or_phreeqc_import_dependency() -> None:
    module_path = (
        Path(__file__).parents[1]
        / "hydrosheaf"
        / "validation"
        / "reaction_competent_baseline.py"
    )
    tree = ast.parse(module_path.read_text(encoding="utf-8"))
    imports = [
        node.module or "" for node in ast.walk(tree) if isinstance(node, ast.ImportFrom)
    ] + [
        alias.name
        for node in ast.walk(tree)
        if isinstance(node, ast.Import)
        for alias in node.names
    ]

    assert all(not name.startswith("hydrosheaf") for name in imports)
    assert all("phreeqc" not in name.lower() for name in imports)


