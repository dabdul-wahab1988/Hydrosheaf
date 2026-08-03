from __future__ import annotations

import pytest

from hydrosheaf.validation import (
    AGE_OUTPUT,
    REACTION_OUTPUT,
    default_specialist_baseline_registry,
    score_age_baseline_outputs,
    score_reaction_baseline_outputs,
)


def _observations() -> dict[str, object]:
    return {
        "tracer_age": {
            "nodes": {
                "young": {"tritium_TU": 8.0, "argon39_pmc": 100.0},
                "old": {"tritium_TU": 0.4, "argon39_pmc": 50.0},
                "missing": {},
            }
        },
        "reaction_chemistry": {
            "edges": {
                "young->old": {
                    "upstream": {"NO3": 1.0, "SO4": 1.0, "Fe": 0.01},
                    "downstream": {"NO3": 0.7, "SO4": 1.0, "Fe": 0.01},
                },
                "same->same": {
                    "upstream": {"Ca": 1.0, "HCO3": 2.0, "SO4": 1.0},
                    "downstream": {"Ca": 1.0, "HCO3": 2.0, "SO4": 1.0},
                },
                "old->missing": {"upstream": {}, "downstream": {}},
            }
        },
    }


def test_specialist_registry_is_explicit_and_fingerprinted() -> None:
    registry = default_specialist_baseline_registry()

    assert registry.names() == (
        "independent_competence_matched_age_specialist",
        "local_thermodynamic_reaction_rules",
        "local_tracer_decay_age",
        "multitracer_atmospheric_history_age",
        "reaction_competent_baseline",
        "stoichiometric_reaction_inverse",
    )
    records = registry.audit_table()
    assert {record["output_kind"] for record in records} == {
        AGE_OUTPUT,
        REACTION_OUTPUT,
    }
    for record in records:
        assert record["uncertainty"]["calibrated"] is False
        assert len(record["fingerprint"]) == 64


def test_local_age_baseline_is_truth_blind_and_abstains_on_missing_tracer() -> None:
    spec = default_specialist_baseline_registry().get("local_tracer_decay_age")
    predictions = spec.predict(_observations())

    assert predictions["young"]["mean_age_years"] == pytest.approx(0.0)
    assert predictions["old"]["mean_age_years"] > predictions["young"]["mean_age_years"]
    assert predictions["missing"]["decision"] == "abstain"
    assert predictions["missing"]["age_95_low"] == 0.0
    assert predictions["missing"]["age_95_high"] == 200.0

    with pytest.raises(ValueError, match="Truth/reference field"):
        spec.predict({"tracer_age": {"nodes": {"x": {"true_age": 4.0}}}})


def test_local_reaction_rules_emit_probabilities_and_explicit_abstention() -> None:
    spec = default_specialist_baseline_registry().get(
        "local_thermodynamic_reaction_rules"
    )
    predictions = spec.predict(_observations())

    denitrification = predictions["young->old"]
    assert denitrification["family"] == "denitrification"
    assert denitrification["decision"] == "select"
    assert sum(denitrification["probabilities"].values()) == pytest.approx(1.0)

    missing = predictions["old->missing"]
    assert missing["family"] == "none"
    assert missing["decision"] == "abstain"


def test_multitracer_age_and_stoichiometric_reaction_are_explicit_upgrades() -> None:
    registry = default_specialist_baseline_registry()
    observations = _observations()
    observations["tracer_age_history"] = {
        "nodes": {
            "young": {
                **observations["tracer_age"]["nodes"]["young"],
                "cfc12_pptv": 520.0,
                "sf6_pptv": 10.0,
                "sample_year": 2025.0,
            },
            "old": {
                **observations["tracer_age"]["nodes"]["old"],
                "cfc12_pptv": 140.0,
                "sf6_pptv": 3.0,
                "sample_year": 2025.0,
            },
        }
    }

    age = registry.get("multitracer_atmospheric_history_age").predict(observations)
    reaction = registry.get("stoichiometric_reaction_inverse").predict(observations)

    assert age["young"]["tracers_used"]
    assert age["young"]["age_95_low"] <= age["young"]["mean_age_years"]
    assert age["young"]["mean_age_years"] <= age["young"]["age_95_high"]
    assert reaction["young->old"]["logits"]
    assert sum(reaction["young->old"]["probabilities"].values()) == pytest.approx(1.0)
    assert reaction["same->same"]["family"] == "none"


def test_specialist_scores_are_truth_release_operations() -> None:
    registry = default_specialist_baseline_registry()
    observations = _observations()
    age_predictions = registry.get("local_tracer_decay_age").predict(observations)
    reaction_predictions = registry.get(
        "local_thermodynamic_reaction_rules"
    ).predict(observations)

    age_score = score_age_baseline_outputs(
        {"young": 0.0, "old": 100.0, "missing": 100.0},
        age_predictions,
    )
    reaction_score = score_reaction_baseline_outputs(
        {"young->old": "denitrification", "old->missing": "carbonate"},
        reaction_predictions,
    )

    assert age_score["status"] == "scored"
    assert age_score["n"] == 3
    assert reaction_score["status"] == "scored"
    assert reaction_score["n"] == 2
    assert reaction_score["metrics"]["accuracy"] >= 0.5


def test_specialist_scores_do_not_accept_missing_truth_as_success() -> None:
    assert score_age_baseline_outputs({}, {})["status"] == "not_available"
    reaction = score_reaction_baseline_outputs(
        {"edge": "carbonate"},
        {},
    )
    assert reaction["status"] == "not_available"
    assert reaction["n_missing_outputs"] == 1
