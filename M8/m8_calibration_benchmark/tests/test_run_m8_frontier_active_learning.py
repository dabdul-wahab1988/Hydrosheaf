"""Tests for the prospective frontier active-learning benchmark."""

from __future__ import annotations

import importlib.util
from pathlib import Path
import sys

import numpy as np


SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "scripts"
    / "run_m8_frontier_active_learning.py"
)
SPEC = importlib.util.spec_from_file_location("run_m8_frontier_active_learning", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def _case(seed: int = 1):
    edge_ids = ("A->B", "A->C", "B->D", "C->D")
    truth = np.asarray([1, 0, 1, 0], dtype=int)
    features = np.asarray(
        [
            [1.0, 1.0, 2.0, 0.7],
            [1.0, 1.5, 1.0, 0.5],
            [1.0, 0.8, 2.5, 0.8],
            [1.0, 1.8, 0.7, 0.4],
        ]
    )
    return MODULE.CaseRecord(
        seed=seed,
        edge_ids=edge_ids,
        sources=("A", "A", "B", "C"),
        targets=("B", "C", "D", "D"),
        truth=truth,
        hydraulic_prior=np.asarray([0.7, 0.5, 0.8, 0.4]),
        features=features,
        outcomes={
            "chemistry_panel": np.asarray([0.1, 0.7, 0.2, 0.8]),
            "age_tracer": np.asarray([4.0, -1.0, 5.0, 0.0]),
            "connectivity_tracer": np.asarray([0.9, 0.1, 0.8, 0.2]),
        },
        candidate_recall=1.0,
        provenance={"imports_hydrosheaf": False},
    )


def test_fit_model_and_closed_loop_strategy_are_finite():
    cases = [_case(seed) for seed in range(1, 8)]
    model = MODULE.fit_measurement_model(cases)
    metric, actions = MODULE.run_strategy(
        _case(99), model, "robust_information_decision_per_cost"
    )
    assert metric["actionable"] == 1
    assert 0.0 <= metric["brier"] <= 1.0
    assert np.isfinite(metric["joint_entropy_reduction_per_cost"])
    assert metric["cost_spent"] <= MODULE.BUDGET
    assert actions


def test_topology_particles_respect_one_outgoing_edge_per_source():
    _, particles, weights = MODULE._topology_particles(_case(31))
    assert particles.shape == (MODULE.N_PARTICLES, 4)
    assert np.all(particles[:, :2].sum(axis=1) <= 1)
    assert np.allclose(weights.sum(), 1.0)


def test_development_model_contains_shared_prior_calibration():
    model = MODULE.fit_measurement_model([_case(seed) for seed in range(1, 8)])
    calibration = model["topology_prior_calibration"]
    assert calibration["class_weight"] == "balanced"
    assert len(calibration["coefficients"]) == 3
    probabilities = MODULE._calibrated_edge_prior(_case(22), model)
    assert np.all((probabilities > 0.0) & (probabilities < 1.0))


def test_options_are_explicit_and_scenario_aligned():
    cases = [_case(seed) for seed in range(1, 8)]
    model = MODULE.fit_measurement_model(cases)
    _, particles, _ = MODULE._topology_particles(_case(11))
    options, lookup = MODULE.build_measurement_options(_case(11), particles, model)
    assert len(options) == 4 * len(MODULE.MEASUREMENT_TYPES)
    assert all(option.cost > 0 for option in options)
    assert all(len(option.scenarios) == 3 for option in options)
    assert set(lookup) == {option.option_id for option in options}


def test_gate_requires_every_predeclared_condition():
    import pandas as pd

    metrics = pd.DataFrame(
        [
            {
                "seed": seed,
                "strategy": "robust_information_decision_per_cost",
                "candidate_recall": 1.0,
                "actionable": 1,
            }
            for seed in range(4)
        ]
    )
    rows = []
    for comparator in ("random_feasible", "legacy_uncertainty_chemistry"):
        rows.extend(
            [
                {"comparator": comparator, "metric": "brier", "ci95_high": -0.01, "ci95_low": -0.03},
                {"comparator": comparator, "metric": "joint_entropy_reduction_per_cost", "ci95_high": 0.03, "ci95_low": 0.01},
                {"comparator": comparator, "metric": "pr_auc", "ci95_high": 0.03, "ci95_low": -0.005},
            ]
        )
    decision = MODULE._gate_decision(metrics, pd.DataFrame(rows))
    assert decision["frontier_active_learning_claim_supported"] is True
    assert decision["gates"][
        "information_efficiency_noninferiority_vs_legacy_uncertainty_chemistry"
    ] is True
    decision["gates"]["candidate_recall"] = False
    assert all(decision["gates"].values()) is False


def test_gate_rejects_material_information_harm_against_legacy():
    import pandas as pd

    metrics = pd.DataFrame(
        [
            {
                "seed": seed,
                "strategy": "robust_information_decision_per_cost",
                "candidate_recall": 1.0,
                "actionable": 1,
            }
            for seed in range(4)
        ]
    )
    rows = []
    for comparator in ("random_feasible", "legacy_uncertainty_chemistry"):
        information_low = -0.011 if comparator == "legacy_uncertainty_chemistry" else 0.01
        rows.extend(
            [
                {"comparator": comparator, "metric": "brier", "ci95_high": -0.01, "ci95_low": -0.03},
                {"comparator": comparator, "metric": "joint_entropy_reduction_per_cost", "ci95_high": 0.03, "ci95_low": information_low},
                {"comparator": comparator, "metric": "pr_auc", "ci95_high": 0.03, "ci95_low": -0.005},
            ]
        )

    decision = MODULE._gate_decision(metrics, pd.DataFrame(rows))

    assert decision["frontier_active_learning_claim_supported"] is False
    assert decision["gates"][
        "information_efficiency_noninferiority_vs_legacy_uncertainty_chemistry"
    ] is False
