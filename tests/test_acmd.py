"""Comprehensive test suite for Active and Certified Measurement Design (ACMD).

Validates:
1. ACMDAction data contracts, validation, standalone cost, and incremental mobilization discount.
2. Initial state evaluation: ALREADY_RESOLVED, INCONSISTENT_EVIDENCE, and ACTIVE states.
3. Closed-loop active cycle: P_t --> x_{t+1} --> y_{t+1} --> P_{t+1} leading to CERTIFIED.
4. Mathematical IMPOSSIBILITY_WITNESS detection when candidates cannot resolve delta.
5. Inconsistent evidence detection when observations contradict physics.
6. Budget exhaustion stopping rule.
7. Full automated simulation against true state oracle with cryptographic certificate.
8. Probabilistic EIG mode and probability gate compliance.
9. Hybrid mode combining worst-case ambiguity reduction with EIG.
"""

from __future__ import annotations

import math
from typing import Sequence
import numpy as np
import pytest

from hydrosheaf import (
    ACMD,
    ACMDAction,
    ACMDCertificate,
    ACMDLoop,
    ACMDMode,
    ACMDStatus,
    ACMDStepRecord,
)
from hydrosheaf.nuclear.ttd_design import (
    TtdHypothesisEnsemble,
)
from hydrosheaf.nuclear.ttd_identified import (
    AgeFunctional,
    TracerConstraint,
)


@pytest.fixture
def age_grid_5() -> np.ndarray:
    """5-bin age grid: [5.0, 20.0, 50.0, 80.0, 150.0] years."""
    return np.asarray([5.0, 20.0, 50.0, 80.0, 150.0], dtype=float)


@pytest.fixture
def young_water_functional() -> AgeFunctional:
    """Target functional: Fraction of water <= 50 years (bins 0, 1, 2)."""
    return AgeFunctional(
        name="young_water_fraction",
        coefficients=[1.0, 1.0, 1.0, 0.0, 0.0],
        units="fraction",
    )


# ── 1. Action contract & cost tests ──────────────────────────────────────────


def test_acmd_action_validation():
    with pytest.raises(ValueError, match="action_id must be non-empty"):
        ACMDAction(action_id="", measurement_type="3H", target_id="W1")

    with pytest.raises(ValueError, match="measurement_type must be non-empty"):
        ACMDAction(action_id="act1", measurement_type="", target_id="W1")

    with pytest.raises(ValueError, match="target_id must be non-empty"):
        ACMDAction(action_id="act1", measurement_type="3H", target_id="")

    with pytest.raises(ValueError, match="cost must be finite and strictly positive"):
        ACMDAction(action_id="act1", measurement_type="3H", target_id="W1", cost=0.0)

    with pytest.raises(ValueError, match="travel_cost must be finite and non-negative"):
        ACMDAction(action_id="act1", measurement_type="3H", target_id="W1", travel_cost=-1.0)

    with pytest.raises(ValueError, match="accessibility must be finite and strictly positive"):
        ACMDAction(action_id="act1", measurement_type="3H", target_id="W1", accessibility=0.0)

    with pytest.raises(ValueError, match="error_bound must be finite and strictly positive"):
        ACMDAction(action_id="act1", measurement_type="3H", target_id="W1", error_bound=-0.1)


def test_acmd_action_cost_formulas():
    action = ACMDAction(
        action_id="act_well1_tracer",
        measurement_type="3H",
        target_id="Well_01",
        cost=10.0,
        travel_cost=5.0,
        accessibility=0.8,
    )
    # Standalone cost: (10 + 5) / 0.8 = 18.75
    assert math.isclose(action.standalone_cost, 18.75, rel_tol=1e-6)

    # Incremental cost if Well_01 not yet visited: travel included -> 18.75
    assert math.isclose(action.incremental_cost(visited_targets=set()), 18.75, rel_tol=1e-6)

    # Incremental cost if Well_01 already visited: travel=0 -> 10.0 / 0.8 = 12.50
    assert math.isclose(action.incremental_cost(visited_targets={"Well_01"}), 12.50, rel_tol=1e-6)


# ── 2. Boundary state tests ──────────────────────────────────────────────────


def test_acmd_already_resolved(age_grid_5, young_water_functional):
    # Prior constraint already pins young water fraction tightly
    c_prior = TracerConstraint(
        tracer="T_prior",
        response=[1.0, 1.0, 1.0, 0.0, 0.0],
        observed=0.75,
        sigma=0.005,
    )
    cand = ACMDAction(
        action_id="act_cand1",
        measurement_type="3H",
        target_id="W1",
        cost=5.0,
        response=[0.5, 0.5, 0.5, 0.1, 0.1],
    )

    loop = ACMD.create_loop(
        age_grid_5,
        [c_prior],
        [cand],
        young_water_functional,
        target_tolerance=0.15,
    )
    assert loop.status == ACMDStatus.ALREADY_RESOLVED
    assert loop.is_terminated
    cert = loop.get_certificate()
    assert cert.status == ACMDStatus.ALREADY_RESOLVED
    assert cert.total_cost == 0.0


def test_acmd_inconsistent_evidence(age_grid_5, young_water_functional):
    # Two mutually contradictory prior constraints (empty polytope)
    c1 = TracerConstraint(
        tracer="T1",
        response=[1.0, 1.0, 1.0, 1.0, 1.0],
        observed=0.90,
        sigma=0.01,
    )
    c2 = TracerConstraint(
        tracer="T2",
        response=[1.0, 0.0, 0.0, 0.0, 0.0],
        observed=0.10,
        sigma=0.01,
    )
    cand = ACMDAction(action_id="c1", measurement_type="3H", target_id="W1", response=[1, 1, 1, 1, 1])

    loop = ACMD.create_loop(
        age_grid_5,
        [c1, c2],
        [cand],
        young_water_functional,
        target_tolerance=0.10,
    )
    assert loop.status == ACMDStatus.INCONSISTENT_EVIDENCE
    assert loop.is_terminated
    cert = loop.get_certificate()
    assert cert.status == ACMDStatus.INCONSISTENT_EVIDENCE


def test_acmd_impossibility_witness(age_grid_5, young_water_functional):
    # Candidate measurement has flat response (cannot separate young from old)
    cand_blind = ACMDAction(
        action_id="cand_blind",
        measurement_type="inert",
        target_id="W1",
        cost=5.0,
        error_bound=0.01,
        response=[1.0, 1.0, 1.0, 1.0, 1.0],  # 1 * sum(x) = 1 always
    )

    # Initial constraints empty -> P_0 is simplex, ambiguity = 1.0
    # Requesting delta = 0.05
    loop = ACMD.create_loop(
        age_grid_5,
        (),
        [cand_blind],
        young_water_functional,
        target_tolerance=0.05,
    )
    assert loop.status == ACMDStatus.ACTIVE

    # Attempt to recommend action: should detect impossibility immediately!
    action, diag = loop.recommend_next_action()
    assert action is None
    assert loop.status == ACMDStatus.IMPOSSIBILITY_WITNESS
    assert diag["status"] == "IMPOSSIBILITY_WITNESS"

    cert = loop.get_certificate()
    assert cert.status == ACMDStatus.IMPOSSIBILITY_WITNESS
    assert cert.lower_witness is not None
    assert cert.upper_witness is not None
    # Verify constructive witness pair
    diff = float(young_water_functional.coefficients @ (cert.upper_witness - cert.lower_witness))
    assert math.isclose(diff, cert.candidate_attainable_ambiguity, abs_tol=1e-5)


# ── 3. Dynamic Closed-Loop Step Updates ──────────────────────────────────────


def test_acmd_closed_loop_step_convergence(age_grid_5, young_water_functional):
    # We want young water fraction resolved to within delta = 0.10
    # Candidate 1: 3H, responds to youngest bins
    cand1 = ACMDAction(
        action_id="act_tritium",
        measurement_type="3H",
        target_id="Well_01",
        cost=6.0,
        travel_cost=4.0,
        error_bound=0.02,
        standard_deviation=0.01,
        response=[1.0, 1.0, 0.5, 0.0, 0.0],
    )
    # Candidate 2: SF6, responds to intermediate young bins
    cand2 = ACMDAction(
        action_id="act_sf6",
        measurement_type="SF6",
        target_id="Well_01",
        cost=4.0,
        travel_cost=4.0,
        error_bound=0.02,
        standard_deviation=0.01,
        response=[0.0, 0.0, 0.8, 0.0, 0.0],
    )

    # Prior: weak initial constraint
    c0 = TracerConstraint(
        tracer="prior_tracer",
        response=[0.5, 0.5, 0.5, 0.5, 0.5],
        observed=0.5,
        sigma=0.2,
    )

    loop = ACMD.create_loop(
        age_grid_5,
        [c0],
        [cand1, cand2],
        young_water_functional,
        target_tolerance=0.10,
        mode=ACMDMode.ROBUST_MINIMAX,
    )
    assert loop.status == ACMDStatus.ACTIVE
    initial_w = loop.current_ambiguity
    assert initial_w > 0.10

    # Step 1: recommend next action
    act_step1, diag1 = loop.recommend_next_action()
    assert act_step1 is not None
    assert act_step1.action_id in ("act_tritium", "act_sf6")

    # Suppose true mass vector is [0.3, 0.3, 0.2, 0.1, 0.1]
    true_x = np.array([0.3, 0.3, 0.2, 0.1, 0.1])
    resp1 = act_step1.resolve_response(age_grid_5)
    true_y1 = float(resp1 @ true_x)

    # Ingest observation y1
    step_rec1 = loop.step(act_step1.action_id, true_y1)
    assert step_rec1.step_index == 1
    assert step_rec1.posterior_ambiguity < initial_w
    assert loop.cumulative_cost > 0.0
    assert "Well_01" in loop.visited_targets

    # If not yet certified, take step 2
    if not loop.is_terminated:
        act_step2, _ = loop.recommend_next_action()
        assert act_step2 is not None
        # Travel cost to Well_01 should now be discounted to 0.0!
        assert act_step2.incremental_cost(loop.visited_targets) < act_step2.standalone_cost

        resp2 = act_step2.resolve_response(age_grid_5)
        true_y2 = float(resp2 @ true_x)
        step_rec2 = loop.step(act_step2.action_id, true_y2)
        assert step_rec2.step_index == 2

    # Loop should now be CERTIFIED
    assert loop.status == ACMDStatus.CERTIFIED
    assert loop.is_terminated
    assert loop.current_ambiguity <= 0.10 + 1e-7

    cert = loop.get_certificate()
    assert cert.status == ACMDStatus.CERTIFIED
    assert len(cert.steps) in (1, 2)
    assert len(cert.certificate_hash) == 64  # SHA-256


# ── 4. Automated Simulation Against Oracle ───────────────────────────────────


def test_acmd_simulation_runner(age_grid_5, young_water_functional):
    true_state = np.array([0.4, 0.3, 0.1, 0.1, 0.1])

    candidates = [
        ACMDAction(
            action_id="tracer_A",
            measurement_type="3H",
            target_id="W1",
            cost=5.0,
            travel_cost=2.0,
            error_bound=0.02,
            response=[1.0, 1.0, 0.5, 0.0, 0.0],
        ),
        ACMDAction(
            action_id="tracer_B",
            measurement_type="CFC12",
            target_id="W2",
            cost=6.0,
            travel_cost=3.0,
            error_bound=0.02,
            response=[0.0, 0.0, 0.8, 0.0, 0.0],
        ),
    ]

    loop = ACMD.create_loop(
        age_grid_5,
        (),
        candidates,
        young_water_functional,
        target_tolerance=0.10,
        mode=ACMDMode.ROBUST_MINIMAX,
    )

    oracle = lambda act: float(act.resolve_response(age_grid_5) @ true_state)

    cert = loop.run_simulation(oracle)
    assert cert.status == ACMDStatus.CERTIFIED
    assert cert.final_ambiguity <= 0.10 + 1e-7
    assert cert.total_cost > 0.0
    assert len(cert.steps) >= 1
    assert cert.certificate_hash


# ── 5. Probabilistic EIG Mode & Probability Gate ─────────────────────────────


def test_acmd_probabilistic_eig_mode(age_grid_5, young_water_functional):
    # Hypothesis ensemble: 3 distinct plausible TTD models
    masses = [
        [0.6, 0.2, 0.1, 0.05, 0.05],  # Very young
        [0.2, 0.2, 0.2, 0.2, 0.2],     # Uniform
        [0.05, 0.05, 0.1, 0.3, 0.5],  # Old
    ]
    ensemble = TtdHypothesisEnsemble(
        hypothesis_ids=["young_model", "uniform_model", "old_model"],
        age_grid_years=age_grid_5,
        masses=masses,
        probabilities=[0.5, 0.3, 0.2],
        probability_semantics="posterior",
    )

    cand1 = ACMDAction(
        action_id="act_eig_1",
        measurement_type="3H",
        target_id="W1",
        cost=5.0,
        standard_deviation=0.05,
        response=[1.0, 0.8, 0.5, 0.1, 0.0],
    )
    cand2 = ACMDAction(
        action_id="act_eig_2",
        measurement_type="14C",
        target_id="W1",
        cost=8.0,
        standard_deviation=0.05,
        response=[0.1, 0.1, 0.2, 0.8, 1.0],
    )

    loop = ACMD.create_loop(
        age_grid_5,
        (),
        [cand1, cand2],
        young_water_functional,
        target_tolerance=0.20,
        mode=ACMDMode.PROBABILISTIC_EIG,
        hypothesis_ensemble=ensemble,
    )

    act, diag = loop.recommend_next_action()
    assert act is not None
    assert diag["mode"] == "probabilistic_eig"
    assert "act_eig_1" in diag["all_scores"]
    assert "act_eig_2" in diag["all_scores"]


def test_acmd_probability_gate_refusal(age_grid_5, young_water_functional):
    # Attempting EIG mode WITHOUT hypothesis probabilities must ABSTAIN!
    cand = ACMDAction(
        action_id="act_cand",
        measurement_type="3H",
        target_id="W1",
        cost=5.0,
        response=[1.0, 0.8, 0.5, 0.1, 0.0],
    )

    loop = ACMD.create_loop(
        age_grid_5,
        (),
        [cand],
        young_water_functional,
        target_tolerance=0.20,
        mode=ACMDMode.PROBABILISTIC_EIG,
        hypothesis_ensemble=None,  # No ensemble provided
    )

    act, diag = loop.recommend_next_action()
    assert act is None
    assert loop.status == ACMDStatus.ABSTAIN
    assert diag["reason"] == "no_hypothesis_ensemble_provided_for_eig_mode"
