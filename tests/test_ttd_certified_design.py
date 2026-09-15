"""Tests for certified measurement design for partially identified groundwater systems."""

from __future__ import annotations

import math
import numpy as np
import pytest

from hydrosheaf.nuclear.joint_lpm import tracer_response_kernel
from hydrosheaf.nuclear.ttd_certified_design import (
    CertifiedCandidateTracer,
    evaluate_worst_case_ambiguity,
    solve_budgeted_minimax_design,
    solve_certified_measurement_design,
)
from hydrosheaf.nuclear.ttd_design import (
    TtdHypothesisEnsemble,
    rank_ttd_measurements,
)
from hydrosheaf.nuclear.ttd_identified import (
    AgeFunctional,
    TracerConstraint,
    solve_ttd_identified_set,
)


@pytest.fixture
def age_grid_5():
    """Simple 5-bin age grid: [5, 20, 50, 80, 150] years."""
    return np.asarray([5.0, 20.0, 50.0, 80.0, 150.0], dtype=float)


@pytest.fixture
def young_water_functional():
    """Target functional: Fraction <= 50 years (bins 0, 1, 2)."""
    return AgeFunctional(
        name="young_water_fraction",
        coefficients=[1.0, 1.0, 1.0, 0.0, 0.0],
        units="fraction",
    )


def test_certified_candidate_input_validation():
    with pytest.raises(ValueError, match="option_id must be non-empty"):
        CertifiedCandidateTracer("", "3H", 2020.0, error_bound=0.5)

    with pytest.raises(ValueError, match="tracer must be non-empty"):
        CertifiedCandidateTracer("opt1", "", 2020.0, error_bound=0.5)

    with pytest.raises(ValueError, match="sample_year must be finite"):
        CertifiedCandidateTracer("opt1", "3H", float("nan"), error_bound=0.5)

    with pytest.raises(ValueError, match="error_bound must be finite and strictly positive"):
        CertifiedCandidateTracer("opt1", "3H", 2020.0, error_bound=-0.1)

    with pytest.raises(ValueError, match="cost must be finite and strictly positive"):
        CertifiedCandidateTracer("opt1", "3H", 2020.0, error_bound=0.5, cost=0.0)

    cand = CertifiedCandidateTracer(
        "opt1", "3H", 2020.0, error_bound=0.5, cost=2.0, response=[1.0, 0.5, 0.2]
    )
    with pytest.raises(ValueError, match="does not match age grid size"):
        cand.resolve_response(np.array([10.0, 20.0]))


def test_empty_measurement_boundary(age_grid_5, young_water_functional):
    # Set up prior constraint: an initial tracer observation bounding mass
    # h = [1.0, 0.8, 0.5, 0.2, 0.05]
    c1 = TracerConstraint(
        tracer="3H",
        response=[1.0, 0.8, 0.5, 0.2, 0.05],
        observed=0.50,
        sigma=0.05,
    )
    # Solve standard identified set to get sharp upper and lower bounds on young_water_fraction
    report = solve_ttd_identified_set(
        age_grid_5,
        [c1],
        [young_water_functional],
        sigma_multiplier=1.96,
    )
    assert "young_water_fraction" in report.bounds
    bound = report.bounds["young_water_fraction"]
    expected_width = bound.upper - bound.lower

    # Now evaluate worst-case ambiguity W(empty)
    eval_empty = evaluate_worst_case_ambiguity(
        age_grid_5,
        [c1],
        (),
        young_water_functional,
        sigma_multiplier=1.96,
    )
    assert eval_empty.status == "FEASIBLE"
    assert math.isclose(eval_empty.worst_case_ambiguity, expected_width, rel_tol=1e-5, abs_tol=1e-6)


def test_ambiguity_monotonicity(age_grid_5, young_water_functional):
    c_prior = TracerConstraint(
        tracer="3H",
        response=[1.0, 0.8, 0.5, 0.2, 0.05],
        observed=0.45,
        sigma=0.04,
    )

    cand1 = CertifiedCandidateTracer(
        option_id="SF6_2020",
        tracer="SF6",
        sample_year=2020.0,
        error_bound=0.02,
        cost=3.0,
        response=[0.9, 0.6, 0.2, 0.0, 0.0],
    )
    cand2 = CertifiedCandidateTracer(
        option_id="14C_2020",
        tracer="14C",
        sample_year=2020.0,
        error_bound=0.03,
        cost=5.0,
        response=[1.0, 0.98, 0.95, 0.80, 0.30],
    )

    w_empty = evaluate_worst_case_ambiguity(
        age_grid_5, [c_prior], (), young_water_functional
    ).worst_case_ambiguity

    w_1 = evaluate_worst_case_ambiguity(
        age_grid_5, [c_prior], [cand1], young_water_functional
    ).worst_case_ambiguity

    w_2 = evaluate_worst_case_ambiguity(
        age_grid_5, [c_prior], [cand2], young_water_functional
    ).worst_case_ambiguity

    w_12 = evaluate_worst_case_ambiguity(
        age_grid_5, [c_prior], [cand1, cand2], young_water_functional
    ).worst_case_ambiguity

    assert w_1 <= w_empty + 1e-7
    assert w_2 <= w_empty + 1e-7
    assert w_12 <= w_1 + 1e-7
    assert w_12 <= w_2 + 1e-7


def test_inconsistent_evidence_detection(age_grid_5, young_water_functional):
    # Construct mutually impossible prior constraints
    c1 = TracerConstraint(
        tracer="T1",
        response=[1.0, 0.0, 0.0, 0.0, 0.0],
        observed=0.90,
        sigma=0.01,
    )
    c2 = TracerConstraint(
        tracer="T2",
        response=[1.0, 0.0, 0.0, 0.0, 0.0],
        observed=0.10,
        sigma=0.01,
    )

    cert = solve_certified_measurement_design(
        age_grid_5,
        [c1, c2],
        [CertifiedCandidateTracer("C1", "3H", 2020.0, 0.1, response=[1, 1, 1, 1, 1])],
        young_water_functional,
        target_tolerance=0.10,
    )
    assert cert.status == "INCONSISTENT_EVIDENCE"
    assert math.isnan(cert.achieved_ambiguity)


def test_already_resolved_detection(age_grid_5, young_water_functional):
    # Prior constraint already tightly pins young water
    c1 = TracerConstraint(
        tracer="T1",
        response=[1.0, 1.0, 1.0, 0.0, 0.0],
        observed=0.80,
        sigma=0.005,
    )

    cand = CertifiedCandidateTracer(
        option_id="cand1",
        tracer="3H",
        sample_year=2020.0,
        error_bound=0.1,
        cost=5.0,
        response=[0.5, 0.5, 0.5, 0.1, 0.1],
    )

    # Tolerance delta = 0.20 is wider than initial ambiguity (which is ~ 2 * 1.96 * 0.005 = 0.0196)
    cert = solve_certified_measurement_design(
        age_grid_5,
        [c1],
        [cand],
        young_water_functional,
        target_tolerance=0.20,
    )
    assert cert.status == "ALREADY_RESOLVED"
    assert cert.total_cost == 0.0
    assert len(cert.selected_option_ids) == 0
    assert cert.achieved_ambiguity <= 0.20


def test_impossibility_witness_detection(age_grid_5, young_water_functional):
    # Set up candidate tracers whose response is constant across the entire age grid
    # A constant response cannot separate young from old water!
    cand_blind = CertifiedCandidateTracer(
        option_id="blind_tracer",
        tracer="blind",
        sample_year=2020.0,
        error_bound=0.01,
        cost=5.0,
        response=[1.0, 1.0, 1.0, 1.0, 1.0],  # exactly 1 * sum(x) = 1 always
    )

    cert = solve_certified_measurement_design(
        age_grid_5,
        (),  # no prior constraints, so P is full simplex
        [cand_blind],
        young_water_functional,
        target_tolerance=0.05,  # request tight 5% resolution
    )

    assert cert.status == "IMPOSSIBILITY_WITNESS"
    assert cert.achieved_ambiguity > 0.05
    assert cert.lower_witness is not None
    assert cert.upper_witness is not None
    # Witness test: young water functional difference between witnesses should equal achieved ambiguity
    diff = young_water_functional.coefficients @ cert.upper_witness - young_water_functional.coefficients @ cert.lower_witness
    assert math.isclose(diff, cert.achieved_ambiguity, abs_tol=1e-5)


def test_cost_optimal_sufficiency_certificate(age_grid_5, young_water_functional):
    # We want young water fraction <= 50 yr resolved to within delta = 0.15
    # Initially without measurements, P is simplex, ambiguity = 1.0
    # Cand_expensive (cost 15): error bound 0.05, highly sensitive to young water
    cand_expensive = CertifiedCandidateTracer(
        option_id="tracer_expensive",
        tracer="T_exp",
        sample_year=2020.0,
        error_bound=0.05,
        cost=15.0,
        response=[1.0, 1.0, 1.0, 0.0, 0.0],
    )
    # Cand_cheap (cost 5): error bound 0.05, also sensitive to young water
    cand_cheap = CertifiedCandidateTracer(
        option_id="tracer_cheap",
        tracer="T_cheap",
        sample_year=2020.0,
        error_bound=0.05,
        cost=5.0,
        response=[1.0, 1.0, 1.0, 0.0, 0.0],
    )
    # Cand_uninformative (cost 2): response insensitive to young water
    cand_uninformative = CertifiedCandidateTracer(
        option_id="tracer_uninf",
        tracer="T_uninf",
        sample_year=2020.0,
        error_bound=0.05,
        cost=2.0,
        response=[0.2, 0.2, 0.2, 0.2, 0.2],
    )

    cert = solve_certified_measurement_design(
        age_grid_5,
        (),
        [cand_expensive, cand_cheap, cand_uninformative],
        young_water_functional,
        target_tolerance=0.15,
    )

    assert cert.status == "CERTIFIED_SUFFICIENT"
    # It should choose tracer_cheap alone!
    assert cert.selected_option_ids == ("tracer_cheap",)
    assert cert.total_cost == 5.0
    assert cert.achieved_ambiguity <= 0.15
    assert len(cert.certificate_hash) == 64


def test_budget_constrained_minimax_design(age_grid_5, young_water_functional):
    # Precision tradeoff candidates
    cand_a = CertifiedCandidateTracer(
        option_id="A",
        tracer="A",
        sample_year=2020.0,
        error_bound=0.10,
        cost=4.0,
        response=[1.0, 1.0, 1.0, 0.0, 0.0],
    )
    cand_b = CertifiedCandidateTracer(
        option_id="B",
        tracer="B",
        sample_year=2020.0,
        error_bound=0.02,
        cost=6.0,
        response=[1.0, 1.0, 1.0, 0.0, 0.0],
    )

    # Budget = 3.0: cannot afford A (4) or B (6), returns ALREADY_RESOLVED (no action within budget)
    cert_b3 = solve_budgeted_minimax_design(
        age_grid_5, (), [cand_a, cand_b], young_water_functional, budget=3.0
    )
    assert cert_b3.selected_option_ids == ()
    assert cert_b3.total_cost == 0.0

    # Budget = 5.0: can afford A (4), achieving ambiguity 0.20
    cert_b5 = solve_budgeted_minimax_design(
        age_grid_5, (), [cand_a, cand_b], young_water_functional, budget=5.0
    )
    assert cert_b5.selected_option_ids == ("A",)
    assert cert_b5.total_cost == 4.0
    assert cert_b5.achieved_ambiguity <= 0.20 + 1e-6

    # Budget = 7.0: can afford B (6), achieving lower ambiguity 0.04
    cert_b7 = solve_budgeted_minimax_design(
        age_grid_5, (), [cand_a, cand_b], young_water_functional, budget=7.0
    )
    assert cert_b7.selected_option_ids == ("B",)
    assert cert_b7.total_cost == 6.0
    assert cert_b7.achieved_ambiguity <= 0.04 + 1e-6
    assert cert_b7.achieved_ambiguity < cert_b5.achieved_ambiguity

    # Complementary synergy: X covers [0, 1], Y covers [2]
    cand_x = CertifiedCandidateTracer(
        option_id="X",
        tracer="X",
        sample_year=2020.0,
        error_bound=0.02,
        cost=4.0,
        response=[1.0, 1.0, 0.0, 0.0, 0.0],
    )
    cand_y = CertifiedCandidateTracer(
        option_id="Y",
        tracer="Y",
        sample_year=2020.0,
        error_bound=0.02,
        cost=5.0,
        response=[0.0, 0.0, 1.0, 0.0, 0.0],
    )
    # With budget 6, neither alone can reduce worst-case ambiguity on young water (bins 0,1,2)
    cert_single = solve_budgeted_minimax_design(
        age_grid_5, (), [cand_x, cand_y], young_water_functional, budget=6.0
    )
    assert cert_single.selected_option_ids == ()

    # With budget 10, both can be purchased, certifying a drop in ambiguity
    cert_synergy = solve_budgeted_minimax_design(
        age_grid_5, (), [cand_x, cand_y], young_water_functional, budget=10.0
    )
    assert set(cert_synergy.selected_option_ids) == {"X", "Y"}
    assert cert_synergy.total_cost == 9.0
    assert cert_synergy.achieved_ambiguity <= 0.10


def test_real_atmospheric_tracer_certified_design():
    """Integration test using actual atmospheric response kernels (3H, SF6, 14C)."""
    age_grid = np.linspace(0.5, 90.5, 30)

    # Anthropocene target: fraction <= 50 yr
    target = AgeFunctional(
        name="anthropocene_fraction",
        coefficients=(age_grid <= 50.0).astype(float),
    )

    # Prior measurement: 3H measured in 2010 with observed = 4.2 TU, sigma = 0.4 TU
    kernel_3h_2010 = tracer_response_kernel("3H", age_grid, 2010.0)
    c_prior = TracerConstraint(
        tracer="3H",
        response=kernel_3h_2010,
        observed=4.2,
        sigma=0.4,
    )

    # Candidate future measurements in 2024:
    # 1. 3H in 2024 (cost 100)
    cand_3h_2024 = CertifiedCandidateTracer(
        option_id="3H_2024",
        tracer="3H",
        sample_year=2024.0,
        error_bound=0.4,
        cost=100.0,
    )
    # 2. SF6 in 2024 (cost 150)
    cand_sf6_2024 = CertifiedCandidateTracer(
        option_id="SF6_2024",
        tracer="SF6",
        sample_year=2024.0,
        error_bound=0.2,
        cost=150.0,
    )
    # 3. 14C in 2024 (cost 350)
    cand_14c_2024 = CertifiedCandidateTracer(
        option_id="14C_2024",
        tracer="14C",
        sample_year=2024.0,
        error_bound=2.0,
        cost=350.0,
    )

    cert = solve_certified_measurement_design(
        age_grid,
        [c_prior],
        [cand_3h_2024, cand_sf6_2024, cand_14c_2024],
        target,
        target_tolerance=0.40,
    )

    assert cert.feasibility_status == "FEASIBLE"
    assert cert.status in {"CERTIFIED_SUFFICIENT", "IMPOSSIBILITY_WITNESS"}
    assert cert.initial_ambiguity > 0.0
    assert cert.achieved_ambiguity <= cert.initial_ambiguity
    assert cert.total_cost >= 0.0


def test_probability_gate_links_to_certified_design():
    """Verify that ttd_design's _probability_gate points user to certified design."""
    ensemble = TtdHypothesisEnsemble(
        hypothesis_ids=("H1", "H2"),
        age_grid_years=(10.0, 50.0, 100.0),
        masses=((0.5, 0.3, 0.2), (0.2, 0.5, 0.3)),
        probabilities=None,
        probability_semantics=None,
    )
    result = rank_ttd_measurements(ensemble, [])
    assert result["status"] == "ABSTAIN"
    assert result["reason"] == "no_probability_model"
    assert result["certified_design_available"] is True
    assert "solve_certified_measurement_design" in result["certified_design_recommendation"]


def test_budgeted_minimax_target_tolerance_semantics(age_grid_5, young_water_functional):
    """Verify that solve_budgeted_minimax_design properly validates against declared target tolerance."""
    cand_a = CertifiedCandidateTracer(
        option_id="A",
        tracer="A",
        sample_year=2020.0,
        error_bound=0.10,
        cost=4.0,
        response=[1.0, 1.0, 1.0, 0.0, 0.0],
    )
    cand_b = CertifiedCandidateTracer(
        option_id="B",
        tracer="B",
        sample_year=2020.0,
        error_bound=0.02,
        cost=6.0,
        response=[1.0, 1.0, 1.0, 0.0, 0.0],
    )

    # Declare target tolerance = 0.05
    target_tol = 0.05

    # Budget 0: cannot afford any candidate. Ambiguity is 1.0 > 0.05. Status must be NO_ACTION_AFFORDABLE.
    cert_0 = solve_budgeted_minimax_design(
        age_grid_5, (), [cand_a, cand_b], young_water_functional,
        budget=0.0, target_tolerance=target_tol,
    )
    assert cert_0.status == "NO_ACTION_AFFORDABLE"
    assert cert_0.target_tolerance == target_tol
    assert cert_0.selected_option_ids == ()
    
    # Budget 5.0: can afford A (cost 4). Achieves ambiguity 0.20.
    # 0.20 is an improvement over 1.0, but still exceeds target 0.05!
    # Status must be BUDGET_EXHAUSTED_INSUFFICIENT, not CERTIFIED_SUFFICIENT!
    cert_5 = solve_budgeted_minimax_design(
        age_grid_5, (), [cand_a, cand_b], young_water_functional,
        budget=5.0, target_tolerance=target_tol,
    )
    assert cert_5.status == "BUDGET_EXHAUSTED_INSUFFICIENT"
    assert cert_5.target_tolerance == target_tol
    assert cert_5.selected_option_ids == ("A",)
    assert cert_5.achieved_ambiguity <= 0.20 + 1e-6

    # Budget 7.0: can afford B (cost 6). Achieves ambiguity 0.04 <= 0.05.
    # Status must be CERTIFIED_SUFFICIENT!
    cert_7 = solve_budgeted_minimax_design(
        age_grid_5, (), [cand_a, cand_b], young_water_functional,
        budget=7.0, target_tolerance=target_tol,
    )
    assert cert_7.status == "CERTIFIED_SUFFICIENT"
    assert cert_7.target_tolerance == target_tol
    assert cert_7.selected_option_ids == ("B",)
    assert cert_7.achieved_ambiguity <= target_tol


def test_solver_numerical_error_not_misclassified_as_infeasible(age_grid_5, young_water_functional, monkeypatch):
    from unittest.mock import MagicMock
    import hydrosheaf.nuclear.ttd_certified_design as tcd

    # Mock linprog to simulate an iteration limit failure (status=1)
    mock_res = MagicMock()
    mock_res.success = False
    mock_res.status = 1  # iteration limit reached, NOT infeasible (status=2)
    mock_res.message = "Iteration limit reached."

    monkeypatch.setattr(tcd, "linprog", lambda *args, **kwargs: mock_res)

    c_prior = TracerConstraint(tracer="3H", response=[1, 0, 0, 0, 0], observed=0.5, sigma=0.1)
    cand = CertifiedCandidateTracer("C1", "SF6", 2020.0, 0.1, response=[1, 1, 1, 1, 1])

    eval_res = evaluate_worst_case_ambiguity(age_grid_5, [c_prior], (), young_water_functional)
    assert eval_res.status == "NUMERICAL_ERROR_1", f"Expected NUMERICAL_ERROR_1, got {eval_res.status}"

    cert = solve_certified_measurement_design(
        age_grid_5, [c_prior], [cand], young_water_functional, target_tolerance=0.1
    )
    assert cert.status == "NUMERICAL_ERROR", f"Expected NUMERICAL_ERROR, got {cert.status}"
    assert cert.status != "INCONSISTENT_EVIDENCE", "Numerical error must not be misclassified as inconsistent evidence!"


def test_all_candidate_numerical_error_abstains_from_certificate(
    age_grid_5, young_water_functional, monkeypatch
):
    """A failed all-candidate LP must not become an impossibility certificate."""
    import hydrosheaf.nuclear.ttd_certified_design as tcd

    candidate = CertifiedCandidateTracer(
        "C1", "SF6", 2020.0, 0.1, response=[1, 1, 1, 1, 1]
    )
    feasible = tcd.AmbiguityEvaluation(
        subset_ids=(),
        worst_case_ambiguity=1.0,
        status="FEASIBLE",
        lower_witness=[0, 0, 0, 0, 1],
        upper_witness=[1, 0, 0, 0, 0],
    )
    failed = tcd.AmbiguityEvaluation(
        subset_ids=("C1",),
        worst_case_ambiguity=float("nan"),
        status="NUMERICAL_ERROR_1",
    )
    evaluations = iter((feasible, failed))
    monkeypatch.setattr(tcd, "evaluate_worst_case_ambiguity", lambda *args, **kwargs: next(evaluations))

    cert = tcd.solve_certified_measurement_design(
        age_grid_5,
        (),
        [candidate],
        young_water_functional,
        target_tolerance=0.1,
    )

    assert cert.status == "NUMERICAL_ERROR"
    assert cert.feasibility_status == "NUMERICAL_ERROR_1"
    assert cert.selected_option_ids == ()
    assert cert.metadata["all_candidates_evaluation_status"] == "NUMERICAL_ERROR_1"


def test_budgeted_all_candidate_numerical_error_abstains_from_budget_result(
    age_grid_5, young_water_functional, monkeypatch
):
    """Budgeted design must preserve an unverified all-candidate reference solve."""
    import hydrosheaf.nuclear.ttd_certified_design as tcd

    candidate = CertifiedCandidateTracer(
        "C1", "SF6", 2020.0, 0.1, response=[1, 1, 1, 1, 1]
    )
    feasible = tcd.AmbiguityEvaluation(
        subset_ids=(),
        worst_case_ambiguity=1.0,
        status="FEASIBLE",
        lower_witness=[0, 0, 0, 0, 1],
        upper_witness=[1, 0, 0, 0, 0],
    )
    failed = tcd.AmbiguityEvaluation(
        subset_ids=("C1",),
        worst_case_ambiguity=float("nan"),
        status="NUMERICAL_ERROR_1",
    )
    evaluations = iter((feasible, failed))
    monkeypatch.setattr(tcd, "evaluate_worst_case_ambiguity", lambda *args, **kwargs: next(evaluations))

    cert = tcd.solve_budgeted_minimax_design(
        age_grid_5,
        (),
        [candidate],
        young_water_functional,
        budget=10.0,
        target_tolerance=0.1,
    )

    assert cert.status == "NUMERICAL_ERROR"
    assert cert.feasibility_status == "NUMERICAL_ERROR_1"
    assert cert.selected_option_ids == ()
    assert cert.metadata["all_candidates_evaluation_status"] == "NUMERICAL_ERROR_1"


def test_budgeted_design_accepts_shared_campaign_cost_function(age_grid_5, young_water_functional):
    """Shared field/shipping charges must be counted once per selected campaign."""
    cand_a = CertifiedCandidateTracer(
        option_id="A_shared",
        tracer="A",
        sample_year=2020.0,
        error_bound=0.10,
        cost=4.0,
        response=[1.0, 1.0, 1.0, 0.0, 0.0],
    )
    cand_b = CertifiedCandidateTracer(
        option_id="B_shared",
        tracer="B",
        sample_year=2020.0,
        error_bound=0.10,
        cost=5.0,
        response=[1.0, 1.0, 1.0, 0.0, 0.0],
    )

    def campaign_cost(selected):
        # A fixed visit/shipment cost is paid once whenever any tracer is
        # selected, rather than once per tracer.
        return (7.0 if selected else 0.0) + sum(c.cost for c in selected)

    cert = solve_budgeted_minimax_design(
        age_grid_5,
        (),
        [cand_a, cand_b],
        young_water_functional,
        budget=12.0,
        target_tolerance=0.05,
        cost_function=campaign_cost,
    )

    assert set(cert.selected_option_ids) != {"A_shared", "B_shared"}
    expected = campaign_cost(tuple(c for c in (cand_a, cand_b) if c.option_id in cert.selected_option_ids))
    assert cert.total_cost == expected
    assert cert.total_cost <= 12.0
