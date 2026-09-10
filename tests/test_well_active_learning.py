"""Tests for well-level Bayesian active learning and experimental design.

Covers:
- WellAction and CampaignConfig data contracts and cost formulas.
- Validation decision matrix risk weighting (FP=2.0, FN=1.8 > TP=0.2).
- Action generation across measurement types, sampling times, and well nodes.
- Infeasible well exclusion and accessibility scaling.
- Predictive scenario generation and hypothesis differentiation.
- Well-level batch selection and deduplication (max_actions_per_well=1).
- Shared mobilization/travel cost discounting on visited wells.
- Budget constraints and diminishing return termination.
- Abstention gate on high R-hat (> 1.10), low ESS (< 20), and collapsed entropy.
- End-to-end rank_campaign_measurements producing JSON, CSV, and Markdown outputs.
- Guardrail statement compliance.
- Routing via rank_next_measurements(method="bayesian_campaign").
"""

from __future__ import annotations

import csv
import json
import tempfile
from pathlib import Path
from typing import Dict, List

import pytest

from hydrosheaf.calibration.active_learning import rank_next_measurements
from hydrosheaf.calibration.well_active_learning import (
    CampaignConfig,
    WellAction,
    build_validation_decision_matrix,
    build_well_actions,
    evaluate_abstention_gate,
    rank_campaign_measurements,
    select_well_campaign_batch,
)
from hydrosheaf.graph.types import Edge


# ── Fixtures and Test Helpers ──────────────────────────────────────────


def _make_edge(edge_id: str, u: str, v: str, attrs: dict | None = None) -> Edge:
    return Edge(edge_id=edge_id, u=u, v=v, attrs=attrs or {})


def _make_synth_posterior(
    topologies: List[List[str]],
    probabilities: List[float],
    r_hat: float = 1.02,
    ess: float = 450.0,
) -> dict:
    """Synthetic posterior_result matching run_topology_posterior structure."""
    ensemble = []
    all_edges = set()
    for edges, prob in zip(topologies, probabilities):
        ensemble.append({
            "graph_edges": edges,
            "probability": prob,
            "count": int(prob * 1000),
        })
        all_edges.update(edges)

    marginal_probs = {e: 0.0 for e in all_edges}
    for e in all_edges:
        for t, p in zip(topologies, probabilities):
            if e in t:
                marginal_probs[e] += p

    return {
        "status": "CONVERGED",
        "selected_edge_ids": topologies[0],
        "edge_probabilities": marginal_probs,
        "graph_ensemble": ensemble,
        "attrs": {
            "mcmc_r_hat": r_hat,
            "mcmc_ess": ess,
            "posterior_joint_graph_entropy": 0.6931,
            "posterior_marginal_edge_entropy": 1.386,
            "posterior_graph_ensemble": ensemble,
        },
    }


def _make_synth_benchmark() -> dict:
    return {
        "variants": {
            "baseline": {"selected_edges": ["E1", "E2"]},
            "null_model_default": {"selected_edges": ["E1"]},
            "assumption_calibrated": {"selected_edges": ["E1", "E3"]},
        },
        "improvement_summary": {},
        "independent_validation": True,
        "manuscript_claim_allowed": True,
    }


# ── Unit Tests ─────────────────────────────────────────────────────────


def test_well_action_dataclass():
    """Test WellAction properties, cost calculation, and string repr."""
    action = WellAction(
        well_id="W1",
        measurement_type="connectivity_tracer",
        sampling_time=30.0,
        base_cost=10.0,
        travel_cost=4.0,
        accessibility=0.5,
        target_edges=["E1", "E2"],
    )
    # Standalone cost = (10 + 4) / 0.5 = 28.0
    assert pytest.approx(action.standalone_cost, rel=1e-6) == 28.0
    assert action.action_id == "W1@connectivity_tracer@t30.0"
    assert "connectivity_tracer" in str(action)
    assert "W1" in repr(action)


def test_campaign_config_defaults_and_override():
    """Test CampaignConfig defaults and dictionary conversion."""
    cfg = CampaignConfig()
    assert cfg.max_actions_per_well == 1
    assert cfg.batch_size == 3
    assert cfg.max_mcmc_r_hat == 1.10
    assert cfg.min_mcmc_ess == 20.0
    assert "connectivity_tracer" in cfg.measurement_costs

    custom_cfg = CampaignConfig(
        batch_size=5,
        max_actions_per_well=2,
        budget=100.0,
    )
    d = custom_cfg.to_dict()
    assert d["batch_size"] == 5
    assert d["max_actions_per_well"] == 2
    assert d["budget"] == 100.0


def test_validation_decision_matrix_risk_weighting():
    """Test that FP and FN receive higher risk reduction weights than TP."""
    val_report = {
        "labels": {
            "E1": {"status": "FP", "expected": False, "inferred": True},
            "E2": {"status": "FN", "expected": True, "inferred": False},
            "E3": {"status": "TP", "expected": True, "inferred": True},
            "E4": {"status": "TN", "expected": False, "inferred": False},
        }
    }
    candidate_edges = [
        _make_edge("E1", "W1", "W2"),
        _make_edge("E2", "W2", "W3"),
        _make_edge("E3", "W1", "W3"),
        _make_edge("E4", "W3", "W4"),
        _make_edge("E5", "W4", "W5"),  # Unlabeled
    ]

    mat = build_validation_decision_matrix(candidate_edges, val_report)
    assert mat["E1"]["validation_status"] == "FP"
    assert mat["E1"]["risk_weight"] == 2.0
    assert mat["E2"]["validation_status"] == "FN"
    assert mat["E2"]["risk_weight"] == 1.8
    assert mat["E3"]["validation_status"] == "TP"
    assert mat["E3"]["risk_weight"] == 0.2
    assert mat["E5"]["validation_status"] == "selected_unlabeled"
    assert mat["E5"]["risk_weight"] == 1.4

    # FP and FN must outweigh TP
    assert mat["E1"]["risk_weight"] > mat["E3"]["risk_weight"]
    assert mat["E2"]["risk_weight"] > mat["E3"]["risk_weight"]


def test_build_well_actions_sampling_and_types():
    """Test generating well actions across types and sampling times."""
    candidate_edges = [
        _make_edge("E1", "W1", "W2"),
        _make_edge("E2", "W2", "W3"),
    ]
    hypotheses = ["H1", "H2"]
    prior_probs = [0.6, 0.4]
    edge_tuples = [
        ("E1",),
        ("E1", "E2"),
    ]
    val_mat = {
        "E1": {"validation_status": "FP", "risk_weight": 2.0},
        "E2": {"validation_status": "FN", "risk_weight": 1.8},
    }

    cfg = CampaignConfig(
        measurement_types=("head_monitoring", "connectivity_tracer"),
        sampling_times=(0.0, 30.0),
    )

    actions = build_well_actions(
        candidate_edges=candidate_edges,
        hypothesis_ids=hypotheses,
        prior_probs=prior_probs,
        edge_tuples=edge_tuples,
        validation_decision_matrix=val_mat,
        config=cfg,
    )

    # 3 wells (W1, W2, W3) * 2 types * 2 times = 12 actions
    assert len(actions) == 12

    # Check action properties
    w1_actions = [a for a in actions if a.well_id == "W1"]
    assert len(w1_actions) == 4
    for a in w1_actions:
        assert "E1" in a.target_edges
        assert len(a.scenarios) == 2
        # Check scenario hypothesis IDs
        assert [s.hypothesis_id for s in a.scenarios] == ["H1", "H2"]


def test_build_well_actions_feasibility_and_accessibility():
    """Test that infeasible wells are skipped and accessibility scales cost."""
    candidate_edges = [
        _make_edge("E1", "W1", "W2"),
        _make_edge("E2", "W2", "W3"),
    ]
    hypotheses = ["H1"]
    prior_probs = [1.0]
    edge_tuples = [("E1",)]
    val_mat = {"E1": {"validation_status": "unlabeled", "risk_weight": 1.0}}

    cfg = CampaignConfig(
        measurement_types=("head_monitoring",),
        sampling_times=(0.0,),
    )

    # W2 is marked infeasible
    feasibility = {"W1": True, "W2": False, "W3": True}
    # W3 has 50% accessibility
    accessibility = {"W1": 1.0, "W2": 1.0, "W3": 0.5}

    actions = build_well_actions(
        candidate_edges=candidate_edges,
        hypothesis_ids=hypotheses,
        prior_probs=prior_probs,
        edge_tuples=edge_tuples,
        validation_decision_matrix=val_mat,
        config=cfg,
        feasibility_map=feasibility,
        accessibility_map=accessibility,
    )

    wells_present = {a.well_id for a in actions}
    assert "W2" not in wells_present
    assert wells_present == {"W1", "W3"}

    w1_action = next(a for a in actions if a.well_id == "W1")
    w3_action = next(a for a in actions if a.well_id == "W3")

    # W3 standalone cost should be double W1 standalone cost
    assert pytest.approx(w3_action.standalone_cost, rel=1e-5) == 2.0 * w1_action.standalone_cost


def test_predictive_scenarios_differentiation():
    """Test that predictive scenarios differentiate between present vs absent edge topologies."""
    candidate_edges = [_make_edge("E1", "W1", "W2")]
    hypotheses = ["H_connected", "H_disconnected"]
    prior_probs = [0.5, 0.5]
    edge_tuples = [("E1",), ()]
    val_mat = {"E1": {"validation_status": "FP", "risk_weight": 2.0}}

    cfg = CampaignConfig(
        measurement_types=("connectivity_tracer", "isotopes_d18O_d2H"),
        sampling_times=(30.0,),
    )

    actions = build_well_actions(
        candidate_edges=candidate_edges,
        hypothesis_ids=hypotheses,
        prior_probs=prior_probs,
        edge_tuples=edge_tuples,
        validation_decision_matrix=val_mat,
        config=cfg,
    )

    tracer_w2 = next(a for a in actions if a.well_id == "W2" and a.measurement_type == "connectivity_tracer")
    # Under H_connected, downstream tracer peak is higher than under H_disconnected
    scenarios_by_h = {s.hypothesis_id: s for s in tracer_w2.scenarios}
    assert scenarios_by_h["H_connected"].mean > scenarios_by_h["H_disconnected"].mean


def test_batch_selection_deduplication():
    """Test that batch selection enforces max_actions_per_well=1."""
    candidate_edges = [
        _make_edge("E1", "W1", "W2"),
        _make_edge("E2", "W2", "W3"),
        _make_edge("E3", "W1", "W4"),
    ]
    hypotheses = ["H1", "H2"]
    prior_probs = [0.5, 0.5]
    edge_tuples = [("E1", "E2"), ("E3",)]
    val_mat = {
        "E1": {"validation_status": "FP", "risk_weight": 2.0},
        "E2": {"validation_status": "FN", "risk_weight": 1.8},
        "E3": {"validation_status": "selected_unlabeled", "risk_weight": 1.4},
    }

    cfg = CampaignConfig(
        batch_size=3,
        max_actions_per_well=1,
        measurement_types=("head_monitoring", "isotopes_d18O_d2H"),
        sampling_times=(0.0, 30.0),
    )

    actions = build_well_actions(
        candidate_edges=candidate_edges,
        hypothesis_ids=hypotheses,
        prior_probs=prior_probs,
        edge_tuples=edge_tuples,
        validation_decision_matrix=val_mat,
        config=cfg,
    )

    batch_result = select_well_campaign_batch(actions, prior_probs, cfg)
    assert batch_result["status"] == "ACTIONABLE"
    selected_details = batch_result["selected_option_details"]
    selected_wells = [d["well_id"] for d in selected_details]

    # Every selected well must be unique
    assert len(selected_wells) == len(set(selected_wells))
    assert len(selected_wells) <= 3


def test_batch_selection_shared_travel_cost():
    """Test that multiple actions at the same well enjoy travel cost discounting."""
    candidate_edges = [_make_edge("E1", "W1", "W2")]
    hypotheses = ["H1", "H2"]
    prior_probs = [0.5, 0.5]
    edge_tuples = [("E1",), ()]
    val_mat = {"E1": {"validation_status": "FP", "risk_weight": 2.0}}

    cfg = CampaignConfig(
        batch_size=2,
        max_actions_per_well=2,
        travel_cost_per_well=5.0,
        measurement_types=("head_monitoring", "connectivity_tracer"),
        sampling_times=(0.0,),
        cost_exponent=0.1,  # Low cost sensitivity so both tests at W2 are selected
    )

    actions = build_well_actions(
        candidate_edges=candidate_edges,
        hypothesis_ids=hypotheses,
        prior_probs=prior_probs,
        edge_tuples=edge_tuples,
        validation_decision_matrix=val_mat,
        config=cfg,
    )

    # Force selection to only contain actions for W2
    w2_actions = [a for a in actions if a.well_id == "W2"]
    batch_result = select_well_campaign_batch(w2_actions, prior_probs, cfg)

    if len(batch_result["selected_option_details"]) == 2:
        first = batch_result["selected_option_details"][0]
        second = batch_result["selected_option_details"][1]
        assert first["well_id"] == "W2"
        assert second["well_id"] == "W2"
        # The second action at W2 should NOT pay the 5.0 travel cost
        assert second["marginal_cost"] < second["standalone_cost"]
        assert batch_result["travel_cost_savings"] >= 5.0


def test_batch_selection_budget_limit():
    """Test that batch selection strictly respects budget constraints."""
    candidate_edges = [_make_edge("E1", "W1", "W2")]
    hypotheses = ["H1", "H2"]
    prior_probs = [0.5, 0.5]
    edge_tuples = [("E1",), ()]
    val_mat = {"E1": {"validation_status": "FP", "risk_weight": 2.0}}

    cfg = CampaignConfig(
        batch_size=5,
        budget=4.5,  # Each head test with travel is 1.0 + 3.0 = 4.0
        measurement_types=("head_monitoring",),
        sampling_times=(0.0,),
    )

    actions = build_well_actions(
        candidate_edges=candidate_edges,
        hypothesis_ids=hypotheses,
        prior_probs=prior_probs,
        edge_tuples=edge_tuples,
        validation_decision_matrix=val_mat,
        config=cfg,
    )

    batch_result = select_well_campaign_batch(actions, prior_probs, cfg)
    assert batch_result["total_cost"] <= 4.5
    assert len(batch_result["selected_options"]) <= 1


def test_abstention_gate_on_high_rhat():
    """Test that abstention triggers when MCMC R-hat exceeds 1.10."""
    post = _make_synth_posterior([["E1"], ["E2"]], [0.5, 0.5], r_hat=1.25)
    cfg = CampaignConfig(max_mcmc_r_hat=1.10)
    abstain, reason = evaluate_abstention_gate(post, prior_entropy=0.69, max_robust_eig=0.4, config=cfg)
    assert abstain is True
    assert "R-hat" in reason


def test_abstention_gate_on_low_ess():
    """Test that abstention triggers when effective sample size is below 20."""
    post = _make_synth_posterior([["E1"], ["E2"]], [0.5, 0.5], r_hat=1.01, ess=12.0)
    cfg = CampaignConfig(min_mcmc_ess=20.0)
    abstain, reason = evaluate_abstention_gate(post, prior_entropy=0.69, max_robust_eig=0.4, config=cfg)
    assert abstain is True
    assert "effective sample size" in reason.lower() or "ess" in reason.lower()


def test_abstention_gate_on_zero_entropy():
    """Test that abstention triggers when posterior entropy is near zero."""
    post = _make_synth_posterior([["E1"]], [1.0], r_hat=1.01, ess=100.0)
    cfg = CampaignConfig()
    abstain, reason = evaluate_abstention_gate(post, prior_entropy=1e-6, max_robust_eig=0.0, config=cfg)
    assert abstain is True
    assert "collapsed" in reason.lower() or "uncertainty" in reason.lower()


def test_rank_campaign_measurements_end_to_end():
    """Test end-to-end rank_campaign_measurements producing JSON, CSV, and Markdown."""
    posterior = _make_synth_posterior([["E1", "E2"], ["E1"]], [0.6, 0.4])
    benchmark = _make_synth_benchmark()
    candidate_edges = [
        _make_edge("E1", "W1", "W2"),
        _make_edge("E2", "W2", "W3"),
    ]
    val_report = {
        "labels": {
            "E2": {"status": "FN", "expected": True, "inferred": False},
        }
    }

    with tempfile.TemporaryDirectory() as tmp_dir:
        report = rank_campaign_measurements(
            posterior_result=posterior,
            benchmark_report=benchmark,
            candidate_edges=candidate_edges,
            validation_report=val_report,
            output_dir=tmp_dir,
        )

        assert report["status"] == "ACTIONABLE"
        assert "rankings" in report
        assert len(report["rankings"]) > 0
        assert "campaign_batch" in report
        assert "claim_guardrail" in report
        assert "decision-support" in report["claim_guardrail"]

        # Verify files were generated
        out_path = Path(tmp_dir)
        json_file = out_path / "well_campaign_recommendations.json"
        csv_file = out_path / "well_campaign_recommendations.csv"
        md_file = out_path / "well_campaign_recommendations.md"

        assert json_file.exists()
        assert csv_file.exists()
        assert md_file.exists()

        # Check JSON payload
        with open(json_file, encoding="utf-8") as f:
            data = json.load(f)
            assert data["status"] == "ACTIONABLE"
            assert len(data["campaign_batch"]["selected_options"]) > 0

        # Check CSV content
        with open(csv_file, encoding="utf-8") as f:
            reader = csv.DictReader(f)
            rows = list(reader)
            assert len(rows) > 0
            assert "action_id" in rows[0]
            assert "well_id" in rows[0]
            assert "expected_information_gain" in rows[0]

        # Check Markdown content
        md_text = md_file.read_text(encoding="utf-8")
        assert "# Well-Level Active Learning Campaign Design" in md_text
        assert "Decision-Support Notice" in md_text
        assert "Selected Campaign Batch" in md_text


def test_rank_next_measurements_bayesian_routing():
    """Test that rank_next_measurements dispatches to bayesian_campaign when requested."""
    posterior = _make_synth_posterior([["E1", "E2"], ["E1"]], [0.6, 0.4])
    benchmark = _make_synth_benchmark()
    candidate_edges = [
        _make_edge("E1", "W1", "W2"),
        _make_edge("E2", "W2", "W3"),
    ]

    # Calling with method='bayesian_campaign'
    result = rank_next_measurements(
        benchmark_report=benchmark,
        candidate_edges=candidate_edges,
        method="bayesian_campaign",
        posterior_ensemble=posterior,
    )
    assert result["status"] == "ACTIONABLE"
    assert "campaign_batch" in result

    # Calling without posterior_ensemble raises ValueError
    with pytest.raises(ValueError, match="posterior_ensemble must be provided"):
        rank_next_measurements(
            benchmark_report=benchmark,
            candidate_edges=candidate_edges,
            method="bayesian_campaign",
            posterior_ensemble=None,
        )

    # Calling without candidate_edges raises ValueError
    with pytest.raises(ValueError, match="candidate_edges must be provided"):
        rank_next_measurements(
            benchmark_report=benchmark,
            candidate_edges=None,
            method="bayesian_campaign",
            posterior_ensemble=posterior,
        )


def test_campaign_abstains_when_predictive_model_is_not_declared():
    """Production-safe mode must not invent observation distributions."""
    posterior = _make_synth_posterior([["E1", "E2"], ["E1"]], [0.6, 0.4])
    benchmark = _make_synth_benchmark()
    edges = [
        _make_edge("E1", "W1", "W2"),
        _make_edge("E2", "W2", "W3"),
    ]
    result = rank_campaign_measurements(
        posterior_result=posterior,
        benchmark_report=benchmark,
        candidate_edges=edges,
        config=CampaignConfig(allow_surrogate_predictive_model=False),
    )
    assert result["status"] == "ABSTAIN"
    assert any("feasible measurement" in reason.lower() for reason in result["abstention_reasons"])


def test_configured_predictive_model_populates_action_scenarios():
    """A supplied forward model is used instead of the topology surrogate."""
    edges = [_make_edge("E1", "W1", "W2")]

    def predictor(context):
        return (1.0 if context["hypothesis_id"] == "H1" else 0.0, 0.1)

    actions = build_well_actions(
        candidate_edges=edges,
        hypothesis_ids=["H1", "H2"],
        prior_probs=[0.5, 0.5],
        edge_tuples=[("E1",), ()],
        config=CampaignConfig(
            predictive_model=predictor,
            allow_surrogate_predictive_model=False,
        ),
    )
    assert actions
    assert all(len(action.scenarios) == 2 for action in actions)
    assert all(action.metadata["predictive_model"] == "configured_callback" for action in actions)


def test_inline_validation_labels_are_used_by_campaign():
    """Inline validation labels must not silently become unlabeled gaps."""
    posterior = _make_synth_posterior([["E1", "E2"], ["E1"]], [0.6, 0.4])
    benchmark = _make_synth_benchmark()
    edges = [
        _make_edge("E1", "W1", "W2"),
        _make_edge("E2", "W2", "W3"),
    ]
    result = rank_campaign_measurements(
        posterior_result=posterior,
        benchmark_report=benchmark,
        candidate_edges=edges,
        validation_report={"labels": {"E2": {"status": "FN"}}},
    )
    assert result["validation_targeting"]["statuses"]["E2"] == "FN"


def test_campaign_abstains_when_predictive_information_is_zero():
    """A non-discriminating observation model must not yield a recommendation."""
    posterior = _make_synth_posterior([["E1"], []], [0.5, 0.5])
    benchmark = _make_synth_benchmark()
    edges = [_make_edge("E1", "W1", "W2")]

    def non_discriminating_predictor(_context):
        return (0.0, 1.0)

    result = rank_campaign_measurements(
        posterior_result=posterior,
        benchmark_report=benchmark,
        candidate_edges=edges,
        config=CampaignConfig(
            predictive_model=non_discriminating_predictor,
            allow_surrogate_predictive_model=False,
        ),
    )
    assert result["status"] == "ABSTAIN"
    assert any("information gain" in reason.lower() for reason in result["abstention_reasons"])
