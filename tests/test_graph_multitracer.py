"""Tests for integrated multi-tracer graph inversion and tracer dropout."""

from __future__ import annotations

import numpy as np
import pytest

from hydrosheaf.benchmarks.ttd_graph_multitracer import (
    MultiTracerBenchmarkConfig,
    evaluate_multitracer_recovery,
    generate_multitracer_case,
)
from hydrosheaf.nuclear.dynamic_kernel_inversion import (
    DynamicTTDInversionConfig,
    solve_dynamic_node_inversion,
)
from hydrosheaf.nuclear.graph_ttd_inversion import (
    GraphTTDInversionConfig,
    solve_node_ttd_inversion,
)
from hydrosheaf.nuclear.multi_tracer_graph_inversion import (
    MultiTracerGraphConfig,
    solve_joint_multitracer_node_inversion,
)


def test_dynamic_node_inversion_stationary_compatibility():
    """Verify that dynamic inversion in stationary mode matches the stationary solver."""
    np.random.seed(42)
    T = 80
    times = np.arange(T)
    lags = (0.0, 7.0, 14.0, 28.0)

    # Input signals
    loc_input = np.sin(2 * np.pi * times / 30.0)
    p_input = np.cos(2 * np.pi * times / 40.0)

    # True model: 40% local, 60% parent with lag 7
    target = np.zeros(T)
    for t in range(T):
        p_val = p_input[t - 7] if t >= 7 else p_input[0]
        target[t] = 0.40 * loc_input[t] + 0.60 * p_val

    obs_idx = np.arange(10, T, 2)
    obs_vals = target[obs_idx]

    # 1. Existing stationary solver
    stat_cfg = GraphTTDInversionConfig(
        lag_grid_days=tuple(int(x) for x in lags),
        young_water_cutoff_days=14,
        regularization_smoothness=0.01,
        regularization_ridge=1e-4,
        regularization_edge_sparsity=0.01,
    )
    stat_res = solve_node_ttd_inversion(
        node_id="test",
        target_times=obs_idx,
        target_values=obs_vals,
        local_input_series=loc_input,
        parent_series_map={"parent_a": p_input},
        config=stat_cfg,
    )

    # 2. New dynamic solver in stationary mode
    dyn_cfg = DynamicTTDInversionConfig(
        lag_grid_days=lags,
        young_water_cutoff_days=14.0,
        kernel_mode="stationary",
        lambda_lag=0.01,
        lambda_edge_sparsity=0.01,
    )
    dyn_res = solve_dynamic_node_inversion(
        target_node="test",
        target_times=obs_idx,
        target_values=obs_vals,
        candidate_parents={"parent_a": p_input},
        local_input=loc_input,
        config=dyn_cfg,
    )

    assert stat_res.status == "ESTIMATED"
    rec = dyn_res["parent_a->test"]
    assert rec.status == "RECOVERED"
    # Young water fractions should be within 0.05
    assert abs(stat_res.young_water_fraction - rec.young_water_fraction) < 0.05


def test_joint_multitracer_inversion_and_loto():
    """Verify joint multi-tracer inversion, LOTO sensitivity, and scoring."""
    cfg = MultiTracerBenchmarkConfig(panel="all", seed=42)
    truth, obs = generate_multitracer_case(cfg)

    # Target node B
    cal_times = obs.time_steps[: obs.training_end_step]
    ho_times = obs.time_steps[obs.training_end_step :]

    target_cal_obs = {t: obs.node_tracer_observations["B"][t][: obs.training_end_step] for t in obs.available_tracers}
    target_ho_obs = {t: obs.node_tracer_observations["B"][t][obs.training_end_step :] for t in obs.available_tracers}

    parent_series = {"A": obs.node_tracer_observations["A"]}
    local_series = {"R": obs.node_tracer_observations["R"]}

    res = solve_joint_multitracer_node_inversion(
        target_node="B",
        target_times=cal_times,
        tracer_observations=target_cal_obs,
        candidate_parents=parent_series,
        local_tracer_inputs=local_series,
        config=MultiTracerGraphConfig(),
        holdout_times=ho_times,
        holdout_tracer_observations=target_ho_obs,
    )

    assert res.status == "RECOVERED"
    assert len(res.active_tracers) == 5
    assert "d18O" in res.per_tracer_rmse
    assert "3H" in res.per_tracer_rmse
    assert len(res.loto_sensitivity) == 5
    assert res.conflict_detected is False

    score = evaluate_multitracer_recovery(truth, obs, res)
    assert score["decision"] == "JUSTIFIED_RECOVERY"


def test_conflicting_tracers_refusal_gate():
    """Verify that contradictory tracer evidence (e.g. modern SF6 + dead 14C) triggers ABSTAIN."""
    cfg = MultiTracerBenchmarkConfig(panel="conflicting_tracers", seed=42)
    truth, obs = generate_multitracer_case(cfg)

    cal_times = obs.time_steps[: obs.training_end_step]
    target_cal_obs = {t: obs.node_tracer_observations["B"][t][: obs.training_end_step] for t in obs.available_tracers}
    parent_series = {"A": obs.node_tracer_observations["A"]}

    res = solve_joint_multitracer_node_inversion(
        target_node="B",
        target_times=cal_times,
        tracer_observations=target_cal_obs,
        candidate_parents=parent_series,
        config=MultiTracerGraphConfig(),
    )

    assert res.status == "ABSTAIN"
    assert res.conflict_detected is True
    assert "tracer_conflict_detected" in res.reason_codes

    score = evaluate_multitracer_recovery(truth, obs, res)
    assert score["decision"] == "CORRECT_ABSTENTION"
