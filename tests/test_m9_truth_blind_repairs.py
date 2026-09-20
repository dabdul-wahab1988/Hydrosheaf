"""Regression tests for the repaired M9 extension contracts."""

from __future__ import annotations

import json

import numpy as np
import pytest

from hydrosheaf.benchmarks.ttd_graph_dynamic_kernel import (
    DynamicKernelBenchmarkConfig,
    generate_dynamic_kernel_case,
)
from hydrosheaf.benchmarks.ttd_graph_multitracer import (
    MultiTracerBenchmarkConfig,
    generate_multitracer_case,
)
from hydrosheaf.nuclear.dynamic_kernel_inversion import (
    DynamicTTDInversionConfig,
    solve_dynamic_node_inversion,
)
from hydrosheaf.nuclear.multi_tracer_graph_inversion import (
    MultiTracerGraphConfig,
    solve_joint_multitracer_node_inversion,
)


def test_dynamic_estimator_rejects_scenario_oracle_and_serializes_dynamic_arrays():
    times = np.arange(80)
    parent = np.sin(times / 8.0)
    local = np.cos(times / 11.0)
    target = 0.25 * local + 0.75 * np.roll(parent, 2)
    target[:2] = target[2]
    cfg = DynamicTTDInversionConfig(
        lag_grid_days=(0.0, 1.0, 2.0, 3.0),
        step_days=1.0,
        kernel_mode="phase",
        season_period_steps=20,
        n_phase_bins=2,
        min_pairs_per_phase=4,
        min_identified_phases=2,
        min_r2=-1.0,
    )
    with pytest.raises(ValueError, match="scenario_hint"):
        solve_dynamic_node_inversion(
            "A",
            times,
            target,
            {"R": parent},
            local_input=local,
            config=cfg,
            scenario_hint="stationary",
        )

    result = solve_dynamic_node_inversion(
        "A", times, target, {"R": parent}, local_input=local, config=cfg
    )["R->A"]
    assert result.status == "RECOVERED"
    assert result.estimated_kernel is not None
    assert result.estimated_mixing_fractions is not None
    assert result.local_recharge_fractions is not None
    assert result.estimated_kernel.shape == (len(times), 4)
    np.testing.assert_allclose(result.estimated_kernel.sum(axis=1), 1.0)
    assert np.ptp(result.estimated_mixing_fractions) >= 0.0
    payload = result.to_dict()
    assert payload["estimated_kernel"] is not None
    assert payload["estimated_mixing_fractions"] is not None
    json.dumps(payload)


def test_truth_and_observation_commitments_cover_visible_contracts():
    truth, observations = generate_dynamic_kernel_case(
        DynamicKernelBenchmarkConfig(seed=11, n_steps=80)
    )
    assert truth.sealed_commitment == observations.sealed_commitment
    assert len(observations.observation_commitment) == 64
    assert truth.local_inputs
    assert observations.local_inputs

    truth_mt, observations_mt = generate_multitracer_case(
        MultiTracerBenchmarkConfig(seed=11, n_steps=80)
    )
    assert truth_mt.sealed_commitment == observations_mt.sealed_commitment
    assert len(observations_mt.observation_commitment) == 64
    assert observations_mt.local_tracer_inputs["B"]


def test_joint_multitracer_fit_reports_refit_loto_and_holdout_metrics():
    truth, observations = generate_multitracer_case(
        MultiTracerBenchmarkConfig(seed=17, panel="isotopes_plus_SF6", n_steps=104)
    )
    train = observations.training_end_step
    target_times = observations.time_steps[:train]
    holdout_times = observations.time_steps[train:]
    calibration = {
        tracer: observations.node_tracer_observations["B"][tracer][:train]
        for tracer in observations.available_tracers
    }
    holdout = {
        tracer: observations.node_tracer_observations["B"][tracer][train:]
        for tracer in observations.available_tracers
    }
    cfg = MultiTracerGraphConfig(
        inversion_config=DynamicTTDInversionConfig(
            lag_grid_days=tuple(float(v) for v in observations.metadata["lag_grid_days"]),
            step_days=observations.metadata["step_days"],
            season_period_steps=26,
            n_phase_bins=4,
            min_r2=-1.0,
        ),
        enforce_holdout_gate=True,
    )
    result = solve_joint_multitracer_node_inversion(
        "B",
        target_times,
        calibration,
        {"A": observations.node_tracer_observations["A"]},
        observations.local_tracer_inputs["B"],
        config=cfg,
        holdout_times=holdout_times,
        holdout_tracer_observations=holdout,
    )
    assert result.diagnostics.get("joint_fit") is True
    assert result.diagnostics.get("loto_refit_status")
    assert set(result.loto_sensitivity) == set(observations.available_tracers)
    assert result.held_out_per_tracer_rmse
    assert result.held_out_per_tracer_normalized_rmse
