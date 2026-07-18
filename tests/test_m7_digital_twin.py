from __future__ import annotations

import sys
from dataclasses import replace
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SCRIPT_DIR = ROOT / "M7" / "m7_integration_benchmark" / "scripts"
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from m7_digital_twin import (  # noqa: E402
    TwinConfig,
    build_models,
    calibrate_forecast_spread,
    evaluate_filter_forecasts,
    run_filter,
    simulate_experiment,
)


def quick_config() -> TwinConfig:
    return TwinConfig(
        n_steps=54,
        forecast_start=36,
        ensemble_size=32,
        n_replicates=2,
        horizons=(1, 3, 6),
        seed=1947,
    )


def test_truth_and_operational_model_are_not_identical() -> None:
    config = quick_config()
    truth, nominal, wrong = build_models(config)
    assert not np.allclose(truth.upstream_weights, nominal.upstream_weights)
    assert not np.allclose(truth.retention, nominal.retention)
    assert truth.nonlinear_recharge > 0.0
    assert truth.hidden_pulse
    assert not nominal.hidden_pulse
    assert not np.allclose(nominal.upstream_weights, wrong.upstream_weights)


def test_future_observations_cannot_change_an_issued_posterior() -> None:
    config = quick_config()
    data = simulate_experiment(config, seed=811)
    _, nominal, _ = build_models(config)
    origin = 39
    altered_batches = list(data.observations)
    for time in range(origin + 1, config.n_steps):
        batch = altered_batches[time]
        altered_batches[time] = replace(batch, values=batch.values + 10_000.0)
    altered = replace(data, observations=tuple(altered_batches))

    original_run = run_filter(
        data,
        nominal,
        config,
        method="updated_twin",
        seed=99,
        stop_assimilation_at=origin,
    )
    altered_run = run_filter(
        altered,
        nominal,
        config,
        method="updated_twin",
        seed=99,
        stop_assimilation_at=origin,
    )
    np.testing.assert_allclose(
        original_run.posterior_ensembles[origin],
        altered_run.posterior_ensembles[origin],
        rtol=0.0,
        atol=0.0,
    )

    original_calibration = calibrate_forecast_spread(
        data,
        original_run,
        nominal,
        config,
        seed=404,
    )
    altered_calibration = calibrate_forecast_spread(
        altered,
        altered_run,
        nominal,
        config,
        seed=404,
    )
    for key, factor in original_calibration.items():
        issue_time = key[0]
        if issue_time <= origin:
            assert factor == altered_calibration[key]


def test_forecast_outputs_are_finite_and_include_unmonitored_nodes() -> None:
    config = quick_config()
    data = simulate_experiment(config, seed=812)
    _, nominal, _ = build_models(config)
    run = run_filter(
        data,
        nominal,
        config,
        method="updated_twin",
        seed=101,
    )
    metrics = evaluate_filter_forecasts(
        data,
        run,
        nominal,
        config,
        replicate=0,
        seed=202,
    )
    assert {"all", "monitored", "unmonitored"} == set(metrics["domain"])
    assert {1, 3, 6} == set(metrics["horizon"])
    assert np.isfinite(metrics[["rmse", "mae", "bias", "coverage90", "crps"]]).all().all()
    assert ((metrics["coverage90"] >= 0.0) & (metrics["coverage90"] <= 1.0)).all()


def test_sequential_updates_reduce_state_error_after_spinup() -> None:
    config = quick_config()
    data = simulate_experiment(config, seed=813)
    _, nominal, _ = build_models(config)
    updated = run_filter(
        data,
        nominal,
        config,
        method="updated_twin",
        seed=303,
    )
    open_loop = run_filter(
        data,
        nominal,
        config,
        method="open_loop",
        seed=303,
        observation_mode="none",
    )
    evaluation = updated.diagnostics.query("time >= 12")
    open_evaluation = open_loop.diagnostics.query("time >= 12")
    assert evaluation["posterior_rmse"].mean() < open_evaluation["posterior_rmse"].mean()
