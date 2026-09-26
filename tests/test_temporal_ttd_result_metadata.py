"""Focused contracts for temporal residence-time result semantics."""

from datetime import datetime, timedelta

import numpy as np

from hydrosheaf.temporal import TemporalNode, TimeSeriesSample
from hydrosheaf.temporal.residence_time import estimate_residence_time_with_details


def _nodes(
    upstream: np.ndarray,
    downstream: np.ndarray,
    ion_order: list[str],
) -> tuple[TemporalNode, TemporalNode]:
    """Build two aligned temporal nodes from one- or two-dimensional arrays."""
    upstream = np.asarray(upstream, dtype=float)
    downstream = np.asarray(downstream, dtype=float)
    if upstream.ndim == 1:
        upstream = upstream[:, None]
    if downstream.ndim == 1:
        downstream = downstream[:, None]

    t0 = datetime(2020, 1, 1)
    upstream_samples = [
        TimeSeriesSample(
            sample_id=f"u{i}",
            node_id="U",
            timestamp=t0 + timedelta(days=i),
            concentrations=row.tolist(),
        )
        for i, row in enumerate(upstream)
    ]
    downstream_samples = [
        TimeSeriesSample(
            sample_id=f"v{i}",
            node_id="V",
            timestamp=t0 + timedelta(days=i),
            concentrations=row.tolist(),
        )
        for i, row in enumerate(downstream)
    ]
    assert upstream.shape[1] == len(ion_order)
    assert downstream.shape[1] == len(ion_order)
    return (
        TemporalNode(node_id="U", samples=upstream_samples),
        TemporalNode(node_id="V", samples=downstream_samples),
    )


def _ttd_signals(n: int = 180) -> tuple[np.ndarray, np.ndarray]:
    t = np.arange(n, dtype=float)
    upstream = np.sin(t / 12.0) + 0.3 * np.cos(t / 6.0)
    kernel = np.array([0.05, 0.15, 0.4, 0.25, 0.15], dtype=float)
    downstream = 0.2 + np.convolve(upstream, kernel, mode="full")[:n]
    return upstream, downstream


def test_ttd_metadata_identifies_finite_grid_per_tracer_fit() -> None:
    upstream, downstream = _ttd_signals()
    node_u, node_v = _nodes(upstream, downstream, ["Cl"])

    tau, uncertainty, used, details, flags = estimate_residence_time_with_details(
        node_u,
        node_v,
        method="ttd",
        tracer_ion="Cl",
        ion_order=["Cl"],
        hydraulic_params={
            "grid_dt_days": 1.0,
            "max_lag_days": 30.0,
            "smoothness_lambda": 0.01,
            "ttd_min_r2": 0.1,
            "attenuation_k_max": 0.02,
            "attenuation_k_steps": 4,
        },
    )

    assert tau > 0.0
    assert uncertainty >= 0.0
    assert "ttd" in used
    assert "ttd_failed_all_tracers" not in flags
    assert "candidates" in details

    metadata = details["result_metadata"]
    assert metadata["inference_family"] == "finite_grid_nonnegative_ttd"
    assert metadata["aggregation"] == "per_tracer_fit"
    assert metadata["tracers_requested"] == ["Cl"]
    assert metadata["tracers_accepted"] == ["Cl"]
    lag_grid = metadata["lag_grid"]
    assert lag_grid["dt_days"] == 1.0
    assert lag_grid["max_lag_days"] == 30.0
    assert lag_grid["smoothness_lambda"] == 0.01
    np.testing.assert_allclose(
        lag_grid["attenuation_k_grid"], np.linspace(0.0, 0.02, 4)
    )


def test_ttd_metadata_identifies_multi_tracer_consensus() -> None:
    upstream, downstream = _ttd_signals()
    node_u, node_v = _nodes(
        np.column_stack([upstream, upstream]),
        np.column_stack([downstream, downstream]),
        ["Cl", "Na"],
    )

    _, _, used, details, flags = estimate_residence_time_with_details(
        node_u,
        node_v,
        method="ttd",
        tracer_ion="Cl,Na",
        ion_order=["Cl", "Na"],
        hydraulic_params={
            "grid_dt_days": 1.0,
            "max_lag_days": 30.0,
            "ttd_min_r2": 0.1,
            "attenuation_k_steps": 2,
        },
    )

    assert "ttd" in used
    assert "ttd_failed_all_tracers" not in flags
    metadata = details["result_metadata"]
    assert metadata["inference_family"] == "finite_grid_nonnegative_ttd"
    assert metadata["aggregation"] == "multi_tracer_consensus"
    assert metadata["tracers_requested"] == ["Cl", "Na"]
    assert set(metadata["tracers_accepted"]) == {"Cl", "Na"}
    assert "consensus" in details


def test_cross_correlation_metadata_preserves_consensus_details() -> None:
    t = np.arange(140, dtype=float)
    upstream = np.sin(t / 11.0) + 0.25 * np.cos(t / 5.0)
    lag_days = 5.0
    downstream = np.sin((t - lag_days) / 11.0) + 0.25 * np.cos(
        (t - lag_days) / 5.0
    )
    node_u, node_v = _nodes(
        np.column_stack([upstream, upstream]),
        np.column_stack([downstream, downstream]),
        ["Cl", "Na"],
    )

    _, _, used, details, flags = estimate_residence_time_with_details(
        node_u,
        node_v,
        method="cross_correlation",
        tracer_ion="Cl,Na",
        ion_order=["Cl", "Na"],
        hydraulic_params={
            "min_peak_corr": 0.0,
            "max_relative_uncertainty": 10.0,
            "max_uncertainty_days": 365.0,
        },
    )

    assert "cross_correlation" in used
    assert "tau_failed_all_tracers" not in flags
    metadata = details["result_metadata"]
    assert metadata["inference_family"] == "cross_correlation_lag"
    assert metadata["aggregation"] == "multi_tracer_consensus"
    assert metadata["tracers_requested"] == ["Cl", "Na"]
    assert set(metadata["tracers_accepted"]) == {"Cl", "Na"}
    assert "candidates" in details
    assert "consensus" in details
