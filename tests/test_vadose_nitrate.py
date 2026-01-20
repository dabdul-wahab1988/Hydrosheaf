from __future__ import annotations

from datetime import datetime, timedelta

import math

import pytest

from hydrosheaf.vadose.nitrate import (
    NitrateLoadingSample,
    predict_no3_breakthrough,
)


def _days(start: datetime, n: int) -> list[datetime]:
    return [start + timedelta(days=i) for i in range(n)]


def test_no3_breakthrough_delta_kernel_shift_and_attenuation() -> None:
    start = datetime(2020, 1, 1)
    timestamps = _days(start, 30)
    recharge = [0.01] * len(timestamps)

    loading = [NitrateLoadingSample(timestamp=start + timedelta(days=0), c_in=100.0)]

    # Delta kernel at lag=10 days on a 1-day grid.
    tau = list(range(0, 21))
    g = [0.0] * len(tau)
    g[10] = 1.0

    k = 0.05
    points, summary = predict_no3_breakthrough(
        edge_id="A->B",
        ttd_tau_days=tau,
        ttd_g=g,
        timestamps=timestamps,
        recharge_m_day=recharge,
        loading=loading,
        k_1_day=k,
        c_in_units="mg/L",
        grid_dt_days=1.0,
    )

    c = {p.timestamp: p.c_wt for p in points}
    assert c[start + timedelta(days=9)] == 0.0
    assert c[start + timedelta(days=10)] == pytest.approx(
        100.0 * math.exp(-k * 10.0), rel=1e-12, abs=1e-12
    )
    assert summary.edge_id == "A->B"
    assert summary.attenuated_fraction > 0.0


def test_no3_breakthrough_constant_input_reaches_constant_output() -> None:
    start = datetime(2020, 1, 1)
    timestamps = _days(start, 20)
    recharge = [0.02] * len(timestamps)
    loading = [
        NitrateLoadingSample(timestamp=start + timedelta(days=i), c_in=50.0)
        for i in range(0, 30)
    ]

    # A simple kernel supported on [0,4] days, normalized by the implementation.
    tau = [0, 1, 2, 3, 4]
    g = [1.0, 1.0, 1.0, 1.0, 1.0]

    points, _ = predict_no3_breakthrough(
        edge_id="U->V",
        ttd_tau_days=tau,
        ttd_g=g,
        timestamps=timestamps,
        recharge_m_day=recharge,
        loading=loading,
        k_1_day=0.0,
        grid_dt_days=1.0,
    )
    # After the kernel support has fully "filled", output equals the constant input.
    for p in points:
        if (p.timestamp - start).days >= 4:
            assert p.c_wt == pytest.approx(50.0, rel=1e-12, abs=1e-12)


def test_no3_breakthrough_attenuation_reduces_total_mass() -> None:
    start = datetime(2020, 1, 1)
    timestamps = _days(start, 60)
    recharge = [0.01] * len(timestamps)
    loading = [
        NitrateLoadingSample(timestamp=start + timedelta(days=i), c_in=10.0)
        for i in range(0, 60)
    ]

    tau = list(range(0, 31))
    g = [1.0] * len(tau)

    _, s0 = predict_no3_breakthrough(
        edge_id="U->V",
        ttd_tau_days=tau,
        ttd_g=g,
        timestamps=timestamps,
        recharge_m_day=recharge,
        loading=loading,
        k_1_day=0.0,
        grid_dt_days=1.0,
    )
    _, s1 = predict_no3_breakthrough(
        edge_id="U->V",
        ttd_tau_days=tau,
        ttd_g=g,
        timestamps=timestamps,
        recharge_m_day=recharge,
        loading=loading,
        k_1_day=0.02,
        grid_dt_days=1.0,
    )
    assert s1.total_mass_delivered < s0.total_mass_delivered
    assert s1.attenuated_fraction > 0.0
