from __future__ import annotations

from dataclasses import dataclass
from typing import List, Sequence, Tuple
import math

import numpy as np



@dataclass(frozen=True)
class TravelTimeSummary:
    mean_days: float
    std_days: float
    p10_days: float
    p50_days: float
    p90_days: float


def _weighted_quantile(values: np.ndarray, weights: np.ndarray, q: float) -> float:
    if values.size == 0:
        return float("nan")
    order = np.argsort(values)
    v = values[order]
    w = weights[order]
    cum = np.cumsum(w)
    if cum[-1] <= 0:
        return float("nan")
    target = q * cum[-1]
    idx = int(np.searchsorted(cum, target, side="left"))
    idx = max(0, min(idx, v.size - 1))
    return float(v[idx])


def advective_travel_time_days(
    *,
    z_m: Sequence[float],
    theta: Sequence[float],
    q_faces_m_day: Sequence[float],
    from_depth_m: float,
    theta_min: float = 1e-6,
    q_min_m_day: float = 1e-9,
) -> float:
    """Compute advective travel time from depth to bottom based on pore velocity.

    Uses dt ≈ dz * theta / q on each cell, with q taken at the lower face of the cell.
    """
    z = np.asarray(z_m, dtype=float)
    th = np.asarray(theta, dtype=float)
    qf = np.asarray(q_faces_m_day, dtype=float)
    if z.size < 2:
        return float("nan")
    dz = float(z[1] - z[0])
    if dz <= 0:
        return float("nan")

    start = float(from_depth_m)
    if start < 0:
        start = 0.0
    if start >= float(z[-1]):
        return 0.0

    # Determine starting cell index.
    i0 = int(np.floor(start / dz))
    i0 = max(0, min(i0, z.size - 2))

    tau = 0.0
    # First partial cell
    z0 = float(z[i0])
    frac = (z0 + dz - start) / dz
    q_cell = float(qf[i0 + 1])
    q_cell = max(q_min_m_day, q_cell)
    theta_cell = max(theta_min, float(0.5 * (th[i0] + th[i0 + 1])))
    tau += frac * dz * theta_cell / q_cell

    # Full cells below
    for i in range(i0 + 1, z.size - 1):
        q_cell = float(qf[i + 1])
        q_cell = max(q_min_m_day, q_cell)
        theta_cell = max(theta_min, float(0.5 * (th[i] + th[i + 1])))
        tau += dz * theta_cell / q_cell

    return float(tau)


def gamma_kernel_grid(
    *,
    mean_days: float,
    cv: float,
    grid_dt_days: float,
    max_lag_days: float,
) -> Tuple[List[float], List[float]]:
    """Construct a discrete gamma PDF on a regular τ grid with given mean and CV."""
    mean = float(mean_days)
    if not (mean > 0):
        return [], []
    cv = float(max(1e-6, cv))
    # Gamma: mean = k*theta, var = k*theta^2, so CV^2 = 1/k -> k = 1/CV^2
    k = 1.0 / (cv * cv)
    theta = mean / k
    # grid
    n = int(np.floor(float(max_lag_days) / float(grid_dt_days))) + 1
    tau = np.linspace(0.0, float(max_lag_days), n)
    # gamma pdf for tau>0
    pdf = np.zeros_like(tau)
    positive = tau > 0
    if np.any(positive):
        t = tau[positive]
        # pdf = t^(k-1) * exp(-t/theta) / (Gamma(k) * theta^k)
        # compute in log-space for stability
        log_pdf = (
            (k - 1.0) * np.log(t)
            - (t / theta)
            - (np.log(math.gamma(k)) + k * np.log(theta))
        )

        pdf[positive] = np.exp(log_pdf)
    # normalize discrete
    area = float(np.sum(pdf) * float(grid_dt_days))
    if area > 0:
        pdf = pdf / area
    return [float(x) for x in tau.tolist()], [float(x) for x in pdf.tolist()]


def mixture_ttd_from_series(
    tau_days: Sequence[float],
    weights: Sequence[float],
    *,
    grid_dt_days: float,
    max_lag_days: float,
    cv: float,
    preferential_flow_fraction: float = 0.0,
    preferential_velocity_factor: float = 10.0,
    preferential_cv: float = 1.0,
) -> Tuple[List[float], List[float], TravelTimeSummary]:
    taus = np.asarray(tau_days, dtype=float)
    w = np.asarray(weights, dtype=float)
    w = np.where(np.isfinite(w) & (w > 0), w, 0.0)
    taus = np.where(np.isfinite(taus) & (taus > 0), taus, np.nan)
    mask = np.isfinite(taus) & (w > 0)
    if not np.any(mask):
        summary = TravelTimeSummary(
            mean_days=float("nan"),
            std_days=float("nan"),
            p10_days=float("nan"),
            p50_days=float("nan"),
            p90_days=float("nan"),
        )
        return [], [], summary
    taus = taus[mask]
    w = w[mask]
    w = w / np.sum(w)

    # Base statistics (matrix/advective only)
    mean = float(np.sum(w * taus))
    var = float(np.sum(w * (taus - mean) ** 2))
    std = float(np.sqrt(max(0.0, var)))
    p10 = _weighted_quantile(taus, w, 0.1)
    p50 = _weighted_quantile(taus, w, 0.5)
    p90 = _weighted_quantile(taus, w, 0.9)
    
    # If preferential flow is active, we should probably update the summary stats
    # to reflect the effective dual-domain behavior.
    # For now, we report the MATRIX (advective) statistics as "baseline",
    # but the kernel below will include the fast path.
    # (Future TODO: Calculate effective moments of the mixture)
    
    summary = TravelTimeSummary(
        mean_days=mean, std_days=std, p10_days=p10, p50_days=p50, p90_days=p90
    )

    # Build mixture kernel as weighted sum of per-time gamma kernels.
    # Now supports Dual-Domain (Matrix + Preferential)
    n = int(np.floor(float(max_lag_days) / float(grid_dt_days))) + 1
    grid = np.linspace(0.0, float(max_lag_days), n)
    pdf_mix = np.zeros_like(grid)
    
    w_matrix = 1.0 - max(0.0, min(1.0, preferential_flow_fraction))
    w_fast = 1.0 - w_matrix
    fast_factor = max(1.0, preferential_velocity_factor)
    
    for tau_i, wi in zip(taus, w):
        # Matrix component
        if w_matrix > 0:
            g_tau, g_pdf = gamma_kernel_grid(
                mean_days=float(tau_i),
                cv=cv,
                grid_dt_days=grid_dt_days,
                max_lag_days=max_lag_days,
            )
            if g_tau:
                pdf_mix += wi * w_matrix * np.asarray(g_pdf, dtype=float)
        
        # Preferential component
        if w_fast > 0:
            tau_fast = float(tau_i) / fast_factor
            g_tau_fast, g_pdf_fast = gamma_kernel_grid(
                mean_days=tau_fast,
                cv=preferential_cv,
                grid_dt_days=grid_dt_days,
                max_lag_days=max_lag_days,
            )
            if g_tau_fast:
                pdf_mix += wi * w_fast * np.asarray(g_pdf_fast, dtype=float)

    area = float(np.sum(pdf_mix) * float(grid_dt_days))
    if area > 0:
        pdf_mix /= area
    return (
        [float(x) for x in grid.tolist()],
        [float(x) for x in pdf_mix.tolist()],
        summary,
    )

