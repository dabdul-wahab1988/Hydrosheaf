from __future__ import annotations

import csv
from dataclasses import dataclass
from datetime import datetime, timedelta
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import numpy as np


@dataclass(frozen=True)
class NitrateLoadingSample:
    timestamp: datetime
    c_in: float  # concentration in infiltrating water (arbitrary units; preserved through convolution)


@dataclass(frozen=True)
class BreakthroughPoint:
    timestamp: datetime
    c_wt: float
    recharge_m_day: float
    flux_mass_per_m2_day: float


@dataclass(frozen=True)
class BreakthroughSummary:
    edge_id: str
    k_1_day: float
    c_in_units: str
    total_mass_delivered: float
    total_mass_no_attenuation: float
    attenuated_fraction: float
    peak_c_wt: float
    peak_time: Optional[str]


def load_no3_loading_csv(
    path: str,
    *,
    time_column: str = "timestamp",
    time_format: str = "%Y-%m-%d",
    concentration_column: str = "NO3_mg_L",
) -> List[NitrateLoadingSample]:
    """Load nitrate input concentration time series at the surface/root-zone.

    The default column name is NO3_mg_L, but any numeric column can be used by overriding
    `concentration_column`.
    """
    samples: List[NitrateLoadingSample] = []
    with open(path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            ts_raw = (row.get(time_column) or "").strip()
            if not ts_raw:
                continue
            try:
                ts = datetime.strptime(ts_raw, time_format)
            except ValueError:
                # common alternatives
                for fmt in ["%Y-%m-%d %H:%M:%S", "%Y/%m/%d", "%m/%d/%Y"]:
                    try:
                        ts = datetime.strptime(ts_raw, fmt)
                        break
                    except ValueError:
                        continue
                else:
                    continue
            raw = row.get(concentration_column)
            if raw in (None, ""):
                continue
            try:
                c = float(raw)
            except ValueError:
                continue
            samples.append(NitrateLoadingSample(timestamp=ts, c_in=float(c)))
    samples.sort(key=lambda s: s.timestamp)
    return samples


def _build_uniform_grid(
    start: datetime,
    end: datetime,
    dt_days: float,
) -> List[datetime]:
    if dt_days <= 0:
        raise ValueError("dt_days must be > 0")
    n = int(np.floor((end - start).total_seconds() / 86400.0 / dt_days)) + 1
    grid: List[datetime] = []
    for i in range(n):
        grid.append(start + timedelta(days=float(i) * float(dt_days)))
    return grid


def _resample_to_grid(
    series: Sequence[Tuple[datetime, float]],
    grid: Sequence[datetime],
    *,
    fill: float = 0.0,
) -> np.ndarray:
    if not grid:
        return np.array([], dtype=float)
    if not series:
        return np.full(len(grid), float(fill), dtype=float)
    times = np.array([t.timestamp() / 86400.0 for t, _ in series], dtype=float)
    vals = np.array([float(v) for _, v in series], dtype=float)
    grid_t = np.array([t.timestamp() / 86400.0 for t in grid], dtype=float)
    # Clamp to range, fill outside with fill
    out = np.interp(grid_t, times, vals, left=float(fill), right=float(fill))
    return out


def predict_no3_breakthrough(
    *,
    edge_id: str,
    ttd_tau_days: Sequence[float],
    ttd_g: Sequence[float],
    timestamps: Sequence[datetime],
    recharge_m_day: Sequence[float],
    loading: Sequence[NitrateLoadingSample],
    k_1_day: float = 0.0,
    c_in_units: str = "mg/L",
    grid_dt_days: Optional[float] = None,
) -> Tuple[List[BreakthroughPoint], BreakthroughSummary]:
    """Predict a nitrate breakthrough curve at the water table using a TTD kernel.

    Model
    -----
      C_wt(t) = ∫ C_in(t-τ) g(τ) exp(-k τ) dτ

    Flux (mass per area per time) is computed as:
      F(t) = R(t) * C_wt(t) * 1000  (assuming C in mg/L and R in m/day)

    Notes
    -----
    - Units for concentration are preserved through convolution (the model is linear).
    - The flux conversion assumes 1 m^3 = 1000 L; if using different units, interpret flux accordingly.
    """
    if len(ttd_tau_days) != len(ttd_g) or not ttd_tau_days:
        raise ValueError("TTD kernel is empty or has mismatched lengths.")
    if not timestamps or len(timestamps) != len(recharge_m_day):
        raise ValueError("timestamps and recharge_m_day must be the same length and non-empty.")
    if not loading:
        raise ValueError("loading is empty.")

    tau = np.asarray(ttd_tau_days, dtype=float)
    g = np.asarray(ttd_g, dtype=float)
    # Ensure nonnegative and normalize to a discrete PDF.
    g = np.where(np.isfinite(g) & (g > 0), g, 0.0)
    if not np.any(g > 0):
        raise ValueError("TTD kernel has no positive mass.")

    # Determine dt used for discrete convolution integral approximation.
    if grid_dt_days is None:
        # Prefer the kernel's grid spacing if available.
        if len(tau) >= 2:
            grid_dt_days = float(max(1e-9, tau[1] - tau[0]))
        else:
            grid_dt_days = 1.0
    dt = float(grid_dt_days)

    area = float(np.sum(g) * dt)
    if area <= 0:
        raise ValueError("TTD kernel normalization failed.")
    g = g / area

    # Build a uniform time grid over the overlap of forcing and loading.
    start = min(min(timestamps), loading[0].timestamp)
    end = max(max(timestamps), loading[-1].timestamp)
    grid = _build_uniform_grid(start, end, dt)

    # Resample loading and recharge onto the same grid.
    cin_series = _resample_to_grid([(s.timestamp, s.c_in) for s in loading], grid, fill=0.0)
    r_series = _resample_to_grid(list(zip(timestamps, recharge_m_day)), grid, fill=0.0)

    # Discrete convolution: C_wt[t] = sum_j cin[t-j] * g[j] * exp(-k*tau[j]) * dt
    att = np.exp(-float(k_1_day) * tau)
    kernel = g * att
    # Convert integral approximation: sum(kernel*cin_shift)*dt
    kernel_dt = kernel * dt

    cwt = np.zeros_like(cin_series)
    for j in range(len(kernel_dt)):
        lag = j
        if lag == 0:
            cwt += cin_series * kernel_dt[j]
        else:
            cwt[lag:] += cin_series[:-lag] * kernel_dt[j]

    # Assemble outputs at the original forcing timestamps (nearest grid points).
    grid_t_days = np.array([t.timestamp() / 86400.0 for t in grid], dtype=float)
    out_points: List[BreakthroughPoint] = []
    for ts, r in zip(timestamps, recharge_m_day):
        t0 = float(ts.timestamp() / 86400.0)
        idx = int(np.argmin(np.abs(grid_t_days - t0)))
        c_val = float(cwt[idx])
        r_val = float(r)
        flux = float(r_val) * float(c_val) * 1000.0
        out_points.append(
            BreakthroughPoint(
                timestamp=ts,
                c_wt=c_val,
                recharge_m_day=r_val,
                flux_mass_per_m2_day=flux,
            )
        )

    # Summary + attenuation fraction relative to k=0 baseline.
    flux_arr = np.array([p.flux_mass_per_m2_day for p in out_points], dtype=float)
    total_mass = float(np.sum(np.maximum(0.0, flux_arr)) * dt)

    # Baseline k=0
    kernel0 = g
    kernel0_dt = kernel0 * dt
    cwt0 = np.zeros_like(cin_series)
    for j in range(len(kernel0_dt)):
        lag = j
        if lag == 0:
            cwt0 += cin_series * kernel0_dt[j]
        else:
            cwt0[lag:] += cin_series[:-lag] * kernel0_dt[j]
    flux0_points = []
    for ts, r in zip(timestamps, recharge_m_day):
        t0 = float(ts.timestamp() / 86400.0)
        idx = int(np.argmin(np.abs(grid_t_days - t0)))
        flux0_points.append(float(r) * float(cwt0[idx]) * 1000.0)
    total_mass0 = float(np.sum(np.maximum(0.0, np.asarray(flux0_points))) * dt)
    atten_frac = 0.0
    if total_mass0 > 0:
        atten_frac = float(np.clip(1.0 - (total_mass / total_mass0), 0.0, 1.0))

    peak_idx = int(np.argmax([p.c_wt for p in out_points])) if out_points else 0
    peak = out_points[peak_idx] if out_points else None

    summary = BreakthroughSummary(
        edge_id=str(edge_id),
        k_1_day=float(k_1_day),
        c_in_units=str(c_in_units),
        total_mass_delivered=float(total_mass),
        total_mass_no_attenuation=float(total_mass0),
        attenuated_fraction=float(atten_frac),
        peak_c_wt=float(peak.c_wt) if peak else float("nan"),
        peak_time=peak.timestamp.strftime("%Y-%m-%d") if peak else None,
    )
    return out_points, summary

