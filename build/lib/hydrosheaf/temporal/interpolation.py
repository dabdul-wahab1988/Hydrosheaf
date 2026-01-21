"""
Time-series interpolation methods.
"""

from datetime import datetime, timedelta
from typing import Dict, List, Optional

import numpy as np

from . import TemporalNode, TimeSeriesSample


def interpolate_to_common_times(
    nodes: Dict[str, TemporalNode],
    target_times: Optional[List[datetime]] = None,
    method: str = "linear",
    frequency_days: int = 30,
) -> Dict[str, TemporalNode]:
    """
    Interpolate all nodes to common time grid.

    Parameters
    ----------
    nodes : Dict[str, TemporalNode]
        Nodes with irregular time-series
    target_times : Optional[List[datetime]]
        Explicit target times. If None, auto-generate from data range.
    method : str
        "linear", "spline", or "nearest"
    frequency_days : int
        If auto-generating times, spacing in days

    Returns
    -------
    Dict[str, TemporalNode]
        Nodes with interpolated samples at common times

    Mathematical Implementation
    ---------------------------
    For each ion j and node i:

    Linear:
        C_ij(t) = C_ij(t_k) + (C_ij(t_{k+1}) - C_ij(t_k)) * (t - t_k) / (t_{k+1} - t_k)
        where t_k <= t < t_{k+1}

    Spline:
        Use cubic spline interpolation

    Nearest:
        C_ij(t) = C_ij(t_k) where k = argmin_k |t - t_k|
    """
    if not nodes:
        return {}

    # Determine target times if not provided
    if target_times is None:
        # Find global time range
        all_times = []
        for node in nodes.values():
            for sample in node.samples:
                all_times.append(sample.timestamp)

        if not all_times:
            return nodes

        t_min = min(all_times)
        t_max = max(all_times)

        # Generate regular grid
        target_times = []
        current = t_min
        while current <= t_max:
            target_times.append(current)
            current += timedelta(days=frequency_days)

    # Interpolate each node
    interpolated_nodes = {}

    for node_id, node in nodes.items():
        if len(node.samples) < 2:
            # Not enough data to interpolate
            interpolated_nodes[node_id] = node
            continue

        # Extract timestamps and concentrations
        original_times = [s.timestamp for s in node.samples]
        n_ions = len(node.samples[0].concentrations)

        # Convert datetimes to floats (days since epoch)
        def to_days(dt):
            return dt.timestamp() / 86400.0

        original_times_float = [to_days(t) for t in original_times]
        target_times_float = [to_days(t) for t in target_times]

        # Interpolate each ion
        interpolated_samples = []

        for target_time, target_time_float in zip(target_times, target_times_float):
            # Skip if outside range
            if (
                target_time_float < original_times_float[0]
                or target_time_float > original_times_float[-1]
            ):
                continue

            concentrations = []

            for j in range(n_ions):
                values = [s.concentrations[j] for s in node.samples]

                if method == "linear":
                    interp_value = _linear_interp(
                        original_times_float, values, target_time_float
                    )
                elif method == "nearest":
                    interp_value = _nearest_interp(
                        original_times_float, values, target_time_float
                    )
                elif method == "spline":
                    interp_value = _spline_interp(
                        original_times_float, values, target_time_float
                    )
                else:
                    interp_value = _linear_interp(
                        original_times_float, values, target_time_float
                    )

                concentrations.append(interp_value)

            # Create interpolated sample
            sample = TimeSeriesSample(
                sample_id=f"{node_id}_{target_time.strftime('%Y%m%d')}",
                node_id=node_id,
                timestamp=target_time,
                concentrations=concentrations,
            )
            interpolated_samples.append(sample)

        # Create new node
        interpolated_nodes[node_id] = TemporalNode(
            node_id=node_id, samples=interpolated_samples
        )

    return interpolated_nodes


def _linear_interp(
    times: List[float], values: List[float], target_time: float
) -> float:
    """
    Linear interpolation.

    Parameters
    ----------
    times : List[float]
        Sorted time points
    values : List[float]
        Values at time points
    target_time : float
        Time to interpolate to

    Returns
    -------
    float
        Interpolated value
    """
    # Find bracketing indices
    if target_time <= times[0]:
        return values[0]
    if target_time >= times[-1]:
        return values[-1]

    for i in range(len(times) - 1):
        if times[i] <= target_time <= times[i + 1]:
            # Linear interpolation
            t0, t1 = times[i], times[i + 1]
            v0, v1 = values[i], values[i + 1]

            if t1 - t0 < 1e-10:
                return v0

            alpha = (target_time - t0) / (t1 - t0)
            return v0 + alpha * (v1 - v0)

    return values[-1]


def _nearest_interp(
    times: List[float], values: List[float], target_time: float
) -> float:
    """
    Nearest neighbor interpolation.
    """
    if target_time <= times[0]:
        return values[0]
    if target_time >= times[-1]:
        return values[-1]

    # Find nearest
    min_dist = float("inf")
    nearest_idx = 0

    for i, t in enumerate(times):
        dist = abs(t - target_time)
        if dist < min_dist:
            min_dist = dist
            nearest_idx = i

    return values[nearest_idx]


def _spline_interp(
    times: List[float], values: List[float], target_time: float
) -> float:
    """
    Cubic spline interpolation using scipy if available, otherwise fallback to linear.
    """
    try:
        from scipy.interpolate import CubicSpline

        cs = CubicSpline(times, values, bc_type="natural")
        return float(cs(target_time))
    except ImportError:
        # Fallback to linear
        return _linear_interp(times, values, target_time)


def align_time_series_recharge(
    node_u: TemporalNode,
    node_v: TemporalNode,
    rain_dates: List[datetime],
    rain_mm: List[float],
    storage_mm: float,
    recharge_fraction: float = 1.0,
) -> tuple[List[TimeSeriesSample], List[TimeSeriesSample]]:
    """
    Align time series using recharge-driven piston flow.
    For each v sample at t_v, finds t_u such that cumulative recharge between t_u and t_v equals storage_mm.
    """
    if not node_u.samples or not node_v.samples or not rain_dates:
        return [], []

    # Sort rainfall
    rain_data = sorted(zip(rain_dates, rain_mm), key=lambda x: x[0])
    r_dates = [x[0] for x in rain_data]
    r_values = np.array([x[1] for x in rain_data], dtype=float) * recharge_fraction

    # Pre-calculate cumulative recharge for fast lookup
    # cum_recharge[i] is recharge from start up to r_dates[i]
    cum_recharge = np.cumsum(r_values)

    # Helper to get cumulative recharge at any time t
    # We map time t to index in r_dates
    r_timestamps = np.array([d.timestamp() for d in r_dates])

    def get_cum_recharge(t: datetime) -> float:
        ts = t.timestamp()
        if ts < r_timestamps[0]:
            return 0.0
        if ts >= r_timestamps[-1]:
            return float(cum_recharge[-1])
        # Linear interp on cumulative curve
        idx = np.searchsorted(r_timestamps, ts)
        if idx == 0:
            return 0.0
        t0, t1 = r_timestamps[idx - 1], r_timestamps[idx]
        c0, c1 = cum_recharge[idx - 1], cum_recharge[idx]
        frac = (ts - t0) / (t1 - t0)
        return float(c0 + frac * (c1 - c0))

    def find_start_time(t_end: datetime, target_vol: float) -> datetime:
        # We want t_start such that Recharge(t_start, t_end) = target_vol
        # cum(t_end) - cum(t_start) = target_vol
        # cum(t_start) = cum(t_end) - target_vol
        c_end = get_cum_recharge(t_end)
        c_target = c_end - target_vol
        if c_target <= 0:
            # Reached beginning of record
            return r_dates[0]

        # Inverse lookup
        idx = np.searchsorted(cum_recharge, c_target)
        if idx == 0:
            return r_dates[0]
        if idx >= len(cum_recharge):
            return r_dates[-1]

        c0, c1 = cum_recharge[idx - 1], cum_recharge[idx]
        t0, t1 = r_timestamps[idx - 1], r_timestamps[idx]

        if c1 - c0 < 1e-9:
            return datetime.fromtimestamp(t0)

        frac = (c_target - c0) / (c1 - c0)
        ts_res = t0 + frac * (t1 - t0)
        return datetime.fromtimestamp(ts_res)

    aligned_u = []
    aligned_v = []

    # Sort upstream samples
    u_samples_sorted = sorted(node_u.samples, key=lambda s: s.timestamp)
    u_times = [s.timestamp for s in u_samples_sorted]
    u_times_float = [t.timestamp() / 86400.0 for t in u_times]
    n_ions = len(u_samples_sorted[0].concentrations)
    u_values_by_ion = [
        [s.concentrations[j] for s in u_samples_sorted] for j in range(n_ions)
    ]

    for v_sample in node_v.samples:
        t_v = v_sample.timestamp
        t_u = find_start_time(t_v, storage_mm)
        t_u_float = t_u.timestamp() / 86400.0

        # Check bounds
        if t_u < u_times[0] or t_u > u_times[-1]:
            continue

        interpolated_concs = []
        for j in range(n_ions):
            val = _linear_interp(u_times_float, u_values_by_ion[j], t_u_float)
            interpolated_concs.append(val)

        synth_u_sample = TimeSeriesSample(
            sample_id=f"piston_{node_u.node_id}_{t_u.strftime('%Y%m%d')}",
            node_id=node_u.node_id,
            timestamp=t_u,
            concentrations=interpolated_concs,
        )
        aligned_u.append(synth_u_sample)
        aligned_v.append(v_sample)

    return aligned_u, aligned_v


def align_time_series(
    node_u: TemporalNode, node_v: TemporalNode, lag_days: float = 0.0
) -> tuple[List[TimeSeriesSample], List[TimeSeriesSample]]:
    """
    Align two time series with a specified lag using linear interpolation.

    Parameters
    ----------
    node_u : TemporalNode
        Upstream node
    node_v : TemporalNode
        Downstream node
    lag_days : float
        Residence time lag in days (v lags behind u by this amount)

    Returns
    -------
    tuple[List[TimeSeriesSample], List[TimeSeriesSample]]
        (aligned_u_samples, aligned_v_samples) with matching timestamps

    Notes
    -----
    For each v sample at time t, estimates u sample at time (t - lag)
    using linear interpolation of the upstream time series.
    """
    if not node_u.samples or not node_v.samples:
        return [], []

    aligned_u = []
    aligned_v = []

    # Sort upstream samples by time just in case
    u_samples_sorted = sorted(node_u.samples, key=lambda s: s.timestamp)
    u_times = [s.timestamp for s in u_samples_sorted]
    u_times_float = [t.timestamp() / 86400.0 for t in u_times]

    # Pre-extract values for each ion to avoid repeated loops
    n_ions = len(u_samples_sorted[0].concentrations)
    u_values_by_ion = []
    for j in range(n_ions):
        u_values_by_ion.append([s.concentrations[j] for s in u_samples_sorted])

    for v_sample in node_v.samples:
        # Target time for upstream
        target_time_u = v_sample.timestamp - timedelta(days=lag_days)
        target_time_u_float = target_time_u.timestamp() / 86400.0

        # Check if we can interpolate (must be within range of upstream data)
        # We allow a small tolerance for extrapolation at boundaries if needed,
        # but strict containment is safer for transport.
        if target_time_u < u_times[0] or target_time_u > u_times[-1]:
            continue

        # Interpolate upstream concentrations
        interpolated_concs = []
        for j in range(n_ions):
            val = _linear_interp(u_times_float, u_values_by_ion[j], target_time_u_float)
            interpolated_concs.append(val)

        # Create synthetic upstream sample
        synth_u_sample = TimeSeriesSample(
            sample_id=f"interp_{node_u.node_id}_{target_time_u.strftime('%Y%m%d')}",
            node_id=node_u.node_id,
            timestamp=target_time_u,  # Note: this is the *lagged* time
            concentrations=interpolated_concs,
        )

        aligned_u.append(synth_u_sample)
        aligned_v.append(v_sample)

    return aligned_u, aligned_v
