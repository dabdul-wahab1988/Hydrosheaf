"""Dynamic, state-dependent graph edge-kernel representations and regularizers.

This module provides immutable, typed representations of time-varying transit-time
distribution (TTD) kernels on directed aquifer flow networks, including:
- DynamicEdgeKernel: 3D tensor of edge kernels h_uv(tau, t)
- Second-order lag curvature operator D_lag^2
- Temporal smoothness operator D_time and Total Variation
- Basis expansions (phase-bin, harmonic Fourier modes, B-splines)
- Exact stationary compatibility conversions
"""

from __future__ import annotations

from dataclasses import dataclass, field
import hashlib
import json
import math
from typing import Any, Mapping, Optional, Sequence

import numpy as np


@dataclass(frozen=True)
class DynamicEdgeKernel:
    """Immutable, typed container for dynamic graph edge-kernel tensors.

    Parameters
    ----------
    edge_ids : tuple[str, ...]
        Tuple of directed edge identifiers (e.g. ("R->A", "A->B")).
    time_grid : np.ndarray
        Strictly increasing array of output observation/simulation time coordinates.
    lag_grid : np.ndarray
        Strictly increasing array of non-negative transit-time lag coordinates.
    values : np.ndarray
        Kernel tensor of shape (n_edges, n_output_times, n_lags).
        values[e, t, tau] represents h_{uv}(tau, t).
    normalization_tolerance : float
        Tolerance for checking simplex normalization sum_tau h(tau, t) = 1.0.
    kernel_mode : str
        One of "stationary", "phase", or "time".
    source_commitment : str
        Optional cryptographic hash/commitment to prevent data leakage.
    metadata : dict[str, Any]
        Arbitrary provenance and parameter metadata.
    """

    edge_ids: tuple[str, ...]
    time_grid: np.ndarray
    lag_grid: np.ndarray
    values: np.ndarray
    normalization_tolerance: float = 1e-5
    kernel_mode: str = "time"
    source_commitment: str = ""
    metadata: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.kernel_mode not in {"stationary", "phase", "time"}:
            raise ValueError(
                f"kernel_mode must be 'stationary', 'phase', or 'time', got {self.kernel_mode!r}"
            )

        if not isinstance(self.edge_ids, tuple):
            object.__setattr__(self, "edge_ids", tuple(self.edge_ids))

        time_grid = np.asarray(self.time_grid, dtype=float)
        lag_grid = np.asarray(self.lag_grid, dtype=float)
        values = np.asarray(self.values, dtype=float)

        # The dataclass is frozen, but NumPy arrays are mutable unless their
        # write flag is disabled.  Copy and freeze the arrays so a benchmark
        # cannot mutate an already-committed kernel after construction.
        time_grid = np.array(time_grid, dtype=float, copy=True)
        lag_grid = np.array(lag_grid, dtype=float, copy=True)
        values = np.array(values, dtype=float, copy=True)

        object.__setattr__(self, "time_grid", time_grid)
        object.__setattr__(self, "lag_grid", lag_grid)
        object.__setattr__(self, "values", values)

        if "step_days" in self.metadata:
            step_days = float(self.metadata["step_days"])
            if not math.isfinite(step_days) or step_days <= 0.0:
                raise ValueError("metadata['step_days'] must be a positive finite value")

        n_edges = len(self.edge_ids)
        n_times = len(time_grid)
        n_lags = len(lag_grid)

        if values.shape != (n_edges, n_times, n_lags):
            raise ValueError(
                f"values shape {values.shape} must match (n_edges={n_edges}, "
                f"n_times={n_times}, n_lags={n_lags})"
            )

        if not np.all(np.isfinite(values)):
            raise ValueError("All dynamic kernel values must be finite numbers.")

        if np.any(values < -1e-9):
            min_val = float(np.min(values))
            raise ValueError(f"Dynamic kernel values must be non-negative; found {min_val}.")

        # Causal support check
        if np.any(lag_grid < -1e-9):
            raise ValueError("Lag grid must contain non-negative lags (causal support).")

        # Monotonicity check
        if n_times > 1 and np.any(np.diff(time_grid) <= 0.0):
            raise ValueError("time_grid must be strictly monotonically increasing.")
        if n_lags > 1 and np.any(np.diff(lag_grid) <= 0.0):
            raise ValueError("lag_grid must be strictly monotonically increasing.")

        # Simplex normalization check for every edge and output time
        row_sums = np.sum(values, axis=-1)
        discrepancies = np.abs(row_sums - 1.0)
        max_disc = float(np.max(discrepancies)) if discrepancies.size > 0 else 0.0
        if max_disc > self.normalization_tolerance:
            raise ValueError(
                f"Dynamic kernel normalization violation: max discrepancy from 1.0 is "
                f"{max_disc:.2e} (tolerance={self.normalization_tolerance:.2e})."
            )

        # In stationary mode, all time slices must be identical
        if self.kernel_mode == "stationary" and n_times > 1:
            first_slice = values[:, 0:1, :]
            if not np.allclose(values, first_slice, atol=1e-7):
                raise ValueError("kernel_mode='stationary' requires identical slices across time.")

        time_grid.flags.writeable = False
        lag_grid.flags.writeable = False
        values.flags.writeable = False

    @property
    def n_edges(self) -> int:
        return len(self.edge_ids)

    @property
    def n_output_times(self) -> int:
        return len(self.time_grid)

    @property
    def n_lags(self) -> int:
        return len(self.lag_grid)

    def edge_index(self, edge_id: str) -> int:
        try:
            return self.edge_ids.index(edge_id)
        except ValueError:
            raise KeyError(f"Edge {edge_id!r} not found in dynamic kernel edges {self.edge_ids}")

    def kernel_at(self, edge_id: str, time_step: int | float) -> np.ndarray:
        """Extract 1D lag distribution h_{uv}(tau, t) at a given time or time step index."""
        e_idx = self.edge_index(edge_id)
        if isinstance(time_step, (int, np.integer)):
            t_idx = int(time_step)
            if not 0 <= t_idx < self.n_output_times:
                raise IndexError(f"Time index {t_idx} out of range [0, {self.n_output_times})")
            return np.array(self.values[e_idx, t_idx, :], copy=True)
        # Interpolate if floating time
        t_val = float(time_step)
        if t_val <= self.time_grid[0]:
            return np.array(self.values[e_idx, 0, :], copy=True)
        if t_val >= self.time_grid[-1]:
            return np.array(self.values[e_idx, -1, :], copy=True)
        # Linear interpolation between adjacent time points
        idx = int(np.searchsorted(self.time_grid, t_val))
        t0, t1 = self.time_grid[idx - 1], self.time_grid[idx]
        alpha = (t_val - t0) / (t1 - t0)
        interp = (1.0 - alpha) * self.values[e_idx, idx - 1, :] + alpha * self.values[e_idx, idx, :]
        interp = np.maximum(interp, 0.0)
        s = float(np.sum(interp))
        return interp / s if s > 0 else interp

    def mean_transit_time(self, edge_id: Optional[str] = None) -> np.ndarray:
        """Compute mean transit time mu(t) = sum_tau tau * h(tau, t)."""
        lags = self.lag_grid
        if edge_id is not None:
            e_idx = self.edge_index(edge_id)
            return self.values[e_idx] @ lags
        return self.values @ lags

    def young_water_fraction(
        self, cutoff_days: float = 90.0, edge_id: Optional[str] = None
    ) -> np.ndarray:
        """Compute young water fraction F_y(t) = sum_{tau <= cutoff} h(tau, t)."""
        mask = self.lag_grid <= cutoff_days
        if not np.any(mask):
            shape = (self.n_output_times,) if edge_id is not None else (self.n_edges, self.n_output_times)
            return np.zeros(shape, dtype=float)
        if edge_id is not None:
            e_idx = self.edge_index(edge_id)
            return np.sum(self.values[e_idx][..., mask], axis=-1)
        return np.sum(self.values[..., mask], axis=-1)

    @classmethod
    def from_stationary(
        cls,
        edge_ids: Sequence[str],
        lag_grid: Sequence[float],
        stationary_kernels: Mapping[str, Sequence[float]] | np.ndarray,
        time_grid: Sequence[float],
        metadata: Optional[Mapping[str, Any]] = None,
    ) -> DynamicEdgeKernel:
        """Embed stationary edge kernels into an identical-across-time DynamicEdgeKernel."""
        edges = tuple(edge_ids)
        t_grid = np.asarray(time_grid, dtype=float)
        l_grid = np.asarray(lag_grid, dtype=float)
        n_edges = len(edges)
        n_times = len(t_grid)
        n_lags = len(l_grid)

        if isinstance(stationary_kernels, np.ndarray):
            if stationary_kernels.shape == (n_edges, n_lags):
                base_mat = stationary_kernels
            else:
                raise ValueError(
                    f"stationary_kernels ndarray shape {stationary_kernels.shape} "
                    f"must be ({n_edges}, {n_lags})"
                )
        else:
            base_mat = np.zeros((n_edges, n_lags), dtype=float)
            for i, edge in enumerate(edges):
                if edge not in stationary_kernels:
                    raise KeyError(f"Missing stationary kernel for edge {edge!r}")
                arr = np.asarray(stationary_kernels[edge], dtype=float)
                if len(arr) != n_lags:
                    raise ValueError(f"Kernel for edge {edge} has length {len(arr)} != {n_lags}")
                base_mat[i] = arr

        # Normalize rows if needed
        row_sums = np.sum(base_mat, axis=1, keepdims=True)
        row_sums = np.where(row_sums > 0, row_sums, 1.0)
        base_mat = np.maximum(base_mat / row_sums, 0.0)

        # Broadcast across time: shape (n_edges, n_times, n_lags)
        values = np.repeat(base_mat[:, np.newaxis, :], n_times, axis=1)

        meta = dict(metadata or {})
        meta["created_from"] = "stationary_embedding"

        return cls(
            edge_ids=edges,
            time_grid=t_grid,
            lag_grid=l_grid,
            values=values,
            kernel_mode="stationary",
            metadata=meta,
        )


def build_lag_curvature_matrix(lags: Sequence[int | float]) -> np.ndarray:
    """Construct non-uniform 2nd-order finite difference curvature operator D_2 on lags.

    For interior lag j:
        D_2[j, j-1] = 2 / (dt1 * (dt1 + dt2))
        D_2[j, j]   = -2 / (dt1 * dt2)
        D_2[j, j+1] = 2 / (dt2 * (dt1 + dt2))
    where dt1 = lag[j] - lag[j-1], dt2 = lag[j+1] - lag[j].
    """
    n = len(lags)
    if n < 3:
        return np.zeros((0, n), dtype=float)

    rows: list[np.ndarray] = []
    for j in range(1, n - 1):
        dt1 = float(lags[j] - lags[j - 1])
        dt2 = float(lags[j + 1] - lags[j])
        if dt1 <= 0.0 or dt2 <= 0.0:
            raise ValueError("Lags must be strictly monotonically increasing.")
        row = np.zeros(n, dtype=float)
        coeff = 2.0 / (dt1 + dt2)
        row[j - 1] = coeff / dt1
        row[j] = -coeff * (1.0 / dt1 + 1.0 / dt2)
        row[j + 1] = coeff / dt2
        rows.append(row)

    return np.vstack(rows)


def build_temporal_smoothness_matrix(time_grid: Sequence[int | float]) -> np.ndarray:
    """Construct 1st-order forward difference operator D_time across output times.

    D_time[t, t] = -1 / dt, D_time[t, t+1] = 1 / dt.
    """
    t_arr = np.asarray(time_grid, dtype=float)
    n = len(t_arr)
    if n < 2:
        return np.zeros((0, n), dtype=float)

    dt = np.diff(t_arr)
    if np.any(dt <= 0.0):
        raise ValueError("time_grid must be strictly monotonically increasing.")

    d_mat = np.zeros((n - 1, n), dtype=float)
    for i in range(n - 1):
        d_mat[i, i] = -1.0 / dt[i]
        d_mat[i, i + 1] = 1.0 / dt[i]
    return d_mat


def build_phase_basis_matrix(
    time_grid: Sequence[int | float],
    season_period: float,
    n_phases: int,
) -> np.ndarray:
    """Construct phase-bin indicator basis matrix B of shape (n_times, n_phases)."""
    t_arr = np.asarray(time_grid, dtype=float)
    n_times = len(t_arr)
    if n_phases < 1:
        raise ValueError("n_phases must be >= 1")
    if season_period <= 0.0:
        raise ValueError("season_period must be positive")

    phases = np.floor(np.mod(t_arr, season_period) / season_period * n_phases).astype(int)
    phases = np.clip(phases, 0, n_phases - 1)

    basis = np.zeros((n_times, n_phases), dtype=float)
    for t in range(n_times):
        basis[t, phases[t]] = 1.0
    return basis


def build_harmonic_basis_matrix(
    time_grid: Sequence[int | float],
    season_period: float,
    n_harmonics: int = 2,
) -> np.ndarray:
    """Construct smooth Fourier harmonic basis [1, cos(w*t), sin(w*t), cos(2w*t), ...].

    Returns matrix of shape (n_times, 2*n_harmonics + 1).
    """
    t_arr = np.asarray(time_grid, dtype=float)
    n_times = len(t_arr)
    if season_period <= 0.0:
        raise ValueError("season_period must be positive")

    cols = [np.ones((n_times, 1), dtype=float)]
    w0 = 2.0 * np.pi / season_period
    for k in range(1, n_harmonics + 1):
        cols.append(np.cos(k * w0 * t_arr)[:, np.newaxis])
        cols.append(np.sin(k * w0 * t_arr)[:, np.newaxis])
    return np.hstack(cols)
