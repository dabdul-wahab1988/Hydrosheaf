"""Causal mass-conserving edge transport operators and mixing algebra for network TTDs.

This module formalizes advective-dispersive age shifting along flow edges:
    g_{u -> v} = T_{uv} @ g_u
and node mixing conservation:
    g_v = rho_v * r_v + sum_{u in Pa(v)} pi_{uv} * T_{uv} @ g_u
where rho_v + sum pi_{uv} = 1.0.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import math
from typing import Any, Mapping, Optional, Tuple

import numpy as np

from .ttd_grid import TTDGrid


@dataclass(frozen=True)
class EdgeTransportOperator:
    """Causal, column-stochastic transition matrix mapping upstream TTD to downstream arrival.

    Parameters
    ----------
    edge_id : Tuple[str, str]
        Directed edge identifier (source_u, target_v).
    matrix : np.ndarray
        Shape (K, K) transition matrix T_{uv} where T_{jk} is the probability of age
        tau_j at node v given age tau_k at node u.
    delta_tau_years : float
        Mean advective travel time along the edge (years).
    dispersion : float
        Dispersion parameter DP = alpha_L / L = D / (v L).
    length_m : float
        Physical flowpath distance (meters).
    pore_velocity_m_y : float
        Effective pore velocity (meters/year).
    metadata : Mapping[str, Any]
        Additional physical parameters (porosity, hydraulic conductivity, etc.).
    """

    edge_id: Tuple[str, str]
    matrix: np.ndarray
    delta_tau_years: float
    dispersion: float = 0.05
    length_m: float = 100.0
    pore_velocity_m_y: float = 10.0
    metadata: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        mat = np.asarray(self.matrix, dtype=float).copy()
        if mat.ndim != 2 or mat.shape[0] != mat.shape[1]:
            raise ValueError(f"Transport matrix must be square (K x K), got shape {mat.shape}.")
        if not np.all(np.isfinite(mat)):
            raise ValueError("Transport matrix contains non-finite values.")
        if np.any(mat < -1e-12):
            raise ValueError("Transport matrix must be non-negative.")
        # Rows are downstream arrival ages j and columns are upstream source
        # ages k.  Causality requires tau_j >= tau_k, so entries strictly
        # above the main diagonal (j < k) are forbidden.
        if np.any(np.triu(mat, k=1) > 1e-12):
            raise ValueError(
                "Transport matrix violates temporal causality: mass appears below "
                "the source-age diagonal."
            )

        delta_tau = float(self.delta_tau_years)
        dispersion = float(self.dispersion)
        length_m = float(self.length_m)
        velocity = float(self.pore_velocity_m_y)
        if not math.isfinite(delta_tau) or delta_tau < 0.0:
            raise ValueError("delta_tau_years must be finite and non-negative.")
        if not math.isfinite(dispersion) or dispersion < 0.0:
            raise ValueError("dispersion must be finite and non-negative.")
        if not math.isfinite(length_m) or length_m <= 0.0:
            raise ValueError("length_m must be finite and strictly positive.")
        if not math.isfinite(velocity) or velocity <= 0.0:
            raise ValueError("pore_velocity_m_y must be finite and strictly positive.")

        # Ensure non-negativity and exact column stochasticity: sum_j T_{jk} = 1.0
        mat[mat < 0.0] = 0.0
        col_sums = np.sum(mat, axis=0)
        zero_cols = col_sums <= 0.0
        if np.any(zero_cols):
            # Fallback for empty columns: identity (no change)
            mat[:, zero_cols] = np.eye(mat.shape[0])[:, zero_cols]
            col_sums = np.sum(mat, axis=0)
        mat = mat / col_sums[None, :]

        mat.setflags(write=False)
        object.__setattr__(self, "matrix", mat)
        object.__setattr__(self, "delta_tau_years", delta_tau)
        object.__setattr__(self, "dispersion", dispersion)
        object.__setattr__(self, "length_m", length_m)
        object.__setattr__(self, "pore_velocity_m_y", velocity)
        object.__setattr__(self, "metadata", dict(self.metadata))

    @property
    def n_bins(self) -> int:
        return self.matrix.shape[0]

    def transport(self, g_upstream: np.ndarray) -> np.ndarray:
        """Propagate upstream mass distribution to downstream arrival: g_down = T @ g_up."""
        arr = np.asarray(g_upstream, dtype=float)
        if arr.shape != (self.n_bins,):
            raise ValueError(f"Input distribution shape {arr.shape} does not match operator dimension {self.n_bins}.")
        return self.matrix @ arr


def build_advection_dispersion_operator(
    edge_id: Tuple[str, str],
    grid: TTDGrid,
    *,
    delta_tau_years: float,
    dispersion: float = 0.05,
    length_m: float = 100.0,
    pore_velocity_m_y: Optional[float] = None,
    metadata: Optional[Mapping[str, Any]] = None,
) -> EdgeTransportOperator:
    """Construct an advection-dispersion edge transport matrix T_{uv}.

    For each upstream age tau_k:
    The arrival age is tau_arrival = tau_k + Delta tau.
    If dispersion > 0:
        h(tau - tau_k) is modeled as a 1D advective-dispersive kernel:
        h(s) = 1 / sqrt(4 pi DP Delta tau s) * exp(-(s - Delta tau)^2 / (4 DP Delta tau s))
        evaluated for s = tau_j - tau_k > 0.
    If dispersion == 0:
        pure piston shift to tau_k + Delta tau.

    Parameters
    ----------
    edge_id : Tuple[str, str]
        (source_node, target_node).
    grid : TTDGrid
        Transit-time grid.
    delta_tau_years : float
        Mean travel time in years (Delta tau = L / v).
    dispersion : float
        Dimensionless dispersion parameter (1/Pe = D / (v L)). Default 0.05.
    length_m : float
        Flow distance in meters.
    pore_velocity_m_y : Optional[float]
        Pore velocity in m/yr (if None, derived as length_m / delta_tau_years).
    """
    dt_travel = float(delta_tau_years)
    dp = float(dispersion)
    if not math.isfinite(dt_travel) or dt_travel < 0.0:
        raise ValueError("delta_tau_years must be finite and non-negative.")
    if not math.isfinite(dp) or dp < 0.0:
        raise ValueError("dispersion must be finite and non-negative.")
    if not math.isfinite(float(length_m)) or float(length_m) <= 0.0:
        raise ValueError("length_m must be finite and strictly positive.")
    if pore_velocity_m_y is not None and (
        not math.isfinite(float(pore_velocity_m_y)) or float(pore_velocity_m_y) <= 0.0
    ):
        raise ValueError("pore_velocity_m_y must be finite and strictly positive when supplied.")
    taus = grid.taus
    k_bins = grid.n_bins

    matrix = np.zeros((k_bins, k_bins), dtype=float)

    if dt_travel < 1e-6:
        # Zero travel time -> identity operator
        matrix = np.eye(k_bins, dtype=float)
    else:
        for k in range(k_bins):
            tau_src = taus[k]
            target_mean = tau_src + dt_travel

            if dp <= 1e-4:
                # Pure piston flow shift: find nearest bin or interpolate between adjacent bins
                if target_mean >= taus[-1]:
                    matrix[-1, k] = 1.0
                else:
                    idx = int(np.searchsorted(taus, target_mean))
                    if idx == 0:
                        matrix[0, k] = 1.0
                    else:
                        t0 = taus[idx - 1]
                        t1 = taus[idx]
                        frac = (target_mean - t0) / max(1e-9, t1 - t0)
                        matrix[idx - 1, k] = 1.0 - frac
                        matrix[idx, k] = frac
            else:
                # Dispersive pulse starting at tau_src.  A small positive
                # evaluation lag avoids the singularity at s = 0 while the
                # exact causal support is retained.
                # Valid arrival ages must satisfy tau_j >= tau_src (causality)
                active_j = np.where(taus >= tau_src)[0]
                if len(active_j) == 0:
                    matrix[-1, k] = 1.0
                    continue

                lags = taus[active_j] - tau_src
                # Avoid division by zero at lag = 0
                lags_safe = np.maximum(lags, 1e-4)

                # 1D dispersion kernel (Ogata-Banks / Kreft-Zuber flux concentration)
                denom = 4.0 * math.pi * dp * dt_travel * lags_safe
                exponent = -((lags_safe - dt_travel) ** 2) / (4.0 * dp * dt_travel * lags_safe)
                pulse = np.exp(np.clip(exponent, -80.0, 0.0)) / np.sqrt(denom)

                # Apply quadrature weights to convert continuous density to discrete mass
                mass_col = pulse * grid.weights[active_j]
                tot_mass = float(np.sum(mass_col))
                if tot_mass > 0.0:
                    mass_col = mass_col / tot_mass
                    matrix[active_j, k] = mass_col
                else:
                    # If pulse lies beyond taus[-1], accumulate in final bin
                    matrix[-1, k] = 1.0

    velocity = pore_velocity_m_y if pore_velocity_m_y is not None else (length_m / max(dt_travel, 1e-3))
    operator_metadata = dict(metadata or {})
    operator_metadata.setdefault(
        "age_grid_truncation_possible",
        bool(dt_travel > 0.0 and np.any(taus + dt_travel > taus[-1])),
    )
    operator_metadata.setdefault(
        "causality_checked",
        True,
    )

    return EdgeTransportOperator(
        edge_id=edge_id,
        matrix=matrix,
        delta_tau_years=dt_travel,
        dispersion=dp,
        length_m=float(length_m),
        pore_velocity_m_y=float(velocity),
        metadata=operator_metadata,
    )


def build_local_recharge_distribution(
    grid: TTDGrid,
    *,
    mean_recharge_age_years: float = 1.0,
    model: str = "exponential",
) -> np.ndarray:
    """Construct an entry transit-time distribution for local modern recharge r_v.

    Parameters
    ----------
    grid : TTDGrid
        The transit-time grid.
    mean_recharge_age_years : float
        Mean age of newly recharged groundwater (default 1.0 year).
    model : str
        "exponential" or "uniform".

    Returns
    -------
    np.ndarray
        Normalized probability mass vector r_v in Delta^(K-1).
    """
    taus = grid.taus
    tau_r = float(mean_recharge_age_years)
    if not math.isfinite(tau_r) or tau_r <= 0.0:
        raise ValueError("mean_recharge_age_years must be finite and strictly positive.")
    tau_r = max(0.1, tau_r)

    if model == "exponential":
        pdf = (1.0 / tau_r) * np.exp(-taus / tau_r)
    elif model == "uniform":
        cutoff = 2.0 * tau_r
        pdf = np.where(taus <= cutoff, 1.0 / cutoff, 0.0)
    else:
        raise ValueError(f"Unknown recharge model: {model}")

    # Convert density to mass vector on discrete bins
    mass = pdf * grid.weights
    total = float(np.sum(mass))
    if not math.isfinite(total) or total <= 0.0:
        raise ValueError("Recharge model produced no finite positive probability mass.")
    return mass / total


@dataclass(frozen=True)
class NodeMixingSpecification:
    """Declared mixing fractions and local recharge model for a network node."""

    node_id: str
    local_fraction: float
    upstream_weights: Mapping[str, float]
    recharge_distribution: np.ndarray

    def __post_init__(self) -> None:
        node_id = str(self.node_id)
        rho = float(self.local_fraction)
        if not math.isfinite(rho) or rho < 0.0:
            raise ValueError("local_fraction must be finite and non-negative.")

        weights: dict[str, float] = {}
        for upstream_id, raw_weight in self.upstream_weights.items():
            weight = float(raw_weight)
            if not math.isfinite(weight) or weight < 0.0:
                raise ValueError("upstream mixing weights must be finite and non-negative.")
            weights[str(upstream_id)] = weight

        r_arr = np.asarray(self.recharge_distribution, dtype=float).copy()
        if r_arr.ndim != 1 or r_arr.size == 0:
            raise ValueError("recharge_distribution must be a non-empty one-dimensional vector.")
        if not np.all(np.isfinite(r_arr)) or np.any(r_arr < 0.0):
            raise ValueError("recharge_distribution must be finite and non-negative.")

        total = rho + sum(weights.values())
        if not math.isfinite(total) or total <= 0.0:
            raise ValueError(
                "local_fraction plus upstream mixing weights must have positive mass."
            )
        rho /= total
        weights = {u: w / total for u, w in weights.items()}

        r_sum = float(np.sum(r_arr))
        if not math.isfinite(r_sum) or r_sum <= 0.0:
            raise ValueError("recharge_distribution must contain positive total mass.")
        r_arr /= r_sum

        r_arr.setflags(write=False)
        object.__setattr__(self, "node_id", node_id)
        object.__setattr__(self, "local_fraction", rho)
        object.__setattr__(self, "upstream_weights", weights)
        object.__setattr__(self, "recharge_distribution", r_arr)

    def composite_distribution(
        self,
        upstream_distributions: Mapping[str, np.ndarray],
        transport_operators: Mapping[Tuple[str, str], EdgeTransportOperator],
    ) -> np.ndarray:
        """Compute the expected downstream arrival distribution g_v under declared mixing."""
        g_comp = self.local_fraction * self.recharge_distribution.copy()
        for parent_id, weight in self.upstream_weights.items():
            if weight <= 0.0:
                continue
            g_up = upstream_distributions.get(parent_id)
            if g_up is None:
                continue
            op = transport_operators.get((parent_id, self.node_id))
            if op is not None:
                g_routed = op.transport(g_up)
            else:
                g_routed = g_up
            g_comp += weight * g_routed
        return g_comp
