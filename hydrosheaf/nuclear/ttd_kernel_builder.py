"""Multi-tracer forward kernel matrix construction for non-parametric TTD inversion.

This module maps multi-tracer observation panels (3H, 14C, 39Ar, 85Kr, SF6, CFCs, 4He)
to linear forward operator matrices A_i on discrete TTD grids, grounded in atmospheric
input histories and geochemical corrections.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import math
from typing import Any, Dict, Mapping, Optional, Sequence, Tuple

import numpy as np

from .input_history import InputHistory
from .joint_lpm import normalize_tracer_key, tracer_response_kernel
from .ttd_grid import TTDGrid


@dataclass(frozen=True)
class TracerObservation:
    """A single tracer observation with its reported measurement uncertainty."""

    tracer: str
    value: float
    sigma: float
    units: str = ""
    weight: float = 1.0

    def __post_init__(self) -> None:
        key = normalize_tracer_key(self.tracer)
        object.__setattr__(self, "tracer", key)
        value = float(self.value)
        sigma = float(self.sigma)
        weight = float(self.weight)
        if not math.isfinite(value):
            raise ValueError(f"Observed value for {self.tracer} must be finite.")
        if not math.isfinite(sigma) or sigma <= 0.0:
            raise ValueError(f"Uncertainty sigma for {self.tracer} must be strictly positive.")
        if not math.isfinite(weight) or weight < 0.0:
            raise ValueError(f"Weight for {self.tracer} cannot be negative.")
        object.__setattr__(self, "value", value)
        object.__setattr__(self, "sigma", sigma)
        object.__setattr__(self, "weight", weight)


@dataclass(frozen=True)
class NodeTracerPanel:
    """Complete observation panel for a network node (well/spring/piezometer)."""

    node_id: str
    sample_year: float
    observations: Tuple[TracerObservation, ...]
    head_m: Optional[float] = None
    chemistry: Mapping[str, float] = field(default_factory=dict)

    def __post_init__(self) -> None:
        node_id = str(self.node_id).strip()
        if not node_id:
            raise ValueError("node_id must be non-empty.")
        sample_year = float(self.sample_year)
        if not math.isfinite(sample_year):
            raise ValueError("sample_year must be finite.")
        head_m = None if self.head_m is None else float(self.head_m)
        if head_m is not None and not math.isfinite(head_m):
            raise ValueError("head_m must be finite when supplied.")
        object.__setattr__(self, "node_id", node_id)
        object.__setattr__(self, "sample_year", sample_year)
        object.__setattr__(self, "head_m", head_m)
        object.__setattr__(self, "observations", tuple(self.observations))
        object.__setattr__(self, "chemistry", dict(self.chemistry))

    @property
    def tracer_names(self) -> Tuple[str, ...]:
        return tuple(obs.tracer for obs in self.observations)


@dataclass(frozen=True)
class MultiTracerForwardSystem:
    """Discrete linear forward system A g = c for a single network node.

    Parameters
    ----------
    node_id : str
        Identifier of the network node.
    sample_year : float
        Sampling year (e.g. 2024.5).
    grid : TTDGrid
        The discrete transit-time grid.
    tracers : Tuple[str, ...]
        Normalized tracer symbols in row order.
    matrix : np.ndarray
        Shape (M, K) forward matrix A where row m is the tracer response kernel.
    observations : np.ndarray
        Shape (M,) observed tracer concentration vector c_obs.
    sigmas : np.ndarray
        Shape (M,) observation uncertainties sigma.
    weights : np.ndarray
        Shape (M,) user/reliability weights w_m.
    head_m : Optional[float]
        Hydraulic head elevation in meters, if measured.
    chemistry : Mapping[str, float]
        Ancillary hydrochemical concentrations (mg/L).
    """

    node_id: str
    sample_year: float
    grid: TTDGrid
    tracers: Tuple[str, ...]
    matrix: np.ndarray
    observations: np.ndarray
    sigmas: np.ndarray
    weights: np.ndarray
    head_m: Optional[float] = None
    chemistry: Mapping[str, float] = field(default_factory=dict)

    def __post_init__(self) -> None:
        m_arr = np.asarray(self.matrix, dtype=float).copy()
        c_arr = np.asarray(self.observations, dtype=float).copy()
        s_arr = np.asarray(self.sigmas, dtype=float).copy()
        w_arr = np.asarray(self.weights, dtype=float).copy()
        tracer_names = tuple(normalize_tracer_key(name) for name in self.tracers)
        node_id = str(self.node_id).strip()
        sample_year = float(self.sample_year)
        head_m = None if self.head_m is None else float(self.head_m)

        n_tracers = len(tracer_names)
        n_bins = self.grid.n_bins

        if not node_id:
            raise ValueError("node_id must be non-empty.")
        if not math.isfinite(sample_year):
            raise ValueError("sample_year must be finite.")
        if head_m is not None and not math.isfinite(head_m):
            raise ValueError("head_m must be finite when supplied.")

        if m_arr.shape != (n_tracers, n_bins):
            raise ValueError(f"Forward matrix shape {m_arr.shape} must be ({n_tracers}, {n_bins}).")
        if c_arr.shape != (n_tracers,):
            raise ValueError(f"Observations shape {c_arr.shape} must match tracer count {n_tracers}.")
        if s_arr.shape != (n_tracers,):
            raise ValueError(f"Sigmas shape {s_arr.shape} must match tracer count {n_tracers}.")
        if w_arr.shape != (n_tracers,):
            raise ValueError(f"Weights shape {w_arr.shape} must match tracer count {n_tracers}.")
        if len(set(tracer_names)) != len(tracer_names):
            raise ValueError("Tracer names must be unique within a forward system.")
        if not np.all(np.isfinite(m_arr)) or not np.all(np.isfinite(c_arr)):
            raise ValueError("Forward matrix and observations must contain only finite values.")
        if not np.all(np.isfinite(s_arr)) or np.any(s_arr <= 0.0):
            raise ValueError("Forward-system sigmas must be finite and strictly positive.")
        if not np.all(np.isfinite(w_arr)) or np.any(w_arr < 0.0):
            raise ValueError("Forward-system weights must be finite and non-negative.")

        object.__setattr__(self, "node_id", node_id)
        object.__setattr__(self, "sample_year", sample_year)
        object.__setattr__(self, "head_m", head_m)
        object.__setattr__(self, "tracers", tracer_names)
        m_arr.setflags(write=False)
        c_arr.setflags(write=False)
        s_arr.setflags(write=False)
        w_arr.setflags(write=False)

        object.__setattr__(self, "matrix", m_arr)
        object.__setattr__(self, "observations", c_arr)
        object.__setattr__(self, "sigmas", s_arr)
        object.__setattr__(self, "weights", w_arr)
        object.__setattr__(self, "chemistry", dict(self.chemistry))

    @property
    def n_tracers(self) -> int:
        return len(self.tracers)

    def predict(self, g: np.ndarray) -> np.ndarray:
        """Compute predicted tracer concentrations: c_pred = A @ g."""
        return self.matrix @ np.asarray(g, dtype=float)

    def residual(self, g: np.ndarray) -> np.ndarray:
        """Compute raw residual vector: r = A @ g - c_obs."""
        return self.predict(g) - self.observations

    def standardized_residual(self, g: np.ndarray) -> np.ndarray:
        """Compute standardized residual vector: (A @ g - c_obs) / sigma."""
        return (self.predict(g) - self.observations) / self.sigmas

    def chi_squared(self, g: np.ndarray) -> float:
        """Compute weighted chi-squared misfit: sum_m w_m * ((A g - c)_m / sigma_m)^2."""
        std_res = self.standardized_residual(g)
        return float(np.sum(self.weights * (std_res ** 2)))

    def condition_number(self) -> float:
        """Compute matrix 2-norm condition number of the row-standardized kernel."""
        if self.n_tracers == 0:
            return 1.0
        row_normed = self.matrix / self.sigmas[:, None]
        try:
            return float(np.linalg.cond(row_normed))
        except Exception:
            return float("inf")

    def effective_rank(self, tolerance: float = 1e-4) -> float:
        """Compute numerical rank of the row-standardized forward kernel."""
        if self.n_tracers == 0:
            return 0.0
        row_normed = self.matrix / self.sigmas[:, None]
        s = np.linalg.svd(row_normed, compute_uv=False)
        if len(s) == 0 or s[0] <= 0:
            return 0.0
        rel_s = s / s[0]
        return float(np.sum(rel_s > tolerance))


def build_forward_system(
    panel: NodeTracerPanel,
    grid: TTDGrid,
    *,
    histories: Optional[Mapping[str, InputHistory]] = None,
    initial_c14_pmc: float = 100.0,
    q_carbon_correction: float = 1.0,
    helium4_background_ccpg: float = 4.6e-8,
    helium4_accumulation_rate_ccpg_per_year: float = 2.0e-11,
    prediction_scale_factors: Optional[Mapping[str, float]] = None,
) -> MultiTracerForwardSystem:
    """Construct the MultiTracerForwardSystem for a given node observation panel.

    Parameters
    ----------
    panel : NodeTracerPanel
        The observed tracers and uncertainties at the node.
    grid : TTDGrid
        The discrete transit-time grid.
    histories : Optional[Mapping[str, InputHistory]]
        Custom tracer recharge input histories (defaults to global if None).
    initial_c14_pmc : float
        Atmospheric initial 14C activity (default 100.0 pmc).
    q_carbon_correction : float
        Carbonate dissolution dilution factor q in (0, 1] (default 1.0).
        Effective initial 14C is q * initial_c14_pmc.
    helium4_background_ccpg : float
        Atmospheric solubility equilibrium baseline for 4He (ccSTP/g).
    helium4_accumulation_rate_ccpg_per_year : float
        Crustal/radiogenic 4He in-growth accumulation rate.
    prediction_scale_factors : Optional[Mapping[str, float]]
        Optional calibration multiplier per tracer.

    Returns
    -------
    MultiTracerForwardSystem
        Immutable linear system ready for joint network inversion.
    """
    if not math.isfinite(float(q_carbon_correction)) or q_carbon_correction <= 0.0 or q_carbon_correction > 1.0:
        raise ValueError(f"q_carbon_correction must be in (0, 1], got {q_carbon_correction}.")

    effective_c14_pmc = float(initial_c14_pmc) * float(q_carbon_correction)

    tracers: list[str] = []
    rows: list[np.ndarray] = []
    observed_vals: list[float] = []
    sigmas: list[float] = []
    weights: list[float] = []

    for obs in panel.observations:
        row = tracer_response_kernel(
            obs.tracer,
            grid.taus,
            panel.sample_year,
            histories=histories,
            initial_c14_pmc=effective_c14_pmc,
            helium4_background_ccpg=helium4_background_ccpg,
            helium4_accumulation_rate_ccpg_per_year=helium4_accumulation_rate_ccpg_per_year,
            prediction_scale_factors=prediction_scale_factors,
        )
        tracers.append(obs.tracer)
        rows.append(row)
        observed_vals.append(float(obs.value))
        sigmas.append(float(obs.sigma))
        weights.append(float(obs.weight))

    if rows:
        matrix = np.vstack(rows)
    else:
        matrix = np.empty((0, grid.n_bins), dtype=float)

    return MultiTracerForwardSystem(
        node_id=panel.node_id,
        sample_year=panel.sample_year,
        grid=grid,
        tracers=tuple(tracers),
        matrix=matrix,
        observations=np.asarray(observed_vals, dtype=float),
        sigmas=np.asarray(sigmas, dtype=float),
        weights=np.asarray(weights, dtype=float),
        head_m=panel.head_m,
        chemistry=panel.chemistry,
    )


def build_network_forward_systems(
    panels: Sequence[NodeTracerPanel],
    grid: TTDGrid,
    **kwargs: Any,
) -> Dict[str, MultiTracerForwardSystem]:
    """Construct forward systems for all observation panels in a network."""
    systems: Dict[str, MultiTracerForwardSystem] = {}
    for panel in panels:
        node_id = str(panel.node_id)
        if node_id in systems:
            raise ValueError(f"Duplicate observation panel for node {node_id!r}.")
        systems[node_id] = build_forward_system(panel, grid, **kwargs)
    return systems
