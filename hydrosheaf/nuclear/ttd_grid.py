"""Multi-scale age discretization, quadrature integration, and curvature operators.

This module provides the discrete mathematical infrastructure for non-parametric,
shape-free transit-time distribution (TTD) inversion across groundwater networks.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np


@dataclass(frozen=True)
class TTDGrid:
    """Discretized groundwater transit-time grid with quadrature integration weights.

    The grid represents discrete age coordinates tau_k for k = 1, ..., K.
    Inferential distributions can be represented either as:
    1. A probability mass vector g in Delta^(K-1) where sum_k g_k = 1.0; or
    2. A probability density function f(tau) where sum_k f(tau_k) * w_k = 1.0.

    Parameters
    ----------
    taus : np.ndarray
        Strictly monotonically increasing age coordinates in years, shape (K,).
    weights : np.ndarray
        Trapezoidal quadrature integration weights, shape (K,).
    young_cutoff_years : float
        Age threshold defining young/modern water (Anthropocene), default 70.0 yr.
    holocene_cutoff_years : float
        Age threshold defining Holocene water, default 11,700.0 yr.
    """

    taus: np.ndarray
    weights: np.ndarray
    young_cutoff_years: float = 70.0
    holocene_cutoff_years: float = 11700.0

    def __post_init__(self) -> None:
        t_arr = np.asarray(self.taus, dtype=float).copy()
        w_arr = np.asarray(self.weights, dtype=float).copy()

        if t_arr.ndim != 1 or t_arr.size < 3:
            raise ValueError("taus must be a 1D array with at least 3 points.")
        if not np.all(np.isfinite(t_arr)):
            raise ValueError("taus must contain only finite numbers.")
        if float(t_arr[0]) < 0.0:
            raise ValueError("Age coordinates must be non-negative.")
        if np.any(np.diff(t_arr) <= 0.0):
            raise ValueError("taus must be strictly monotonically increasing.")

        if w_arr.shape != t_arr.shape:
            raise ValueError("weights must match taus shape.")
        if not np.all(np.isfinite(w_arr)) or np.any(w_arr <= 0.0):
            raise ValueError("Quadrature weights must be strictly positive.")

        if self.young_cutoff_years <= 0.0:
            raise ValueError("young_cutoff_years must be positive.")
        if self.holocene_cutoff_years <= self.young_cutoff_years:
            raise ValueError("holocene_cutoff_years must exceed young_cutoff_years.")

        t_arr.setflags(write=False)
        w_arr.setflags(write=False)
        object.__setattr__(self, "taus", t_arr)
        object.__setattr__(self, "weights", w_arr)

    @property
    def n_bins(self) -> int:
        """Number of discrete age bins K."""
        return int(self.taus.size)

    @property
    def max_age(self) -> float:
        """Maximum age represented in the grid (years)."""
        return float(self.taus[-1])

    @property
    def min_age(self) -> float:
        """Minimum age represented in the grid (years)."""
        return float(self.taus[0])

    def mean_transit_time(self, mass_or_pdf: np.ndarray, *, is_pdf: bool = False) -> float:
        """Calculate the expected transit time E[tau] in years."""
        arr = np.asarray(mass_or_pdf, dtype=float)
        if arr.shape != self.taus.shape:
            raise ValueError(f"Array shape {arr.shape} does not match grid size {self.n_bins}.")
        if is_pdf:
            total_mass = float(np.sum(arr * self.weights))
            if total_mass <= 0.0:
                return float("nan")
            return float(np.sum(self.taus * arr * self.weights) / total_mass)
        else:
            total_mass = float(np.sum(arr))
            if total_mass <= 0.0:
                return float("nan")
            return float(np.sum(self.taus * arr) / total_mass)

    def young_water_fraction(self, mass_or_pdf: np.ndarray, *, is_pdf: bool = False) -> float:
        """Fraction of water younger than young_cutoff_years."""
        arr = np.asarray(mass_or_pdf, dtype=float)
        mask = self.taus <= self.young_cutoff_years
        if is_pdf:
            total_mass = float(np.sum(arr * self.weights))
            if total_mass <= 0.0:
                return 0.0
            return float(np.sum((arr * self.weights)[mask]) / total_mass)
        else:
            total_mass = float(np.sum(arr))
            if total_mass <= 0.0:
                return 0.0
            return float(np.sum(arr[mask]) / total_mass)

    def holocene_fraction(self, mass_or_pdf: np.ndarray, *, is_pdf: bool = False) -> float:
        """Fraction of water between young_cutoff_years and holocene_cutoff_years."""
        arr = np.asarray(mass_or_pdf, dtype=float)
        mask = (self.taus > self.young_cutoff_years) & (self.taus <= self.holocene_cutoff_years)
        if is_pdf:
            total_mass = float(np.sum(arr * self.weights))
            if total_mass <= 0.0:
                return 0.0
            return float(np.sum((arr * self.weights)[mask]) / total_mass)
        else:
            total_mass = float(np.sum(arr))
            if total_mass <= 0.0:
                return 0.0
            return float(np.sum(arr[mask]) / total_mass)

    def pleistocene_fraction(self, mass_or_pdf: np.ndarray, *, is_pdf: bool = False) -> float:
        """Fraction of paleowater older than holocene_cutoff_years."""
        arr = np.asarray(mass_or_pdf, dtype=float)
        mask = self.taus > self.holocene_cutoff_years
        if is_pdf:
            total_mass = float(np.sum(arr * self.weights))
            if total_mass <= 0.0:
                return 0.0
            return float(np.sum((arr * self.weights)[mask]) / total_mass)
        else:
            total_mass = float(np.sum(arr))
            if total_mass <= 0.0:
                return 0.0
            return float(np.sum(arr[mask]) / total_mass)

    def age_fractions(self, mass_or_pdf: np.ndarray, *, is_pdf: bool = False) -> dict[str, float]:
        """Return Anthropocene, Holocene, and Pleistocene fraction breakdown."""
        return {
            "anthropocene": self.young_water_fraction(mass_or_pdf, is_pdf=is_pdf),
            "holocene": self.holocene_fraction(mass_or_pdf, is_pdf=is_pdf),
            "pleistocene": self.pleistocene_fraction(mass_or_pdf, is_pdf=is_pdf),
        }

    def quantiles(
        self, mass_or_pdf: np.ndarray, qs: Sequence[float] = (0.1, 0.25, 0.5, 0.75, 0.9), *, is_pdf: bool = False
    ) -> np.ndarray:
        """Calculate age quantiles from the discrete distribution."""
        arr = np.asarray(mass_or_pdf, dtype=float)
        if is_pdf:
            masses = arr * self.weights
        else:
            masses = arr.copy()
        total = float(np.sum(masses))
        if total <= 0.0:
            return np.full(len(qs), np.nan)
        cdf = np.cumsum(masses) / total
        return np.interp(np.asarray(qs, dtype=float), cdf, self.taus)


def compute_trapezoidal_weights(taus: np.ndarray) -> np.ndarray:
    """Compute non-uniform trapezoidal quadrature integration weights for an age grid."""
    n = len(taus)
    if n < 2:
        raise ValueError("At least two points required to compute quadrature weights.")
    weights = np.empty(n, dtype=float)
    weights[0] = 0.5 * (taus[1] - taus[0])
    weights[-1] = 0.5 * (taus[-1] - taus[-2])
    weights[1:-1] = 0.5 * (taus[2:] - taus[:-2])
    return weights


def build_multiscale_ttd_grid(
    *,
    max_age_years: float = 50000.0,
    dt_young: float = 0.5,
    dt_holocene: float = 50.0,
    dt_pleistocene: float = 500.0,
    young_cutoff_years: float = 70.0,
    holocene_cutoff_years: float = 11700.0,
) -> TTDGrid:
    """Construct a three-era multi-scale discrete age grid.

    Parameters
    ----------
    max_age_years : float
        Maximum transit time represented in the grid (e.g. 50,000 years).
    dt_young : float
        Step size in modern Anthropocene era [0, young_cutoff] (sub-annual, e.g. 0.5 yr).
    dt_holocene : float
        Step size in Holocene era [young_cutoff, holocene_cutoff] (decadal, e.g. 50 yr).
    dt_pleistocene : float
        Step size in Pleistocene era [holocene_cutoff, max_age] (century, e.g. 500 yr).
    young_cutoff_years : float
        Boundary of modern tracer sensitivity (default 70 yr).
    holocene_cutoff_years : float
        Boundary of Holocene epoch (default 11,700 yr).

    Returns
    -------
    TTDGrid
        The configured immutable grid and weights.
    """
    if max_age_years <= holocene_cutoff_years:
        raise ValueError(f"max_age_years ({max_age_years}) must exceed holocene_cutoff_years ({holocene_cutoff_years}).")
    if dt_young <= 0.0 or dt_holocene <= 0.0 or dt_pleistocene <= 0.0:
        raise ValueError("Grid step sizes must be strictly positive.")

    # Era 1: Modern window (captures seasonal 3H, SF6, CFC transients)
    era1 = np.arange(0.0, young_cutoff_years + 0.5 * dt_young, dt_young)

    # Era 2: Holocene window (captures 39Ar and modern 14C)
    era2 = np.arange(young_cutoff_years, holocene_cutoff_years + 0.5 * dt_holocene, dt_holocene)

    # Era 3: Pleistocene window (captures deep 14C and 4He accumulation)
    era3 = np.arange(holocene_cutoff_years, max_age_years + 0.5 * dt_pleistocene, dt_pleistocene)

    taus = np.unique(np.round(np.concatenate([era1, era2, era3, [max_age_years]]), 6))
    taus = taus[taus <= max_age_years]
    weights = compute_trapezoidal_weights(taus)

    return TTDGrid(
        taus=taus,
        weights=weights,
        young_cutoff_years=young_cutoff_years,
        holocene_cutoff_years=holocene_cutoff_years,
    )


def build_uniform_ttd_grid(
    *,
    max_age_years: float = 100.0,
    dt_years: float = 0.5,
    young_cutoff_years: float = 70.0,
) -> TTDGrid:
    """Construct a uniform discrete age grid (useful for shallow, young systems)."""
    if max_age_years <= 0.0 or dt_years <= 0.0:
        raise ValueError("max_age_years and dt_years must be positive.")
    taus = np.arange(0.0, max_age_years + 0.5 * dt_years, dt_years, dtype=float)
    weights = compute_trapezoidal_weights(taus)
    holocene = max(young_cutoff_years * 2.0, max_age_years + 10.0)
    return TTDGrid(
        taus=taus,
        weights=weights,
        young_cutoff_years=young_cutoff_years,
        holocene_cutoff_years=holocene,
    )


def build_d1_difference_matrix(grid: TTDGrid) -> np.ndarray:
    """Construct non-uniform 1st-order finite difference operator D_1 in R^((K-1) x K)."""
    taus = grid.taus
    n = len(taus)
    dmat = np.zeros((n - 1, n), dtype=float)
    dt = np.diff(taus)
    for j in range(n - 1):
        inv_dt = 1.0 / dt[j]
        dmat[j, j] = -inv_dt
        dmat[j, j + 1] = inv_dt
    return dmat


def build_d2_curvature_matrix(grid: TTDGrid) -> np.ndarray:
    """Construct non-uniform 2nd-order finite difference curvature operator D_2 in R^((K-2) x K).

    Approximates (d^2 g / d tau^2) on non-uniform discretization grids:
    (D_2 g)_j = (2 / (dt_1 + dt_2)) * [(g_{j+1} - g_j)/dt_2 - (g_j - g_{j-1})/dt_1]
    """
    taus = grid.taus
    n = len(taus)
    if n < 3:
        return np.zeros((0, n), dtype=float)

    rows: list[np.ndarray] = []
    for j in range(1, n - 1):
        dt1 = float(taus[j] - taus[j - 1])
        dt2 = float(taus[j + 1] - taus[j])
        if dt1 <= 0.0 or dt2 <= 0.0:
            raise ValueError("Lags must be strictly monotonically increasing.")
        row = np.zeros(n, dtype=float)
        coeff = 2.0 / (dt1 + dt2)
        row[j - 1] = coeff / dt1
        row[j] = -coeff * (1.0 / dt1 + 1.0 / dt2)
        row[j + 1] = coeff / dt2
        rows.append(row)

    return np.vstack(rows)


def build_mass_aware_curvature_matrix(grid: TTDGrid) -> np.ndarray:
    """Construct a curvature penalty for a probability-mass vector.

    ``TTDGrid`` uses two related representations: a probability mass vector
    ``g`` with ``sum(g) == 1`` and a density ``f`` with
    ``sum(f * weights) == 1``.  The finite-difference operator returned by
    :func:`build_d2_curvature_matrix` acts on values of a function sampled at
    the age coordinates, so it is a density operator.  On a non-uniform grid,
    applying it directly to ``g`` makes the regularizer depend on the local
    bin width.  This helper transforms mass to density before applying the
    curvature operator and weights the rows by the corresponding integration
    measure.

    The returned matrix ``L`` can be used directly in a quadratic penalty:

    ``||L @ g||_2**2``.
    """
    d2_density = build_d2_curvature_matrix(grid)
    if d2_density.shape[0] == 0:
        return d2_density

    # Interior rows are located at tau[1:-1].  Weighting by the square root of
    # the local quadrature measure approximates the integral of curvature^2.
    density_from_mass = np.diag(1.0 / grid.weights)
    row_measure = np.sqrt(grid.weights[1:-1])
    return row_measure[:, None] * (d2_density @ density_from_mass)


def shannon_entropy(mass: np.ndarray, *, eps: float = 1e-12) -> float:
    """Compute discrete Shannon entropy -sum_k p_k ln(p_k) of a mass vector."""
    p = np.asarray(mass, dtype=float)
    total = float(np.sum(p))
    if total <= 0.0:
        return 0.0
    p_norm = np.maximum(p / total, eps)
    p_norm = p_norm / np.sum(p_norm)
    return float(-np.sum(p_norm * np.log(p_norm)))


def wasserstein_1d(
    mass1: np.ndarray,
    mass2: np.ndarray,
    grid: TTDGrid,
) -> float:
    """Compute 1D Wasserstein-1 (Earth Mover's Distance) between two age distributions.

    For 1D distributions on ordered grid tau_k:
    W_1(g_1, g_2) = integral_0^infty |F_1(tau) - F_2(tau)| d tau
                  ~= sum_{k=0}^{K-2} |CDF_1(tau_k) - CDF_2(tau_k)|
                     * (tau_{k+1} - tau_k)
    """
    m1 = np.asarray(mass1, dtype=float)
    m2 = np.asarray(mass2, dtype=float)
    if m1.shape != (grid.n_bins,) or m2.shape != (grid.n_bins,):
        raise ValueError("Mass vectors must have the same length as grid.taus.")
    if (
        not np.all(np.isfinite(m1))
        or not np.all(np.isfinite(m2))
        or np.any(m1 < 0.0)
        or np.any(m2 < 0.0)
    ):
        raise ValueError("Mass vectors must be finite and non-negative.")
    tot1 = float(np.sum(m1))
    tot2 = float(np.sum(m2))
    if tot1 <= 0.0 or tot2 <= 0.0:
        return float("nan")
    cdf1 = np.cumsum(m1) / tot1
    cdf2 = np.cumsum(m2) / tot2

    # For discrete masses located at the grid coordinates, the CDF gap is
    # constant on each interval [tau_i, tau_{i+1}).  Integrating over those
    # intervals avoids the half-weight endpoint bias of applying trapezoidal
    # node weights directly to CDF values.
    gap = np.abs(cdf1 - cdf2)
    return float(np.sum(gap[:-1] * np.diff(grid.taus)))
