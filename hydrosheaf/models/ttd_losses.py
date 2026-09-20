"""Wasserstein and time-frequency loss functions for groundwater TTD inversion.

This module provides specialized loss metrics and regularization objectives:
1. Distribution-space 1D Wasserstein distances (W1, W2^2, unbalanced OT)
2. Time-series Wasserstein transport for signed/shifted tracer series
3. Robust time-domain losses (Huber, weighted RMSE, likelihood NLL)
4. Spectral/frequency-domain errors (spectral power, seasonal phase, coherence)
5. Wavelet-band multiresolution residual diagnostics
6. Composite loss evaluation with complete machine-readable decomposition
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Sequence, Tuple

import numpy as np
from scipy import signal


@dataclass(frozen=True)
class LossConfig:
    """Configuration for composite TTD observation and discrepancy losses.

    Parameters
    ----------
    time_loss : str
        Base time-domain loss: 'huber', 'rmse', 'mae', 'nll'.
    wasserstein_ttd_weight : float
        Weight lambda_w_ttd on normalized TTD Wasserstein distance.
    wasserstein_series_weight : float
        Weight lambda_w_series on time-series transport distance.
    spectral_weight : float
        Weight lambda_freq on seasonal spectral/phase error.
    wavelet_weight : float
        Weight lambda_wave on wavelet-band multiresolution error.
    huber_delta : float
        Threshold parameter delta for Huber loss.
    frequency_bands : tuple[tuple[float, float], ...]
        Target frequency bands [f_min, f_max] for band-pass power errors.
    wavelet_name : str
        Wavelet family name (e.g. 'db4', 'haar').
    wavelet_levels : int
        Number of dyadic decomposition levels.
    missing_data_policy : str
        'ignore_mask', 'drop_nan', or 'error'.
    signed_series_policy : str
        'shift_baseline', 'huber_residual', or 'positive_negative'.
    """

    time_loss: str = "huber"
    wasserstein_ttd_weight: float = 0.0
    wasserstein_series_weight: float = 0.0
    spectral_weight: float = 0.0
    coherence_weight: float = 0.0
    wavelet_weight: float = 0.0
    huber_delta: float = 1.0
    frequency_bands: tuple[tuple[float, float], ...] = ()
    wavelet_name: str = "db4"
    wavelet_levels: int = 3
    missing_data_policy: str = "ignore_mask"
    signed_series_policy: str = "shift_baseline"

    def __post_init__(self) -> None:
        if self.time_loss not in {"huber", "rmse", "mae", "nll"}:
            raise ValueError(f"Invalid time_loss {self.time_loss!r}")
        if self.huber_delta <= 0.0:
            raise ValueError("huber_delta must be positive.")
        if self.missing_data_policy not in {"ignore_mask", "drop_nan", "error"}:
            raise ValueError(f"Invalid missing_data_policy {self.missing_data_policy!r}")
        if self.signed_series_policy not in {"shift_baseline", "huber_residual", "positive_negative"}:
            raise ValueError(f"Invalid signed_series_policy {self.signed_series_policy!r}")
        for name in (
            "wasserstein_ttd_weight",
            "wasserstein_series_weight",
            "spectral_weight",
            "coherence_weight",
            "wavelet_weight",
        ):
            if float(getattr(self, name)) < 0.0:
                raise ValueError(f"{name} cannot be negative")
        if self.wavelet_levels < 1:
            raise ValueError("wavelet_levels must be at least one")
        for band in self.frequency_bands:
            if len(band) != 2 or band[0] < 0.0 or band[1] <= band[0]:
                raise ValueError("frequency_bands must contain increasing non-negative pairs")


# ---------------------------------------------------------------------------
# 1. Distribution-space 1D Wasserstein distances
# ---------------------------------------------------------------------------

def _cumulative_distribution(pdf: np.ndarray, grid: np.ndarray) -> np.ndarray:
    """Compute the right-continuous CDF of point masses on ``grid``.

    TTD arrays in HydroSheaf are probability masses at lag coordinates, not
    samples of a density.  Treating them as a trapezoidal density distorted
    non-uniform grids and made the W2 interpolation ambiguous.
    """
    total = float(np.sum(pdf))
    if total <= 0.0:
        raise ValueError("distribution must have positive mass")
    return np.cumsum(pdf, dtype=float) / total


def _validate_distribution(
    values: Sequence[float] | np.ndarray,
    lag_grid: Sequence[float] | np.ndarray,
    name: str,
) -> tuple[np.ndarray, np.ndarray]:
    arr = np.asarray(values, dtype=float).reshape(-1)
    grid = np.asarray(lag_grid, dtype=float).reshape(-1)
    if len(arr) != len(grid) or len(arr) < 2:
        raise ValueError(f"{name} and lag_grid must have equal length >= 2")
    if not np.all(np.isfinite(arr)) or not np.all(np.isfinite(grid)):
        raise ValueError(f"{name} and lag_grid must be finite")
    if np.any(arr < -1e-9):
        raise ValueError(f"{name} must be non-negative for Wasserstein calculation")
    if np.any(np.diff(grid) <= 0.0):
        raise ValueError("lag_grid must be strictly increasing (strictly monotonically increasing)")
    if float(np.sum(arr)) <= 0.0:
        raise ValueError(f"{name} must have strictly positive mass")
    return np.maximum(arr, 0.0), grid


def _quantile_from_point_masses(mass: np.ndarray, grid: np.ndarray, u: np.ndarray) -> np.ndarray:
    """Inverse CDF using piecewise-linear interpolation between mass nodes."""
    cdf = _cumulative_distribution(mass, grid)
    # Include a zero-probability anchor at the first lag.  Search on the
    # right-continuous CDF gives the physical support point for a point mass.
    out = np.empty_like(u, dtype=float)
    for i, prob in enumerate(u):
        idx = int(np.searchsorted(cdf, prob, side="left"))
        out[i] = grid[min(idx, len(grid) - 1)]
    return out


def wasserstein_ttd_w1(
    predicted: Sequence[float] | np.ndarray,
    observed: Sequence[float] | np.ndarray,
    lag_grid: Sequence[float] | np.ndarray,
) -> float:
    """Compute 1-Wasserstein (Earth Mover's) distance between two 1D TTD distributions.

    For 1D distributions with CDFs F_P and F_Q:
        W_1(P, Q) = int_0^infty |F_P(tau) - F_Q(tau)| d tau

    Parameters
    ----------
    predicted : array_like
        Predicted non-negative distribution values h_pred(tau).
    observed : array_like
        Target/truth non-negative distribution values h_obs(tau).
    lag_grid : array_like
        Strictly increasing lag coordinates.

    Returns
    -------
    float
        Exact 1-Wasserstein distance.
    """
    p, lags = _validate_distribution(predicted, lag_grid, "predicted")
    q, _ = _validate_distribution(observed, lags, "observed")
    cdf_p = _cumulative_distribution(p, lags)
    cdf_q = _cumulative_distribution(q, lags)
    # For point masses, integrate the CDF difference on each half-open interval.
    return float(np.sum(np.abs(cdf_p[:-1] - cdf_q[:-1]) * np.diff(lags)))


def wasserstein_ttd_w2(
    predicted: Sequence[float] | np.ndarray,
    observed: Sequence[float] | np.ndarray,
    lag_grid: Sequence[float] | np.ndarray,
    n_quantiles: int = 100,
) -> float:
    """Compute squared 2-Wasserstein distance W_2^2(P, Q) = int_0^1 (F_P^-1(u) - F_Q^-1(u))^2 du."""
    if n_quantiles < 2:
        raise ValueError("n_quantiles must be at least two")
    p, lags = _validate_distribution(predicted, lag_grid, "predicted")
    q, _ = _validate_distribution(observed, lags, "observed")
    u_grid = np.linspace(0.0, 1.0, n_quantiles)
    inv_p = _quantile_from_point_masses(p, lags, u_grid)
    inv_q = _quantile_from_point_masses(q, lags, u_grid)
    return float(np.mean((inv_p - inv_q) ** 2))


def unbalanced_wasserstein_ttd(
    predicted: Sequence[float] | np.ndarray,
    observed: Sequence[float] | np.ndarray,
    lag_grid: Sequence[float] | np.ndarray,
    mass_penalty: float = 10.0,
) -> float:
    """Unbalanced Wasserstein distance allowing mass difference with linear penalty."""
    p = np.asarray(predicted, dtype=float)
    q = np.asarray(observed, dtype=float)
    lags = np.asarray(lag_grid, dtype=float)

    sp = float(np.sum(p))
    sq = float(np.sum(q))
    mass_diff = abs(sp - sq)

    if sp > 0 and sq > 0:
        w1_norm = wasserstein_ttd_w1(p / sp, q / sq, lags)
    else:
        w1_norm = float(np.ptp(lags))

    # Min transport scale times common mass + penalty on difference
    common_mass = min(sp, sq)
    return float(common_mass * w1_norm + mass_penalty * mass_diff)


# ---------------------------------------------------------------------------
# 2. Time-series transport for signed/shifted series
# ---------------------------------------------------------------------------

def wasserstein_time_series(
    predicted: Sequence[float] | np.ndarray,
    observed: Sequence[float] | np.ndarray,
    time_grid: Sequence[float] | np.ndarray,
    policy: str = "shift_baseline",
) -> float:
    """Compute time-series optimal transport distance for tracer series."""
    y_p = np.asarray(predicted, dtype=float)
    y_o = np.asarray(observed, dtype=float)
    t = np.asarray(time_grid, dtype=float)

    if y_p.shape != y_o.shape or y_p.shape != t.shape:
        raise ValueError("predicted, observed, and time_grid must have the same shape")
    if len(t) < 2 or np.any(~np.isfinite(t)) or np.any(np.diff(t) <= 0.0):
        raise ValueError("time_grid must be finite and strictly increasing")

    finite = np.isfinite(y_p) & np.isfinite(y_o)
    if not np.any(finite):
        return 0.0
    y_p, y_o, t = y_p[finite], y_o[finite], t[finite]

    if policy == "shift_baseline":
        min_val = min(float(np.min(y_p)), float(np.min(y_o)))
        offset = abs(min_val) + 1.0 if min_val <= 0.0 else 0.0
        pos_p = y_p + offset
        pos_o = y_o + offset
        return wasserstein_ttd_w1(pos_p, pos_o, t)
    elif policy == "huber_residual":
        return huber_loss(y_p, y_o)
    elif policy == "positive_negative":
        p_plus, p_minus = np.maximum(y_p, 0.0), np.maximum(-y_p, 0.0)
        o_plus, o_minus = np.maximum(y_o, 0.0), np.maximum(-y_o, 0.0)
        d_plus = wasserstein_ttd_w1(p_plus + 1e-4, o_plus + 1e-4, t)
        d_minus = wasserstein_ttd_w1(p_minus + 1e-4, o_minus + 1e-4, t)
        return float(0.5 * (d_plus + d_minus))
    else:
        raise ValueError(f"Unknown policy {policy!r}")


# ---------------------------------------------------------------------------
# 3. Time-domain observation losses
# ---------------------------------------------------------------------------

def _apply_mask(
    y_pred: np.ndarray,
    y_obs: np.ndarray,
    weights: Optional[np.ndarray] = None,
    mask: Optional[np.ndarray] = None,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    p = np.asarray(y_pred, dtype=float)
    o = np.asarray(y_obs, dtype=float)
    if p.shape != o.shape:
        raise ValueError("predicted and observed must have the same shape")
    w = np.ones_like(p) if weights is None else np.asarray(weights, dtype=float)
    if w.shape != p.shape:
        raise ValueError("weights must have the same shape as predicted")
    if mask is not None and np.asarray(mask).shape != p.shape:
        raise ValueError("mask must have the same shape as predicted")
    if not np.all(np.isfinite(w)) or np.any(w < 0.0):
        raise ValueError("weights must be finite and non-negative")

    valid = np.isfinite(p) & np.isfinite(o) & (w > 0)
    if mask is not None:
        valid = valid & np.asarray(mask, dtype=bool)

    if not np.any(valid):
        return np.array([]), np.array([]), np.array([])
    return p[valid], o[valid], w[valid]


def weighted_rmse(
    predicted: Sequence[float] | np.ndarray,
    observed: Sequence[float] | np.ndarray,
    weights: Optional[Sequence[float] | np.ndarray] = None,
    mask: Optional[Sequence[bool] | np.ndarray] = None,
) -> float:
    p, o, w = _apply_mask(np.asarray(predicted), np.asarray(observed), weights, mask)
    if len(p) == 0:
        return 0.0
    return float(np.sqrt(np.sum(w * (p - o) ** 2) / np.sum(w)))


def weighted_mae(
    predicted: Sequence[float] | np.ndarray,
    observed: Sequence[float] | np.ndarray,
    weights: Optional[Sequence[float] | np.ndarray] = None,
    mask: Optional[Sequence[bool] | np.ndarray] = None,
) -> float:
    p, o, w = _apply_mask(np.asarray(predicted), np.asarray(observed), weights, mask)
    if len(p) == 0:
        return 0.0
    return float(np.sum(w * np.abs(p - o)) / np.sum(w))


def huber_loss(
    predicted: Sequence[float] | np.ndarray,
    observed: Sequence[float] | np.ndarray,
    delta: float = 1.0,
    weights: Optional[Sequence[float] | np.ndarray] = None,
    mask: Optional[Sequence[bool] | np.ndarray] = None,
) -> float:
    """Huber robust loss: 0.5 * e^2 for |e| <= delta, else delta * (|e| - 0.5 * delta)."""
    p, o, w = _apply_mask(np.asarray(predicted), np.asarray(observed), weights, mask)
    if len(p) == 0:
        return 0.0
    errors = np.abs(p - o)
    quadratic = np.minimum(errors, delta)
    linear = errors - quadratic
    loss = 0.5 * quadratic**2 + delta * linear
    return float(np.sum(w * loss) / np.sum(w))


def likelihood_weighted_loss(
    predicted: Sequence[float] | np.ndarray,
    observed: Sequence[float] | np.ndarray,
    sigmas: Optional[Sequence[float] | np.ndarray] = None,
    mask: Optional[Sequence[bool] | np.ndarray] = None,
) -> float:
    """Gaussian negative log-likelihood loss weighted by 1 / sigma^2."""
    p = np.asarray(predicted, dtype=float)
    o = np.asarray(observed, dtype=float)
    s = np.ones_like(p) if sigmas is None else np.maximum(np.asarray(sigmas, dtype=float), 1e-4)
    w = 1.0 / (s**2)
    p_val, o_val, w_val = _apply_mask(p, o, w, mask)
    if len(p_val) == 0:
        return 0.0
    return float(0.5 * np.sum(w_val * (p_val - o_val) ** 2) / len(p_val))


# ---------------------------------------------------------------------------
# 4. Frequency-domain losses
# ---------------------------------------------------------------------------

def spectral_energy_loss(
    predicted: Sequence[float] | np.ndarray,
    observed: Sequence[float] | np.ndarray,
    dt: float = 1.0,
    period: Optional[float] = None,
    n_harmonics: int = 3,
) -> float:
    """Compute spectral power error between predicted and observed signals.

    If period is supplied, focuses on the fundamental seasonal frequency w0 = 2*pi / period
    and its integer harmonics. Otherwise computes normalized total spectral power discrepancy.
    """
    p = np.asarray(predicted, dtype=float)
    o = np.asarray(observed, dtype=float)
    valid = np.isfinite(p) & np.isfinite(o)
    if np.sum(valid) < 8:
        return 0.0

    p_clean = p[valid] - np.mean(p[valid])
    o_clean = o[valid] - np.mean(o[valid])
    n = len(p_clean)

    # FFT power spectrum
    fft_p = np.fft.rfft(p_clean)
    fft_o = np.fft.rfft(o_clean)
    freqs = np.fft.rfftfreq(n, d=dt)

    power_p = np.abs(fft_p) ** 2 / n
    power_o = np.abs(fft_o) ** 2 / n

    if period is not None and period > 0.0:
        f0 = 1.0 / period
        err = 0.0
        for k in range(1, n_harmonics + 1):
            target_f = k * f0
            idx = int(np.argmin(np.abs(freqs - target_f)))
            err += abs(float(power_p[idx] - power_o[idx]))
        denom = float(np.sum(power_o)) + 1e-6
        return float(err / denom)

    # Total relative power error across all frequencies
    return float(np.sum(np.abs(power_p - power_o)) / (np.sum(power_o) + 1e-6))


def phase_error_at_frequency(
    predicted: Sequence[float] | np.ndarray,
    observed: Sequence[float] | np.ndarray,
    target_frequency: float,
    dt: float = 1.0,
) -> float:
    """Compute phase discrepancy |phi_pred - phi_obs| in [0, pi] at target frequency."""
    p = np.asarray(predicted, dtype=float)
    o = np.asarray(observed, dtype=float)
    valid = np.isfinite(p) & np.isfinite(o)
    if np.sum(valid) < 8:
        return 0.0

    p_c = p[valid] - np.mean(p[valid])
    o_c = o[valid] - np.mean(o[valid])
    n = len(p_c)

    freqs = np.fft.rfftfreq(n, d=dt)
    idx = int(np.argmin(np.abs(freqs - target_frequency)))

    fft_p = np.fft.rfft(p_c)[idx]
    fft_o = np.fft.rfft(o_c)[idx]

    angle_p = np.angle(fft_p)
    angle_o = np.angle(fft_o)

    diff = abs(angle_p - angle_o)
    return float(min(diff, 2.0 * np.pi - diff))


def cross_spectral_coherence_error(
    predicted: Sequence[float] | np.ndarray,
    observed: Sequence[float] | np.ndarray,
    dt: float = 1.0,
) -> float:
    """Compute 1 - mean(magnitude_squared_coherence) via Welch's method."""
    p = np.asarray(predicted, dtype=float)
    o = np.asarray(observed, dtype=float)
    valid = np.isfinite(p) & np.isfinite(o)
    if np.sum(valid) < 16:
        return 0.0

    p_c = p[valid]
    o_c = o[valid]
    nperseg = min(len(p_c), 32)
    f, coh = signal.coherence(p_c, o_c, fs=1.0 / dt, nperseg=nperseg)
    if not np.any(np.isfinite(coh)):
        raise ValueError("coherence calculation produced no finite values")
    return float(1.0 - np.nanmean(coh))


# ---------------------------------------------------------------------------
# 5. Wavelet losses
# ---------------------------------------------------------------------------

def wavelet_band_loss(
    predicted: Sequence[float] | np.ndarray,
    observed: Sequence[float] | np.ndarray,
    wavelet: str = "db4",
    levels: int = 3,
) -> float:
    """Multiresolution dyadic sub-band energy error.

    Uses PyWavelets if available; falls back to an exact dyadic Haar filter bank.
    """
    p = np.asarray(predicted, dtype=float)
    o = np.asarray(observed, dtype=float)
    valid = np.isfinite(p) & np.isfinite(o)
    if np.sum(valid) < 8:
        return 0.0

    p_c = p[valid]
    o_c = o[valid]

    try:
        import pywt
        # Decompose both signals
        max_level = pywt.dwt_max_level(len(p_c), wavelet)
        lev = min(levels, max_level)
        if lev < 1:
            return float(np.mean((p_c - o_c) ** 2))

        coeffs_p = pywt.wavedec(p_c, wavelet, level=lev)
        coeffs_o = pywt.wavedec(o_c, wavelet, level=lev)

        total_err = 0.0
        total_energy = 0.0
        for cp, co in zip(coeffs_p, coeffs_o):
            energy_p = np.sum(cp**2)
            energy_o = np.sum(co**2)
            total_err += abs(energy_p - energy_o)
            total_energy += energy_o
        return float(total_err / (total_energy + 1e-6))
    except ImportError:
        # Self-contained dyadic Haar filter bank fallback.  Compare the same
        # sub-band energies as the PyWavelets branch rather than the energy of
        # the residual alone; this keeps the metric a discrepancy between the
        # two signals under either implementation.
        def haar_energies(values: np.ndarray) -> list[float]:
            energies: list[float] = []
            current = values
            for _ in range(min(levels, 3)):
                if len(current) < 2:
                    break
                if len(current) % 2:
                    current = current[:-1]
                approx = (current[0::2] + current[1::2]) / np.sqrt(2.0)
                detail = (current[0::2] - current[1::2]) / np.sqrt(2.0)
                energies.append(float(np.sum(detail**2)))
                current = approx
            energies.append(float(np.sum(current**2)))
            return energies

        ep, eo = haar_energies(p_c), haar_energies(o_c)
        denom = float(np.sum(eo)) + 1e-12
        return float(np.sum(np.abs(np.asarray(ep) - np.asarray(eo))) / denom)


# ---------------------------------------------------------------------------
# 6. Composite Loss Function
# ---------------------------------------------------------------------------

@dataclass
class CompositeLossResult:
    total_loss: float
    decomposition: dict[str, float]


def evaluate_composite_ttd_loss(
    predicted_signal: np.ndarray,
    observed_signal: np.ndarray,
    time_grid: np.ndarray,
    predicted_ttd: Optional[np.ndarray] = None,
    observed_ttd: Optional[np.ndarray] = None,
    lag_grid: Optional[np.ndarray] = None,
    weights: Optional[np.ndarray] = None,
    mask: Optional[np.ndarray] = None,
    season_period: Optional[float] = None,
    config: Optional[LossConfig] = None,
) -> CompositeLossResult:
    """Evaluate compound TTD loss with full machine-readable decomposition.

    L = L_time
        + lambda_w_ttd * L_Wasserstein_TTD
        + lambda_w_series * L_Wasserstein_series
        + lambda_freq * L_frequency
        + lambda_wave * L_wavelet
    """
    cfg = config or LossConfig()
    decomp: dict[str, float] = {}

    # 1. Base time-domain loss
    if cfg.time_loss == "huber":
        l_time = huber_loss(predicted_signal, observed_signal, delta=cfg.huber_delta, weights=weights, mask=mask)
    elif cfg.time_loss == "rmse":
        l_time = weighted_rmse(predicted_signal, observed_signal, weights=weights, mask=mask)
    elif cfg.time_loss == "mae":
        l_time = weighted_mae(predicted_signal, observed_signal, weights=weights, mask=mask)
    elif cfg.time_loss == "nll":
        l_time = likelihood_weighted_loss(predicted_signal, observed_signal, sigmas=weights, mask=mask)
    else:
        l_time = huber_loss(predicted_signal, observed_signal, delta=cfg.huber_delta, weights=weights, mask=mask)
    decomp["time_loss"] = float(l_time)

    # 2. Wasserstein TTD loss
    l_w_ttd = 0.0
    if cfg.wasserstein_ttd_weight > 0.0:
        if predicted_ttd is None or observed_ttd is None or lag_grid is None:
            raise ValueError(
                "wasserstein_ttd_weight is positive but predicted_ttd, observed_ttd, "
                "or lag_grid was not supplied"
            )
        l_w_ttd = wasserstein_ttd_w1(predicted_ttd, observed_ttd, lag_grid)
    decomp["wasserstein_ttd"] = float(l_w_ttd)

    # 3. Wasserstein time-series loss
    l_w_series = 0.0
    if cfg.wasserstein_series_weight > 0.0:
        l_w_series = wasserstein_time_series(
            predicted_signal, observed_signal, time_grid, policy=cfg.signed_series_policy
        )
    decomp["wasserstein_series"] = float(l_w_series)

    # 4. Spectral/frequency loss
    l_freq = 0.0
    if cfg.spectral_weight > 0.0:
        dt = float(time_grid[1] - time_grid[0]) if len(time_grid) > 1 else 1.0
        l_freq = spectral_energy_loss(predicted_signal, observed_signal, dt=dt, period=season_period)
    decomp["spectral_loss"] = float(l_freq)

    l_coherence = 0.0
    if cfg.coherence_weight > 0.0:
        dt = float(time_grid[1] - time_grid[0]) if len(time_grid) > 1 else 1.0
        l_coherence = cross_spectral_coherence_error(predicted_signal, observed_signal, dt=dt)
    decomp["coherence_loss"] = float(l_coherence)

    # 5. Wavelet loss
    l_wave = 0.0
    if cfg.wavelet_weight > 0.0:
        l_wave = wavelet_band_loss(predicted_signal, observed_signal, wavelet=cfg.wavelet_name, levels=cfg.wavelet_levels)
    decomp["wavelet_loss"] = float(l_wave)

    total = (
        l_time
        + cfg.wasserstein_ttd_weight * l_w_ttd
        + cfg.wasserstein_series_weight * l_w_series
        + cfg.spectral_weight * l_freq
        + cfg.coherence_weight * l_coherence
        + cfg.wavelet_weight * l_wave
    )
    decomp["total_loss"] = float(total)

    return CompositeLossResult(total_loss=float(total), decomposition=decomp)
