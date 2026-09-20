"""Hydrosheaf models package."""

from .ttd_losses import (
    LossConfig,
    wasserstein_ttd_w1,
    wasserstein_ttd_w2,
    unbalanced_wasserstein_ttd,
    wasserstein_time_series,
    weighted_rmse,
    weighted_mae,
    huber_loss,
    likelihood_weighted_loss,
    spectral_energy_loss,
    phase_error_at_frequency,
    cross_spectral_coherence_error,
    wavelet_band_loss,
    evaluate_composite_ttd_loss,
)

__all__ = [
    "LossConfig",
    "wasserstein_ttd_w1",
    "wasserstein_ttd_w2",
    "unbalanced_wasserstein_ttd",
    "wasserstein_time_series",
    "weighted_rmse",
    "weighted_mae",
    "huber_loss",
    "likelihood_weighted_loss",
    "spectral_energy_loss",
    "phase_error_at_frequency",
    "cross_spectral_coherence_error",
    "wavelet_band_loss",
    "evaluate_composite_ttd_loss",
]
