"""Tests for Wasserstein, time-domain, frequency, and wavelet TTD loss functions."""

from __future__ import annotations

import numpy as np
import pytest

from hydrosheaf.models.ttd_losses import (
    LossConfig,
    evaluate_composite_ttd_loss,
    huber_loss,
    phase_error_at_frequency,
    spectral_energy_loss,
    unbalanced_wasserstein_ttd,
    wasserstein_time_series,
    wasserstein_ttd_w1,
    wasserstein_ttd_w2,
    wavelet_band_loss,
    weighted_mae,
    weighted_rmse,
)


def test_wasserstein_w1_identity_and_translation():
    lags = np.linspace(0.0, 50.0, 51)
    p = np.zeros(51)
    p[10] = 1.0  # Dirac at lag 10

    # W1(p, p) = 0
    assert pytest.approx(wasserstein_ttd_w1(p, p, lags), abs=1e-7) == 0.0

    # Translate by 5 steps (from lag 10 to lag 15)
    q = np.zeros(51)
    q[15] = 1.0  # Dirac at lag 15
    w1_dist = wasserstein_ttd_w1(p, q, lags)
    # Expected distance is exactly 5.0
    assert pytest.approx(w1_dist, abs=1e-5) == 5.0

    # Non-negativity
    assert wasserstein_ttd_w1(p, q, lags) >= 0.0


def test_wasserstein_w2_translation():
    lags = np.linspace(0.0, 50.0, 51)
    p = np.zeros(51)
    p[10] = 1.0
    q = np.zeros(51)
    q[14] = 1.0

    # For pure translation by delta=4, W2^2 = delta^2 = 16.0
    w2_dist = wasserstein_ttd_w2(p, q, lags)
    assert pytest.approx(w2_dist, rel=0.10) == 16.0


def test_wasserstein_w2_rejects_invalid_inputs_and_composite_does_not_hide_missing_ttd():
    with pytest.raises(ValueError, match="strictly increasing"):
        wasserstein_ttd_w2([1.0, 0.0], [0.5, 0.5], [0.0, 0.0])
    with pytest.raises(ValueError, match="predicted_ttd"):
        evaluate_composite_ttd_loss(
            predicted_signal=np.arange(10.0),
            observed_signal=np.arange(10.0),
            time_grid=np.arange(10.0),
            config=LossConfig(wasserstein_ttd_weight=1.0),
        )


def test_unbalanced_wasserstein():
    lags = np.linspace(0.0, 20.0, 21)
    p = np.zeros(21)
    p[5] = 1.0
    q = np.zeros(21)
    q[5] = 0.8  # 20% mass missing

    uw = unbalanced_wasserstein_ttd(p, q, lags, mass_penalty=10.0)
    # Transport cost between normalized is 0, penalty on 0.2 missing mass is 0.2 * 10 = 2.0
    assert pytest.approx(uw, abs=1e-5) == 2.0


def test_time_domain_losses():
    y_pred = np.array([1.0, 2.0, 3.0, 4.0])
    y_obs = np.array([1.0, 2.5, 3.0, 5.0])
    # Errors: [0, 0.5, 0, 1.0]

    # RMSE = sqrt(mean([0, 0.25, 0, 1.0])) = sqrt(1.25 / 4) = sqrt(0.3125) = 0.5590
    assert pytest.approx(weighted_rmse(y_pred, y_obs), abs=1e-4) == 0.5590
    # MAE = mean([0, 0.5, 0, 1.0]) = 1.5 / 4 = 0.375
    assert pytest.approx(weighted_mae(y_pred, y_obs), abs=1e-4) == 0.375

    # Huber loss with delta=1.0: for e <= 1.0, 0.5*e^2 -> 0.5*(0.25 + 1.0)/4 = 0.15625
    assert pytest.approx(huber_loss(y_pred, y_obs, delta=1.0), abs=1e-4) == 0.15625

    # Missing value handling
    y_obs_nan = np.array([1.0, 2.5, np.nan, 5.0])
    rmse_nan = weighted_rmse(y_pred, y_obs_nan)
    assert np.isfinite(rmse_nan)


def test_frequency_domain_losses():
    t = np.arange(100, dtype=float)
    # Fundamental frequency f0 = 1/20 = 0.05
    sig1 = np.sin(2 * np.pi * t / 20.0)
    sig2 = np.sin(2 * np.pi * t / 20.0)

    # Identical signals have zero spectral energy error
    err_same = spectral_energy_loss(sig1, sig2, period=20.0)
    assert pytest.approx(err_same, abs=1e-5) == 0.0

    # Frequency shifted signal
    sig_shifted = np.sin(2 * np.pi * t / 15.0)
    err_diff = spectral_energy_loss(sig1, sig_shifted, period=20.0)
    assert err_diff > 0.1

    # Phase error: pi/2 shift gives phase error ~pi/2
    sig_cos = np.cos(2 * np.pi * t / 20.0)
    phase_err = phase_error_at_frequency(sig1, sig_cos, target_frequency=0.05)
    assert pytest.approx(phase_err, abs=0.2) == np.pi / 2.0


def test_wavelet_band_loss():
    t = np.arange(64, dtype=float)
    sig = np.sin(2 * np.pi * t / 16.0)
    assert pytest.approx(wavelet_band_loss(sig, sig), abs=1e-5) == 0.0

    noisy = sig + 0.5 * np.random.randn(64)
    loss = wavelet_band_loss(sig, noisy)
    assert loss > 0.0


def test_composite_loss_decomposition():
    t = np.arange(50, dtype=float)
    y_pred = np.sin(t)
    y_obs = np.sin(t) + 0.1

    lags = np.linspace(0, 10, 11)
    p_ttd = np.ones(11) / 11
    q_ttd = np.ones(11) / 11

    cfg = LossConfig(
        time_loss="huber",
        wasserstein_ttd_weight=1.0,
        spectral_weight=0.5,
        wavelet_weight=0.2,
    )

    res = evaluate_composite_ttd_loss(
        predicted_signal=y_pred,
        observed_signal=y_obs,
        time_grid=t,
        predicted_ttd=p_ttd,
        observed_ttd=q_ttd,
        lag_grid=lags,
        config=cfg,
    )

    assert "time_loss" in res.decomposition
    assert "wasserstein_ttd" in res.decomposition
    assert "spectral_loss" in res.decomposition
    assert "wavelet_loss" in res.decomposition
    assert "total_loss" in res.decomposition
    assert res.total_loss > 0.0
