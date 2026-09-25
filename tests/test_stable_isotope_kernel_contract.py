"""Focused contracts for stable-isotope capabilities and scalar TTD kernels."""

import numpy as np
import pytest

from hydrosheaf.nuclear.joint_lpm import tracer_response_kernel
from hydrosheaf.nuclear.tracer_registry import build_default_tracer_registry


@pytest.mark.parametrize("tracer", ["d18O", "18O", "d2H", "2H", "δ18O", "δ2H"])
def test_scalar_age_grid_kernel_rejects_stable_isotopes(tracer: str) -> None:
    with pytest.raises(ValueError, match="recharge-history-aware"):
        tracer_response_kernel(tracer, [0.0, 1.0, 5.0], 2024.0)


def test_registry_distinguishes_time_history_and_scalar_kernel_capabilities() -> None:
    registry = build_default_tracer_registry()

    assert registry["d18O"].response_family == "stable_isotope_time_history"
    assert registry["d2H"].response_family == "stable_isotope_time_history"
    assert not registry["d18O"].supports_scalar_age_grid_kernel
    assert not registry["d2H"].supports_scalar_age_grid_kernel

    assert registry["3H"].response_family == "radioactive_decay"
    assert registry["3H"].supports_scalar_age_grid_kernel
    assert registry["SF6"].response_family == "gas_input_history"
    assert registry["SF6"].supports_scalar_age_grid_kernel


def test_existing_c14_scalar_kernel_remains_supported() -> None:
    ages = np.array([0.0, 10.0, 100.0], dtype=float)

    kernel = tracer_response_kernel("14C", ages, 2024.0)

    assert kernel.shape == ages.shape
    assert np.all(np.isfinite(kernel))
    assert kernel[0] == pytest.approx(100.0)
    assert kernel[0] > kernel[1] > kernel[2]
