"""
Vadose-Saturated Zone Coupling Module.

This module provides functions for coupling vadose zone recharge output
with saturated zone transport models.
"""

from dataclasses import dataclass, field
from typing import List, Optional

import numpy as np

from ..config import Config
from .flopy_1d import (
    TransportResult,
    build_1d_transport_model,
    run_1d_transport,
    check_flopy_available,
)


@dataclass
class VadoseCouplingResult:
    """Result of coupled vadose-saturated simulation."""

    # Vadose zone output (input to saturated zone)
    vadose_times: np.ndarray
    vadose_recharge_m_day: np.ndarray
    vadose_concentration: np.ndarray

    # Saturated zone output
    saturated_result: Optional[TransportResult] = None

    # Combined output at observation point
    combined_times: Optional[np.ndarray] = None
    combined_concentration: Optional[np.ndarray] = None

    # Attenuation metrics
    attenuation_factor: float = 1.0
    mean_travel_time_days: float = 0.0
    peak_concentration: float = 0.0
    peak_time_days: float = 0.0

    # Diagnostics
    success: bool = True
    warnings: List[str] = field(default_factory=list)


def couple_vadose_saturated(
    vadose_times: np.ndarray,
    vadose_concentration: np.ndarray,
    vadose_recharge_m_day: Optional[np.ndarray] = None,
    config: Optional[Config] = None,
    use_analytical: bool = True,
    # Aquifer properties (can be overridden from config)
    aquifer_length_m: float = 100.0,
    aquifer_thickness_m: float = 10.0,
    hydraulic_k_m_day: float = 1.0,
    porosity: float = 0.25,
    dispersivity_m: float = 1.0,
    denitrification_k_1_day: float = 0.001,
    head_gradient: float = 0.01,
) -> VadoseCouplingResult:
    """
    Couple vadose zone recharge with saturated zone transport.

    This function takes the output from a vadose zone simulation (e.g.,
    nitrate breakthrough curve) and uses it as input to a 1D saturated
    zone transport model.

    Parameters
    ----------
    vadose_times : np.ndarray
        Time points from vadose simulation (days)
    vadose_concentration : np.ndarray
        Concentration at water table from vadose zone (mmol/L)
    vadose_recharge_m_day : np.ndarray, optional
        Recharge rate from vadose zone (m/day). If None, assumes constant.
    config : Config, optional
        Hydrosheaf configuration for transport parameters
    use_analytical : bool
        If True, use analytical solution (no FloPy required)
    aquifer_length_m : float
        Length of saturated flow path (m)
    aquifer_thickness_m : float
        Aquifer thickness (m)
    hydraulic_k_m_day : float
        Hydraulic conductivity (m/day)
    porosity : float
        Effective porosity
    dispersivity_m : float
        Longitudinal dispersivity (m)
    denitrification_k_1_day : float
        First-order decay rate for denitrification (1/day)
    head_gradient : float
        Hydraulic gradient (dimensionless)

    Returns
    -------
    VadoseCouplingResult
        Coupled simulation results
    """
    warnings_list = []

    # Override with config if provided
    if config is not None:
        aquifer_thickness_m = config.aquifer_thickness_m
        porosity = config.aquifer_porosity
        hydraulic_k_m_day = config.aquifer_hydraulic_k_m_day
        dispersivity_m = config.dispersivity_m
        denitrification_k_1_day = config.denitrification_k_1_day

    # Calculate pore velocity
    velocity_m_day = hydraulic_k_m_day * head_gradient / porosity

    # Default recharge if not provided
    if vadose_recharge_m_day is None:
        vadose_recharge_m_day = np.ones_like(vadose_times) * 0.001  # Default 1 mm/day

    # Ensure arrays are numpy arrays
    vadose_times = np.asarray(vadose_times)
    vadose_concentration = np.asarray(vadose_concentration)
    vadose_recharge_m_day = np.asarray(vadose_recharge_m_day)

    # Estimate mean travel time through saturated zone
    mean_travel_time = aquifer_length_m / velocity_m_day

    # Create output time array (extend to capture full breakthrough)
    max_vadose_time = vadose_times[-1] if len(vadose_times) > 0 else 365.0
    total_time = max_vadose_time + 3 * mean_travel_time
    n_output_times = max(100, int(total_time / 10))
    output_times = np.linspace(0, total_time, n_output_times)

    if use_analytical:
        # Use convolution approach for time-varying input
        # The output concentration is the convolution of input with impulse response
        combined_concentration = _convolve_transport(
            input_times=vadose_times,
            input_concentration=vadose_concentration,
            output_times=output_times,
            distance_m=aquifer_length_m,
            velocity_m_day=velocity_m_day,
            dispersivity_m=dispersivity_m,
            decay_rate_1_day=denitrification_k_1_day,
        )

        saturated_result = TransportResult(
            times=output_times,
            concentrations=combined_concentration,
            velocity_m_day=velocity_m_day,
            dispersivity_m=dispersivity_m,
            decay_rate_1_day=denitrification_k_1_day,
            success=True,
        )

    else:
        # Use FloPy for numerical solution
        if not check_flopy_available():
            warnings_list.append("FloPy not available, falling back to analytical")
            return couple_vadose_saturated(
                vadose_times=vadose_times,
                vadose_concentration=vadose_concentration,
                vadose_recharge_m_day=vadose_recharge_m_day,
                config=config,
                use_analytical=True,
                aquifer_length_m=aquifer_length_m,
                aquifer_thickness_m=aquifer_thickness_m,
                hydraulic_k_m_day=hydraulic_k_m_day,
                porosity=porosity,
                dispersivity_m=dispersivity_m,
                denitrification_k_1_day=denitrification_k_1_day,
                head_gradient=head_gradient,
            )

        # Build and run FloPy model
        # Use mean input concentration as source
        mean_input_conc = np.mean(vadose_concentration)

        try:
            mf, mt, params = build_1d_transport_model(
                aquifer_length_m=aquifer_length_m,
                aquifer_thickness_m=aquifer_thickness_m,
                n_cells=50,
                hydraulic_k_m_day=hydraulic_k_m_day,
                porosity=porosity,
                dispersivity_m=dispersivity_m,
                decay_rate_1_day=denitrification_k_1_day,
                head_upstream_m=10.0,
                head_downstream_m=10.0 - head_gradient * aquifer_length_m,
                n_stress_periods=10,
                perlen_days=total_time / 10,
                source_concentration=mean_input_conc,
            )

            saturated_result = run_1d_transport(mf, mt, params, cleanup=True)

            # Check if FloPy run actually succeeded (executables might be missing)
            if (
                not saturated_result.success
                or len(saturated_result.concentrations) <= 1
            ):
                raise RuntimeError(
                    f"FloPy model failed: {saturated_result.warnings}. "
                    "MODFLOW/MT3DMS executables may not be installed."
                )

            combined_concentration = saturated_result.concentrations
            output_times = saturated_result.times

        except Exception as e:
            warnings_list.append(
                f"FloPy error: {str(e)}. Falling back to analytical solution."
            )
            # Fallback to analytical
            combined_concentration = _convolve_transport(
                input_times=vadose_times,
                input_concentration=vadose_concentration,
                output_times=output_times,
                distance_m=aquifer_length_m,
                velocity_m_day=velocity_m_day,
                dispersivity_m=dispersivity_m,
                decay_rate_1_day=denitrification_k_1_day,
            )
            saturated_result = TransportResult(
                times=output_times,
                concentrations=combined_concentration,
                success=True,
            )

    # Calculate attenuation metrics
    peak_input = np.max(vadose_concentration) if len(vadose_concentration) > 0 else 0.0
    peak_output = (
        np.max(combined_concentration) if len(combined_concentration) > 0 else 0.0
    )
    attenuation_factor = peak_output / peak_input if peak_input > 0 else 0.0

    peak_idx = (
        np.argmax(combined_concentration) if len(combined_concentration) > 0 else 0
    )
    peak_time = output_times[peak_idx] if len(output_times) > peak_idx else 0.0

    return VadoseCouplingResult(
        vadose_times=vadose_times,
        vadose_recharge_m_day=vadose_recharge_m_day,
        vadose_concentration=vadose_concentration,
        saturated_result=saturated_result,
        combined_times=output_times,
        combined_concentration=combined_concentration,
        attenuation_factor=attenuation_factor,
        mean_travel_time_days=mean_travel_time,
        peak_concentration=peak_output,
        peak_time_days=peak_time,
        success=True,
        warnings=warnings_list,
    )


def _convolve_transport(
    input_times: np.ndarray,
    input_concentration: np.ndarray,
    output_times: np.ndarray,
    distance_m: float,
    velocity_m_day: float,
    dispersivity_m: float,
    decay_rate_1_day: float = 0.0,
) -> np.ndarray:
    """
    Convolve input concentration with transport impulse response.

    Uses the analytical 1D advection-dispersion-decay solution as the
    impulse response function.
    """
    from scipy.interpolate import interp1d

    # Dispersion coefficient
    D = dispersivity_m * velocity_m_day

    # Interpolate input to finer grid for convolution
    dt = min(1.0, (output_times[-1] - output_times[0]) / 1000)
    fine_times = np.arange(0, output_times[-1] + dt, dt)

    # Interpolate input concentration
    if len(input_times) > 1:
        interp_func = interp1d(
            input_times,
            input_concentration,
            kind="linear",
            fill_value=0.0,
            bounds_error=False,
        )
        fine_input = interp_func(fine_times)
    else:
        fine_input = np.zeros_like(fine_times)
        if len(input_times) == 1:
            idx = np.argmin(np.abs(fine_times - input_times[0]))
            fine_input[idx] = input_concentration[0]

    # Compute impulse response (Green's function)
    impulse_response = np.zeros_like(fine_times)
    x = distance_m

    for i, t in enumerate(fine_times):
        if t > 0:
            # Gaussian-like impulse response
            # beta = np.sqrt(velocity_m_day**2 + 4 * decay_rate_1_day * D)
            sigma = np.sqrt(2 * D * t)
            if sigma > 0:
                # Peak arrival
                # t_peak = x / velocity_m_day
                # Impulse response approximation
                exp_decay = np.exp(-decay_rate_1_day * t)
                gauss = np.exp(
                    -((x - velocity_m_day * t) ** 2) / (4 * D * t)
                ) / np.sqrt(4 * np.pi * D * t)
                impulse_response[i] = gauss * exp_decay * velocity_m_day

    # Normalize impulse response
    if np.sum(impulse_response) > 0:
        impulse_response = impulse_response / (np.sum(impulse_response) * dt)

    # Convolve
    convolved = (
        np.convolve(fine_input, impulse_response, mode="full")[: len(fine_times)] * dt
    )

    # Interpolate back to output times
    if len(fine_times) > 1:
        output_func = interp1d(
            fine_times, convolved, kind="linear", fill_value=0.0, bounds_error=False
        )
        result = output_func(output_times)
    else:
        result = np.zeros_like(output_times)

    return np.maximum(result, 0.0)


def estimate_saturated_travel_time(
    distance_m: float,
    hydraulic_k_m_day: float,
    head_gradient: float,
    porosity: float,
) -> float:
    """
    Estimate mean travel time through saturated zone.

    Parameters
    ----------
    distance_m : float
        Flow path length (m)
    hydraulic_k_m_day : float
        Hydraulic conductivity (m/day)
    head_gradient : float
        Hydraulic gradient (dimensionless)
    porosity : float
        Effective porosity

    Returns
    -------
    float
        Mean travel time (days)
    """
    velocity = hydraulic_k_m_day * head_gradient / porosity
    if velocity > 0:
        return distance_m / velocity
    return float("inf")


def estimate_denitrification_attenuation(
    travel_time_days: float,
    decay_rate_1_day: float,
) -> float:
    """
    Estimate nitrate attenuation due to denitrification.

    Parameters
    ----------
    travel_time_days : float
        Travel time through reactive zone (days)
    decay_rate_1_day : float
        First-order decay rate (1/day)

    Returns
    -------
    float
        Attenuation factor (0-1, where 1 = no attenuation)
    """
    return np.exp(-decay_rate_1_day * travel_time_days)
