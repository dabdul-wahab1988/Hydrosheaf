"""
Transport validation and calibration helpers.
"""

from typing import Dict, List

import numpy as np
import pandas as pd
from scipy import signal

from ..calibration.problem import CalibrationProblem
from ..calibration.definitions import AdjustableParameter, Observation
from .analytical import ade_1d_green


class TransportValidationProblem(CalibrationProblem):
    """
    Calibrates 1D ADE parameters (v, alpha, k, lag, tau, dilution, baseflow)
    to match observed nitrate concentrations.
    """

    def __init__(
        self,
        t_grid: np.ndarray,
        c_in: np.ndarray,
        obs_df: pd.DataFrame,
        distance_x: float,
        use_log: bool = True,
        fit_lag: bool = True,
        use_second_stage: bool = False,
        use_seasonal: bool = False,
        use_hybrid: bool = False,
    ):
        self.t_grid = t_grid
        self.c_in = c_in
        self.distance_x = distance_x
        self.obs_df = obs_df
        self.use_log = use_log
        self.fit_lag = fit_lag
        self.use_second_stage = use_second_stage
        self.use_seasonal = use_seasonal
        self.use_hybrid = use_hybrid

        self.default_params = {
            "v": 5.0,
            "alpha": 10.0,
            "k": 0.001,
            "lag": 0.0,
            "tau": 10.0,
            "dilution": 0.5,
            "baseflow": 10.0,
            "v2": 5.0,
            "alpha2": 10.0,
            "k2": 0.001,
            "dilution_amp": 0.0,
            "dilution_phase": 0.0,
            "baseflow_amp": 0.0,
            "baseflow_phase": 0.0,
        }

        self.parameters = [
            AdjustableParameter("v", self.default_params["v"], 0.1, 100.0),
            AdjustableParameter(
                "alpha", self.default_params["alpha"], 0.1, 500.0, log_transform=True
            ),
            AdjustableParameter("k", self.default_params["k"], 0.0, 0.5),
            AdjustableParameter(
                "tau", self.default_params["tau"], 0.1, 200.0, log_transform=True
            ),
            AdjustableParameter("dilution", self.default_params["dilution"], 0.01, 1.5),
            AdjustableParameter("baseflow", self.default_params["baseflow"], 0.0, 100.0),
        ]

        if self.fit_lag:
            self.parameters.append(
                AdjustableParameter("lag", self.default_params["lag"], 0.0, 200.0)
            )

        if self.use_second_stage:
            self.parameters.extend([
                AdjustableParameter("v2", self.default_params["v2"], 0.1, 100.0),
                AdjustableParameter(
                    "alpha2",
                    self.default_params["alpha2"],
                    0.1,
                    500.0,
                    log_transform=True,
                ),
                AdjustableParameter("k2", self.default_params["k2"], 0.0, 0.5),
            ])

        if self.use_seasonal:
            self.parameters.extend([
                AdjustableParameter(
                    "dilution_amp", self.default_params["dilution_amp"], 0.0, 1.0
                ),
                AdjustableParameter(
                    "dilution_phase", self.default_params["dilution_phase"], 0.0, 365.0
                ),
                AdjustableParameter(
                    "baseflow_amp", self.default_params["baseflow_amp"], 0.0, 1.0
                ),
                AdjustableParameter(
                    "baseflow_phase", self.default_params["baseflow_phase"], 0.0, 365.0
                ),
            ])

        self.param_limits = {p.name: (p.lower_bound, p.upper_bound) for p in self.parameters}

        self.observations: List[Observation] = []
        self.observation_values_raw: Dict[str, float] = {}
        linear_weight = 0.5 if self.use_hybrid else 1.0
        log_weight = 0.5
        for _, row in obs_df.iterrows():
            obs_name = f"obs_{int(row['days'])}"
            raw_value = float(row["NO3_mg_L"])
            self.observation_values_raw[obs_name] = raw_value
            if self.use_hybrid:
                self.observations.append(Observation(obs_name, raw_value, weight=linear_weight))
                self.observations.append(
                    Observation(f"log_{obs_name}", self._log_transform(raw_value), weight=log_weight)
                )
            else:
                self.observations.append(Observation(obs_name, self._transform(raw_value)))

    def _transform(self, value: float) -> float:
        if self.use_log:
            return self._log_transform(value)
        return float(value)

    @staticmethod
    def _log_transform(value: float) -> float:
        return float(np.log1p(max(value, 0.0)))

    def _sanitize_params(self, param_values: Dict[str, float]) -> Dict[str, float]:
        safe = dict(self.default_params)
        for name, value in param_values.items():
            if value is None:
                continue
            try:
                value_float = float(value)
            except (TypeError, ValueError):
                continue
            if not np.isfinite(value_float):
                continue
            safe[name] = value_float

        for name, (low, high) in self.param_limits.items():
            value = safe.get(name, self.default_params.get(name, low))
            if not np.isfinite(value):
                value = self.default_params.get(name, low)
            if value < low:
                value = low
            if value > high:
                value = high
            safe[name] = value

        return safe

    def _apply_lag(self, series: np.ndarray, lag: float) -> np.ndarray:
        return np.interp(self.t_grid - lag, self.t_grid, series, left=0.0, right=0.0)

    def _apply_reservoir(self, series: np.ndarray, tau: float) -> np.ndarray:
        if tau <= 0:
            return series
        dt = self.t_grid[1] - self.t_grid[0] if len(self.t_grid) > 1 else 1.0
        alpha = dt / max(tau, 1e-9)
        smoothed = np.zeros_like(series)
        for i in range(len(series)):
            if i == 0:
                smoothed[i] = series[i]
            else:
                smoothed[i] = smoothed[i - 1] + alpha * (series[i] - smoothed[i - 1])
        return smoothed

    def _simulate(self, param_values: Dict[str, float]) -> np.ndarray:
        safe = self._sanitize_params(param_values)
        v = safe["v"]
        alpha = safe["alpha"]
        k = safe["k"]
        lag = safe["lag"]
        tau = safe["tau"]
        dilution = safe["dilution"]
        baseflow = safe["baseflow"]

        D = max(alpha * v, 1e-9)
        c_in_lagged = self._apply_lag(self.c_in, lag)
        g = ade_1d_green(self.t_grid, self.distance_x, v, D, k)
        c_out_full = signal.fftconvolve(c_in_lagged, g, mode="full")
        stage_output = c_out_full[:len(self.t_grid)]

        if self.use_second_stage:
            v2 = safe["v2"]
            alpha2 = safe["alpha2"]
            k2 = safe["k2"]
            D2 = alpha2 * v2
            g2 = ade_1d_green(self.t_grid, self.distance_x, v2, D2, k2)
            stage_output = signal.fftconvolve(stage_output, g2, mode="full")[:len(self.t_grid)]

        dilution_series = dilution
        baseflow_series = baseflow
        if self.use_seasonal:
            season = 2 * np.pi * (self.t_grid / 365.0)
            dilution_amp = safe["dilution_amp"]
            dilution_phase = 2 * np.pi * (
                safe["dilution_phase"] / 365.0
            )
            baseflow_amp = safe["baseflow_amp"]
            baseflow_phase = 2 * np.pi * (
                safe["baseflow_phase"] / 365.0
            )
            dilution_series = dilution * (1.0 + dilution_amp * np.sin(season + dilution_phase))
            baseflow_series = baseflow * (1.0 + baseflow_amp * np.sin(season + baseflow_phase))
            dilution_series = np.maximum(dilution_series, 0.0)
            baseflow_series = np.maximum(baseflow_series, 0.0)

        mixed = stage_output * dilution_series
        smoothed = self._apply_reservoir(mixed, tau)
        return smoothed + baseflow_series

    def get_parameters(self) -> List[AdjustableParameter]:
        return self.parameters

    def get_observations(self) -> List[Observation]:
        return self.observations

    def get_observation_targets(self) -> Dict[str, float]:
        return dict(self.observation_values_raw)

    def simulate_series(self, param_values: Dict[str, float]) -> np.ndarray:
        return self._simulate(param_values)

    def run_model(self, param_values: Dict[str, float]) -> Dict[str, float]:
        try:
            c_out = self._simulate(param_values)
        except Exception:
            return self._penalty_results()

        results: Dict[str, float] = {}
        for obs_name in self.observation_values_raw.keys():
            try:
                day = int(obs_name.split("_")[1])
                if 0 <= day < len(c_out):
                    raw_value = float(c_out[day])
                    if not np.isfinite(raw_value):
                        raw_value = -9999.0
                    if self.use_hybrid:
                        results[obs_name] = raw_value
                        results[f"log_{obs_name}"] = self._log_transform(raw_value)
                    else:
                        results[obs_name] = self._transform(raw_value)
                else:
                    results[obs_name] = -9999.0
                    if self.use_hybrid:
                        results[f"log_{obs_name}"] = -9999.0
            except Exception:
                results[obs_name] = -9999.0
                if self.use_hybrid:
                    results[f"log_{obs_name}"] = -9999.0
        return results

    def run_model_raw(self, param_values: Dict[str, float]) -> Dict[str, float]:
        try:
            c_out = self._simulate(param_values)
        except Exception:
            return {name: -9999.0 for name in self.observation_values_raw.keys()}
        results: Dict[str, float] = {}
        for obs_name in self.observation_values_raw.keys():
            try:
                day = int(obs_name.split("_")[1])
                if 0 <= day < len(c_out):
                    raw_value = float(c_out[day])
                    if not np.isfinite(raw_value):
                        raw_value = -9999.0
                    results[obs_name] = raw_value
                else:
                    results[obs_name] = -9999.0
            except Exception:
                results[obs_name] = -9999.0
        return results

    def _penalty_results(self) -> Dict[str, float]:
        results: Dict[str, float] = {}
        for obs_name in self.observation_values_raw.keys():
            results[obs_name] = -9999.0
            if self.use_hybrid:
                results[f"log_{obs_name}"] = -9999.0
        return results
