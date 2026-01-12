from __future__ import annotations

import csv
from dataclasses import dataclass
from datetime import datetime
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import math
import numpy as np

try:
    from scipy.optimize import least_squares  # type: ignore
except Exception:  # pragma: no cover
    least_squares = None

from .contracts import VadoseForcingSample, VadoseLayer, VadoseProfile, VadoseRunConfig
from .richards1d import run_richards_column


@dataclass(frozen=True)
class ThetaObservation:
    timestamp: datetime
    depth_m: float
    theta: float


@dataclass(frozen=True)
class VadoseCalibrationResult:
    ks_multiplier: float
    kc: float
    cost: float
    nfev: int
    status: int
    message: str


def load_theta_observations_csv(
    path: str,
    *,
    time_column: str = "timestamp",
    time_format: str = "%Y-%m-%d",
    depth_column: str = "depth_m",
    theta_column: str = "theta",
) -> List[ThetaObservation]:
    obs: List[ThetaObservation] = []
    with open(path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            ts_raw = (row.get(time_column) or "").strip()
            if not ts_raw:
                continue
            try:
                ts = datetime.strptime(ts_raw, time_format)
            except ValueError:
                continue
            depth_raw = row.get(depth_column)
            theta_raw = row.get(theta_column)
            if depth_raw in (None, "") or theta_raw in (None, ""):
                continue
            try:
                depth = float(depth_raw)
                theta = float(theta_raw)
            except ValueError:
                continue
            if not (0.0 <= theta <= 1.0):
                continue
            obs.append(ThetaObservation(timestamp=ts, depth_m=depth, theta=theta))
    obs.sort(key=lambda o: o.timestamp)
    return obs


def scale_profile_ks(profile: VadoseProfile, ks_multiplier: float) -> VadoseProfile:
    layers: List[VadoseLayer] = []
    for layer in profile.layers:
        if layer.ks_m_day is None:
            layers.append(layer)
            continue
        layers.append(
            VadoseLayer(
                thickness_m=layer.thickness_m,
                texture=layer.texture,
                theta_r=layer.theta_r,
                theta_s=layer.theta_s,
                alpha_1_m=layer.alpha_1_m,
                n=layer.n,
                ks_m_day=float(layer.ks_m_day) * float(ks_multiplier),
                l=layer.l,
            )
        )
    return VadoseProfile(
        profile_id=profile.profile_id,
        depth_m=profile.depth_m,
        layers=layers,
        root_depth_m=profile.root_depth_m,
    )


def calibrate_ks_and_kc(
    profile: VadoseProfile,
    forcing: Sequence[VadoseForcingSample],
    observations: Sequence[ThetaObservation],
    *,
    config: Optional[VadoseRunConfig] = None,
    water_table_depth_m: Optional[Sequence[Tuple[datetime, float]]] = None,
    ks_bounds: Tuple[float, float] = (0.1, 10.0),
    kc_bounds: Tuple[float, float] = (0.3, 1.5),
    max_nfev: int = 25,
) -> VadoseCalibrationResult:
    """Calibrate a global K_s multiplier and Kc against observed θ(z,t).

    This is intended as a minimal, field-realistic calibration step (Tier 2).
    """
    if least_squares is None:
        raise RuntimeError("SciPy optimize is required for vadose calibration (scipy.optimize.least_squares).")
    cfg0 = config or VadoseRunConfig()
    if not observations:
        raise ValueError("No theta observations provided for calibration.")

    # Map observation depths to nearest node index once (grid fixed).
    dz = float(cfg0.dz_m)
    n_nodes = int(np.floor(float(profile.depth_m) / dz)) + 1
    z = np.linspace(0.0, float(profile.depth_m), n_nodes)

    obs_depth_idx = []
    obs_time = []
    obs_theta = []
    for o in observations:
        idx = int(np.argmin(np.abs(z - float(o.depth_m))))
        obs_depth_idx.append(idx)
        obs_time.append(o.timestamp)
        obs_theta.append(float(o.theta))
    obs_theta_arr = np.asarray(obs_theta, dtype=float)

    # Helper: map timestamps to nearest forcing index.
    forcing_times = [f.timestamp for f in forcing]
    def _time_to_index(t: datetime) -> int:
        # nearest time
        diffs = [abs((tt - t).total_seconds()) for tt in forcing_times]
        return int(np.argmin(diffs))

    obs_time_idx = [ _time_to_index(t) for t in obs_time ]

    # Param vector: [log_ks_mult, kc]
    x0 = np.array([0.0, float(cfg0.kc)], dtype=float)
    lower = np.array([math.log10(float(ks_bounds[0])), float(kc_bounds[0])], dtype=float)
    upper = np.array([math.log10(float(ks_bounds[1])), float(kc_bounds[1])], dtype=float)

    def residuals(x: np.ndarray) -> np.ndarray:
        log_ks_mult = float(x[0])
        kc = float(x[1])
        ks_mult = 10.0 ** log_ks_mult
        prof = scale_profile_ks(profile, ks_mult)
        cfg = VadoseRunConfig(
            dz_m=cfg0.dz_m,
            max_picard_iter=cfg0.max_picard_iter,
            picard_tol_psi_m=cfg0.picard_tol_psi_m,
            theta_min=cfg0.theta_min,
            q_min_m_day=cfg0.q_min_m_day,
            mass_balance_tol_m=cfg0.mass_balance_tol_m,
            min_converged_fraction=cfg0.min_converged_fraction,
            kc=kc,
            evaporation_fraction=cfg0.evaporation_fraction,
            feddes_h_wilt_m=cfg0.feddes_h_wilt_m,
            feddes_h_opt_m=cfg0.feddes_h_opt_m,
            feddes_h_anoxic_m=cfg0.feddes_h_anoxic_m,
            lower_boundary=cfg0.lower_boundary,
            require_wt_deeper_than_profile=cfg0.require_wt_deeper_than_profile,
            ttd_grid_dt_days=cfg0.ttd_grid_dt_days,
            ttd_max_lag_days=cfg0.ttd_max_lag_days,
            ttd_default_cv=cfg0.ttd_default_cv,
        )
        sim = run_richards_column(prof, forcing, config=cfg, water_table_depth_m=water_table_depth_m)
        pred = []
        for ti, zi in zip(obs_time_idx, obs_depth_idx):
            ti = max(0, min(ti, len(sim.theta) - 1))
            zi = max(0, min(zi, len(sim.theta[ti]) - 1))
            pred.append(float(sim.theta[ti][zi]))
        pred_arr = np.asarray(pred, dtype=float)
        return pred_arr - obs_theta_arr

    res = least_squares(residuals, x0=x0, bounds=(lower, upper), max_nfev=int(max_nfev))
    ks_mult = float(10.0 ** float(res.x[0]))
    kc = float(res.x[1])
    return VadoseCalibrationResult(
        ks_multiplier=ks_mult,
        kc=kc,
        cost=float(res.cost),
        nfev=int(res.nfev),
        status=int(res.status),
        message=str(res.message),
    )
