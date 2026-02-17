from __future__ import annotations

import csv
from dataclasses import dataclass
from datetime import datetime
from typing import List, Optional, Sequence, Tuple

import math
import numpy as np

try:
    from scipy.optimize import least_squares  # type: ignore
except Exception:  # pragma: no cover
    least_squares = None

from .contracts import VadoseForcingSample, VadoseLayer, VadoseProfile, VadoseRunConfig
from .richards1d import run_richards_column
from .nitrate import NitrateLoadingSample, predict_no3_breakthrough
from .ttd import advective_travel_time_days, mixture_ttd_from_series


@dataclass(frozen=True)
class ThetaObservation:
    timestamp: datetime
    depth_m: float
    theta: float


@dataclass(frozen=True)
class TracerObservation:
    timestamp: datetime
    concentration: float
    # Optional: observation type/location? For now assume it's at the bottom (recharge).


@dataclass(frozen=True)
class VadoseCalibrationResult:
    ks_multiplier: float
    kc: float
    preferential_flow_fraction: float
    ttd_cv: float
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
                pore_connectivity=layer.pore_connectivity,
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
    tracer_observations: Optional[Sequence[TracerObservation]] = None,
    nitrate_loading: Optional[Sequence[NitrateLoadingSample]] = None,
    config: Optional[VadoseRunConfig] = None,
    water_table_depth_m: Optional[Sequence[Tuple[datetime, float]]] = None,
    ks_bounds: Tuple[float, float] = (0.1, 10.0),
    kc_bounds: Tuple[float, float] = (0.3, 1.5),
    preferential_flow_bounds: Tuple[float, float] = (0.0, 0.4),
    ttd_cv_bounds: Tuple[float, float] = (0.1, 2.0),
    max_nfev: int = 25,
) -> VadoseCalibrationResult:
    """Calibrate Ks, Kc, and optionally preferential flow parameters against θ and tracer data."""
    if least_squares is None:
        raise RuntimeError(
            "SciPy optimize is required for vadose calibration (scipy.optimize.least_squares)."
        )
    cfg0 = config or VadoseRunConfig()
    if not observations and not tracer_observations:
        raise ValueError("No observations provided for calibration.")

    # --- Setup Theta Observations ---
    obs_theta_arr = np.array([], dtype=float)
    obs_depth_idx = []
    obs_time_idx = []
    
    # Map observation depths to nearest node index once (grid fixed).
    dz = float(cfg0.dz_m)
    n_nodes = int(np.floor(float(profile.depth_m) / dz)) + 1
    z = np.linspace(0.0, float(profile.depth_m), n_nodes)
    
    # Helper: map timestamps to nearest forcing index.
    forcing_times = [f.timestamp for f in forcing]
    def _time_to_index(t: datetime) -> int:
        diffs = [abs((tt - t).total_seconds()) for tt in forcing_times]
        return int(np.argmin(diffs))

    if observations:
        obs_depth_idx_list = []
        obs_time_list = []
        obs_theta_list = []
        for o in observations:
            idx = int(np.argmin(np.abs(z - float(o.depth_m))))
            obs_depth_idx_list.append(idx)
            obs_time_list.append(o.timestamp)
            obs_theta_list.append(float(o.theta))
        obs_theta_arr = np.asarray(obs_theta_list, dtype=float)
        obs_time_idx = [_time_to_index(t) for t in obs_time_list]
        obs_depth_idx = obs_depth_idx_list

    # --- Setup Tracer Observations ---
    tracer_obs_arr = np.array([], dtype=float)
    tracer_time_idx = []
    tracer_scale = 1.0
    
    if tracer_observations:
        if not nitrate_loading:
             raise ValueError("nitrate_loading required if tracer_observations are used")
        
        tr_obs_list = []
        tr_time_list = []
        for o in tracer_observations:
            tr_obs_list.append(float(o.concentration))
            tr_time_list.append(o.timestamp)
        tracer_obs_arr = np.asarray(tr_obs_list, dtype=float)
        # Normalize scale for optimization balance
        mean_obs = np.mean(tracer_obs_arr)
        if mean_obs > 1e-6:
             tracer_scale = 1.0 / mean_obs
        
        # We need to map tracer obs times to the output of predict_no3_breakthrough
        # But predict_no3_breakthrough returns points at forcing timestamps (usually).
        # We'll just map to nearest forcing index like theta.
        tracer_time_idx = [_time_to_index(t) for t in tr_time_list]


    # Param vector: [log_ks_mult, kc, pref_frac, ttd_cv]
    # We include pref_frac and ttd_cv even if no tracer obs? 
    # If no tracer obs, transport is unconstrained. We should lock them or ignore them.
    # For simplicity, we always include them but if no tracer data, they won't affect cost much 
    # (except maybe regularizing or if we add priors).
    # Better: only optimize them if tracer data exists.
    
    use_transport = bool(tracer_observations and nitrate_loading)
    
    # Initial values
    x0_list = [0.0, float(cfg0.kc)]
    lower_list = [math.log10(float(ks_bounds[0])), float(kc_bounds[0])]
    upper_list = [math.log10(float(ks_bounds[1])), float(kc_bounds[1])]
    
    if use_transport:
        x0_list.extend([float(cfg0.preferential_flow_fraction), float(cfg0.ttd_default_cv)])
        lower_list.extend([float(preferential_flow_bounds[0]), float(ttd_cv_bounds[0])])
        upper_list.extend([float(preferential_flow_bounds[1]), float(ttd_cv_bounds[1])])

    x0 = np.array(x0_list, dtype=float)
    lower = np.array(lower_list, dtype=float)
    upper = np.array(upper_list, dtype=float)

    def residuals(x: np.ndarray) -> np.ndarray:
        log_ks_mult = float(x[0])
        kc = float(x[1])
        
        pref_frac = float(cfg0.preferential_flow_fraction)
        ttd_cv = float(cfg0.ttd_default_cv)
        
        if use_transport:
            pref_frac = float(x[2])
            ttd_cv = float(x[3])

        ks_mult = 10.0**log_ks_mult
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
            ttd_default_cv=ttd_cv,
            preferential_flow_fraction=pref_frac,
            preferential_velocity_factor=cfg0.preferential_velocity_factor,
            preferential_cv=cfg0.preferential_cv,
        )
        
        sim = run_richards_column(
            prof, forcing, config=cfg, water_table_depth_m=water_table_depth_m
        )
        
        # Theta residuals
        res_list = []
        if len(obs_time_idx) > 0:
            pred = []
            for ti, zi in zip(obs_time_idx, obs_depth_idx):
                ti = max(0, min(ti, len(sim.theta) - 1))
                zi = max(0, min(zi, len(sim.theta[ti]) - 1))
                pred.append(float(sim.theta[ti][zi]))
            pred_arr = np.asarray(pred, dtype=float)
            res_list.append(pred_arr - obs_theta_arr)
            
        # Tracer residuals
        if use_transport and nitrate_loading:
            # Need to compute TTD and run convolution
            # We need tau series.
            # Reuse logic from run.py roughly
            z_arr = np.asarray(sim.z_m, dtype=float)
            weights = np.maximum(0.0, np.asarray(sim.recharge_m_day, dtype=float))
            
            tau_series = []
            # We assume tracer breakthrough at bottom (depth_m)
            # advective_travel_time_days needs u_depth_m=0 (surface loading)
            for t_idx in range(len(sim.timestamps)):
                theta_t = np.asarray(sim.theta[t_idx], dtype=float)
                q_faces = np.asarray(sim.q_faces_m_day[t_idx], dtype=float)
                tau = advective_travel_time_days(
                    z_m=z_arr.tolist(),
                    theta=theta_t.tolist(),
                    q_faces_m_day=q_faces.tolist(),
                    from_depth_m=0.0,
                    theta_min=float(cfg.theta_min),
                    q_min_m_day=float(cfg.q_min_m_day),
                )
                tau_series.append(float(tau))

            tau_grid, g, _ = mixture_ttd_from_series(
                tau_series,
                weights.tolist(),
                grid_dt_days=float(cfg.ttd_grid_dt_days),
                max_lag_days=float(cfg.ttd_max_lag_days),
                cv=float(cfg.ttd_default_cv),
                preferential_flow_fraction=float(cfg.preferential_flow_fraction),
                preferential_velocity_factor=float(cfg.preferential_velocity_factor),
                preferential_cv=float(cfg.preferential_cv),
            )
            
            # Run convolution
            # predict_no3_breakthrough returns list of BreakthroughPoint
            out_points, _ = predict_no3_breakthrough(
                edge_id="calibration",
                ttd_tau_days=tau_grid,
                ttd_g=g,
                timestamps=sim.timestamps,
                recharge_m_day=sim.recharge_m_day,
                loading=nitrate_loading,
                grid_dt_days=float(cfg.ttd_grid_dt_days)
            )
            
            # Map predictions to obs times
            # out_points are aligned with sim.timestamps (forcing)
            # So we can use tracer_time_idx directly
            pred_tr = []
            for ti in tracer_time_idx:
                 ti = max(0, min(ti, len(out_points) - 1))
                 pred_tr.append(float(out_points[ti].c_wt))
            
            pred_tr_arr = np.asarray(pred_tr, dtype=float)
            res_tr = (pred_tr_arr - tracer_obs_arr) * tracer_scale
            res_list.append(res_tr)

        if not res_list:
            return np.array([0.0]) # Should not happen
            
        return np.concatenate(res_list)

    res = least_squares(residuals, x0=x0, bounds=(lower, upper), max_nfev=int(max_nfev))
    
    ks_mult = float(10.0 ** float(res.x[0]))
    kc_final = float(res.x[1])
    pref_frac_final = float(res.x[2]) if use_transport else float(cfg0.preferential_flow_fraction)
    ttd_cv_final = float(res.x[3]) if use_transport else float(cfg0.ttd_default_cv)
    
    return VadoseCalibrationResult(
        ks_multiplier=ks_mult,
        kc=kc_final,
        preferential_flow_fraction=pref_frac_final,
        ttd_cv=ttd_cv_final,
        cost=float(res.cost),
        nfev=int(res.nfev),
        status=int(res.status),
        message=str(res.message),
    )

