from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from ..physics.priors import PhysicsPrior
from .contracts import (
    VadoseForcingSample,
    VadoseLinksRow,
    VadoseProfile,
    VadoseRunConfig,
)
from .richards1d import run_richards_column
from .ttd import advective_travel_time_days, mixture_ttd_from_series


@dataclass(frozen=True)
class VadoseRunResult:
    profile_id: str
    timestamps: List[datetime]
    recharge_m_day: List[float]
    surface_flux_m_day: List[float]
    # Optional TTD kernel (mixture gamma approximation) from a specified depth
    ttd_tau_days: Optional[List[float]] = None
    ttd_g: Optional[List[float]] = None


@dataclass(frozen=True)
class VadoseEdgePrior:
    u: str
    v: str
    p_uv: float
    tt_mean_days: float
    tt_std_days: float
    tt_p10_days: float
    tt_p90_days: float
    source: str = "vadose_richards"

    def to_physics_prior(self) -> PhysicsPrior:
        return PhysicsPrior(
            u=self.u,
            v=self.v,
            p_uv=self.p_uv,
            tt_mean_days=self.tt_mean_days,
            tt_std_days=self.tt_std_days,
            tt_p10_days=self.tt_p10_days,
            tt_p90_days=self.tt_p90_days,
            source=self.source,
        )


def run_vadose_profile(
    profile: VadoseProfile,
    forcing: Sequence[VadoseForcingSample],
    *,
    config: Optional[VadoseRunConfig] = None,
    water_table_depth_m: Optional[Sequence[Tuple[datetime, float]]] = None,
) -> VadoseRunResult:
    cfg = config or VadoseRunConfig()
    sim = run_richards_column(
        profile, forcing, config=cfg, water_table_depth_m=water_table_depth_m
    )
    return VadoseRunResult(
        profile_id=profile.profile_id,
        timestamps=list(sim.timestamps),
        recharge_m_day=list(sim.recharge_m_day),
        surface_flux_m_day=list(sim.surface_flux_m_day),
    )


def build_vadose_edge_priors(
    profile: VadoseProfile,
    forcing: Sequence[VadoseForcingSample],
    links: Sequence[VadoseLinksRow],
    *,
    config: Optional[VadoseRunConfig] = None,
    water_table_depth_m: Optional[Sequence[Tuple[datetime, float]]] = None,
) -> Tuple[List[VadoseEdgePrior], Dict[str, Dict[str, object]]]:
    """Run vadose physics and build PhysicsPrior-compatible edge summaries.

    Returns:
      - per-link edge priors (mean/std/quantiles of advective travel time, recharge-weighted)
      - diagnostics per link, including optional mixture TTD kernel for reproducibility
    """
    cfg = config or VadoseRunConfig()
    sim = run_richards_column(
        profile, forcing, config=cfg, water_table_depth_m=water_table_depth_m
    )

    # dz = float(sim.z_m[1] - sim.z_m[0]) if len(sim.z_m) >= 2 else float(cfg.dz_m)
    z = np.asarray(sim.z_m, dtype=float)

    # Precompute daily q_faces and daily taus for each link.
    priors: List[VadoseEdgePrior] = []
    details: Dict[str, Dict[str, object]] = {}
    details["_vadose_run"] = {
        "profile_id": profile.profile_id,
        "timestamps": [t.strftime("%Y-%m-%d") for t in sim.timestamps],
        "recharge_m_day": [float(v) for v in sim.recharge_m_day],
        "surface_flux_m_day": [float(v) for v in sim.surface_flux_m_day],
    }

    # Use recharge as weight (proxy for which times contribute to water reaching the water table).
    weights = np.maximum(0.0, np.asarray(sim.recharge_m_day, dtype=float))
    # Avoid all-zero weights.
    if float(np.sum(weights)) <= 0:
        weights = np.ones_like(weights)

    # Run-level diagnostics for transparency.
    conv = [bool(d.converged) for d in sim.diagnostics]
    conv_frac = float(np.mean(conv)) if conv else float("nan")
    max_abs_mb = (
        float(np.max([abs(float(d.mass_balance_error_m)) for d in sim.diagnostics]))
        if sim.diagnostics
        else float("nan")
    )
    run_flags: List[str] = []
    if np.isfinite(conv_frac) and conv_frac < float(cfg.min_converged_fraction):
        run_flags.append(f"low_convergence_fraction:{conv_frac:.3f}")
    if np.isfinite(max_abs_mb) and max_abs_mb > float(cfg.mass_balance_tol_m):
        run_flags.append(f"mass_balance_error>{float(cfg.mass_balance_tol_m):g}")

    for link in links:
        tau_series: List[float] = []
        for t_idx in range(len(sim.timestamps)):
            theta_t = np.asarray(sim.theta[t_idx], dtype=float)
            q_faces = np.asarray(sim.q_faces_m_day[t_idx], dtype=float)
            tau = advective_travel_time_days(
                z_m=z.tolist(),
                theta=theta_t.tolist(),
                q_faces_m_day=q_faces.tolist(),
                from_depth_m=float(link.u_depth_m),
                theta_min=float(cfg.theta_min),
                q_min_m_day=float(cfg.q_min_m_day),
            )
            tau_series.append(float(tau))

        tau_grid, g, summary = mixture_ttd_from_series(
            tau_series,
            weights.tolist(),
            grid_dt_days=float(cfg.ttd_grid_dt_days),
            max_lag_days=float(cfg.ttd_max_lag_days),
            cv=float(cfg.ttd_default_cv),
            preferential_flow_fraction=float(cfg.preferential_flow_fraction),
            preferential_velocity_factor=float(cfg.preferential_velocity_factor),
            preferential_cv=float(cfg.preferential_cv),
        )


        prior = VadoseEdgePrior(
            u=link.u,
            v=link.v,
            p_uv=float(link.p_uv),
            tt_mean_days=float(summary.mean_days),
            tt_std_days=float(summary.std_days),
            tt_p10_days=float(summary.p10_days),
            tt_p90_days=float(summary.p90_days),
        )
        priors.append(prior)
        details[f"{link.u}->{link.v}"] = {
            "profile_id": profile.profile_id,
            "profile_depth_m": float(profile.depth_m),
            "solver_converged_fraction": conv_frac,
            "solver_max_abs_mass_balance_error_m": max_abs_mb,
            "solver_flags": run_flags,
            "u_depth_m": float(link.u_depth_m),
            "tau_series_days": tau_series,
            "tau_weights": weights.tolist(),
            "ttd_tau_days": tau_grid,
            "ttd_g": g,
        }

    return priors, details
