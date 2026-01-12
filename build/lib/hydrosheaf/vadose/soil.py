from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional

import math
import numpy as np


@dataclass(frozen=True)
class VGParams:
    theta_r: float
    theta_s: float
    alpha_1_m: float
    n: float
    ks_m_day: float
    l: float = 0.5

    @property
    def m(self) -> float:
        return 1.0 - 1.0 / self.n


# Typical van Genuchten parameters by USDA texture class (approximate).
# Sources: commonly-used tabulated values (e.g., Carsel & Parrish 1988 style tables).
_TEXTURE_TABLE: Dict[str, VGParams] = {
    "sand": VGParams(theta_r=0.045, theta_s=0.43, alpha_1_m=14.5, n=2.68, ks_m_day=7.1),
    "loamy_sand": VGParams(theta_r=0.057, theta_s=0.41, alpha_1_m=12.4, n=2.28, ks_m_day=3.5),
    "sandy_loam": VGParams(theta_r=0.065, theta_s=0.41, alpha_1_m=7.5, n=1.89, ks_m_day=1.1),
    "loam": VGParams(theta_r=0.078, theta_s=0.43, alpha_1_m=3.6, n=1.56, ks_m_day=0.25),
    "silt_loam": VGParams(theta_r=0.067, theta_s=0.45, alpha_1_m=2.0, n=1.41, ks_m_day=0.10),
    "silt": VGParams(theta_r=0.034, theta_s=0.46, alpha_1_m=1.6, n=1.37, ks_m_day=0.06),
    "sandy_clay_loam": VGParams(theta_r=0.100, theta_s=0.39, alpha_1_m=5.9, n=1.48, ks_m_day=0.32),
    "clay_loam": VGParams(theta_r=0.095, theta_s=0.41, alpha_1_m=1.9, n=1.31, ks_m_day=0.05),
    "silty_clay_loam": VGParams(theta_r=0.089, theta_s=0.43, alpha_1_m=1.0, n=1.23, ks_m_day=0.02),
    "sandy_clay": VGParams(theta_r=0.100, theta_s=0.38, alpha_1_m=2.7, n=1.23, ks_m_day=0.12),
    "silty_clay": VGParams(theta_r=0.070, theta_s=0.36, alpha_1_m=1.0, n=1.09, ks_m_day=0.01),
    "clay": VGParams(theta_r=0.068, theta_s=0.38, alpha_1_m=0.8, n=1.09, ks_m_day=0.01),
}


def vg_from_texture(texture: str) -> Optional[VGParams]:
    key = texture.strip().lower().replace(" ", "_")
    return _TEXTURE_TABLE.get(key)


def theta_from_psi(psi_m: np.ndarray, p: VGParams) -> np.ndarray:
    """Water content θ(ψ) for pressure head psi (m)."""
    psi = np.asarray(psi_m, dtype=float)
    theta = np.empty_like(psi)
    theta_sat = float(p.theta_s)
    theta_res = float(p.theta_r)

    # Saturated when psi >= 0 (pressure head non-negative).
    saturated = psi >= 0.0
    theta[saturated] = theta_sat

    unsat = ~saturated
    if np.any(unsat):
        h = np.abs(psi[unsat])
        m = p.m
        se = 1.0 / np.power(1.0 + np.power(p.alpha_1_m * h, p.n), m)
        theta[unsat] = theta_res + se * (theta_sat - theta_res)
    return theta


def C_from_psi(psi_m: np.ndarray, p: VGParams) -> np.ndarray:
    """Specific moisture capacity C(ψ) = dθ/dψ (1/m)."""
    psi = np.asarray(psi_m, dtype=float)
    C = np.zeros_like(psi)
    unsat = psi < 0.0
    if not np.any(unsat):
        return C

    h = np.abs(psi[unsat])
    m = p.m
    alpha = p.alpha_1_m
    n = p.n
    theta_sat = p.theta_s
    theta_res = p.theta_r
    # dSe/dh for van Genuchten:
    # Se = (1 + (a h)^n)^(-m)
    # dSe/dh = -m * (1 + (a h)^n)^(-m-1) * n * (a^n) * h^(n-1)
    ahn = np.power(alpha * h, n)
    denom = np.power(1.0 + ahn, m + 1.0)
    dse_dh = -m * n * np.power(alpha, n) * np.power(h, n - 1.0) / denom
    dtheta_dh = (theta_sat - theta_res) * dse_dh
    # dθ/dψ: for psi<0, h = -psi, dh/dψ = -1.
    C[unsat] = -dtheta_dh
    return C


def K_from_psi(psi_m: np.ndarray, p: VGParams) -> np.ndarray:
    """Hydraulic conductivity K(ψ) (m/day)."""
    psi = np.asarray(psi_m, dtype=float)
    K = np.empty_like(psi)
    Ks = float(p.ks_m_day)
    saturated = psi >= 0.0
    K[saturated] = Ks
    unsat = ~saturated
    if np.any(unsat):
        h = np.abs(psi[unsat])
        m = p.m
        se = 1.0 / np.power(1.0 + np.power(p.alpha_1_m * h, p.n), m)
        # Mualem-van Genuchten relative conductivity
        # Kr = Se^l * [1 - (1 - Se^(1/m))^m]^2
        se_clip = np.clip(se, 0.0, 1.0)
        term = 1.0 - np.power(1.0 - np.power(se_clip, 1.0 / max(1e-12, m)), m)
        Kr = np.power(se_clip, p.l) * (term * term)
        K[unsat] = Ks * Kr
    return K


def feddes_alpha(psi_m: np.ndarray, *, h_anoxic_m: float, h_opt_m: float, h_wilt_m: float) -> np.ndarray:
    """Simplified Feddes stress response α(ψ) in [0,1].

    - Too wet (anoxic): psi >= h_anoxic -> alpha = 0
    - Optimal: h_opt <= psi < h_anoxic -> alpha = 1
    - Drying: h_wilt < psi < h_opt -> alpha ramps linearly 0..1
    - Wilting: psi <= h_wilt -> alpha = 0
    """
    psi = np.asarray(psi_m, dtype=float)
    alpha = np.zeros_like(psi)
    # Optimal (not too wet)
    optimal = (psi < h_anoxic_m) & (psi >= h_opt_m)
    alpha[optimal] = 1.0
    # Drying ramp
    ramp = (psi < h_opt_m) & (psi > h_wilt_m)
    if np.any(ramp):
        alpha[ramp] = (psi[ramp] - h_wilt_m) / (h_opt_m - h_wilt_m)
    return np.clip(alpha, 0.0, 1.0)


