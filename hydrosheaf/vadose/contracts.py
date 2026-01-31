from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
from typing import List, Literal, Optional


@dataclass(frozen=True)
class VadoseForcingSample:
    """Daily (or sub-daily) forcing used by the vadose column model."""

    timestamp: datetime
    precipitation_mm_day: float
    et0_mm_day: float
    irrigation_mm_day: float = 0.0


@dataclass(frozen=True)
class VadoseLayer:
    """A soil layer parameterization.

    Parameters use van Genuchten-Mualem form:
      theta_r, theta_s: residual/saturated water content (vol/vol)
      alpha_1_m: alpha parameter (1/m)
      n: shape parameter (dimensionless)
      ks_m_day: saturated hydraulic conductivity (m/day)
      l: pore-connectivity parameter (dimensionless, default 0.5)

    `texture` is an optional label that can be mapped to typical parameters
    when explicit hydraulic parameters are not provided.
    """

    thickness_m: float
    texture: Optional[str] = None
    theta_r: Optional[float] = None
    theta_s: Optional[float] = None
    alpha_1_m: Optional[float] = None
    n: Optional[float] = None
    ks_m_day: Optional[float] = None
    pore_connectivity: float = 0.5

    # Anisotropy and Geometry (Quasi-2D support via tensor projection)
    # dip_angle_deg: Angle of layering relative to horizontal (0 = flat).
    # anisotropy_ratio: K_parallel / K_perpendicular (default 1.0 = isotropic).
    dip_angle_deg: float = 0.0
    anisotropy_ratio: float = 1.0


@dataclass(frozen=True)
class VadoseProfile:
    """Vertical 1D column definition."""

    profile_id: str
    depth_m: float
    layers: List[VadoseLayer]
    # Optional: a representative rooting depth (m). When None, config root depth is used.
    root_depth_m: Optional[float] = None


@dataclass(frozen=True)
class VadoseRunConfig:
    """Numerical + process configuration for the vadose solver."""

    dz_m: float = 0.1
    max_picard_iter: int = 25
    picard_tol_psi_m: float = 1e-3
    theta_min: float = 1e-6
    q_min_m_day: float = 1e-9
    mass_balance_tol_m: float = 1e-3
    min_converged_fraction: float = 0.9

    # ET controls
    kc: float = 1.0
    evaporation_fraction: float = 0.3  # fraction of ETp treated as surface evaporation

    # Root uptake stress (Feddes-style, simplified piecewise in pressure head, meters)
    feddes_h_wilt_m: float = -150.0
    feddes_h_opt_m: float = -3.0
    feddes_h_anoxic_m: float = -0.1

    # Lower boundary control
    lower_boundary: Literal["free_drainage", "constant_head_from_wt", "seepage_face"] = "free_drainage"
    require_wt_deeper_than_profile: bool = True

    # Travel-time distribution (TTD) kernel construction (gamma mixture approximation)
    ttd_grid_dt_days: float = 5.0
    ttd_max_lag_days: float = 3650.0
    ttd_default_cv: float = (
        0.7  # coefficient of variation for gamma kernel around advective mean
    )

    # Dual-domain / Preferential flow parameters (statistical approximation)
    preferential_flow_fraction: float = 0.0  # Fraction of flow bypassing matrix (0.0 to 1.0)
    preferential_velocity_factor: float = 10.0  # How much faster is the preferential path?
    preferential_cv: float = 1.0  # CV for the preferential component



@dataclass(frozen=True)
class VadoseLinksRow:
    """Mapping row used to produce edge priors (vadose node -> groundwater node)."""

    u: str
    v: str
    u_depth_m: float
    p_uv: float = 1.0
