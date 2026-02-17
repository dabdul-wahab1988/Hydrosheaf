from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
import math
from typing import List, Optional, Sequence, Tuple

import numpy as np

try:
    from scipy.sparse import diags  # type: ignore
    from scipy.sparse.linalg import spsolve  # type: ignore
except Exception:  # pragma: no cover
    diags = None
    spsolve = None

from .contracts import VadoseForcingSample, VadoseProfile, VadoseRunConfig
from .soil import (
    C_from_psi,
    K_from_psi,
    VGParams,
    feddes_alpha,
    theta_from_psi,
    vg_from_texture,
)


@dataclass(frozen=True)
class RichardsStepDiagnostics:
    converged: bool
    n_iter: int
    max_delta_psi: float
    mass_balance_error_m: float


@dataclass
class RichardsSimulation:
    timestamps: List[datetime]
    z_m: List[float]  # depth grid (m, positive downward)
    psi_m: List[List[float]]  # pressure head by time, per node
    theta: List[List[float]]  # water content by time, per node
    q_faces_m_day: List[
        List[float]
    ]  # Darcy flux profile by time, faces length = n_nodes+1
    recharge_m_day: List[float]  # bottom flux (m/day)
    surface_flux_m_day: List[float]  # applied top flux (m/day)
    transpiration_m_day: List[float]  # total transpiration sink (m/day)
    diagnostics: List[RichardsStepDiagnostics]


def _layer_params_by_node(profile: VadoseProfile, z_m: np.ndarray) -> List[VGParams]:
    params: List[VGParams] = []
    layer_edges = []
    depth = 0.0
    for layer in profile.layers:
        depth += float(layer.thickness_m)
        layer_edges.append(depth)

    def _layer_for_depth(d: float) -> int:
        for idx, edge in enumerate(layer_edges):
            if d <= edge + 1e-12:
                return idx
        return len(layer_edges) - 1

    for d in z_m:
        layer = profile.layers[_layer_for_depth(float(d))]
        ks_val = float(layer.ks_m_day) if layer.ks_m_day is not None else 0.0
        
        # Apply quasi-2D tensor anisotropy correction if defined
        # K_zz = K_par * sin^2(alpha) + K_perp * cos^2(alpha)
        # where K_perp = ks, K_par = ks * anisotropy
        if layer.dip_angle_deg != 0.0 or layer.anisotropy_ratio != 1.0:
            rad = math.radians(float(layer.dip_angle_deg))
            s2 = math.sin(rad) ** 2
            c2 = math.cos(rad) ** 2
            # Use provided ks as the base perpendicular conductivity if not texture-derived
            # If texture derived, we handle it below
            pass # calculated inside the block or after texture lookup

        final_ks = 0.0
        
        if (
            layer.theta_r is not None
            and layer.theta_s is not None
            and layer.alpha_1_m is not None
            and layer.n is not None
            and layer.ks_m_day is not None
        ):
            # Explicit parameters
            base_ks = float(layer.ks_m_day)
            rad = math.radians(float(layer.dip_angle_deg))
            s2 = math.sin(rad) ** 2
            c2 = math.cos(rad) ** 2
            ratio = float(layer.anisotropy_ratio)
            # K_eff = (base * ratio) * sin^2 + base * cos^2 = base * (ratio * s2 + c2)
            final_ks = base_ks * (ratio * s2 + c2)
            
            params.append(
                VGParams(
                    theta_r=float(layer.theta_r),
                    theta_s=float(layer.theta_s),
                    alpha_1_m=float(layer.alpha_1_m),
                    n=float(layer.n),
                    ks_m_day=final_ks,
                    pore_connectivity=float(layer.pore_connectivity),
                )
            )
            continue
            
        if layer.texture:
            p = vg_from_texture(layer.texture)
            if p is not None:
                # Texture derived base parameters
                base_ks = float(p.ks_m_day)
                rad = math.radians(float(layer.dip_angle_deg))
                s2 = math.sin(rad) ** 2
                c2 = math.cos(rad) ** 2
                ratio = float(layer.anisotropy_ratio)
                final_ks = base_ks * (ratio * s2 + c2)
                
                # Clone p with new ks
                params.append(VGParams(
                    theta_r=p.theta_r,
                    theta_s=p.theta_s,
                    alpha_1_m=p.alpha_1_m,
                    n=p.n,
                    ks_m_day=final_ks,
                    pore_connectivity=p.pore_connectivity
                ))
                continue
                
        raise ValueError(
            "Vadose layer missing hydraulic params (theta_r/theta_s/alpha/n/ks) and no known texture mapping."
        )
    return params


def _root_weights(z_m: np.ndarray, root_depth_m: float) -> np.ndarray:
    if root_depth_m <= 0:
        return np.zeros_like(z_m)
    w = np.where(z_m <= root_depth_m, 1.0, 0.0)
    return w


def run_richards_column(
    profile: VadoseProfile,
    forcing: Sequence[VadoseForcingSample],
    *,
    config: Optional[VadoseRunConfig] = None,
    water_table_depth_m: Optional[Sequence[Tuple[datetime, float]]] = None,
    initial_psi_m: float = -2.0,
) -> RichardsSimulation:
    if config is None:
        config = VadoseRunConfig()
    if diags is None or spsolve is None:
        raise RuntimeError(
            "SciPy is required for the vadose Richards solver (scipy.sparse)."
        )

    if not forcing:
        raise ValueError("forcing is empty")

    dz = float(config.dz_m)
    if dz <= 0:
        raise ValueError("dz_m must be > 0")
    depth_m = float(profile.depth_m)
    if depth_m <= dz:
        raise ValueError("profile.depth_m must be > dz_m")

    n_nodes = int(np.floor(depth_m / dz)) + 1
    z = np.linspace(0.0, depth_m, n_nodes)
    params_by_node = _layer_params_by_node(profile, z)

    root_depth = (
        float(profile.root_depth_m) if profile.root_depth_m is not None else depth_m
    )
    root_depth = float(min(max(root_depth, 0.0), depth_m))
    root_w = _root_weights(z, root_depth)
    root_norm = float(np.sum(root_w) * dz) if float(np.sum(root_w)) > 0 else 1.0

    psi_n = np.full(n_nodes, float(initial_psi_m), dtype=float)

    timestamps: List[datetime] = []
    psi_hist: List[List[float]] = []
    theta_hist: List[List[float]] = []
    recharge: List[float] = []
    surface_flux: List[float] = []
    transpiration: List[float] = []
    q_faces_hist: List[List[float]] = []
    diags_out: List[RichardsStepDiagnostics] = []

    def _wt_at_time(t: datetime) -> Optional[float]:
        if not water_table_depth_m:
            return None
        # nearest previous value (piecewise constant)
        best = None
        for tt, val in water_table_depth_m:
            if tt <= t:
                best = float(val)
            else:
                break
        return best

    for i, sample in enumerate(forcing):
        if i == 0:
            dt_days = 1.0
        else:
            dt_days = max(
                1e-9,
                (sample.timestamp - forcing[i - 1].timestamp).total_seconds() / 86400.0,
            )
        dt = float(dt_days)

        etp_m_day = float(config.kc) * float(sample.et0_mm_day) / 1000.0
        evap_m_day = float(config.evaporation_fraction) * etp_m_day
        transp_m_day = (1.0 - float(config.evaporation_fraction)) * etp_m_day

        infil_m_day = (
            float(sample.precipitation_mm_day) + float(sample.irrigation_mm_day)
        ) / 1000.0
        # Atmospheric top flux (positive downward)
        q_top = infil_m_day - evap_m_day

        # Lower boundary: optional constant head from water table depth (requires wt deeper than profile)
        psi_bottom_bc = None
        is_seepage_active = False
        
        if config.lower_boundary == "constant_head_from_wt":
            wt = _wt_at_time(sample.timestamp)
            if wt is None:
                raise ValueError(
                    "lower_boundary=constant_head_from_wt requires water_table_depth_m time series"
                )
            if config.require_wt_deeper_than_profile and wt <= depth_m + 1e-12:
                raise ValueError(
                    "water table depth must be deeper than profile depth when require_wt_deeper_than_profile=True"
                )
            # Approximate hydrostatic suction at column bottom: psi ≈ -(wt - depth)
            psi_bottom_bc = float(-(wt - depth_m))
        elif config.lower_boundary == "seepage_face":
            # Initial guess for seepage state based on previous time step
            if psi_n[-1] >= 0.0:
                is_seepage_active = True

        theta_n = np.array(
            [
                theta_from_psi(np.array([psi_n[j]]), params_by_node[j])[0]
                for j in range(n_nodes)
            ],
            dtype=float,
        )

        psi_k = psi_n.copy()
        converged = False
        max_delta = float("inf")
        n_iter = 0

        for n_iter in range(1, int(config.max_picard_iter) + 1):
            # Constitutive terms at iterate
            theta_k = np.array(
                [
                    theta_from_psi(np.array([psi_k[j]]), params_by_node[j])[0]
                    for j in range(n_nodes)
                ],
                dtype=float,
            )
            C_k = np.array(
                [
                    C_from_psi(np.array([psi_k[j]]), params_by_node[j])[0]
                    for j in range(n_nodes)
                ],
                dtype=float,
            )
            # Ensure non-negative capacity
            C_k = np.maximum(C_k, 0.0)

            K_k = np.array(
                [
                    K_from_psi(np.array([psi_k[j]]), params_by_node[j])[0]
                    for j in range(n_nodes)
                ],
                dtype=float,
            )
            # Face conductivities: harmonic mean
            K_face = np.zeros(n_nodes - 1, dtype=float)
            for j in range(n_nodes - 1):
                a = float(K_k[j])
                b = float(K_k[j + 1])
                if a <= 0 or b <= 0:
                    K_face[j] = 0.0
                else:
                    K_face[j] = 2.0 * a * b / (a + b)

            # Root uptake sink (transpiration), distributed and stressed by psi
            alpha = feddes_alpha(
                psi_k,
                h_anoxic_m=float(config.feddes_h_anoxic_m),
                h_opt_m=float(config.feddes_h_opt_m),
                h_wilt_m=float(config.feddes_h_wilt_m),
            )
            sink = np.zeros(n_nodes, dtype=float)
            if transp_m_day > 0 and np.any(root_w > 0):
                sink = alpha * root_w * (transp_m_day / root_norm)

            # Build tri-diagonal system A psi^{k+1} = rhs
            main = np.zeros(n_nodes, dtype=float)
            lower = np.zeros(n_nodes - 1, dtype=float)
            upper = np.zeros(n_nodes - 1, dtype=float)
            rhs = np.zeros(n_nodes, dtype=float)

            # Interior nodes
            for j in range(1, n_nodes - 1):
                Km = float(K_face[j - 1])
                Kp = float(K_face[j])
                main[j] = (C_k[j] / dt) + (Km + Kp) / (dz * dz)
                lower[j - 1] = -Km / (dz * dz)
                upper[j] = -Kp / (dz * dz)
                rhs[j] = (
                    (C_k[j] / dt) * psi_k[j]
                    + (theta_n[j] - theta_k[j]) / dt
                    + (Km - Kp) / dz
                    - sink[j]
                )

            # Top node j=0 with specified flux q_top
            if n_nodes >= 2:
                Kp = float(K_face[0])
                main[0] = (C_k[0] / dt) + (Kp) / (dz * dz)
                upper[0] = -Kp / (dz * dz)
                rhs[0] = (
                    (C_k[0] / dt) * psi_k[0]
                    + (theta_n[0] - theta_k[0]) / dt
                    + (q_top - Kp) / dz
                    - sink[0]
                )

            # Bottom node j=n_nodes-1
            j = n_nodes - 1
            Km = float(K_face[-1]) if n_nodes >= 2 else 0.0
            
            # Determine BC for this iteration
            iter_dirichlet_val: Optional[float] = None
            iter_flux_q: float = float(K_k[-1])  # Default free drainage
            
            if psi_bottom_bc is not None:
                iter_dirichlet_val = float(psi_bottom_bc)
            elif config.lower_boundary == "seepage_face":
                if is_seepage_active:
                    iter_dirichlet_val = 0.0
                else:
                    iter_flux_q = 0.0 # Zero flux when unsaturated/closed
            
            if iter_dirichlet_val is not None:
                # Impose Dirichlet on bottom node by a strong penalty row.
                main[j] = 1.0
                if n_nodes >= 2:
                    lower[-1] = 0.0
                rhs[j] = float(iter_dirichlet_val)
            else:
                main[j] = (C_k[j] / dt) + (Km) / (dz * dz)
                if n_nodes >= 2:
                    lower[-1] = -Km / (dz * dz)
                rhs[j] = (
                    (C_k[j] / dt) * psi_k[j]
                    + (theta_n[j] - theta_k[j]) / dt
                    + (Km - iter_flux_q) / dz
                    - sink[j]
                )

            A = diags(
                diagonals=[lower, main, upper],
                offsets=[-1, 0, 1],
                shape=(n_nodes, n_nodes),
                format="csc",
            )
            psi_next = np.array(spsolve(A, rhs), dtype=float)
            
            # Update Seepage Face state for next iteration
            if config.lower_boundary == "seepage_face":
                if is_seepage_active:
                    # Check flux. If inflow (q < 0), switch to closed.
                    # q_bot = K_face_bot * (1 - (psi_bot - psi_up)/dz)
                    # psi_bot = 0
                    if n_nodes >= 2:
                        # Use K_face from this iter (approx)
                        kf = Km # Km is K_face[-1]
                        psi_up = psi_next[-2]
                        q_check = kf * (1.0 + psi_up / dz)
                        if q_check < -1e-9: # tolerance
                            is_seepage_active = False
                else:
                    # Check pressure. If saturated, switch to open.
                    if psi_next[-1] >= 0.0:
                        is_seepage_active = True

            max_delta = float(np.max(np.abs(psi_next - psi_k)))
            psi_k = psi_next
            if max_delta <= float(config.picard_tol_psi_m):
                converged = True
                break

        psi_np1 = psi_k
        theta_np1 = np.array(
            [
                theta_from_psi(np.array([psi_np1[j]]), params_by_node[j])[0]
                for j in range(n_nodes)
            ],
            dtype=float,
        )

        # Fluxes for diagnostics and recharge
        K_np1 = np.array(
            [
                K_from_psi(np.array([psi_np1[j]]), params_by_node[j])[0]
                for j in range(n_nodes)
            ],
            dtype=float,
        )
        q_faces = np.zeros(n_nodes + 1, dtype=float)
        q_faces[0] = q_top
        for j in range(n_nodes - 1):
            Kf = 2.0 * K_np1[j] * K_np1[j + 1] / max(1e-12, (K_np1[j] + K_np1[j + 1]))
            q_faces[j + 1] = Kf * (1.0 - (psi_np1[j + 1] - psi_np1[j]) / dz)
        # Bottom face
        if psi_bottom_bc is not None:
            # Darcy flux with prescribed head (Dirichlet)
            Kf = float(K_np1[-1])
            # q = -K * (d_psi/dz - 1) = K * (1 - (psi_v - psi_u)/dz)
            # here z is positive down, so flux is K * (1 - (psi_bot - psi_prev)/dz)
            q_faces[-1] = Kf * (1.0 - (float(psi_bottom_bc) - psi_np1[-2]) / dz)
        else:
            # Free drainage (Unit gradient: q = K)
            q_faces[-1] = float(K_np1[-1])


        # Mass balance check over the column depth:
        storage_n = float(np.sum(theta_n) * dz)
        storage_np1 = float(np.sum(theta_np1) * dz)
        sink_total = float(np.sum(sink) * dz)
        inflow = float(q_faces[0])
        outflow = float(q_faces[-1])
        # storage change = (inflow - outflow - sink) * dt
        mb_err = (storage_np1 - storage_n) - (inflow - outflow - sink_total) * dt

        timestamps.append(sample.timestamp)
        psi_hist.append([float(v) for v in psi_np1])
        theta_hist.append([float(v) for v in theta_np1])
        q_faces_hist.append([float(v) for v in q_faces.tolist()])
        recharge.append(float(q_faces[-1]))
        surface_flux.append(float(q_top))
        transpiration.append(float(transp_m_day))
        diags_out.append(
            RichardsStepDiagnostics(
                converged=bool(converged),
                n_iter=int(n_iter),
                max_delta_psi=float(max_delta),
                mass_balance_error_m=float(mb_err),
            )
        )

        psi_n = psi_np1

    return RichardsSimulation(
        timestamps=timestamps,
        z_m=[float(v) for v in z.tolist()],
        psi_m=psi_hist,
        theta=theta_hist,
        q_faces_m_day=q_faces_hist,
        recharge_m_day=recharge,
        surface_flux_m_day=surface_flux,
        transpiration_m_day=transpiration,
        diagnostics=diags_out,
    )
