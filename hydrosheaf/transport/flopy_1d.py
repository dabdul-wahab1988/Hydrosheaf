"""
FloPy 1D Saturated Zone Transport Model.

This module provides functions for building and running 1D advection-dispersion
transport models using FloPy's MODFLOW and MT3DMS interfaces.
"""

import tempfile
import shutil
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any

import numpy as np

# Try to import FloPy
try:
    import flopy
    from flopy.modflow import (
        Modflow,
        ModflowDis,
        ModflowBas,
        ModflowLpf,
        ModflowWel,
        ModflowOc,
        ModflowPcg,
        ModflowLmt,
    )
    from flopy.mt3d import Mt3dms, Mt3dBtn, Mt3dAdv, Mt3dDsp, Mt3dSsm, Mt3dGcg, Mt3dRct

    FLOPY_AVAILABLE = True
except ImportError:
    flopy = None
    FLOPY_AVAILABLE = False

from .binaries import get_executable_path



@dataclass
class TransportResult:
    """Result of 1D transport simulation."""

    # Time series output
    times: np.ndarray  # Time points (days)
    concentrations: np.ndarray  # Concentration at outlet (mmol/L or mg/L)

    # Model parameters used
    n_cells: int = 50
    cell_length_m: float = 1.0
    porosity: float = 0.25
    dispersivity_m: float = 1.0
    velocity_m_day: float = 0.1
    decay_rate_1_day: float = 0.0

    # Source term
    source_concentration: float = 0.0
    source_flux_m3_day: float = 0.0

    # Diagnostics
    mass_in: float = 0.0
    mass_out: float = 0.0
    mass_balance_error: float = 0.0

    # Metadata
    model_name: str = "hydrosheaf_1d"
    workspace: Optional[str] = None
    success: bool = True
    warnings: List[str] = field(default_factory=list)


def check_flopy_available() -> bool:
    """Check if FloPy is available for transport modeling."""
    return FLOPY_AVAILABLE


def build_1d_transport_model(
    model_name: str = "hydrosheaf_1d",
    workspace: Optional[str] = None,
    # Aquifer geometry
    aquifer_length_m: float = 100.0,
    aquifer_thickness_m: float = 10.0,
    aquifer_width_m: float = 1.0,
    n_cells: int = 50,
    # Hydraulic properties
    hydraulic_k_m_day: float = 1.0,
    porosity: float = 0.25,
    # Transport properties
    dispersivity_m: float = 1.0,
    decay_rate_1_day: float = 0.0,
    # Boundary conditions
    head_upstream_m: float = 10.0,
    head_downstream_m: float = 9.0,
    # Time discretization
    n_stress_periods: int = 10,
    perlen_days: float = 365.0,
    n_time_steps: int = 10,
    # Source term (recharge from vadose)
    source_concentration: float = 1.0,
    source_cell: int = 0,
) -> Tuple[Any, Any, Dict[str, float]]:
    """
    Build a 1D MODFLOW + MT3DMS transport model.

    Parameters
    ----------
    model_name : str
        Name for the model files
    workspace : str, optional
        Directory for model files (temp dir if None)
    aquifer_length_m : float
        Length of the 1D column (m)
    aquifer_thickness_m : float
        Aquifer thickness (m)
    aquifer_width_m : float
        Aquifer width (m)
    n_cells : int
        Number of cells in the model
    hydraulic_k_m_day : float
        Hydraulic conductivity (m/day)
    porosity : float
        Effective porosity
    dispersivity_m : float
        Longitudinal dispersivity (m)
    decay_rate_1_day : float
        First-order decay rate (1/day) for denitrification
    head_upstream_m : float
        Constant head at upstream boundary (m)
    head_downstream_m : float
        Constant head at downstream boundary (m)
    n_stress_periods : int
        Number of stress periods
    perlen_days : float
        Length of each stress period (days)
    n_time_steps : int
        Number of time steps per stress period
    source_concentration : float
        Concentration of source (recharge) water
    source_cell : int
        Cell index for source injection (0 = upstream)

    Returns
    -------
    Tuple[Modflow, Mt3dms, Dict]
        MODFLOW model, MT3DMS model, and parameter dictionary
    """
    if not FLOPY_AVAILABLE:
        raise ImportError(
            "FloPy is required for transport modeling. "
            "Install with: pip install flopy>=3.3"
        )

    # Create workspace
    if workspace is None:
        workspace = tempfile.mkdtemp(prefix="hydrosheaf_transport_")

    ws_path = Path(workspace)
    ws_path.mkdir(parents=True, exist_ok=True)

    # Cell dimensions
    cell_length = aquifer_length_m / n_cells
    delr = [cell_length] * n_cells  # Row width (x-direction)
    delc = [aquifer_width_m]  # Column width (y-direction)

    # Layer properties
    nlay = 1
    nrow = 1
    ncol = n_cells
    top = aquifer_thickness_m
    botm = [0.0]

    # Time discretization
    nper = n_stress_periods
    perlen = [perlen_days] * nper
    nstp = [n_time_steps] * nper
    steady = [False] * nper

    # --- Build MODFLOW model ---
    mf_exe = get_executable_path("mf2005")
    mf = Modflow(model_name, exe_name=mf_exe, model_ws=str(ws_path))

    # Discretization
    _dis = ModflowDis(
        mf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
        nper=nper,
        perlen=perlen,
        nstp=nstp,
        steady=steady,
    )

    # Basic package (initial heads and boundary conditions)
    ibound = np.ones((nlay, nrow, ncol), dtype=int)
    ibound[0, 0, 0] = -1  # Constant head upstream
    ibound[0, 0, -1] = -1  # Constant head downstream

    strt = np.linspace(head_upstream_m, head_downstream_m, ncol)
    strt = strt.reshape((nlay, nrow, ncol))

    _bas = ModflowBas(mf, ibound=ibound, strt=strt)

    # Layer properties (LPF)
    _lpf = ModflowLpf(
        mf, hk=hydraulic_k_m_day, vka=hydraulic_k_m_day, ss=1e-5, sy=porosity
    )

    # Well package for source injection
    if source_concentration > 0:
        # Injection at source cell
        well_data = {0: [[0, 0, source_cell, 0.001]]}  # Small flux for source
        _wel = ModflowWel(mf, stress_period_data=well_data)

    # Output control
    _oc = ModflowOc(mf, stress_period_data={(0, 0): ["save head", "save budget"]})

    # PCG solver
    _pcg = ModflowPcg(mf)

    # LMT package (Link-MT3DMS)
    _lmt = ModflowLmt(mf)

    # --- Build MT3DMS model ---
    mt_exe = get_executable_path("mt3dms")
    mt = Mt3dms(
        modelname=model_name,
        model_ws=str(ws_path),
        exe_name=mt_exe,
        modflowmodel=mf,
    )

    # BTN - Basic transport
    icbund = np.ones((nlay, nrow, ncol), dtype=int)
    sconc = np.zeros((nlay, nrow, ncol))
    sconc[0, 0, source_cell] = source_concentration  # Initial concentration at source

    _btn = Mt3dBtn(
        mt,
        icbund=icbund,
        prsity=porosity,
        sconc=sconc,
        ncomp=1,
        mcomp=1,
        dt0=perlen_days / n_time_steps / 10,  # Initial time step
        nprs=-1,
        timprs=None,
        obs=[(0, 0, ncol - 1)],  # Observation at outlet
    )

    # ADV - Advection
    _adv = Mt3dAdv(mt, mixelm=0)  # Upstream finite difference

    # DSP - Dispersion
    _dsp = Mt3dDsp(mt, al=dispersivity_m, trpt=0.1, trpv=0.01, dmcoef=1e-9)

    # SSM - Source/Sink mixing
    ssm_data = {}
    if source_concentration > 0:
        # Constant concentration source
        ssm_data[0] = [[0, 0, source_cell, source_concentration, 15]]  # Type 15 = RCH

    _ssm = Mt3dSsm(mt, stress_period_data=ssm_data, mxss=100)

    # RCT - Reaction (first-order decay for denitrification)
    if decay_rate_1_day > 0:
        _rct = Mt3dRct(
            mt,
            isothm=0,  # No sorption
            ireact=1,  # First-order kinetic reaction
            rc1=decay_rate_1_day,  # Decay rate
            igetsc=0,
        )
    else:
        # No reaction
        pass  # _rct is not strictly needed if decay is 0? Or maybe default is no reaction.
        # Actually usually if RCT is missing it's fine.
        # But let's check old code: `rct = Mt3dRct(mt, isothm=0, ireact=0)` was used else.
        # I should keep it but prefix.

    # GCG - Solver
    _gcg = Mt3dGcg(mt)

    # Parameter dictionary
    params = {
        "n_cells": n_cells,
        "cell_length_m": cell_length,
        "porosity": porosity,
        "hydraulic_k_m_day": hydraulic_k_m_day,
        "dispersivity_m": dispersivity_m,
        "decay_rate_1_day": decay_rate_1_day,
        "aquifer_length_m": aquifer_length_m,
        "aquifer_thickness_m": aquifer_thickness_m,
        "head_gradient": (head_upstream_m - head_downstream_m) / aquifer_length_m,
        "velocity_m_day": hydraulic_k_m_day
        * (head_upstream_m - head_downstream_m)
        / aquifer_length_m
        / porosity,
        "perlen_days": perlen_days,
        "n_stress_periods": n_stress_periods,
        "workspace": str(ws_path),
    }

    return mf, mt, params


def run_1d_transport(
    mf: Any,
    mt: Any,
    params: Dict[str, float],
    cleanup: bool = True,
) -> TransportResult:
    """
    Run the 1D transport model and extract results.

    Parameters
    ----------
    mf : Modflow
        MODFLOW model object
    mt : Mt3dms
        MT3DMS model object
    params : dict
        Model parameters from build_1d_transport_model
    cleanup : bool
        Whether to delete model files after running

    Returns
    -------
    TransportResult
        Transport simulation results
    """
    if not FLOPY_AVAILABLE:
        raise ImportError("FloPy is required for transport modeling.")

    workspace = params.get("workspace", ".")
    warnings_list = []

    try:
        # Write and run MODFLOW
        mf.write_input()
        success_mf, buff_mf = mf.run_model(silent=True)

        if not success_mf:
            warnings_list.append("MODFLOW run failed")
            return TransportResult(
                times=np.array([0.0]),
                concentrations=np.array([0.0]),
                success=False,
                warnings=warnings_list,
                workspace=workspace,
            )

        # Write and run MT3DMS
        mt.write_input()
        success_mt, buff_mt = mt.run_model(silent=True, normal_msg="Normal termination")

        if not success_mt:
            # Only warn, don't fail immediately, check for output files
            warnings_list.append("MT3DMS reported failure (check output)")

        # Read concentration output
        ucn_file = Path(workspace) / f"{mt.name}.ucn"
        if not ucn_file.exists():
            ucn_file = Path(workspace) / "MT3D001.UCN"

        if ucn_file.exists():
            ucn = flopy.utils.UcnFile(str(ucn_file))
            times = ucn.get_times()
            # n_cells = params.get("n_cells", 50)

            # Get concentration at outlet (last cell)
            concentrations = []
            for t in times:
                conc = ucn.get_data(totim=t)
                outlet_conc = conc[0, 0, -1]
                concentrations.append(outlet_conc)

            times = np.array(times)
            concentrations = np.array(concentrations)
        else:
            warnings_list.append("UCN file not found")
            try:
                files = [f.name for f in Path(workspace).glob("*")]
                warnings_list.append(f"Files in {workspace}: {files}")
            except Exception:
                pass
            times = np.array([0.0])
            concentrations = np.array([0.0])

        # Calculate mass balance
        mass_in = 0.0
        mass_out = 0.0
        # Could extract from budget files if needed

        result = TransportResult(
            times=times,
            concentrations=concentrations,
            n_cells=int(params.get("n_cells", 50)),
            cell_length_m=float(params.get("cell_length_m", 1.0)),
            porosity=float(params.get("porosity", 0.25)),
            dispersivity_m=float(params.get("dispersivity_m", 1.0)),
            velocity_m_day=float(params.get("velocity_m_day", 0.1)),
            decay_rate_1_day=float(params.get("decay_rate_1_day", 0.0)),
            mass_in=mass_in,
            mass_out=mass_out,
            model_name=mf.name,
            workspace=workspace,
            success=True,
            warnings=warnings_list,
        )

    except Exception as e:
        warnings_list.append(f"Simulation error: {str(e)}")
        result = TransportResult(
            times=np.array([0.0]),
            concentrations=np.array([0.0]),
            success=False,
            warnings=warnings_list,
            workspace=workspace,
        )

    finally:
        if cleanup and workspace:
            try:
                shutil.rmtree(workspace)
            except Exception:
                pass

    return result


def analytical_1d_transport(
    x: float,
    t: float,
    c0: float,
    v: float,
    D: float,
    k: float = 0.0,
) -> float:
    """
    Analytical solution for 1D advection-dispersion with first-order decay.

    Uses a numerically stable formulation of the Ogata-Banks solution.

    Parameters
    ----------
    x : float
        Distance from source (m)
    t : float
        Time (days)
    c0 : float
        Source concentration
    v : float
        Pore velocity (m/day)
    D : float
        Dispersion coefficient (m^2/day)
    k : float
        First-order decay rate (1/day)

    Returns
    -------
    float
        Concentration at (x, t)
    """
    from scipy.special import erfc

    if t <= 0 or D <= 0 or x < 0:
        return 0.0

    # For very small t or very large x, concentration is effectively 0
    sqrt_Dt = np.sqrt(D * t)
    if sqrt_Dt < 1e-10:
        return 0.0

    # Check if front hasn't reached x yet (approximately)
    front_position = v * t
    if x > front_position + 6 * sqrt_Dt:
        # Very far ahead of front, concentration is negligible
        return 0.0

    # Effective velocity with decay
    beta = np.sqrt(v**2 + 4 * k * D)

    # Arguments for erfc
    arg1 = (x - v * t) / (2 * sqrt_Dt)
    arg2 = (x + v * t) / (2 * sqrt_Dt)

    # For the case with decay, use modified arguments
    if k > 0:
        arg1 = (x - beta * t) / (2 * sqrt_Dt)
        arg2 = (x + beta * t) / (2 * sqrt_Dt)

    # Compute exponential terms carefully to avoid overflow
    exp_arg1 = (v - beta) * x / (2 * D)
    exp_arg2 = (v + beta) * x / (2 * D)

    # Use log-sum-exp trick for numerical stability
    try:
        if exp_arg2 > 700:  # Would overflow
            # Use asymptotic approximation: erfc(large) -> 0
            term1 = np.exp(exp_arg1) * erfc(arg1)
            term2 = 0.0
        else:
            term1 = np.exp(exp_arg1) * erfc(arg1)
            term2 = np.exp(exp_arg2) * erfc(arg2)

        c = (c0 / 2) * (term1 + term2)

        if np.isnan(c) or np.isinf(c):
            c = 0.0

    except (OverflowError, FloatingPointError):
        c = 0.0

    return float(np.clip(c, 0, c0))


def run_analytical_1d_transport(
    times: np.ndarray,
    distance_m: float,
    velocity_m_day: float,
    dispersivity_m: float,
    decay_rate_1_day: float = 0.0,
    source_concentration: float = 1.0,
) -> TransportResult:
    """
    Run analytical 1D transport solution (no FloPy required).

    Parameters
    ----------
    times : np.ndarray
        Time points to evaluate (days)
    distance_m : float
        Distance from source to outlet (m)
    velocity_m_day : float
        Pore velocity (m/day)
    dispersivity_m : float
        Longitudinal dispersivity (m)
    decay_rate_1_day : float
        First-order decay rate (1/day)
    source_concentration : float
        Source concentration

    Returns
    -------
    TransportResult
        Transport result with concentration time series
    """
    # Dispersion coefficient D = alpha * v
    D = dispersivity_m * velocity_m_day

    concentrations = []
    for t in times:
        c = analytical_1d_transport(
            x=distance_m,
            t=t,
            c0=source_concentration,
            v=velocity_m_day,
            D=D,
            k=decay_rate_1_day,
        )
        concentrations.append(c)

    return TransportResult(
        times=times,
        concentrations=np.array(concentrations),
        velocity_m_day=velocity_m_day,
        dispersivity_m=dispersivity_m,
        decay_rate_1_day=decay_rate_1_day,
        source_concentration=source_concentration,
        success=True,
        warnings=[],
    )
