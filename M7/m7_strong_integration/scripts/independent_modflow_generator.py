"""Independent MODFLOW 6 / MODPATH 7 aquifer generator for M7.2.

This module intentionally imports no ``hydrosheaf`` code.  MODFLOW supplies
heads and cell-by-cell flow, MODPATH supplies particle routes and travel times,
and a separate nonlinear aqueous process simulator supplies chemistry and
tracers.  HydroSheaf sees only the returned observations.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
from pathlib import Path
import subprocess
from typing import Dict, List, Tuple

import flopy
import numpy as np
from scipy.ndimage import gaussian_filter


ION_ORDER = [
    "Ca",
    "Mg",
    "Na",
    "K",
    "HCO3",
    "Cl",
    "SO4",
    "NO3",
    "F",
    "Fe",
    "PO4",
    "SiO2",
]


@dataclass(frozen=True)
class IndependentAquifer:
    seed: int
    observations: Tuple[Dict[str, object], ...]
    true_edges: Tuple[Tuple[str, str], ...]
    true_ages_years: Dict[str, float]
    true_processes: Dict[str, str]
    pathline_rows: Tuple[Dict[str, object], ...]
    provenance: Dict[str, object]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _version(executable: Path) -> str:
    completed = subprocess.run(
        [str(executable), "-v"],
        check=False,
        capture_output=True,
        text=True,
    )
    text = (completed.stdout or completed.stderr).strip()
    return text.splitlines()[0] if text else "version_unavailable"


def _build_flow_model(
    seed: int,
    workspace: Path,
    mf6_executable: Path,
) -> Tuple[object, np.ndarray, Dict[str, object]]:
    rng = np.random.default_rng(seed)
    nlay, nrow, ncol = 1, 10, 20
    delr = delc = 150.0
    sim = flopy.mf6.MFSimulation(
        sim_name=f"m7flow_{seed}",
        sim_ws=workspace,
        exe_name=str(mf6_executable),
    )
    flopy.mf6.ModflowTdis(
        sim,
        time_units="DAYS",
        nper=1,
        perioddata=[(1.0, 1, 1.0)],
    )
    ims = flopy.mf6.ModflowIms(
        sim,
        complexity="MODERATE",
        outer_dvclose=1e-9,
        inner_dvclose=1e-9,
        outer_maximum=200,
        inner_maximum=300,
    )
    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname="gwf",
        save_flows=True,
    )
    sim.register_ims_package(ims, [gwf.name])
    flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=25.0,
        botm=0.0,
    )
    flopy.mf6.ModflowGwfic(gwf, strt=95.0)

    log_k = gaussian_filter(
        rng.normal(0.0, 0.65, size=(nrow, ncol)),
        sigma=1.25,
        mode="reflect",
    )
    hydraulic_k = np.clip(8.0 * np.exp(log_k), 1.0, 35.0)[None, :, :]
    # A low-K lens forces pathlines to respond to heterogeneity rather than
    # following a trivial straight analytical gradient.
    lens_row = 3 + seed % 3
    hydraulic_k[:, lens_row : lens_row + 2, 7:13] *= 0.18
    flopy.mf6.ModflowGwfnpf(
        gwf,
        icelltype=0,
        k=hydraulic_k,
        save_specific_discharge=True,
    )

    left_head = 101.0 + 0.15 * (seed % 4)
    right_head = 90.0
    constant_heads = []
    for row in range(nrow):
        constant_heads.append(((0, row, 0), left_head))
        constant_heads.append(((0, row, ncol - 1), right_head))
    flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data=constant_heads,
        save_flows=True,
    )
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord="gwf.hds",
        budget_filerecord="gwf.cbc",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )
    sim.write_simulation(silent=True)
    success, report = sim.run_simulation(silent=True, report=True)
    if not success:
        raise RuntimeError("MODFLOW 6 failed:\n" + "\n".join(report[-20:]))
    heads = flopy.utils.HeadFile(workspace / "gwf.hds").get_data()
    metadata = {
        "nlay": nlay,
        "nrow": nrow,
        "ncol": ncol,
        "delr_m": delr,
        "delc_m": delc,
        "porosity": 0.25,
        "left_head_m": left_head,
        "right_head_m": right_head,
        "hydraulic_k_min_m_day": float(hydraulic_k.min()),
        "hydraulic_k_max_m_day": float(hydraulic_k.max()),
        "heterogeneous_k": True,
        "low_k_lens": True,
    }
    return gwf, heads, metadata


def _run_pathlines(
    seed: int,
    workspace: Path,
    gwf: object,
    mp7_executable: Path,
) -> List[np.ndarray]:
    nrow = int(gwf.modelgrid.nrow)
    ncol = int(gwf.modelgrid.ncol)
    source_rows = [1, nrow // 2, nrow - 2]
    nodes = [row * ncol + 1 for row in source_rows]
    model_name = f"m7mp_{seed}"
    mp = flopy.modpath.Modpath7.create_mp7(
        modelname=model_name,
        trackdir="forward",
        flowmodel=gwf,
        model_ws=workspace,
        rowcelldivisions=1,
        columncelldivisions=1,
        layercelldivisions=1,
        nodes=nodes,
        porosity=0.25,
        exe_name=str(mp7_executable),
    )
    mp.write_input()
    success, report = mp.run_model(silent=True, report=True)
    if not success:
        raise RuntimeError("MODPATH 7 failed:\n" + "\n".join(report[-20:]))
    pathline_file = flopy.utils.PathlineFile(workspace / f"{model_name}.mppth")
    return [np.sort(path, order="time") for path in pathline_file.get_alldata()]


def _chemistry_step(
    chemistry: np.ndarray,
    travel_years: float,
    process: str,
    rng: np.random.Generator,
) -> np.ndarray:
    """Independent nonlinear process simulator in mmol/L."""

    state = chemistry.copy()
    dt = max(0.05, float(travel_years))
    index = {ion: idx for idx, ion in enumerate(ION_ORDER)}
    if process == "carbonate_weathering":
        saturation = state[index["Ca"]] / 3.2
        extent = max(0.0, 0.34 * (1.0 - np.exp(-dt / 8.0)) * (1.0 - saturation))
        state[index["Ca"]] += extent
        state[index["Mg"]] += 0.22 * extent
        state[index["HCO3"]] += 1.82 * extent + 0.05 * extent**2
    elif process == "silicate_weathering":
        extent = 0.075 * np.log1p(dt) * np.exp(-state[index["SiO2"]] / 7.0)
        state[index["Na"]] += 0.72 * extent
        state[index["K"]] += 0.18 * extent
        state[index["Ca"]] += 0.12 * extent
        state[index["HCO3"]] += 0.86 * extent
        state[index["SiO2"]] += 1.55 * extent
    elif process == "sulfate_reduction":
        fraction = min(0.72, 1.0 - np.exp(-0.055 * dt))
        consumed = state[index["SO4"]] * fraction
        state[index["SO4"]] -= consumed
        state[index["HCO3"]] += 1.72 * consumed
        state[index["NO3"]] *= np.exp(-0.11 * dt)
        state[index["Fe"]] += 0.025 * consumed
    elif process == "denitrification":
        fraction = min(0.88, 1.0 - np.exp(-0.090 * dt))
        consumed = state[index["NO3"]] * fraction
        state[index["NO3"]] -= consumed
        state[index["HCO3"]] += 0.74 * consumed
    elif process == "iron_reduction":
        redox_gate = np.exp(-5.0 * state[index["NO3"]])
        extent = 0.018 * np.log1p(dt) * redox_gate
        state[index["Fe"]] += extent
        state[index["HCO3"]] += 0.92 * extent
    elif process == "carbonate_precipitation":
        activity_proxy = state[index["Ca"]] * state[index["HCO3"]]
        extent = min(
            0.45 * state[index["Ca"]],
            max(0.0, 0.018 * dt * (activity_proxy - 0.75)),
        )
        state[index["Ca"]] -= extent
        state[index["HCO3"]] -= min(
            0.90 * state[index["HCO3"]],
            1.85 * extent,
        )
    else:
        raise ValueError(f"Unknown independent process: {process}")

    # Heteroscedastic laboratory/process error and a conservative chloride
    # perturbation prevent exact recovery by a linear reaction dictionary.
    state *= rng.lognormal(0.0, 0.018 + 0.0015 * dt, state.size)
    state[index["Cl"]] *= rng.lognormal(0.0, 0.008)
    return np.clip(state, 1e-6, None)


def _tracer_observations(
    age_years: float,
    rng: np.random.Generator,
) -> Dict[str, float]:
    # These equations are deliberately not HydroSheaf's input-history lookup.
    tritium = (
        6.2 * np.exp(-age_years / 17.0)
        + 0.50 * np.exp(-((age_years - 34.0) / 9.5) ** 2)
        + rng.normal(0.0, 0.10)
    )
    argon39 = (
        100.0
        * np.exp(-np.log(2.0) * age_years / 269.0)
        * (1.0 - 0.05 * (1.0 - np.exp(-age_years / 55.0)))
        + rng.normal(0.0, 1.2)
    )
    return {
        "tritium_TU": float(max(0.02, tritium)),
        "argon39_pmc": float(max(0.05, argon39)),
    }


def generate_independent_aquifer(
    seed: int,
    workspace: Path,
    mf6_executable: Path,
    mp7_executable: Path,
) -> IndependentAquifer:
    """Run the external simulators and return a blind HydroSheaf case."""

    workspace = Path(workspace)
    workspace.mkdir(parents=True, exist_ok=True)
    mf6_executable = Path(mf6_executable).resolve()
    mp7_executable = Path(mp7_executable).resolve()
    if not mf6_executable.exists() or not mp7_executable.exists():
        raise FileNotFoundError(
            "Both official mf6 and mp7 executables are required."
        )

    rng = np.random.default_rng(int(seed))
    gwf, heads, model_metadata = _build_flow_model(
        int(seed), workspace, mf6_executable
    )
    pathlines = _run_pathlines(
        int(seed), workspace, gwf, mp7_executable
    )
    nrow = int(model_metadata["nrow"])
    ncol = int(model_metadata["ncol"])
    delr = float(model_metadata["delr_m"])
    delc = float(model_metadata["delc_m"])
    milestone_columns = (1, 6, 11, 16)
    # Each case exercises every preregistered reaction family.  The reducing
    # sequence follows the chemically defensible electron-acceptor progression
    # nitrate -> sulfate -> ferric iron.
    process_paths = (
        (
            "carbonate_weathering",
            "carbonate_precipitation",
            "silicate_weathering",
        ),
        (
            "denitrification",
            "sulfate_reduction",
            "iron_reduction",
        ),
        (
            "silicate_weathering",
            "carbonate_weathering",
            "denitrification",
        ),
    )

    observations: List[Dict[str, object]] = []
    true_edges: List[Tuple[str, str]] = []
    true_ages: Dict[str, float] = {}
    true_processes: Dict[str, str] = {}
    pathline_rows: List[Dict[str, object]] = []
    base = np.asarray(
        [0.62, 0.27, 0.48, 0.055, 1.62, 0.24, 0.42, 0.30, 0.007, 0.003, 0.004, 0.55],
        dtype=float,
    )

    for particle_index, path in enumerate(pathlines):
        previous_node: str | None = None
        previous_time = 0.0
        chemistry = base * rng.lognormal(0.0, 0.055, base.size)
        for milestone_index, column in enumerate(milestone_columns):
            records = path[
                np.floor(np.maximum(path["x"], 0.0) / delr).astype(int)
                >= column
            ]
            if records.size == 0:
                raise RuntimeError(
                    f"Particle {particle_index} did not reach column {column}."
                )
            record = records[0]
            row = min(
                nrow - 1,
                max(
                    0,
                    nrow - 1 - int(np.floor(max(float(record["y"]), 0.0) / delc)),
                ),
            )
            actual_column = min(
                ncol - 1,
                max(0, int(np.floor(max(float(record["x"]), 0.0) / delr))),
            )
            elapsed_days = float(record["time"])
            age_years = 1.5 + elapsed_days / 365.25
            node_id = f"MF{seed}_P{particle_index}_M{milestone_index}"
            if previous_node is not None:
                process = process_paths[particle_index][milestone_index - 1]
                chemistry = _chemistry_step(
                    chemistry,
                    (elapsed_days - previous_time) / 365.25,
                    process,
                    rng,
                )
                edge_id = f"{previous_node}->{node_id}"
                true_edges.append((previous_node, node_id))
                true_processes[edge_id] = process
            else:
                process = "source"

            x = float(record["x"])
            y = float(record["y"])
            head = float(heads[0, row, actual_column])
            tracer = _tracer_observations(age_years, rng)
            p_h = 6.85 + 0.018 * age_years
            if process == "carbonate_precipitation":
                p_h += 0.55
            sample: Dict[str, object] = {
                "site_id": node_id,
                "sample_id": node_id,
                "x_m": x,
                "y_m": y,
                "lat": 7.0 + y / 111_000.0,
                "lon": -1.5 + x / 111_000.0,
                "head_meas": head + rng.normal(0.0, 0.10),
                "hydraulic_head": head + rng.normal(0.0, 0.10),
                "elevation": head + 12.0,
                "screen_depth": 18.0,
                "well_depth": 40.0,
                "aquifer_unit": "external_modflow_synthetic",
                "aquifer_layer": 1,
                "pH": float(np.clip(p_h + rng.normal(0.0, 0.04), 6.3, 8.8)),
                "temp_c": float(25.0 + rng.normal(0.0, 0.35)),
                "sample_date": 2025.5,
                **tracer,
            }
            for ion_index, ion in enumerate(ION_ORDER):
                sample[ion] = float(chemistry[ion_index])
            observations.append(sample)
            true_ages[node_id] = age_years
            pathline_rows.append(
                {
                    "particle": particle_index,
                    "milestone": milestone_index,
                    "node_id": node_id,
                    "modpath_node_zero_based": int(record["node"]),
                    "row_zero_based": row,
                    "column_zero_based": actual_column,
                    "x_m": x,
                    "y_m": y,
                    "travel_time_days": elapsed_days,
                    "age_years": age_years,
                }
            )
            previous_node = node_id
            previous_time = elapsed_days

    provenance: Dict[str, object] = {
        "generator": "independent_modflow_modpath_nonlinear_reactive_v1",
        "imports_hydrosheaf": False,
        "mf6_executable": str(mf6_executable),
        "mf6_sha256": _sha256(mf6_executable),
        "mf6_version": _version(mf6_executable),
        "mp7_executable": str(mp7_executable),
        "mp7_sha256": _sha256(mp7_executable),
        "mp7_version": _version(mp7_executable),
        "flopy_version": flopy.__version__,
        "chemistry_generator": "independent nonlinear process equations",
        "tracer_generator": "independent analytic 3H/39Ar surrogate",
        **model_metadata,
    }
    return IndependentAquifer(
        seed=int(seed),
        observations=tuple(observations),
        true_edges=tuple(true_edges),
        true_ages_years=true_ages,
        true_processes=true_processes,
        pathline_rows=tuple(pathline_rows),
        provenance=provenance,
    )
