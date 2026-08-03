"""Independent analytic-lattice groundwater generator for programme validation.

This family deliberately does not import HydroSheaf, MODFLOW, MODPATH, or the
M7 generator.  It uses a curvilinear three-path geometry, an analytic head
field, segment-wise mass-balance chemistry, and a separate random stream.  The
returned observations share the public HydroSheaf schema, while all truth is
kept in the generator record and never placed in an observation row.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np


ION_ORDER = (
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
)

NODE_COORDINATES: Tuple[Tuple[str, float, float], ...] = (
    ("R", 0.0, 700.0),
    ("U1", 320.0, 850.0),
    ("L1", 300.0, 550.0),
    ("N1", 450.0, 700.0),
    ("U2", 650.0, 1_000.0),
    ("L2", 650.0, 380.0),
    ("N2", 900.0, 700.0),
    ("U3", 1_000.0, 1_100.0),
    ("L3", 1_000.0, 300.0),
    ("U4", 1_350.0, 980.0),
    ("L4", 1_350.0, 420.0),
    ("M", 1_750.0, 700.0),
)
TRUE_EDGE_PROCESSES: Tuple[Tuple[str, str, str], ...] = (
    ("R", "U1", "carbonate_weathering"),
    ("U1", "U2", "silicate_weathering"),
    ("U2", "U3", "carbonate_precipitation"),
    ("U3", "U4", "denitrification"),
    ("U4", "M", "mixing"),
    ("R", "L1", "sulfate_reduction"),
    ("L1", "L2", "iron_reduction"),
    ("L2", "L3", "silicate_weathering"),
    ("L3", "L4", "carbonate_weathering"),
    ("L4", "M", "mixing"),
    ("R", "N1", "denitrification"),
    ("N1", "N2", "sulfate_reduction"),
    ("N2", "M", "mixing"),
)
AGE_BY_NODE = {
    "R": 3.5,
    "U1": 11.0,
    "L1": 12.0,
    "N1": 14.0,
    "U2": 22.0,
    "L2": 25.0,
    "N2": 31.0,
    "U3": 38.0,
    "L3": 41.0,
    "U4": 57.0,
    "L4": 60.0,
    "M": 84.0,
}
BASE_CHEMISTRY = np.asarray(
    [1.15, 0.31, 0.52, 0.065, 1.75, 0.28, 0.43, 0.48, 0.012, 0.004, 0.006, 0.58],
    dtype=float,
)


@dataclass(frozen=True)
class IndependentLatticeAquifer:
    """Blind observations plus withheld truth for one analytic-lattice case."""

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


def _head(x_m: float, y_m: float, path_index: int, milestone: int) -> float:
    """Analytic head field with a monotone primary gradient and curvature."""

    return float(
        121.0
        - 0.0048 * float(x_m)
        - 0.0024 * float(y_m)
        + 0.25 * np.sin(float(x_m) / 470.0 + path_index * 0.7)
        - 0.08 * milestone
    )


def _advance_chemistry(
    chemistry: np.ndarray,
    process: str,
    rng: np.random.Generator,
) -> np.ndarray:
    """Apply an independent segment mass-balance rule."""

    updated = np.asarray(chemistry, dtype=float).copy()
    effects = {
        "carbonate_weathering": (0.16, 0.04, 0.05, 0.12, 0.0, 0.0),
        "silicate_weathering": (0.04, 0.02, 0.13, 0.10, 0.0, 0.16),
        "carbonate_precipitation": (-0.10, -0.02, 0.0, -0.08, 0.0, -0.02),
        "denitrification": (0.02, 0.0, 0.01, 0.04, -0.20, 0.0),
        "sulfate_reduction": (0.04, 0.01, 0.0, 0.06, -0.16, 0.0),
        "iron_reduction": (0.03, 0.01, 0.0, 0.03, 0.0, 0.02),
    }
    ca, mg, na, hco3, nitrate, sio2 = effects.get(process, (0.0,) * 6)
    updated[0] += ca
    updated[1] += mg
    updated[2] += na
    updated[4] += hco3
    updated[7] = max(0.02, updated[7] + nitrate)
    updated[11] += sio2
    if process == "iron_reduction":
        updated[9] += 0.018
        updated[6] = max(0.02, updated[6] - 0.04)
    if process == "sulfate_reduction":
        updated[6] = max(0.02, updated[6] - 0.18)
    noise = rng.normal(0.0, 0.008, size=updated.size)
    return np.maximum(0.001, updated * (1.0 + noise))


def _tracer_panel(
    age_years: float,
    process_history: Tuple[str, ...],
    rng: np.random.Generator,
) -> Dict[str, float]:
    """Generate smooth analytic age tracers with a different forward model."""

    age = float(age_years)
    reducing_steps = sum(
        process in {"denitrification", "sulfate_reduction", "iron_reduction"}
        for process in process_history
    )
    tritium = 0.15 + 6.4 * np.exp(-age / 17.0) + 0.10 * np.cos(age / 8.0)
    argon39 = 97.0 * np.exp(-age / 330.0)
    c14 = 96.0 * np.exp(-np.log(2.0) * age / 5_730.0)
    c14 *= np.exp(-0.055 * reducing_steps)
    cfc12 = 520.0 * np.exp(-age / 38.0) / (1.0 + 0.006 * age)
    sf6 = 9.0 * np.exp(-age / 42.0)
    he4 = 4.2e-8 + 1.15e-10 * age
    he3 = 0.86 * max(0.0, tritium) * (
        np.exp(np.log(2.0) * min(age, 70.0) / 12.32) - 1.0
    )
    values = {
        "tritium_TU": tritium,
        "argon39_pmc": argon39,
        "c14_pmc": c14,
        "cfc11_pptv": 0.72 * cfc12,
        "cfc12_pptv": cfc12,
        "cfc113_pptv": 0.16 * cfc12,
        "sf6_pptv": sf6,
        "he4_ccpg": he4,
        "h3_he3_TU": he3,
    }
    noisy: Dict[str, float] = {}
    for name, value in values.items():
        scale = max(abs(float(value)) * 0.025, 0.01)
        noisy[name] = float(max(0.0, float(value) + rng.normal(0.0, scale)))
    noisy["tritium_sigma_TU"] = 0.15
    noisy["argon39_sigma_pmc"] = 2.0
    return noisy


def generate_independent_lattice(seed: int, workspace: Path | None = None) -> IndependentLatticeAquifer:
    """Generate one independent analytic-lattice case without external solvers."""

    del workspace  # kept in the signature so ensemble runners can share an adapter
    seed = int(seed)
    rng = np.random.default_rng(seed * 7_919 + 17)
    observations: List[Dict[str, object]] = []
    true_edges = [
        (f"AL{seed}_{source}", f"AL{seed}_{target}")
        for source, target, _ in TRUE_EDGE_PROCESSES
    ]
    true_ages: Dict[str, float] = {}
    true_processes = {
        f"AL{seed}_{source}->AL{seed}_{target}": process
        for source, target, process in TRUE_EDGE_PROCESSES
    }
    pathline_rows: List[Dict[str, object]] = []
    parent_edges: Dict[str, list[tuple[str, str]]] = {}
    for source, target, process in TRUE_EDGE_PROCESSES:
        parent_edges.setdefault(target, []).append((source, process))
    chemistry_by_node: Dict[str, np.ndarray] = {}
    history_by_node: Dict[str, Tuple[str, ...]] = {}

    for node_index, (node_label, x_m, y_m) in enumerate(NODE_COORDINATES):
        node_id = f"AL{seed}_{node_label}"
        parents = parent_edges.get(node_label, [])
        if not parents:
            chemistry = BASE_CHEMISTRY * rng.lognormal(0.0, 0.035, BASE_CHEMISTRY.size)
            process = "source"
            process_history: Tuple[str, ...] = ()
        elif len(parents) == 1:
            parent, process = parents[0]
            chemistry = _advance_chemistry(chemistry_by_node[parent], process, rng)
            process_history = history_by_node[parent] + (process,)
        else:
            chemistry = np.mean(
                [chemistry_by_node[parent] for parent, _ in parents], axis=0
            )
            chemistry = chemistry * (1.0 + rng.normal(0.0, 0.012, chemistry.size))
            process = "mixing"
            longest_history = max(
                (history_by_node[parent] for parent, _ in parents),
                key=len,
            )
            process_history = longest_history + (process,)
        chemistry_by_node[node_label] = np.maximum(0.001, chemistry)
        history_by_node[node_label] = process_history

        age = float(AGE_BY_NODE[node_label] + rng.normal(0.0, 0.18))
        if age <= 0.0:
            raise RuntimeError(f"Analytic-lattice node {node_label} has non-positive age.")
        true_ages[node_id] = age
        branch_index = 0 if node_label.startswith("U") else 1 if node_label.startswith("L") else 2
        head_true = _head(x_m, y_m, branch_index, node_index)
        head = float(head_true + rng.normal(0.0, 0.06))
        hydraulic_head = float(
            head_true
            + 0.32 * np.sin(x_m / 175.0 + y_m / 240.0)
            + rng.normal(0.0, 0.14)
        )
        d18o = -4.8 - 0.0007 * y_m + 0.05 * branch_index
        d2h = 8.0 * d18o + 10.0 + rng.normal(0.0, 0.45)
        reaction_steps = len(process_history)
        d13c = -18.5 + 0.55 * reaction_steps + (0.9 if "carbonate" in process else 0.0)
        sr87 = 0.7090 + 0.00075 * reaction_steps + 0.0002 * branch_index
        p_h = 6.72 + 0.014 * age + (0.28 if process == "carbonate_precipitation" else 0.0)
        aquifer_layer = 2 if node_label.startswith("U") else 1
        sample: Dict[str, object] = {
            "site_id": node_id,
            "sample_id": node_id,
            "x_m": float(x_m),
            "y_m": float(y_m),
            "lat": 7.2 + float(y_m) / 111_000.0,
            "lon": -1.8 + float(x_m) / 111_000.0,
            "head_meas": head,
            "hydraulic_head": hydraulic_head,
            "elevation": head_true + 15.0,
            "screen_depth": 20.0 + 1.5 * (aquifer_layer - 1),
            "well_depth": 48.0 + 0.4 * node_index,
            "aquifer_unit": "analytic_lattice_external",
            "aquifer_layer": aquifer_layer,
            "pH": float(np.clip(p_h + rng.normal(0.0, 0.035), 6.3, 8.8)),
            "temp_c": float(24.2 + rng.normal(0.0, 0.25)),
            "sample_date": 2025.5,
            "d18O_permil": float(d18o + rng.normal(0.0, 0.06)),
            "d2H_permil": float(d2h),
            "d13C_permil": float(d13c + rng.normal(0.0, 0.12)),
            "sr87_sr86": float(sr87 + rng.normal(0.0, 0.00003)),
            **_tracer_panel(age, process_history, rng),
        }
        for ion, value in zip(ION_ORDER, chemistry_by_node[node_label]):
            sample[ion] = float(max(0.001, value))
        observations.append(sample)
        pathline_rows.append(
            {
                "node_id": node_id,
                "node_label": node_label,
                "parent_count": len(parents),
                "travel_age_years": age,
                "x_m": float(x_m),
                "y_m": float(y_m),
                "head": head_true,
            }
        )

    provenance: Dict[str, object] = {
        "generator": "independent_analytic_lattice_reactive_v1",
        "generator_family": "analytic_lattice",
        "imports_hydrosheaf": False,
        "geometry": "three_branching_curvilinear_paths_with_terminal_merge",
        "head_forward_model": "analytic_curved_gradient",
        "chemistry_forward_model": "segment_mass_balance_increment",
        "tracer_forward_model": "smooth_exponential_pulse_surrogates",
        "random_stream": "PCG64 seed_affine_transform_7919_plus_17",
        "head_channel_relation": "same_latent_head_with_local_bias_and_independent_noise",
        "head_channel_covariance_model": {
            "head_meas_sigma_m": 0.06,
            "hydraulic_head_sigma_m": 0.14,
            "measurement_error_correlation": 0.0,
            "hydraulic_head_local_bias_amplitude_m": 0.32,
            "channels_are_independent": False,
            "combination": "primary_only_with_discrepancy",
        },
        "observation_stress_scenarios": ["structured_missingness", "left_censoring", "tracer_permutation"],
        "generator_source_sha256": _sha256(Path(__file__).resolve()),
        "external_simulators": [],
    }
    return IndependentLatticeAquifer(
        seed=seed,
        observations=tuple(observations),
        true_edges=tuple(true_edges),
        true_ages_years=true_ages,
        true_processes=true_processes,
        pathline_rows=tuple(pathline_rows),
        provenance=provenance,
    )


__all__ = ["IndependentLatticeAquifer", "generate_independent_lattice"]
