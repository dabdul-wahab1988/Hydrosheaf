"""Independent multilayer mixing generator for held-out validation.

This family is deliberately separate from the MODFLOW/MODPATH and analytic
lattice generators.  It uses a hand-specified three-layer directed network,
conservative endmember mixing, process increments, and a two-component age
mixture.  It imports no HydroSheaf or M7 generator code.
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

NODE_COORDINATES: Tuple[Tuple[str, float, float, int], ...] = (
    ("R", 0.0, 700.0, 1),
    ("A1", 280.0, 900.0, 2),
    ("B1", 280.0, 500.0, 1),
    ("C1", 280.0, 700.0, 3),
    ("A2", 620.0, 960.0, 2),
    ("B2", 620.0, 440.0, 1),
    ("C2", 620.0, 700.0, 3),
    ("A3", 900.0, 940.0, 2),
    ("B3", 900.0, 460.0, 1),
    ("M", 1_180.0, 700.0, 2),
    ("OUT", 1_520.0, 700.0, 2),
)
TRUE_EDGE_PROCESSES: Tuple[Tuple[str, str, str], ...] = (
    ("R", "A1", "carbonate_weathering"),
    ("A1", "A2", "silicate_weathering"),
    ("A2", "A3", "silicate_weathering"),
    ("A3", "M", "mixing"),
    ("R", "B1", "sulfate_reduction"),
    ("B1", "B2", "iron_reduction"),
    ("B2", "B3", "carbonate_weathering"),
    ("B3", "M", "mixing"),
    ("R", "C1", "denitrification"),
    ("C1", "C2", "sulfate_reduction"),
    ("A2", "C2", "mixing"),
    ("C2", "M", "mixing"),
    ("B2", "M", "mixing"),
    ("M", "OUT", "carbonate_precipitation"),
)
AGE_BY_NODE = {
    "R": 2.5,
    "A1": 8.0,
    "B1": 10.0,
    "C1": 7.0,
    "A2": 17.0,
    "B2": 21.0,
    "C2": 28.0,
    "A3": 32.0,
    "B3": 35.0,
    "M": 48.0,
    "OUT": 66.0,
}
BASE_CHEMISTRY = np.asarray(
    [1.05, 0.28, 0.50, 0.06, 1.55, 0.30, 0.52, 0.62, 0.014, 0.006, 0.008, 0.42],
    dtype=float,
)


@dataclass(frozen=True)
class IndependentMixingAquifer:
    """Blind observations and withheld truth for one mixing case."""

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


def _head(x_m: float, y_m: float, layer: int, node_index: int) -> float:
    return float(
        138.0
        - 0.0042 * float(x_m)
        - 0.65 * float(layer)
        - 0.00025 * float(y_m)
        + 0.04 * np.sin(float(x_m) / 270.0 + node_index * 0.3)
    )


def _advance_chemistry(
    chemistry: np.ndarray,
    process: str,
    rng: np.random.Generator,
) -> np.ndarray:
    updated = np.asarray(chemistry, dtype=float).copy()
    effects = {
        "carbonate_weathering": (0.12, 0.03, 0.03, 0.10, 0.0, 0.0),
        "silicate_weathering": (0.03, 0.02, 0.14, 0.06, 0.0, 0.12),
        "carbonate_precipitation": (-0.09, -0.02, 0.0, -0.07, 0.0, -0.01),
        "denitrification": (0.01, 0.0, 0.01, 0.03, -0.22, 0.0),
        "sulfate_reduction": (0.02, 0.01, 0.0, 0.04, -0.12, 0.0),
        "iron_reduction": (0.02, 0.01, 0.0, 0.02, 0.0, 0.01),
    }
    ca, mg, na, hco3, nitrate, sio2 = effects.get(process, (0.0,) * 6)
    updated[0] += ca
    updated[1] += mg
    updated[2] += na
    updated[4] += hco3
    updated[7] = max(0.02, updated[7] + nitrate)
    updated[11] += sio2
    if process == "iron_reduction":
        updated[9] += 0.016
        updated[6] = max(0.02, updated[6] - 0.03)
    if process == "sulfate_reduction":
        updated[6] = max(0.02, updated[6] - 0.13)
    return np.maximum(0.001, updated * (1.0 + rng.normal(0.0, 0.01, updated.size)))


def _tracer_panel(
    age_years: float,
    process_history: Tuple[str, ...],
    rng: np.random.Generator,
) -> Dict[str, float]:
    """Generate a two-component mixture rather than a single exponential clock."""

    age = float(age_years)
    redox_steps = sum(
        process in {"denitrification", "sulfate_reduction", "iron_reduction"}
        for process in process_history
    )
    young_fraction = float(np.clip(0.72 - 0.004 * age, 0.18, 0.72))
    young = np.exp(-age / 16.0)
    old = np.exp(-age / 270.0)
    mixture = young_fraction * young + (1.0 - young_fraction) * old
    tritium = 0.12 + 6.0 * (young_fraction * young + 0.025 * old)
    argon39 = 96.0 * mixture
    c14 = 98.0 * (0.35 * young + 0.65 * old) * np.exp(-0.045 * redox_steps)
    cfc12 = 520.0 * (0.78 * young + 0.22 * old)
    sf6 = 10.0 * (0.86 * young + 0.14 * old)
    he4 = 4.2e-8 + 1.35e-10 * age + 0.4e-8 * (1.0 - young_fraction)
    he3 = 0.82 * max(0.0, tritium) * (
        np.exp(np.log(2.0) * min(age, 70.0) / 12.32) - 1.0
    )
    values = {
        "tritium_TU": tritium,
        "argon39_pmc": argon39,
        "c14_pmc": c14,
        "cfc11_pptv": 0.68 * cfc12,
        "cfc12_pptv": cfc12,
        "cfc113_pptv": 0.15 * cfc12,
        "sf6_pptv": sf6,
        "he4_ccpg": he4,
        "h3_he3_TU": he3,
    }
    noisy: Dict[str, float] = {}
    for name, value in values.items():
        scale = max(abs(float(value)) * 0.022, 0.01)
        noisy[name] = float(max(0.0, float(value) + rng.normal(0.0, scale)))
    noisy["tritium_sigma_TU"] = 0.15
    noisy["argon39_sigma_pmc"] = 2.0
    return noisy


def generate_independent_mixing(
    seed: int,
    workspace: Path | None = None,
) -> IndependentMixingAquifer:
    """Generate one multilayer mixing case without external simulators."""

    del workspace
    seed = int(seed)
    rng = np.random.default_rng(seed * 12_347 + 53)
    observations: List[Dict[str, object]] = []
    true_edges = [
        (f"IM{seed}_{source}", f"IM{seed}_{target}")
        for source, target, _ in TRUE_EDGE_PROCESSES
    ]
    true_ages = {
        f"IM{seed}_{node}": float(age + rng.normal(0.0, 0.14))
        for node, age in AGE_BY_NODE.items()
    }
    true_processes = {
        f"IM{seed}_{source}->IM{seed}_{target}": process
        for source, target, process in TRUE_EDGE_PROCESSES
    }
    parent_edges: Dict[str, list[tuple[str, str]]] = {}
    for source, target, process in TRUE_EDGE_PROCESSES:
        parent_edges.setdefault(target, []).append((source, process))
    chemistry_by_node: Dict[str, np.ndarray] = {}
    history_by_node: Dict[str, Tuple[str, ...]] = {}
    pathline_rows: List[Dict[str, object]] = []

    for node_index, (node_label, x_m, y_m, layer) in enumerate(NODE_COORDINATES):
        node_id = f"IM{seed}_{node_label}"
        parents = parent_edges.get(node_label, [])
        if not parents:
            chemistry = BASE_CHEMISTRY * rng.lognormal(0.0, 0.04, BASE_CHEMISTRY.size)
            process = "source"
            process_history: Tuple[str, ...] = ()
        elif len(parents) == 1:
            parent, process = parents[0]
            chemistry = _advance_chemistry(chemistry_by_node[parent], process, rng)
            process_history = history_by_node[parent] + (process,)
        else:
            chemistry = np.average(
                [chemistry_by_node[parent] for parent, _ in parents],
                axis=0,
                weights=np.arange(1, len(parents) + 1, dtype=float),
            )
            chemistry = np.maximum(
                0.001,
                chemistry * (1.0 + rng.normal(0.0, 0.014, chemistry.size)),
            )
            process = "mixing"
            longest_history = max(
                (history_by_node[parent] for parent, _ in parents),
                key=len,
            )
            process_history = longest_history + (process,)
        chemistry_by_node[node_label] = np.maximum(0.001, chemistry)
        history_by_node[node_label] = process_history

        age = float(true_ages[node_id])
        head_true = _head(x_m, y_m, layer, node_index)
        head = float(head_true + rng.normal(0.0, 0.08))
        hydraulic_head = float(
            head_true
            + 0.22 * np.sin(x_m / 210.0 + layer)
            + 0.08 * (layer - 2)
            + rng.normal(0.0, 0.12)
        )
        branch_index = 0 if node_label.startswith("A") else 1 if node_label.startswith("B") else 2
        d18o = -4.4 - 0.0005 * y_m - 0.12 * layer + 0.06 * branch_index
        d2h = 8.0 * d18o + 10.0 + rng.normal(0.0, 0.40)
        reaction_steps = len(process_history)
        d13c = -19.0 + 0.48 * reaction_steps + (0.8 if "carbonate" in process else 0.0)
        sr87 = 0.7088 + 0.00065 * reaction_steps + 0.00025 * layer
        p_h = 6.70 + 0.016 * age + (0.22 if process == "carbonate_precipitation" else 0.0)
        sample: Dict[str, object] = {
            "site_id": node_id,
            "sample_id": node_id,
            "x_m": float(x_m),
            "y_m": float(y_m),
            "lat": 7.1 + float(y_m) / 111_000.0,
            "lon": -1.7 + float(x_m) / 111_000.0,
            "head_meas": head,
            "hydraulic_head": hydraulic_head,
            "elevation": head_true + 13.0,
            "screen_depth": 18.0 + 2.0 * layer,
            "well_depth": 42.0 + 0.8 * node_index,
            "aquifer_unit": "independent_mixing_external",
            "aquifer_layer": layer,
            "pH": float(np.clip(p_h + rng.normal(0.0, 0.04), 6.2, 8.9)),
            "temp_c": float(24.6 + rng.normal(0.0, 0.22)),
            "sample_date": 2025.5,
            "d18O_permil": float(d18o + rng.normal(0.0, 0.07)),
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
                "layer": layer,
                "head": head_true,
            }
        )

    provenance: Dict[str, object] = {
        "generator": "independent_multilayer_mixing_ttd_v1",
        "generator_family": "independent_mixing",
        "imports_hydrosheaf": False,
        "geometry": "three_layer_branching_shortcut_and_terminal_merge",
        "head_forward_model": "piecewise_layer_gradient_with_structured_bias",
        "chemistry_forward_model": "weighted_endmember_mixing_plus_process_increments",
        "tracer_forward_model": "two_component_mixture_age_surrogates",
        "random_stream": "PCG64 seed_affine_transform_12347_plus_53",
        "head_channel_relation": "shared_latent_head_with_structured_secondary_bias",
        "head_channel_covariance_model": {
            "head_meas_sigma_m": 0.08,
            "hydraulic_head_sigma_m": 0.12,
            "measurement_error_correlation": 0.0,
            "hydraulic_head_structured_bias_amplitude_m": 0.22,
            "channels_are_independent": False,
            "combination": "primary_only_with_discrepancy",
        },
        "observation_stress_scenarios": [
            "structured_missingness",
            "left_censoring",
            "combined_stress",
        ],
        "generator_source_sha256": _sha256(Path(__file__).resolve()),
        "external_simulators": [],
    }
    return IndependentMixingAquifer(
        seed=seed,
        observations=tuple(observations),
        true_edges=tuple(true_edges),
        true_ages_years=true_ages,
        true_processes=true_processes,
        pathline_rows=tuple(pathline_rows),
        provenance=provenance,
    )


__all__ = ["IndependentMixingAquifer", "generate_independent_mixing"]
