"""Independent synthetic-aquifer truth generator for the blind M7.1 benchmark.

This module intentionally imports no ``hydrosheaf`` code.  The inference module
receives only ``observations``; hidden topology, ages, and process labels remain
in the evaluation layer.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Tuple

import numpy as np


ION_ORDER = ["Ca", "Mg", "Na", "K", "HCO3", "Cl", "SO4", "NO3", "F", "Fe", "PO4"]


@dataclass(frozen=True)
class BlindAquifer:
    seed: int
    archetype: str
    age_regime: str
    observations: Tuple[Dict[str, object], ...]
    true_edges: Tuple[Tuple[str, str], ...]
    true_ages_years: Dict[str, float]
    true_processes: Dict[str, str]


_PROCESS_VECTORS = {
    # Deliberately perturbed, non-stoichiometric response vectors.  These are
    # not copied from HydroSheaf's mineral/reaction dictionary.
    "carbonate_weathering": np.array(
        [0.43, 0.12, 0.01, 0.00, 0.79, 0.01, 0.02, -0.01, 0.000, 0.000, 0.000]
    ),
    "silicate_weathering": np.array(
        [0.05, 0.06, 0.31, 0.07, 0.38, 0.01, 0.01, -0.02, 0.004, 0.003, 0.001]
    ),
    "sulfate_reduction": np.array(
        [0.02, 0.00, 0.01, 0.00, 0.24, 0.00, -0.37, -0.05, 0.000, 0.010, 0.002]
    ),
}


def _truth_topology(archetype: str) -> List[Tuple[int, int]]:
    base = [(0, 1), (1, 2), (2, 3), (1, 4), (4, 5), (5, 6), (2, 7), (7, 8), (8, 9)]
    if archetype == "convergent":
        return base + [(5, 8)]
    if archetype == "leaky":
        return base + [(4, 7)]
    return base


def _stable_process(edge_index: int, archetype: str) -> str:
    names = ("carbonate_weathering", "silicate_weathering", "sulfate_reduction")
    offset = {"branching": 0, "convergent": 1, "leaky": 2}[archetype]
    return names[(edge_index + offset) % len(names)]


def generate_blind_aquifer(seed: int) -> BlindAquifer:
    """Generate one replicated aquifer with hidden joint truth.

    The chemistry contains nonlinear concentration response, unmodelled
    dilution/mixing, heteroscedastic error, and an occasional contamination
    pulse.  This prevents an inverse crime against the linear inverse model.
    """

    rng = np.random.default_rng(int(seed))
    archetypes = ("branching", "convergent", "leaky")
    archetype = archetypes[int(seed) % len(archetypes)]
    edges_idx = _truth_topology(archetype)
    node_ids = [f"A{seed}_N{i:02d}" for i in range(10)]

    # Coordinates are genuine degrees; spacing is roughly 0.8--1.1 km.
    x = np.array([0, 1, 2, 3, 2, 3, 4, 3, 4, 5], dtype=float)
    y = np.array([0, 0, 0, 0, 1, 1, 1, -1, -1, -1], dtype=float)
    x += rng.normal(0.0, 0.06, size=10)
    y += rng.normal(0.0, 0.06, size=10)
    lat = 7.10 + y * 0.0082
    lon = -1.45 + x * 0.0082

    parents: Dict[int, List[int]] = {i: [] for i in range(10)}
    for u, v in edges_idx:
        parents[v].append(u)

    ages = np.zeros(10, dtype=float)
    heads = np.zeros(10, dtype=float)
    chemistry = np.zeros((10, len(ION_ORDER)), dtype=float)
    d18o = np.zeros(10, dtype=float)
    d2h = np.zeros(10, dtype=float)
    chemistry[0] = np.array(
        [0.55, 0.24, 0.48, 0.06, 1.55, 0.22, 0.16, 0.31, 0.006, 0.003, 0.004]
    ) * rng.lognormal(0.0, 0.08, len(ION_ORDER))
    heads[0] = 128.0 + rng.normal(0.0, 1.0)
    ages[0] = 2.5 + rng.uniform(0.0, 2.0)
    d18o[0] = -4.8 + rng.normal(0.0, 0.08)
    d2h[0] = 7.6 * d18o[0] + 8.0 + rng.normal(0.0, 0.5)

    process_by_edge: Dict[str, str] = {}
    edge_counter = 0
    for v in range(1, 10):
        ps = parents[v]
        if not ps:
            ps = [max(0, v - 1)]
        weights = rng.dirichlet(np.ones(len(ps)) * 4.0)
        parent_chem = sum(w * chemistry[p] for w, p in zip(weights, ps))
        parent_age = float(sum(w * ages[p] for w, p in zip(weights, ps)))
        parent_head = float(sum(w * heads[p] for w, p in zip(weights, ps)))
        parent_d18o = float(sum(w * d18o[p] for w, p in zip(weights, ps)))

        travel = rng.uniform(2.0, 9.0)
        ages[v] = parent_age + travel
        heads[v] = parent_head - rng.uniform(2.5, 7.0)
        process = _stable_process(edge_counter, archetype)
        response = _PROCESS_VECTORS[process].copy()
        extent = rng.uniform(0.12, 0.55)

        # Independent nonlinear response plus an unmodelled dilute recharge
        # fraction.  The quadratic term is small but systematic.
        nonlinear = extent * response + 0.08 * extent**2 * np.sign(response) * response**2
        dilution = rng.uniform(0.00, 0.10)
        recharge = chemistry[0] * rng.lognormal(0.0, 0.05, len(ION_ORDER))
        chemistry[v] = (1.0 - dilution) * (parent_chem + nonlinear) + dilution * recharge
        chemistry[v] *= rng.lognormal(0.0, 0.025 + 0.003 * v, len(ION_ORDER))
        chemistry[v] = np.clip(chemistry[v], 1.0e-5, None)

        d18o[v] = parent_d18o + rng.normal(0.025, 0.07)
        d2h[v] = 7.6 * d18o[v] + 8.0 + rng.normal(0.0, 0.7)
        for p in ps:
            edge_id = f"{node_ids[p]}->{node_ids[v]}"
            process_by_edge[edge_id] = process
            edge_counter += 1

    age_regime = "old_tracer_uninformative" if int(seed) % 4 == 0 else "tracer_informative"
    if age_regime == "old_tracer_uninformative":
        ages += 75.0

    # An off-model nitrate pulse challenges chemistry-only ranking in a subset
    # of aquifers without changing the hidden topology.
    if int(seed) % 5 == 0:
        chemistry[7, ION_ORDER.index("NO3")] += rng.uniform(0.25, 0.55)

    observations: List[Dict[str, object]] = []
    for i, node_id in enumerate(node_ids):
        observed_head = heads[i] + rng.normal(0.0, 1.2)
        # A smooth, independent tracer surrogate: informative through ~55 y
        # but close to the detection floor for old water.
        tritium = 7.0 * np.exp(-ages[i] / 18.0) + 0.35 * np.exp(
            -((ages[i] - 32.0) / 8.0) ** 2
        )
        tritium = max(0.02, float(tritium + rng.normal(0.0, 0.10)))
        row: Dict[str, object] = {
            "site_id": node_id,
            "sample_id": node_id,
            "lat": float(lat[i]),
            "lon": float(lon[i]),
            "elevation": float(observed_head + 10.0 + rng.normal(0.0, 0.5)),
            "head_meas": float(observed_head),
            "hydraulic_head": float(observed_head),
            "screen_depth": float(18.0 + 2.0 * y[i] + rng.normal(0.0, 1.5)),
            "well_depth": float(45.0 + 3.0 * y[i] + rng.normal(0.0, 2.0)),
            "aquifer_unit": "synthetic_blind",
            "aquifer_layer": 1,
            "d18O": float(d18o[i] + rng.normal(0.0, 0.08)),
            "d2H": float(d2h[i] + rng.normal(0.0, 0.6)),
            "tritium_TU": tritium,
            "sample_date": 2025.5,
            "pH": float(7.1 + 0.10 * i + rng.normal(0.0, 0.08)),
            "temp_c": float(25.0 + rng.normal(0.0, 0.7)),
        }
        for ion_index, ion in enumerate(ION_ORDER):
            value = float(chemistry[i, ion_index])
            # Blind missingness affects only non-core diagnostic ions.
            if ion in {"F", "Fe", "PO4"} and rng.random() < 0.08:
                row[ion] = None
            else:
                row[ion] = value
        observations.append(row)

    true_edges = tuple((node_ids[u], node_ids[v]) for u, v in edges_idx)
    return BlindAquifer(
        seed=int(seed),
        archetype=archetype,
        age_regime=age_regime,
        observations=tuple(observations),
        true_edges=true_edges,
        true_ages_years={node_ids[i]: float(ages[i]) for i in range(10)},
        true_processes=process_by_edge,
    )


def poison_truth(case: BlindAquifer) -> BlindAquifer:
    """Return the same observations with deliberately false hidden truth.

    Used as a leakage test: inference outputs must remain byte-for-byte equal.
    """

    reversed_edges = tuple((v, u) for u, v in reversed(case.true_edges))
    return BlindAquifer(
        seed=case.seed,
        archetype="poisoned",
        age_regime="poisoned",
        observations=case.observations,
        true_edges=reversed_edges,
        true_ages_years={key: 999.0 - value for key, value in case.true_ages_years.items()},
        true_processes={key: "poisoned" for key in case.true_processes},
    )
