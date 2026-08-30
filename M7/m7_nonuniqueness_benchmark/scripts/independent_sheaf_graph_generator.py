"""Independent scalar-section generator for the M7.4 comparator benchmark.

The generator deliberately imports no HydroSheaf package code.  It creates
directed candidate graphs, noisy and partially missing node observations, and
edge-specific affine maps.  Ground truth is returned separately from the
inference-facing records so the benchmark runner can fit on development cases
and join locked-test truth only during scoring.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Sequence

import numpy as np


SCENARIOS = (
    "identity_limit",
    "heterogeneous_affine",
    "incompatible_cycles",
    "noisy_missing",
)


@dataclass(frozen=True)
class CandidateMap:
    edge_id: str
    u: str
    v: str
    prior_probability: float
    alpha: float
    offset: float
    is_true_edge: int
    is_corrupted_edge: int


@dataclass(frozen=True)
class IndependentSectionCase:
    seed: int
    scenario: str
    node_ids: tuple[str, ...]
    observations: Mapping[str, float | None]
    latent_truth: Mapping[str, float]
    candidates: tuple[CandidateMap, ...]
    true_edge_ids: frozenset[str]
    corrupted_edge_ids: frozenset[str]
    provenance: Mapping[str, object]


def scenario_for_seed(seed: int, first_seed: int) -> str:
    """Assign consecutive seeds to balanced, deterministic scenario strata."""

    return SCENARIOS[(int(seed) - int(first_seed)) % len(SCENARIOS)]


def _true_map(rng: np.random.Generator, scenario: str) -> tuple[float, float]:
    if scenario == "identity_limit":
        return 1.0, 0.0
    alpha = float(rng.choice(np.asarray([0.82, 0.92, 1.08, 1.18])))
    offset = float(rng.choice(np.asarray([1.5, 3.0, 4.5, 6.0])))
    return alpha, offset


def _false_map(rng: np.random.Generator, scenario: str) -> tuple[float, float]:
    if scenario == "identity_limit":
        return 1.0, 0.0
    alpha = float(rng.choice(np.asarray([0.78, 0.88, 1.12, 1.22])))
    offset = float(rng.choice(np.asarray([-5.0, -2.5, 2.5, 5.0, 7.5])))
    return alpha, offset


def _edge_id(u: str, v: str) -> str:
    return f"{u}->{v}"


def _draw_prior(rng: np.random.Generator, is_true: bool) -> float:
    # Deliberately overlapping distributions: the prior is useful but weak.
    value = rng.beta(5.0, 4.0) if is_true else rng.beta(4.0, 5.0)
    return float(np.clip(value, 0.05, 0.95))


def _candidate_pairs(
    rng: np.random.Generator,
    node_ids: Sequence[str],
    parents: Mapping[int, int],
    n_false: int,
    scenario: str,
) -> list[tuple[int, int, int]]:
    """Return false candidate pairs and an incompatible-edge indicator."""

    true_pairs = {(parent, child) for child, parent in parents.items()}
    chosen: list[tuple[int, int, int]] = []
    seen = set(true_pairs)

    if scenario == "incompatible_cycles":
        # Reverse true edges first.  Their maps oppose the planted directed
        # section and therefore create explicit cycle/holonomy conflicts.
        order = rng.permutation(np.asarray(sorted(true_pairs), dtype=int))
        for parent, child in order:
            pair = (int(child), int(parent))
            if pair not in seen:
                chosen.append((pair[0], pair[1], 1))
                seen.add(pair)
            if len(chosen) >= n_false // 2:
                break

    attempts = 0
    while len(chosen) < n_false and attempts < 50_000:
        attempts += 1
        u = int(rng.integers(0, len(node_ids)))
        v = int(rng.integers(0, len(node_ids)))
        if u == v or (u, v) in seen:
            continue
        # Non-cycle scenarios retain the same downstream orientation as the
        # true tree.  The cycle scenario deliberately permits reverse edges.
        if scenario != "incompatible_cycles" and u > v:
            u, v = v, u
        if u == v or (u, v) in seen:
            continue
        corrupted = int(scenario == "incompatible_cycles" and u > v)
        chosen.append((u, v, corrupted))
        seen.add((u, v))
    if len(chosen) != n_false:
        raise RuntimeError("Could not construct the requested candidate graph.")
    return chosen


def generate_independent_section_case(
    seed: int,
    scenario: str,
    *,
    n_nodes: int = 20,
    false_edge_ratio: float = 2.0,
) -> IndependentSectionCase:
    """Generate one controlled case without calling HydroSheaf code."""

    if scenario not in SCENARIOS:
        raise ValueError(f"Unknown scenario: {scenario}")
    if n_nodes < 8:
        raise ValueError("n_nodes must be at least 8")

    rng = np.random.default_rng(int(seed))
    node_ids = tuple(f"N{index:02d}" for index in range(n_nodes))

    parents: dict[int, int] = {}
    true_maps: dict[tuple[int, int], tuple[float, float]] = {}
    latent = {0: float(rng.uniform(12.0, 18.0))}
    process_sigma = 0.0 if scenario == "identity_limit" else 0.08
    for child in range(1, n_nodes):
        low = max(0, child - 4)
        parent = int(rng.integers(low, child))
        parents[child] = parent
        alpha, offset = _true_map(rng, scenario)
        true_maps[(parent, child)] = (alpha, offset)
        latent[child] = float(
            alpha * latent[parent] + offset + rng.normal(0.0, process_sigma)
        )

    noise_sigma = {
        "identity_limit": 0.30,
        "heterogeneous_affine": 0.35,
        "incompatible_cycles": 0.40,
        "noisy_missing": 1.20,
    }[scenario]
    missing_rate = {
        "identity_limit": 0.15,
        "heterogeneous_affine": 0.20,
        "incompatible_cycles": 0.45,
        "noisy_missing": 0.40,
    }[scenario]

    observed_values = np.asarray(
        [latent[index] + rng.normal(0.0, noise_sigma) for index in range(n_nodes)],
        dtype=float,
    )
    missing = rng.random(n_nodes) < missing_rate
    missing[0] = False
    # Retain at least six anchors in every case.
    if int((~missing).sum()) < 6:
        missing[rng.choice(np.arange(1, n_nodes), size=6, replace=False)] = False

    n_false = int(round((n_nodes - 1) * float(false_edge_ratio)))
    false_pairs = _candidate_pairs(rng, node_ids, parents, n_false, scenario)
    if scenario == "incompatible_cycles":
        # Hide at least one endpoint of each planted reverse edge where possible;
        # this makes the global-section test non-reducible to a local residual.
        for u, v, corrupted in false_pairs:
            if corrupted and u != 0 and v != 0 and not (missing[u] or missing[v]):
                missing[u] = True

    candidates: list[CandidateMap] = []
    true_edge_ids: set[str] = set()
    for child, parent in parents.items():
        u, v = node_ids[parent], node_ids[child]
        edge_id = _edge_id(u, v)
        alpha, offset = true_maps[(parent, child)]
        true_edge_ids.add(edge_id)
        candidates.append(
            CandidateMap(
                edge_id=edge_id,
                u=u,
                v=v,
                prior_probability=_draw_prior(rng, True),
                alpha=float(alpha),
                offset=float(offset),
                is_true_edge=1,
                is_corrupted_edge=0,
            )
        )

    corrupted_edge_ids: set[str] = set()
    for u_index, v_index, corrupted in false_pairs:
        u, v = node_ids[u_index], node_ids[v_index]
        edge_id = _edge_id(u, v)
        if corrupted:
            alpha, offset = 1.22, 7.5
            corrupted_edge_ids.add(edge_id)
        else:
            alpha, offset = _false_map(rng, scenario)
        candidates.append(
            CandidateMap(
                edge_id=edge_id,
                u=u,
                v=v,
                prior_probability=_draw_prior(rng, False),
                alpha=float(alpha),
                offset=float(offset),
                is_true_edge=0,
                is_corrupted_edge=int(corrupted),
            )
        )

    rng.shuffle(candidates)
    observations = {
        node_ids[index]: (None if missing[index] else float(observed_values[index]))
        for index in range(n_nodes)
    }
    latent_truth = {node_ids[index]: float(latent[index]) for index in range(n_nodes)}
    provenance = {
        "generator": "independent_sheaf_graph_generator.py",
        "seed": int(seed),
        "scenario": scenario,
        "n_nodes": int(n_nodes),
        "n_true_edges": int(n_nodes - 1),
        "n_false_edges": int(n_false),
        "observation_noise_sigma": float(noise_sigma),
        "missing_rate_realised": float(missing.mean()),
        "imports_hydrosheaf": False,
        "excluded_capabilities": ["temporal_series", "graph_3d", "vadose"],
    }
    return IndependentSectionCase(
        seed=int(seed),
        scenario=scenario,
        node_ids=node_ids,
        observations=observations,
        latent_truth=latent_truth,
        candidates=tuple(candidates),
        true_edge_ids=frozenset(true_edge_ids),
        corrupted_edge_ids=frozenset(corrupted_edge_ids),
        provenance=provenance,
    )

