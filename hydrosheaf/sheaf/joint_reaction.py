"""Joint chemical-state and reaction-extent optimisation on a sheaf graph.

The legacy workflow first fits a reaction vector per edge, incorporates that
vector into an edge offset, and then solves for node states.  This module
instead estimates all node concentrations and all edge reaction extents in one
convex objective, with fixed transport maps.  It deliberately does not claim
to be a joint residence-time, hydraulic, and chemistry inversion: those models
need explicit, validated physical coupling terms before they can share an
objective responsibly.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import numpy as np

from .directed_section import DirectedEdgeMap, SPECIES_CHARGES


@dataclass(frozen=True)
class JointReactionSectionResult:
    """Solution and audit data from a joint state--reaction solve."""

    node_states: Dict[str, List[float]]
    reaction_extents: Dict[str, List[float]]
    reaction_labels: Dict[str, List[str]]
    edge_residuals: Dict[str, float]
    objective: float
    iterations: int
    converged: bool
    optimality: float


@dataclass(frozen=True)
class _EdgeBlock:
    edge_map: DirectedEdgeMap
    source_index: int
    target_index: int
    transport_offset: np.ndarray
    reaction_matrix: np.ndarray
    reaction_slice: slice
    penalty_scales: np.ndarray
    signed_mask: np.ndarray
    labels: Tuple[str, ...]


def _as_vector(
    values: Sequence[float],
    *,
    name: str,
    size: int,
) -> np.ndarray:
    vector = np.asarray(values, dtype=float)
    if vector.shape != (size,) or not np.all(np.isfinite(vector)):
        raise ValueError(f"{name} must contain {size} finite values.")
    return vector


def _reaction_metadata(
    edge_map: DirectedEdgeMap,
    *,
    n_species: int,
) -> Tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
    Tuple[str, ...],
    np.ndarray,
    np.ndarray,
]:
    if edge_map.transport_offset is None:
        raise ValueError(
            "Joint state--reaction solving requires a separate transport_offset "
            f"for edge {edge_map.edge.edge_id!r}. Build maps with build_edge_maps()."
        )

    transport_offset = _as_vector(
        edge_map.transport_offset,
        name=f"transport_offset for edge {edge_map.edge.edge_id!r}",
        size=n_species,
    )

    raw_matrix = edge_map.reaction_matrix
    if raw_matrix is None:
        raise ValueError(
            "Joint state--reaction solving requires reaction_matrix metadata "
            f"for edge {edge_map.edge.edge_id!r}. Build maps with build_edge_maps()."
        )
    if not raw_matrix:
        matrix = np.empty((0, n_species), dtype=float)
    else:
        matrix = np.asarray(raw_matrix, dtype=float)
        if matrix.ndim != 2 or matrix.shape[1] != n_species or not np.all(np.isfinite(matrix)):
            raise ValueError(
                f"reaction_matrix for edge {edge_map.edge.edge_id!r} must have "
                f"shape (n_reactions, {n_species}) and finite values."
            )

    n_reactions = matrix.shape[0]
    raw_labels = edge_map.reaction_labels
    if raw_labels is None:
        labels = tuple(f"reaction_{index}" for index in range(n_reactions))
    else:
        labels = tuple(str(label) for label in raw_labels)
        if len(labels) != n_reactions:
            raise ValueError(
                f"reaction_labels for edge {edge_map.edge.edge_id!r} must match "
                "reaction_matrix rows."
            )

    raw_penalties = edge_map.reaction_penalty_scales
    penalty_scales = np.ones(n_reactions, dtype=float) if raw_penalties is None else np.asarray(raw_penalties, dtype=float)
    if penalty_scales.shape != (n_reactions,) or not np.all(np.isfinite(penalty_scales)) or np.any(penalty_scales <= 0.0):
        raise ValueError(
            f"reaction_penalty_scales for edge {edge_map.edge.edge_id!r} must be "
            "positive finite values matching reaction_matrix rows."
        )

    raw_signed = edge_map.signed_reaction_mask
    signed_mask = np.zeros(n_reactions, dtype=bool) if raw_signed is None else np.asarray(raw_signed, dtype=bool)
    if signed_mask.shape != (n_reactions,):
        raise ValueError(
            f"signed_reaction_mask for edge {edge_map.edge.edge_id!r} must match "
            "reaction_matrix rows."
        )

    raw_initial = edge_map.reaction_extents
    initial = np.zeros(n_reactions, dtype=float) if raw_initial is None else np.asarray(raw_initial, dtype=float)
    if initial.shape != (n_reactions,) or not np.all(np.isfinite(initial)):
        raise ValueError(
            f"reaction_extents for edge {edge_map.edge.edge_id!r} must be finite "
            "values matching reaction_matrix rows."
        )
    initial = initial.copy()
    initial[~signed_mask] = np.maximum(initial[~signed_mask], 0.0)

    return transport_offset, matrix, penalty_scales, labels, initial, signed_mask


def _initial_states(
    node_ids: Sequence[str],
    node_obs: Mapping[str, Optional[Sequence[float]]],
    *,
    n_species: int,
    non_negative: bool,
) -> np.ndarray:
    observations: Dict[Tuple[int, int], float] = {}
    index = {node_id: position for position, node_id in enumerate(node_ids)}
    means: List[List[float]] = [[] for _ in range(n_species)]

    for node_id, values in node_obs.items():
        if node_id not in index:
            raise ValueError(f"node_obs contains unknown node {node_id!r}.")
        if values is None:
            continue
        if len(values) != n_species:
            raise ValueError(
                f"Observation vector for node {node_id!r} must contain {n_species} values."
            )
        for species_index, value in enumerate(values):
            if value is None:
                continue
            numeric = float(value)
            if not np.isfinite(numeric):
                raise ValueError(f"Observation for node {node_id!r} must be finite.")
            observations[(index[node_id], species_index)] = numeric
            means[species_index].append(numeric)

    state = np.zeros((len(node_ids), n_species), dtype=float)
    for species_index, values in enumerate(means):
        if values:
            state[:, species_index] = float(np.mean(values))
    for (node_index, species_index), value in observations.items():
        state[node_index, species_index] = value
    if non_negative:
        state = np.maximum(state, 0.0)
    return state


def solve_joint_reaction_section(
    node_ids: Iterable[str],
    edge_maps: Iterable[DirectedEdgeMap],
    node_obs: Mapping[str, Optional[Sequence[float]]],
    species_names: Sequence[str],
    *,
    species_weights: Optional[Sequence[float]] = None,
    obs_weight: float = 1.0,
    charge_balance_weight: float = 0.0,
    diag_eps: float = 1e-6,
    lambda_l1: float = 0.0,
    lambda_l2: float = 0.0,
    max_iter: int = 500,
    tol: float = 1e-7,
    non_negative: bool = True,
) -> JointReactionSectionResult:
    """Jointly fit node chemistry and edge reaction extents.

    With fixed edge transport maps, this minimises the convex objective

    ``0.5 * sum_e ||alpha_e x_u + t_e + S_e.T z_e - x_v||_W^2``

    plus observed-state fit, optional charge-balance, state ridge, and sparse
    reaction penalties.  A full proximal-gradient step updates every state and
    every edge extent together; the solver never substitutes a pre-fitted
    reaction offset into the state equation.
    """
    nodes = list(node_ids)
    if not nodes:
        return JointReactionSectionResult({}, {}, {}, {}, 0.0, 0, True, 0.0)
    if len(nodes) != len(set(nodes)):
        raise ValueError("node_ids must be unique.")
    if not species_names:
        raise ValueError("species_names must not be empty.")
    if not np.isfinite(obs_weight) or obs_weight < 0.0:
        raise ValueError("obs_weight must be finite and non-negative.")
    if not np.isfinite(charge_balance_weight) or charge_balance_weight < 0.0:
        raise ValueError("charge_balance_weight must be finite and non-negative.")
    if not np.isfinite(diag_eps) or diag_eps < 0.0:
        raise ValueError("diag_eps must be finite and non-negative.")
    if not np.isfinite(lambda_l1) or lambda_l1 < 0.0:
        raise ValueError("lambda_l1 must be finite and non-negative.")
    if not np.isfinite(lambda_l2) or lambda_l2 < 0.0:
        raise ValueError("lambda_l2 must be finite and non-negative.")
    if max_iter < 1:
        raise ValueError("max_iter must be at least 1.")
    if not np.isfinite(tol) or tol < 0.0:
        raise ValueError("tol must be finite and non-negative.")

    n_species = len(species_names)
    if species_weights is None:
        weights = np.ones(n_species, dtype=float)
    else:
        weights = _as_vector(species_weights, name="species_weights", size=n_species)
    if np.any(weights < 0.0):
        raise ValueError("species_weights must be non-negative.")

    index = {node_id: position for position, node_id in enumerate(nodes)}
    edge_maps_list = list(edge_maps)
    if len({edge_map.edge.edge_id for edge_map in edge_maps_list}) != len(edge_maps_list):
        raise ValueError("edge_maps must have unique edge IDs.")

    cursor = 0
    blocks: List[_EdgeBlock] = []
    initial_extents: List[np.ndarray] = []
    for edge_map in edge_maps_list:
        if edge_map.edge.u not in index or edge_map.edge.v not in index:
            raise ValueError(
                f"Edge {edge_map.edge.edge_id!r} refers to a node absent from node_ids."
            )
        if not np.isfinite(edge_map.alpha) or not np.isfinite(edge_map.weight) or edge_map.weight < 0.0:
            raise ValueError(
                f"Edge {edge_map.edge.edge_id!r} must have finite alpha and non-negative weight."
            )
        (
            transport_offset,
            reaction_matrix,
            penalty_scales,
            labels,
            initial,
            signed_mask,
        ) = _reaction_metadata(edge_map, n_species=n_species)
        extent_slice = slice(cursor, cursor + reaction_matrix.shape[0])
        cursor += reaction_matrix.shape[0]
        blocks.append(
            _EdgeBlock(
                edge_map=edge_map,
                source_index=index[edge_map.edge.u],
                target_index=index[edge_map.edge.v],
                transport_offset=transport_offset,
                reaction_matrix=reaction_matrix,
                reaction_slice=extent_slice,
                penalty_scales=penalty_scales,
                signed_mask=signed_mask,
                labels=labels,
            )
        )
        initial_extents.append(initial)

    state = _initial_states(
        nodes,
        node_obs,
        n_species=n_species,
        non_negative=non_negative,
    )
    extents = (
        np.concatenate(initial_extents) if initial_extents else np.empty(0, dtype=float)
    )
    penalties = (
        np.concatenate([block.penalty_scales for block in blocks])
        if blocks
        else np.empty(0, dtype=float)
    )
    signed = (
        np.concatenate([block.signed_mask for block in blocks])
        if blocks
        else np.empty(0, dtype=bool)
    )

    observations: List[Tuple[int, int, float]] = []
    for node_id, values in node_obs.items():
        if values is None:
            continue
        for species_index, value in enumerate(values):
            if value is not None:
                observations.append((index[node_id], species_index, float(value)))
    charges = np.asarray([SPECIES_CHARGES.get(name, 0.0) for name in species_names])

    # A diagonal majorizer of the smooth Hessian lets a single proximal update
    # move concentrations and reaction extents on their natural scales.  A
    # global step would be dominated by strongly weighted observations and
    # would make the reaction block unnecessarily slow.  The row-sum bound is
    # conservative enough to majorize every edge, observation, ridge, and
    # charge-balance quadratic while retaining a simultaneous full-vector step.
    state_majorizer = np.zeros_like(state)
    extent_majorizer = np.zeros_like(extents)
    for block in blocks:
        edge_weight = float(block.edge_map.weight)
        for species_index, species_weight in enumerate(weights):
            if edge_weight == 0.0 or species_weight == 0.0:
                continue
            coefficients = [
                ("state", block.source_index, species_index, float(block.edge_map.alpha)),
                ("state", block.target_index, species_index, -1.0),
            ]
            coefficients.extend(
                (
                    "extent",
                    block.reaction_slice.start + reaction_index,
                    0,
                    float(coefficient),
                )
                for reaction_index, coefficient in enumerate(
                    block.reaction_matrix[:, species_index]
                )
                if coefficient != 0.0
            )
            absolute_sum = sum(abs(coefficient) for *_, coefficient in coefficients)
            for kind, first_index, second_index, coefficient in coefficients:
                contribution = edge_weight * species_weight * abs(coefficient) * absolute_sum
                if kind == "state":
                    state_majorizer[first_index, second_index] += contribution
                else:
                    extent_majorizer[first_index] += contribution
    for node_index, species_index, _ in observations:
        state_majorizer[node_index, species_index] += obs_weight
    state_majorizer += diag_eps
    if charge_balance_weight and np.any(charges):
        state_majorizer += charge_balance_weight * np.abs(charges) * np.sum(np.abs(charges))
    extent_majorizer += lambda_l2
    state_majorizer = np.maximum(state_majorizer, 1e-12)
    extent_majorizer = np.maximum(extent_majorizer, 1e-12)

    def smooth_value_and_gradient(
        current_state: np.ndarray,
        current_extents: np.ndarray,
    ) -> Tuple[float, np.ndarray, np.ndarray]:
        value = 0.0
        grad_state = np.zeros_like(current_state)
        grad_extents = np.zeros_like(current_extents)
        for block in blocks:
            residual = (
                float(block.edge_map.alpha) * current_state[block.source_index]
                + block.transport_offset
                - current_state[block.target_index]
            )
            if block.reaction_matrix.shape[0]:
                residual = residual + block.reaction_matrix.T @ current_extents[block.reaction_slice]
            weighted_residual = float(block.edge_map.weight) * weights * residual
            value += 0.5 * float(np.dot(residual, weighted_residual))
            grad_state[block.source_index] += float(block.edge_map.alpha) * weighted_residual
            grad_state[block.target_index] -= weighted_residual
            if block.reaction_matrix.shape[0]:
                grad_extents[block.reaction_slice] += block.reaction_matrix @ weighted_residual

        for node_index, species_index, observed in observations:
            residual = current_state[node_index, species_index] - observed
            value += 0.5 * obs_weight * residual * residual
            grad_state[node_index, species_index] += obs_weight * residual

        if diag_eps:
            value += 0.5 * diag_eps * float(np.sum(current_state * current_state))
            grad_state += diag_eps * current_state
        if charge_balance_weight and np.any(charges):
            imbalance = current_state @ charges
            value += 0.5 * charge_balance_weight * float(np.dot(imbalance, imbalance))
            grad_state += charge_balance_weight * np.outer(imbalance, charges)
        if lambda_l2 and current_extents.size:
            value += 0.5 * lambda_l2 * float(np.dot(current_extents, current_extents))
            grad_extents += lambda_l2 * current_extents
        return value, grad_state, grad_extents

    def full_objective(current_state: np.ndarray, current_extents: np.ndarray) -> float:
        smooth, _, _ = smooth_value_and_gradient(current_state, current_extents)
        if current_extents.size:
            smooth += lambda_l1 * float(np.dot(penalties, np.abs(current_extents)))
        return smooth

    converged = False
    optimality = float("inf")
    iteration = 0
    for iteration in range(1, max_iter + 1):
        _, grad_state, grad_extents = smooth_value_and_gradient(state, extents)
        next_state = state - grad_state / state_majorizer
        if non_negative:
            next_state = np.maximum(next_state, 0.0)
        next_extents = extents - grad_extents / extent_majorizer
        if next_extents.size:
            thresholds = lambda_l1 * penalties / extent_majorizer
            next_extents[signed] = np.sign(next_extents[signed]) * np.maximum(
                np.abs(next_extents[signed]) - thresholds[signed], 0.0
            )
            next_extents[~signed] = np.maximum(
                next_extents[~signed] - thresholds[~signed], 0.0
            )

        delta_state = next_state - state
        delta_extents = next_extents - extents

        optimality = max(
            float(np.max(np.abs(delta_state) * state_majorizer))
            if delta_state.size
            else 0.0,
            float(np.max(np.abs(delta_extents) * extent_majorizer))
            if delta_extents.size
            else 0.0,
        )
        state, extents = next_state, next_extents
        if optimality <= tol:
            converged = True
            break

    reaction_extents: Dict[str, List[float]] = {}
    reaction_labels: Dict[str, List[str]] = {}
    edge_residuals: Dict[str, float] = {}
    for block in blocks:
        edge_id = block.edge_map.edge.edge_id
        edge_extents = extents[block.reaction_slice]
        reaction_extents[edge_id] = edge_extents.tolist()
        reaction_labels[edge_id] = list(block.labels)
        residual = (
            float(block.edge_map.alpha) * state[block.source_index]
            + block.transport_offset
            - state[block.target_index]
        )
        if block.reaction_matrix.shape[0]:
            residual = residual + block.reaction_matrix.T @ edge_extents
        weighted_sq = float(block.edge_map.weight) * float(np.dot(weights, residual * residual))
        edge_residuals[edge_id] = float(np.sqrt(max(0.0, weighted_sq)))

    return JointReactionSectionResult(
        node_states={node_id: state[position].tolist() for position, node_id in enumerate(nodes)},
        reaction_extents=reaction_extents,
        reaction_labels=reaction_labels,
        edge_residuals=edge_residuals,
        objective=full_objective(state, extents),
        iterations=iteration,
        converged=converged,
        optimality=optimality,
    )
