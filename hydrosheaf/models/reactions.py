"""Reaction dictionary and sparse fitting."""

from dataclasses import dataclass
from typing import Iterable, List, Optional, Sequence, Tuple

from ..config import Config, DEFAULT_ION_ORDER
from ..data.minerals import get_mineral_stoich, MINERAL_LIBRARY



def _vector_from_coeffs(coeffs: dict, ion_order: Iterable[str]) -> List[float]:
    return [float(coeffs.get(ion, 0.0)) for ion in ion_order]


def build_reaction_dictionary(
    config: Config,
) -> Tuple[List[List[float]], List[str], List[bool]]:
    ion_order = config.ion_order or DEFAULT_ION_ORDER
    kappa = config.denit_kappa

    reactions: List[Tuple[str, dict, bool]] = []

    # 1. Add User-Selected Minerals
    for name in config.active_minerals:
        # Check standard library
        try:
            stoich = get_mineral_stoich(name)
            reactions.append((name, stoich, True))
        except ValueError:
            # Could range for custom loaded minerals later
            continue

    # 2. Add Process Reactions (Non-Mineral)
    # These are always available or controlled by specific flags
    reactions.append(("NO3src", {"NO3": 1}, False))
    reactions.append(("denit", {"HCO3": kappa, "NO3": -1}, False))
    
    # 3. Add Exchange Reactions (Controlled by config flag)
    if config.exchange_enabled:
        reactions.append(("CaNa_exch", {"Ca": 1, "Na": -2}, False))
        reactions.append(("MgNa_exch", {"Mg": 1, "Na": -2}, False))


    labels = [label for label, _, _ in reactions]
    mineral_mask = [is_mineral for _, _, is_mineral in reactions]
    matrix = [_vector_from_coeffs(coeffs, ion_order) for _, coeffs, _ in reactions]
    return matrix, labels, mineral_mask


def _dot(a: Iterable[float], b: Iterable[float]) -> float:
    return sum(x * y for x, y in zip(a, b))


def _combine_reactions(matrix: Sequence[Sequence[float]], weights: Sequence[float]) -> List[float]:
    if not matrix:
        return []
    result = [0.0] * len(matrix[0])
    for reaction, weight in zip(matrix, weights):
        for idx, value in enumerate(reaction):
            result[idx] += value * weight
    return result


def _weighted_vectors(
    reaction_matrix: Sequence[Sequence[float]],
    residual: Sequence[float],
    weights: Sequence[float],
) -> Tuple[List[List[float]], List[float]]:
    if not weights:
        weighted_matrix = [list(row) for row in reaction_matrix]
        weighted_residual = list(residual)
    else:
        weighted_matrix = [
            [value * (weight ** 0.5) for value, weight in zip(row, weights)]
            for row in reaction_matrix
        ]
        weighted_residual = [value * (weight ** 0.5) for value, weight in zip(residual, weights)]
    return weighted_matrix, weighted_residual


def _weighted_norm_sq(values: Sequence[float], weights: Sequence[float]) -> float:
    if not weights:
        return sum(v * v for v in values)
    return sum(w * v * v for v, w in zip(values, weights))


@dataclass
class ReactionFit:
    extents: List[float]
    residual: List[float]
    residual_norm: float
    l1_norm: float
    iterations: int
    converged: bool


def _soft_threshold(value: float, threshold: float) -> float:
    if value > threshold:
        return value - threshold
    if value < -threshold:
        return value + threshold
    return 0.0


def fit_reactions(
    residual: List[float],
    reaction_matrix: List[List[float]],
    weights: List[float],
    lambda_l1: float,
    max_iter: int = 200,
    tol: float = 1e-6,
    signed_mask: Optional[List[bool]] = None,
    lb: Optional[List[float]] = None,
    ub: Optional[List[float]] = None,
) -> ReactionFit:
    if not reaction_matrix:
        return ReactionFit([], residual, _weighted_norm_sq(residual, weights), 0.0, 0, True)

    weighted_matrix, weighted_residual = _weighted_vectors(reaction_matrix, residual, weights)
    m = len(weighted_matrix)

    gram = [[_dot(weighted_matrix[i], weighted_matrix[j]) for j in range(m)] for i in range(m)]
    s_r = [_dot(weighted_matrix[i], weighted_residual) for i in range(m)]

    if lb is not None and len(lb) != m:
        raise ValueError("lb length must match reaction matrix size.")
    if ub is not None and len(ub) != m:
        raise ValueError("ub length must match reaction matrix size.")
    if signed_mask is None:
        signed_mask = [False] * m
    if len(signed_mask) != m:
        raise ValueError("signed_mask length must match reaction matrix size.")

    z = [0.0] * m
    converged = False
    for iteration in range(1, max_iter + 1):
        max_delta = 0.0
        for j in range(m):
            rho = s_r[j] - sum(gram[j][k] * z[k] for k in range(m) if k != j)
            denom = gram[j][j] + 1e-12
            updated = _soft_threshold(rho, lambda_l1 / 2.0) / denom
            if not signed_mask[j]:
                updated = max(0.0, updated)
            if lb is not None:
                updated = max(lb[j], updated)
            if ub is not None:
                updated = min(ub[j], updated)
            max_delta = max(max_delta, abs(updated - z[j]))
            z[j] = updated
        if max_delta <= tol:
            converged = True
            break

    fitted = _combine_reactions(reaction_matrix, z)
    post_residual = [r - f for r, f in zip(residual, fitted)]
    residual_norm = _weighted_norm_sq(post_residual, weights)
    l1_norm = sum(abs(value) for value in z)
    return ReactionFit(z, post_residual, residual_norm, l1_norm, iteration, converged)
