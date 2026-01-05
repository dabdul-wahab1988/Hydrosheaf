"""Sheaf-style residual helpers."""

from typing import Iterable, List, Sequence

from .transport import apply_affine


def _combine_reactions(matrix: Sequence[Sequence[float]], weights: Sequence[float]) -> List[float]:
    if not matrix:
        return []
    result = [0.0] * len(matrix[0])
    for reaction, weight in zip(matrix, weights):
        for idx, value in enumerate(reaction):
            result[idx] += value * weight
    return result


def edge_residual(
    x_u: Iterable[float],
    x_v: Iterable[float],
    transport_matrix: Sequence[Sequence[float]],
    transport_offset: Sequence[float],
    reaction_matrix: Sequence[Sequence[float]],
    reaction_extents: Sequence[float],
) -> List[float]:
    transport_pred = apply_affine(list(map(list, transport_matrix)), list(transport_offset), x_u)
    reaction_pred = _combine_reactions(reaction_matrix, reaction_extents)
    return [v - (t + r) for v, t, r in zip(x_v, transport_pred, reaction_pred)]
