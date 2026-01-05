"""Transport model fitting for evaporation and mixing."""

from typing import Iterable, List, Tuple


def _dot(a: Iterable[float], b: Iterable[float]) -> float:
    return sum(x * y for x, y in zip(a, b))


def _weighted_dot(a: Iterable[float], b: Iterable[float], w: Iterable[float]) -> float:
    return sum(wi * x * y for wi, x, y in zip(w, a, b))


def _weighted_norm_sq(a: Iterable[float], w: Iterable[float]) -> float:
    return sum(wi * x * x for wi, x in zip(w, a))


def _residual(x_v: Iterable[float], x_pred: Iterable[float]) -> List[float]:
    return [v - p for v, p in zip(x_v, x_pred)]


def _norm_sq(residual: Iterable[float], w: Iterable[float]) -> float:
    return sum(wi * r * r for wi, r in zip(w, residual))


def fit_evaporation(
    x_u: Iterable[float],
    x_v: Iterable[float],
    weights: Iterable[float],
) -> Tuple[float, List[float], float]:
    denom = _weighted_norm_sq(x_u, weights)
    if denom == 0:
        gamma = 1.0
    else:
        gamma = _weighted_dot(x_u, x_v, weights) / denom
    gamma = max(1.0, gamma)
    x_pred = [gamma * x for x in x_u]
    residual = _residual(x_v, x_pred)
    return gamma, residual, _norm_sq(residual, weights)


def fit_mixing(
    x_u: Iterable[float],
    x_v: Iterable[float],
    x_end: Iterable[float],
    weights: Iterable[float],
) -> Tuple[float, List[float], float]:
    delta = [e - u for e, u in zip(x_end, x_u)]
    num = _weighted_dot(delta, [v - u for v, u in zip(x_v, x_u)], weights)
    denom = _weighted_norm_sq(delta, weights)
    if denom == 0:
        f = 0.0
    else:
        f = num / denom
    f = max(0.0, min(1.0, f))
    x_pred = [u + f * d for u, d in zip(x_u, delta)]
    residual = _residual(x_v, x_pred)
    return f, residual, _norm_sq(residual, weights)


def homogeneous_embed(values: Iterable[float]) -> List[float]:
    return [float(v) for v in values] + [1.0]


def evaporation_affine(gamma: float, size: int) -> Tuple[List[List[float]], List[float]]:
    matrix = [[0.0 for _ in range(size)] for _ in range(size)]
    for idx in range(size):
        matrix[idx][idx] = gamma
    offset = [0.0] * size
    return matrix, offset


def mixing_affine(
    f: float,
    x_end: Iterable[float],
) -> Tuple[List[List[float]], List[float]]:
    end_list = [float(v) for v in x_end]
    size = len(end_list)
    matrix = [[0.0 for _ in range(size)] for _ in range(size)]
    for idx in range(size):
        matrix[idx][idx] = 1.0 - f
    offset = [f * v for v in end_list]
    return matrix, offset


def apply_affine(matrix: List[List[float]], offset: List[float], x_u: Iterable[float]) -> List[float]:
    x_list = list(x_u)
    result = []
    for row, bias in zip(matrix, offset):
        value = sum(r * x for r, x in zip(row, x_list)) + bias
        result.append(value)
    return result
