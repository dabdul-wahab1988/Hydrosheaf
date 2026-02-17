"""Reaction dictionary and sparse fitting."""

from dataclasses import dataclass
from typing import Iterable, List, Mapping, Optional, Sequence, Tuple

from ..config import Config, DEFAULT_ION_ORDER
from ..data.minerals import get_mineral_stoich


def _vector_from_coeffs(coeffs: Mapping[str, float], ion_order: Iterable[str]) -> List[float]:
    return [float(coeffs.get(ion, 0.0)) for ion in ion_order]


def build_reaction_dictionary(
    config: Config,
) -> Tuple[List[List[float]], List[str], List[bool]]:
    ion_order = config.ion_order or DEFAULT_ION_ORDER
    kappa = config.denit_kappa

    reactions: List[Tuple[str, Mapping[str, float], bool]] = []

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


def _combine_reactions(
    matrix: Sequence[Sequence[float]], weights: Sequence[float]
) -> List[float]:
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
            [value * (weight**0.5) for value, weight in zip(row, weights)]
            for row in reaction_matrix
        ]
        weighted_residual = [
            value * (weight**0.5) for value, weight in zip(residual, weights)
        ]
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

    # Uncertainty fields (populated if uncertainty_method != "none")
    extents_std: Optional[List[float]] = None
    extents_ci_low: Optional[List[float]] = None
    extents_ci_high: Optional[List[float]] = None
    uncertainty_result: Optional[object] = (
        None  # UncertaintyResult from uncertainty module
    )


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
        return ReactionFit(
            [], residual, _weighted_norm_sq(residual, weights), 0.0, 0, True
        )

    weighted_matrix, weighted_residual = _weighted_vectors(
        reaction_matrix, residual, weights
    )
    m = len(weighted_matrix)

    gram = [
        [_dot(weighted_matrix[i], weighted_matrix[j]) for j in range(m)]
        for i in range(m)
    ]

    # Compute condition number for multicollinearity warning
    trace_gram = sum(gram[i][i] for i in range(m))
    frob_sq = sum(gram[i][j]**2 for i in range(m) for j in range(m))
    # Approximate condition number using trace/Frobenius ratio
    if trace_gram > 1e-12:
        approx_cond = (frob_sq**0.5) / (trace_gram / m)
        if approx_cond > 30:
            import warnings
            warnings.warn(
                f"Reaction matrix may be ill-conditioned (approx. κ={approx_cond:.1f}). "
                "Consider reducing active minerals or enabling exchange reactions."
            )

    s_r = [_dot(weighted_matrix[i], weighted_residual) for i in range(m)]

    if lb is not None and ub is not None:
        if len(lb) != m or len(ub) != m:
            raise ValueError("lb and ub must match matrix size.")
        for i, (l, u) in enumerate(zip(lb, ub)):
            if l is not None and u is not None and l > u:
                raise ValueError(f"Lower bound exceeds upper bound at index {i}: {l} > {u}")

    if signed_mask is None:
        signed_mask = [False] * m
    if len(signed_mask) != m:
        raise ValueError("signed_mask length must match reaction matrix size.")

    # Adaptive ridge regression parameter for numerical stability
    # Use max diagonal element to scale epsilon, plus a small absolute floor
    max_diag = 0.0
    for i in range(m):
        if gram[i][i] > max_diag:
            max_diag = gram[i][i]
    ridge_epsilon = max_diag * 1e-10 + 1e-20

    z = [0.0] * m
    converged = False
    iteration = 0
    
    # Store objective history to detect cycling
    prev_obj = float('inf')
    
    for iteration in range(1, max_iter + 1):
        max_delta = 0.0
        for j in range(m):
            # Calculate partial residual correlation
            # Use math.fsum if available for precision, but sum is usually fine for small m
            dot_prod = sum(gram[j][k] * z[k] for k in range(m) if k != j)
            rho = s_r[j] - dot_prod
            
            denom = gram[j][j] + ridge_epsilon
            
            # If denominator is effectively zero (zero reaction vector), force coeff to 0
            if denom < 1e-15:
                updated = 0.0
            else:
                updated = _soft_threshold(rho, lambda_l1 / 2.0) / denom
            
            if not signed_mask[j]:
                updated = max(0.0, updated)
            if lb is not None and lb[j] is not None:
                updated = max(lb[j], updated)
            if ub is not None and ub[j] is not None:
                updated = min(ub[j], updated)
            
            delta = abs(updated - z[j])
            if delta > max_delta:
                max_delta = delta

            z[j] = updated
            
        if max_delta <= tol:
            converged = True
            break
            
        # Optional: Check objective function every 10 iterations to detect cycling
        # (omitted for speed unless we want strict guarantees)


    fitted = _combine_reactions(reaction_matrix, z)
    post_residual = [r - f for r, f in zip(residual, fitted)]
    residual_norm = _weighted_norm_sq(post_residual, weights)
    l1_norm = sum(abs(value) for value in z)
    return ReactionFit(z, post_residual, residual_norm, l1_norm, iteration, converged)
