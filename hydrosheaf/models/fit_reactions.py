"""Reaction fitting helpers."""

from typing import List, Optional

from .reactions import ReactionFit, fit_reactions as _fit


def fit_reactions(
    residual: List[float],
    reaction_matrix: List[List[float]],
    weights: List[float],
    lambda_l1: float,
    max_iter: int,
    tol: float,
    signed_mask: Optional[List[bool]] = None,
    lb: Optional[List[float]] = None,
    ub: Optional[List[float]] = None,
) -> ReactionFit:
    return _fit(
        residual,
        reaction_matrix,
        weights=weights,
        lambda_l1=lambda_l1,
        max_iter=max_iter,
        tol=tol,
        signed_mask=signed_mask,
        lb=lb,
        ub=ub,
    )
