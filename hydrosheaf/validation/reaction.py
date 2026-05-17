"""Validation helpers for sparse hydrogeochemical inverse modelling.

These functions support the M5 framing:

Hydrosheaf is a sparse linear inverse reaction model screened and stress-tested
using PHREEQC thermodynamic and forward-validation diagnostics. It is not a
fully coupled nonlinear PHREEQC inverse solver.
"""
from __future__ import annotations

import math
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import numpy as np

from ..config import Config
from ..models.reactions import fit_reactions


def _safe_float(value: Any) -> Optional[float]:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(number):
        return None
    return number


def _weighted_norm_sq(values: Sequence[float], weights: Sequence[float]) -> float:
    if not weights:
        return float(sum(value * value for value in values))
    return float(sum(w * value * value for value, w in zip(values, weights)))


def _active_reactions(
    labels: Sequence[str],
    extents: Sequence[float],
    *,
    threshold: float,
) -> List[Dict[str, Any]]:
    rows = []
    for label, extent in zip(labels, extents):
        value = float(extent)
        if abs(value) >= threshold:
            rows.append({"reaction": label, "extent": value})
    return rows


def _rank_diagnostics(matrix: Sequence[Sequence[float]]) -> Dict[str, Any]:
    if not matrix:
        return {
            "n_ions": 0,
            "n_reactions": 0,
            "rank": 0,
            "underdetermined": False,
            "rank_deficient": False,
            "condition_number": float("nan"),
            "ill_conditioned": False,
        }
    arr = np.array(matrix, dtype=float).T
    n_ions, n_reactions = arr.shape
    rank = int(np.linalg.matrix_rank(arr))
    if n_ions == 0 or n_reactions == 0:
        condition = float("nan")
    else:
        condition = float(np.linalg.cond(arr))
    return {
        "n_ions": int(n_ions),
        "n_reactions": int(n_reactions),
        "rank": rank,
        "underdetermined": bool(n_reactions > n_ions),
        "rank_deficient": bool(rank < min(n_ions, n_reactions)),
        "condition_number": condition,
        "ill_conditioned": bool(math.isfinite(condition) and condition >= 100.0),
    }


def thermodynamic_bound_violations(
    labels: Sequence[str],
    extents: Sequence[float],
    bounds: Optional[Mapping[str, Any]],
    *,
    tolerance: float = 1.0e-8,
) -> List[Dict[str, Any]]:
    """Return extent rows that violate PHREEQC-derived SI bounds."""

    if not bounds:
        return []
    lb = list(bounds.get("lb") or [])
    ub = list(bounds.get("ub") or [])
    active = dict(bounds.get("constraints_active") or {})
    if not lb or not ub:
        return []
    rows: List[Dict[str, Any]] = []
    for idx, (label, extent) in enumerate(zip(labels, extents)):
        value = float(extent)
        lower = lb[idx] if idx < len(lb) else None
        upper = ub[idx] if idx < len(ub) else None
        lower_value = _safe_float(lower)
        upper_value = _safe_float(upper)
        violation = False
        if lower_value is not None and value < lower_value - tolerance:
            violation = True
        if upper_value is not None and value > upper_value + tolerance:
            violation = True
        if violation:
            rows.append(
                {
                    "reaction": label,
                    "extent": value,
                    "lower_bound": lower_value,
                    "upper_bound": upper_value,
                    "constraint": active.get(label, "unknown"),
                }
            )
    return rows


def fit_sparse_reaction_once(
    transport_residual: Sequence[float],
    reaction_matrix: Sequence[Sequence[float]],
    reaction_labels: Sequence[str],
    weights: Sequence[float],
    *,
    lambda_l1: float,
    lambda_l2: float = 0.0,
    penalty_scales: Optional[List[float]] = None,
    signed_mask: Optional[List[bool]] = None,
    lb: Optional[List[float]] = None,
    ub: Optional[List[float]] = None,
    selected_threshold: float = 1.0e-6,
    max_iter: int = 300,
    tol: float = 1.0e-6,
) -> Dict[str, Any]:
    """Fit sparse linear reactions and report transport/reaction separation."""

    before = _weighted_norm_sq(transport_residual, weights)
    fit = fit_reactions(
        list(map(float, transport_residual)),
        [list(map(float, row)) for row in reaction_matrix],
        list(map(float, weights)),
        float(lambda_l1),
        lambda_l2=lambda_l2,
        penalty_scales=penalty_scales,
        max_iter=max_iter,
        tol=tol,
        signed_mask=signed_mask,
        lb=lb,
        ub=ub,
    )
    after = float(fit.residual_norm)
    active = _active_reactions(reaction_labels, fit.extents, threshold=selected_threshold)
    return {
        "lambda_l1": float(lambda_l1),
        "lambda_l2": float(lambda_l2),
        "transport_residual_norm": before,
        "reaction_residual_norm": after,
        "residual_reduction_fraction": (before - after) / before if before > 0 else 0.0,
        "l1_norm": float(fit.l1_norm),
        "n_selected_reactions": len(active),
        "selected_reactions": active,
        "extents": [float(value) for value in fit.extents],
        "post_reaction_residual": [float(value) for value in fit.residual],
        "iterations": int(fit.iterations),
        "converged": bool(fit.converged),
    }


def l1_penalty_sensitivity(
    transport_residual: Sequence[float],
    reaction_matrix: Sequence[Sequence[float]],
    reaction_labels: Sequence[str],
    weights: Sequence[float],
    lambda_grid: Iterable[float],
    *,
    lambda_l2: float = 0.0,
    penalty_scales: Optional[List[float]] = None,
    signed_mask: Optional[List[bool]] = None,
    lb: Optional[List[float]] = None,
    ub: Optional[List[float]] = None,
    selected_threshold: float = 1.0e-6,
) -> List[Dict[str, Any]]:
    """Refit the reaction dictionary across an L1 penalty grid."""

    rows = []
    for lambda_l1 in sorted(float(value) for value in lambda_grid):
        rows.append(
            fit_sparse_reaction_once(
                transport_residual,
                reaction_matrix,
                reaction_labels,
                weights,
                lambda_l1=lambda_l1,
                lambda_l2=lambda_l2,
                penalty_scales=penalty_scales,
                signed_mask=signed_mask,
                lb=lb,
                ub=ub,
                selected_threshold=selected_threshold,
            )
        )
    return rows


def missing_ion_sensitivity(
    transport_residual: Sequence[float],
    reaction_matrix: Sequence[Sequence[float]],
    reaction_labels: Sequence[str],
    ion_order: Sequence[str],
    weights: Sequence[float],
    missing_ion_sets: Iterable[Iterable[str]],
    *,
    lambda_l1: float,
    lambda_l2: float = 0.0,
    penalty_scales: Optional[List[float]] = None,
    selected_threshold: float = 1.0e-6,
) -> List[Dict[str, Any]]:
    """Refit after dropping selected ions by setting their weights to zero."""

    rows: List[Dict[str, Any]] = []
    ion_index = {ion: idx for idx, ion in enumerate(ion_order)}
    for missing in missing_ion_sets:
        missing_set = {str(ion) for ion in missing}
        sensitivity_weights = list(map(float, weights))
        for ion in missing_set:
            idx = ion_index.get(ion)
            if idx is not None and idx < len(sensitivity_weights):
                sensitivity_weights[idx] = 0.0
        fit = fit_sparse_reaction_once(
            transport_residual,
            reaction_matrix,
            reaction_labels,
            sensitivity_weights,
            lambda_l1=lambda_l1,
            lambda_l2=lambda_l2,
            penalty_scales=penalty_scales,
            selected_threshold=selected_threshold,
        )
        fit["missing_ions"] = sorted(missing_set)
        fit["n_missing_ions"] = len(missing_set)
        rows.append(fit)
    return rows


def validate_sparse_inverse_reaction_model(
    upstream: Sequence[float],
    downstream: Sequence[float],
    post_transport: Sequence[float],
    reaction_matrix: Sequence[Sequence[float]],
    reaction_labels: Sequence[str],
    config: Config,
    *,
    lambda_grid: Iterable[float],
    missing_ion_sets: Optional[Iterable[Iterable[str]]] = None,
    phreeqc_bounds: Optional[Mapping[str, Any]] = None,
    selected_threshold: float = 1.0e-6,
    penalty_scales: Optional[List[float]] = None,
) -> Dict[str, Any]:
    """Run the M5 sparse inverse-reaction validation harness."""

    x_v = list(map(float, downstream))
    x_transport = list(map(float, post_transport))
    if len(x_v) != len(x_transport):
        raise ValueError("downstream and post_transport vectors must have the same length.")
    residual = [v - t for v, t in zip(x_v, x_transport)]
    weights = list(map(float, config.weights))
    rank = _rank_diagnostics(reaction_matrix)
    l1_rows = l1_penalty_sensitivity(
        residual,
        reaction_matrix,
        reaction_labels,
        weights,
        lambda_grid,
        lambda_l2=config.lambda_l2,
        penalty_scales=penalty_scales,
        lb=list(phreeqc_bounds.get("lb")) if phreeqc_bounds and phreeqc_bounds.get("lb") else None,
        ub=list(phreeqc_bounds.get("ub")) if phreeqc_bounds and phreeqc_bounds.get("ub") else None,
        selected_threshold=selected_threshold,
    )
    best = min(l1_rows, key=lambda row: (row["reaction_residual_norm"], row["n_selected_reactions"]))
    missing_rows = missing_ion_sensitivity(
        residual,
        reaction_matrix,
        reaction_labels,
        config.ion_order,
        weights,
        missing_ion_sets or [],
        lambda_l1=float(best["lambda_l1"]),
        lambda_l2=config.lambda_l2,
        penalty_scales=penalty_scales,
        selected_threshold=selected_threshold,
    )
    violations = thermodynamic_bound_violations(
        reaction_labels,
        best["extents"],
        phreeqc_bounds,
    )
    flags = {
        "thermodynamic_bound_violation": bool(violations),
        "missing_ion_sensitive": _is_missing_ion_sensitive(best, missing_rows),
        "l1_penalty_sensitive": _is_l1_sensitive(l1_rows),
        "linearized_dictionary_limit": bool(
            rank["underdetermined"] or rank["rank_deficient"] or rank["ill_conditioned"]
        ),
        "not_fully_coupled_phreeqc_inverse_solver": True,
    }
    return {
        "model_framing": (
            "sparse linear inverse reaction model screened and stress-tested using "
            "PHREEQC thermodynamic and forward-validation diagnostics"
        ),
        "not_a_fully_coupled_nonlinear_phreeqc_inverse_solver": True,
        "transport_reaction_separation": {
            "upstream": list(map(float, upstream)),
            "post_transport": x_transport,
            "downstream": x_v,
            "transport_residual": residual,
            "transport_residual_norm": _weighted_norm_sq(residual, weights),
        },
        "reaction_dictionary": {
            "labels": list(reaction_labels),
            "rank_diagnostics": rank,
        },
        "l1_penalty_sensitivity": l1_rows,
        "missing_ion_sensitivity": missing_rows,
        "best_fit": best,
        "thermodynamic_bound_violations": violations,
        "flags": flags,
        "claim_guardrail": (
            "Report Hydrosheaf as a sparse linear inverse reaction model with "
            "PHREEQC screening/forward checks. Do not describe it as a fully "
            "coupled nonlinear PHREEQC inverse solver."
        ),
    }


def _is_l1_sensitive(rows: Sequence[Mapping[str, Any]]) -> bool:
    if len(rows) < 2:
        return False
    selected = {int(row["n_selected_reactions"]) for row in rows}
    residuals = [float(row["reaction_residual_norm"]) for row in rows]
    if len(selected) > 1:
        return True
    min_res = min(residuals)
    max_res = max(residuals)
    return min_res > 0 and max_res / min_res >= 2.0


def _is_missing_ion_sensitive(best: Mapping[str, Any], rows: Sequence[Mapping[str, Any]]) -> bool:
    if not rows:
        return False
    baseline = {row["reaction"] for row in best.get("selected_reactions", [])}
    baseline_norm = float(best.get("reaction_residual_norm", 0.0))
    for row in rows:
        current = {item["reaction"] for item in row.get("selected_reactions", [])}
        if current != baseline:
            return True
        current_norm = float(row.get("reaction_residual_norm", 0.0))
        if baseline_norm > 0 and current_norm / baseline_norm >= 2.0:
            return True
    return False
