"""Systematic Gaussian sigma estimation for head gradient computation.

Provides three independent, MODPATH-free methods to estimate the optimal
Gaussian smoothing sigma for MODFLOW head gradient computation.  Each method
probes a different structural property of the gradient field across a sweep
of sigma values, and the module recommends a consensus sigma.

Methods
-------
1. *Circular Variance Elbow*
   Measures how tightly flow-direction angles cluster as smoothing increases.
   The elbow is the sigma where the reduction in circular variance per unit
   sigma falls below a fraction of its maximum — i.e. the point of
   diminishing returns for directional alignment.

2. *Gradient Magnitude Stability*
   Tracks the median gradient magnitude and flags the first sigma where the
   relative change from the previous sigma drops below 5 %, indicating the
   magnitude field has stabilised.

3. *Projection IQR Crossover*
   Projects the gradient vector at each source node onto the edge direction
   for a set of candidate edges.  Finds the smallest sigma at which the
   inter-quartile range of these projections exceeds a minimum threshold
   (1e-8) with at least four projections — the point where the gradient
   field provides meaningful directional discrimination.

Recommendation
--------------
The recommended sigma is ``max(elbow, stability, crossover)`` rounded to the
nearest integer in the sweep range.  The method that produced the maximum
(the "binding" method) is reported for diagnostic transparency.
"""

from __future__ import annotations

import math
import re
import sys
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from .modflow_head import (
    GridGeometry,
    _cell_to_ijk,
    _gaussian_smooth_layer,
    _head_array_from_map,
    _safe_gradient_1d,
)

# ---------------------------------------------------------------------------
# Default sigma sweep range
# ---------------------------------------------------------------------------

_DEFAULT_SIGMA_RANGE: List[int] = [0, 1, 2, 3, 4, 6, 8, 12, 16, 24, 32]

# Fraction of maximum variance reduction below which we declare the elbow.
_ELBOW_FRACTION: float = 0.10

# Relative-change threshold for magnitude stabilisation.
_STABILITY_THRESHOLD: float = 0.05

# Minimum IQR and projection count for the crossover method.
_MIN_IQR: float = 1e-8
_MIN_PROJECTIONS: int = 4


# ===================================================================
# Public API
# ===================================================================


def compute_gradient_sweep(
    head_map: Dict[int, float],
    grid: GridGeometry,
    sigma_range: Optional[List[int]] = None,
) -> Dict[int, Dict[int, Tuple[float, float, float]]]:
    """Compute head-gradient maps at every sigma in one efficient pass.

    The MODFLOW head array is built **once** and then smoothed at different
    sigma values, avoiding redundant array construction.

    Parameters
    ----------
    head_map:
        Output of :func:`~hydrosheaf.physics.modflow_head.parse_fhd` (or
        ``parse_hds``).  Maps 1-based MODFLOW cell index to head value.
    grid:
        Grid-geometry descriptor.
    sigma_range:
        Sigma values (in cell-width units) to sweep.  Defaults to
        ``[0, 1, 2, 3, 4, 6, 8, 12, 16, 24, 32]``.

    Returns
    -------
    ``{sigma: {cell_index: (gx, gy, gz)}}``
        One gradient map per sigma value.
    """
    if sigma_range is None:
        sigma_range = list(_DEFAULT_SIGMA_RANGE)

    # --- Build head array once (avoid rebuilding per sigma) ---------------
    head_array = _head_array_from_map(head_map, grid.ncol, grid.nrow, grid.nlay)

    # Pre-compute rotation matrix (same for all sigma)
    theta = math.radians(grid.rotation_deg)
    cos_t = math.cos(theta)
    sin_t = math.sin(theta)

    ncol = grid.ncol
    nrow = grid.nrow
    nlay = grid.nlay
    dx = grid.dx
    dy = grid.dy

    result: Dict[int, Dict[int, Tuple[float, float, float]]] = {}

    for sigma in sigma_range:
        # Work on a copy so smoothing at one sigma does not affect the next
        smoothed = head_array.copy()

        if float(sigma) > 0:
            for k in range(nlay):
                smoothed[:, :, k] = _gaussian_smooth_layer(
                    smoothed[:, :, k], float(sigma)
                )

        gradient_map: Dict[int, Tuple[float, float, float]] = {}
        for cell_idx in head_map:
            col, row, lay = _cell_to_ijk(cell_idx, ncol, nrow)
            if not (0 <= col < ncol and 0 <= row < nrow and 0 <= lay < nlay):
                continue

            h_center = smoothed[col, row, lay]
            if np.isnan(h_center):
                continue

            gx_grid = _safe_gradient_1d(smoothed[:, row, lay], col, dx, ncol)
            gy_grid = _safe_gradient_1d(smoothed[col, :, lay], row, dy, nrow)

            if np.isnan(gx_grid) or np.isnan(gy_grid):
                continue

            # Rotate grid gradient to world (projected) coordinates
            gx_world = gx_grid * cos_t - gy_grid * sin_t
            gy_world = gx_grid * sin_t + gy_grid * cos_t
            gz: float = 0.0

            if nlay > 1:
                gz = _safe_gradient_1d(smoothed[col, row, :], lay, 1.0, nlay)

            gradient_map[cell_idx] = (float(gx_world), float(gy_world), float(gz))

        result[int(sigma)] = gradient_map

    return result


def estimate_smoothing_scale(
    head_map: Dict[int, float],
    grid: GridGeometry,
    node_df: pd.DataFrame,
    sigma_range: Optional[List[int]] = None,
) -> Dict[str, Any]:
    """Estimate the optimal Gaussian smoothing sigma for head-gradient computation.

    Runs three independent, MODPATH-free estimation methods across a sweep of
    sigma values and returns a consensus recommendation.

    Parameters
    ----------
    head_map:
        Dict mapping 1-based MODFLOW cell index to head value.
    grid:
        Grid-geometry descriptor.
    node_df:
        DataFrame with at least ``node_id``, ``x``, ``y`` columns.
        Node IDs are expected to follow the ``"cell_{int}"`` convention
        (e.g. ``"cell_130749"``).
    sigma_range:
        Sigma values to sweep.  Defaults to ``[0, 1, 2, 3, 4, 6, 8, 12, 16, 24, 32]``.

    Returns
    -------
    dict with keys:

    * ``"recommended_sigma"`` — int, the recommended sigma.
    * ``"binding_method"`` — str, which method(s) produced the maximum sigma.
    * ``"method_1_circular_variance"`` — dict with ``elbow_sigma``, ``curve``,
      ``variances``.
    * ``"method_2_magnitude_stability"`` — dict with ``stability_sigma``,
      ``curve``, ``magnitudes``.
    * ``"method_3_iqr_crossover"`` — dict with ``crossover_sigma``, ``curve``,
      ``iqr_values``.
    * ``"diagnostics"`` — dict with per-method per-sigma detail.
    """
    if sigma_range is None:
        sigma_range = list(_DEFAULT_SIGMA_RANGE)

    # Ensure sorted for consistent sweep order
    sigma_range_sorted = sorted(sigma_range)

    # --- Step 1: compute gradient sweep ---------------------------------
    gradient_maps = compute_gradient_sweep(head_map, grid, sigma_range_sorted)

    # --- Step 2: extract node IDs and build supporting structures -------
    node_ids = node_df["node_id"].astype(str).tolist()

    sample_map = _build_sample_map_from_node_df(node_df)
    candidate_edges = _build_grid_adjacency_edges(node_df, grid)

    # --- Step 3: run the three estimation methods -----------------------
    cv_result = _circular_variance_sweep(
        gradient_maps, sigma_range_sorted, node_ids
    )
    ms_result = _magnitude_stability_sweep(
        gradient_maps, sigma_range_sorted, node_ids
    )
    iqr_result = _iqr_crossover_sweep(
        gradient_maps, sigma_range_sorted, sample_map, candidate_edges
    )

    # --- Step 4: consensus recommendation --------------------------------
    elbow_sigma = cv_result.get("elbow_sigma", 0)
    stability_sigma = ms_result.get("stability_sigma", 0)
    crossover_sigma = iqr_result.get("crossover_sigma", 0)

    candidates = {
        "circular_variance_elbow": elbow_sigma,
        "magnitude_stability": stability_sigma,
        "iqr_crossover": crossover_sigma,
    }

    max_sigma = max(candidates.values())
    binding = [name for name, val in candidates.items() if val == max_sigma]

    # Round to nearest integer present in the sweep range
    recommended = int(
        min(sigma_range_sorted, key=lambda s: abs(s - max_sigma))
    )

    return {
        "recommended_sigma": recommended,
        "binding_method": binding,
        "binding_detail": {
            "elbow_sigma": elbow_sigma,
            "stability_sigma": stability_sigma,
            "crossover_sigma": crossover_sigma,
        },
        "method_1_circular_variance": cv_result,
        "method_2_magnitude_stability": ms_result,
        "method_3_iqr_crossover": iqr_result,
        "diagnostics": {
            "sigma_range": sigma_range_sorted,
        },
    }


# ===================================================================
# Internal helpers
# ===================================================================


def _extract_cell_index(node_id: str) -> Optional[int]:
    """Extract MODFLOW 1-based cell index from a ``"cell_{int}"`` node ID.

    Returns *None* if the node ID does not follow the expected convention.
    """
    m = re.search(r"cell_(\d+)", node_id)
    if m:
        return int(m.group(1))
    # Fallback: try bare integer string
    try:
        return int(node_id)
    except ValueError:
        return None


def _build_sample_map_from_node_df(
    node_df: pd.DataFrame,
) -> Dict[str, Dict[str, float]]:
    """Convert *node_df* to the sample-map format used by graph methods.

    Each entry is ``{node_id: {"x": x, "y": y}}``.
    Tries ``x`` / ``lon`` and ``y`` / ``lat`` column names.
    """
    sample_map: Dict[str, Dict[str, float]] = {}
    for _, row in node_df.iterrows():
        nid = str(row["node_id"])
        x_val: float = 0.0
        y_val: float = 0.0
        for xk in ("x", "lon", "easting"):
            if xk in row and not pd.isna(row[xk]):
                x_val = float(row[xk])
                break
        for yk in ("y", "lat", "northing"):
            if yk in row and not pd.isna(row[yk]):
                y_val = float(row[yk])
                break
        sample_map[nid] = {"x": x_val, "y": y_val}
    return sample_map


def _build_grid_adjacency_edges(
    node_df: pd.DataFrame,
    grid: GridGeometry,
) -> List[Tuple[str, str]]:
    """Build directed candidate edges from 4-connected grid adjacency.

    For every cell that appears in *node_df*, creates a directed edge to
    each of its cardinal neighbours (N, S, E, W) **within the same layer**
    provided the neighbour cell also exists in *node_df*.

    Returns a list of ``(u_node_id, v_node_id)`` tuples.
    """
    # cell_index (1-based) → node_id
    cell_to_node: Dict[int, str] = {}
    for _, row in node_df.iterrows():
        nid = str(row["node_id"])
        ci = _extract_cell_index(nid)
        if ci is not None:
            cell_to_node[ci] = nid

    cells_per_layer = grid.ncol * grid.nrow
    edges: List[Tuple[str, str]] = []

    for cell_idx, u_nid in cell_to_node.items():
        col, row, lay = _cell_to_ijk(cell_idx, grid.ncol, grid.nrow)
        # 4-connected neighbours in the same layer
        for dc, dr in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            nc = col + dc
            nr = row + dr
            if 0 <= nc < grid.ncol and 0 <= nr < grid.nrow:
                neighbour_cell = lay * cells_per_layer + nr * grid.ncol + nc + 1
                v_nid = cell_to_node.get(neighbour_cell)
                if v_nid is not None:
                    edges.append((u_nid, v_nid))

    return edges


# ===================================================================
# Method 1: Circular Variance Elbow
# ===================================================================


def _circular_variance_sweep(
    gradient_maps: Dict[int, Dict[int, Tuple[float, float, float]]],
    sigma_range: List[int],
    node_ids: List[str],
) -> Dict[str, Any]:
    """Estimate sigma via circular-variance elbow detection.

    For each sigma:

    1. Compute the flow-direction angle θ = atan2(-gy, -gx) at every node.
    2. Compute the circular variance V = 1 - R, where R is the mean
       resultant length of the unit vectors exp(iθ).
       V = 1 means random directions; V = 0 means perfect alignment.
    3. Find the **elbow**: the sigma where ΔV/Δσ (the rate of variance
       reduction) drops below ``_ELBOW_FRACTION`` of its maximum value.

    Parameters
    ----------
    gradient_maps:
        ``{sigma: {cell_index: (gx, gy, gz)}}`` from :func:`compute_gradient_sweep`.
    sigma_range:
        Sorted list of sigma values (same keys as *gradient_maps*).
    node_ids:
        List of node-id strings (e.g. ``"cell_130749"``).

    Returns
    -------
    dict with keys ``"elbow_sigma"``, ``"curve"``, ``"variances"``.
    """
    sigmas_sorted = sorted(sigma_range)
    variances: List[float] = []
    n_nodes_with_gradient: List[int] = []

    for sigma in sigmas_sorted:
        gmap = gradient_maps.get(sigma, {})
        angles_rad: List[float] = []
        count = 0

        for nid in node_ids:
            cell_idx = _extract_cell_index(nid)
            if cell_idx is None:
                continue
            grad = gmap.get(cell_idx)
            if grad is None:
                continue
            gx, gy, _gz = grad
            # Skip nodes with near-zero gradient (no meaningful direction)
            if abs(gx) < 1e-15 and abs(gy) < 1e-15:
                continue
            theta = math.atan2(-gy, -gx)
            angles_rad.append(theta)
            count += 1

        n_nodes_with_gradient.append(count)

        if len(angles_rad) < 2:
            variances.append(float("nan"))
            continue

        # Mean resultant length R
        c = sum(math.cos(a) for a in angles_rad) / len(angles_rad)
        s = sum(math.sin(a) for a in angles_rad) / len(angles_rad)
        R = math.sqrt(c * c + s * s)
        V = 1.0 - R  # circular variance
        variances.append(V)

    # --- Build (sigma, variance) curve (skip NaN entries) ----------------
    curve = [
        (sigmas_sorted[i], variances[i])
        for i in range(len(sigmas_sorted))
        if not math.isnan(variances[i])
    ]

    # --- Elbow detection -------------------------------------------------
    elbow_sigma = _find_variance_elbow(sigmas_sorted, variances)

    return {
        "elbow_sigma": elbow_sigma,
        "curve": curve,
        "variances": variances,
        "n_nodes_with_gradient": n_nodes_with_gradient,
    }


def _find_variance_elbow(
    sigmas: List[int],
    variances: List[float],
) -> int:
    """Find the sigma where the per-step reduction in circular variance
    falls below ``_ELBOW_FRACTION`` of the maximum reduction.

    Works with non-uniform sigma steps by computing ΔV/Δσ.
    """
    valid_indices = [
        i for i in range(len(variances)) if not math.isnan(variances[i])
    ]
    if len(valid_indices) < 3:
        # Not enough data for elbow detection → return smallest non-zero sigma
        return sigmas[1] if len(sigmas) > 1 else sigmas[0]

    # Compute successive reductions per unit sigma
    reductions: List[float] = []
    reduction_sigmas: List[int] = []
    for j in range(1, len(valid_indices)):
        i_prev = valid_indices[j - 1]
        i_curr = valid_indices[j]
        dV = variances[i_prev] - variances[i_curr]  # positive = improvement
        d_sigma = float(sigmas[i_curr] - sigmas[i_prev])
        if d_sigma > 0 and dV > 0:
            reductions.append(dV / d_sigma)
            reduction_sigmas.append(sigmas[i_curr])

    if not reductions:
        return sigmas[1] if len(sigmas) > 1 else sigmas[0]

    max_reduction = max(reductions)
    threshold = _ELBOW_FRACTION * max_reduction

    # Find the first sigma where the reduction rate drops below threshold
    for k, red in enumerate(reductions):
        if red < threshold:
            # The elbow is the sigma *before* the drop-off
            # i.e. the last sigma that still gave good benefit
            if k > 0:
                return reduction_sigmas[k - 1]
            return sigmas[0]

    # No clear elbow → all sigma values still giving good returns
    # Return the highest sigma tested
    return reduction_sigmas[-1] if reduction_sigmas else sigmas[-1]


# ===================================================================
# Method 2: Gradient Magnitude Stability
# ===================================================================


def _magnitude_stability_sweep(
    gradient_maps: Dict[int, Dict[int, Tuple[float, float, float]]],
    sigma_range: List[int],
    node_ids: List[str],
) -> Dict[str, Any]:
    """Estimate sigma via gradient-magnitude stabilisation.

    For each sigma, compute the **median** gradient magnitude across all
    nodes that have a gradient.  The stabilisation sigma is the first sigma
    where the relative change in median magnitude from the previous sigma is
    below ``_STABILITY_THRESHOLD`` (5 %).

    Parameters
    ----------
    gradient_maps:
        ``{sigma: {cell_index: (gx, gy, gz)}}``.
    sigma_range:
        Sorted list of sigma values.
    node_ids:
        List of node-id strings.

    Returns
    -------
    dict with keys ``"stability_sigma"``, ``"curve"``, ``"magnitudes"``.
    """
    sigmas_sorted = sorted(sigma_range)
    median_magnitudes: List[float] = []

    for sigma in sigmas_sorted:
        gmap = gradient_maps.get(sigma, {})
        mags: List[float] = []

        for nid in node_ids:
            cell_idx = _extract_cell_index(nid)
            if cell_idx is None:
                continue
            grad = gmap.get(cell_idx)
            if grad is None:
                continue
            gx, gy, gz = grad
            mags.append(math.sqrt(gx * gx + gy * gy + gz * gz))

        if mags:
            median_magnitudes.append(float(np.median(mags)))
        else:
            median_magnitudes.append(float("nan"))

    # --- Build curve -----------------------------------------------------
    curve = [
        (sigmas_sorted[i], median_magnitudes[i])
        for i in range(len(sigmas_sorted))
        if not math.isnan(median_magnitudes[i])
    ]

    # --- Stabilisation detection -----------------------------------------
    stability_sigma = _find_stability_sigma(sigmas_sorted, median_magnitudes)

    return {
        "stability_sigma": stability_sigma,
        "curve": curve,
        "magnitudes": median_magnitudes,
    }


def _find_stability_sigma(
    sigmas: List[int],
    magnitudes: List[float],
) -> int:
    """Return the first sigma where the relative change in median magnitude
    drops below ``_STABILITY_THRESHOLD``.

    Relative change = |mag(σ) - mag(σ_prev)| / mag(σ_prev).
    """
    valid_indices = [
        i for i in range(len(magnitudes)) if not math.isnan(magnitudes[i])
    ]
    if len(valid_indices) < 2:
        return sigmas[0]

    for j in range(1, len(valid_indices)):
        i_prev = valid_indices[j - 1]
        i_curr = valid_indices[j]
        mag_prev = magnitudes[i_prev]
        mag_curr = magnitudes[i_curr]
        if mag_prev <= 1e-15:
            continue
        rel_change = abs(mag_curr - mag_prev) / mag_prev
        if rel_change < _STABILITY_THRESHOLD:
            return sigmas[i_curr]

    # Never stabilised → return the largest sigma tested
    return sigmas[valid_indices[-1]]


# ===================================================================
# Method 3: Projection IQR Crossover
# ===================================================================


def _iqr_crossover_sweep(
    gradient_maps: Dict[int, Dict[int, Tuple[float, float, float]]],
    sigma_range: List[int],
    sample_map: Dict[str, Dict[str, float]],
    candidate_edges: List[Tuple[str, str]],
) -> Dict[str, Any]:
    """Estimate sigma via projection IQR crossover.

    For each sigma, compute ``proj = dot(-∇h(u), d_uv)`` for every
    candidate edge where the source node *u* has a gradient and both
    endpoints have coordinates.  The crossover sigma is the **smallest**
    sigma at which:

    * ``IQR(projections) > 1e-8``, AND
    * ``n_projections >= 4``.

    This is the point where the gradient field provides sufficient
    directional discrimination for projected-gradient priors.

    Parameters
    ----------
    gradient_maps:
        ``{sigma: {cell_index: (gx, gy, gz)}}``.
    sigma_range:
        Sorted list of sigma values.
    sample_map:
        ``{node_id: {"x": x, "y": y}}``.
    candidate_edges:
        List of ``(u_node_id, v_node_id)`` tuples.

    Returns
    -------
    dict with keys ``"crossover_sigma"``, ``"curve"``, ``"iqr_values"``,
    ``"n_projections"``.
    """
    sigmas_sorted = sorted(sigma_range)
    iqr_values: List[float] = []
    n_projections_list: List[int] = []

    # Pre-compute edge vectors once (they don't depend on sigma)
    edge_vectors: Dict[Tuple[str, str], Optional[Tuple[float, float]]] = {}
    for u, v in candidate_edges:
        su = sample_map.get(u)
        sv = sample_map.get(v)
        if su is None or sv is None:
            edge_vectors[(u, v)] = None
            continue
        dx = sv["x"] - su["x"]
        dy = sv["y"] - su["y"]
        dist = math.sqrt(dx * dx + dy * dy)
        if dist < 1e-9:
            edge_vectors[(u, v)] = None
            continue
        edge_vectors[(u, v)] = (dx / dist, dy / dist)

    for sigma in sigmas_sorted:
        gmap = gradient_maps.get(sigma, {})
        projections: List[float] = []

        for u, v in candidate_edges:
            d_uv = edge_vectors.get((u, v))
            if d_uv is None:
                continue

            cell_idx = _extract_cell_index(u)
            if cell_idx is None:
                continue
            grad = gmap.get(cell_idx)
            if grad is None:
                continue

            gx, gy, _gz = grad
            # Flow-direction vector = -grad(h)
            flow_x = -gx
            flow_y = -gy

            # Ignore near-zero flow (no meaningful direction)
            if abs(flow_x) < 1e-15 and abs(flow_y) < 1e-15:
                continue

            proj = flow_x * d_uv[0] + flow_y * d_uv[1]
            projections.append(proj)

        n_proj = len(projections)
        n_projections_list.append(n_proj)

        if n_proj >= _MIN_PROJECTIONS:
            proj_arr = np.array(projections)
            q1 = float(np.percentile(proj_arr, 25))
            q3 = float(np.percentile(proj_arr, 75))
            iqr = q3 - q1
            iqr_values.append(iqr)
        else:
            iqr_values.append(float("nan"))

    # --- Build curve -----------------------------------------------------
    curve = [
        (sigmas_sorted[i], iqr_values[i])
        for i in range(len(sigmas_sorted))
        if not math.isnan(iqr_values[i])
    ]

    # --- Crossover detection ---------------------------------------------
    crossover_sigma = _find_iqr_crossover(
        sigmas_sorted, iqr_values, n_projections_list
    )

    return {
        "crossover_sigma": crossover_sigma,
        "curve": curve,
        "iqr_values": iqr_values,
        "n_projections": n_projections_list,
    }


def _find_iqr_crossover(
    sigmas: List[int],
    iqr_values: List[float],
    n_projections: List[int],
) -> int:
    """Return the smallest sigma where IQR > ``_MIN_IQR`` AND
    ``n_projections >= _MIN_PROJECTIONS``.

    If no sigma passes both criteria, return the largest sigma tested.
    """
    for i in range(len(sigmas)):
        if n_projections[i] >= _MIN_PROJECTIONS:
            if not math.isnan(iqr_values[i]) and iqr_values[i] > _MIN_IQR:
                return sigmas[i]

    # No crossover found → return the largest sigma (most smoothing gives
    # the widest projection spread)
    return sigmas[-1]


# ===================================================================
# CLI entry point
# ===================================================================


def _cli() -> None:
    """Command-line entry point for the sigma estimation module.

    Usage::

        python -m hydrosheaf.physics.head_gradient_scale \\
            path/to/model.fhd \\
            --ncol 183 --nrow 202 --nlay 8 \\
            --dx 110.0 --dy 73.333 \\
            --rotation -12.0 \\
            --nodes path/to/node_mapping.csv

    The node-mapping CSV must contain ``node_id``, ``x``, ``y`` columns.
    """
    import argparse

    parser = argparse.ArgumentParser(
        description="Estimate optimal Gaussian smoothing sigma for head gradients",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Example:\n"
            "  python -m hydrosheaf.physics.head_gradient_scale base.fhd \\\n"
            "      --ncol 183 --nrow 202 --nlay 8 --dx 110.0 --dy 73.333 \\\n"
            "      --rotation -12.0 --nodes modpath_node_mapping.csv"
        ),
    )
    parser.add_argument("fhd", type=str, help="Path to MODFLOW .fhd file")
    parser.add_argument("--ncol", type=int, required=True, help="Number of columns")
    parser.add_argument("--nrow", type=int, required=True, help="Number of rows")
    parser.add_argument("--nlay", type=int, default=1, help="Number of layers")
    parser.add_argument("--dx", type=float, required=True, help="Cell width (model units)")
    parser.add_argument("--dy", type=float, required=True, help="Cell height (model units)")
    parser.add_argument(
        "--rotation", type=float, default=0.0, help="Grid rotation ANGROT (degrees CCW)"
    )
    parser.add_argument("--origin-x", type=float, default=0.0, help="Grid origin X")
    parser.add_argument("--origin-y", type=float, default=0.0, help="Grid origin Y")
    parser.add_argument(
        "--nodes", type=str, required=True, help="Path to node-mapping CSV"
    )
    parser.add_argument(
        "--sigma-range",
        type=str,
        default=None,
        help="Comma-separated sigma values (default: 0,1,2,3,4,6,8,12,16,24,32)",
    )
    args = parser.parse_args()

    # --- Parse inputs ----------------------------------------------------
    from .modflow_head import parse_fhd

    fhd_path = Path(args.fhd)
    if not fhd_path.exists():
        print(f"Error: .fhd file not found: {fhd_path}", file=sys.stderr)
        sys.exit(1)

    print(f"Parsing {fhd_path} ...")
    head_map = parse_fhd(fhd_path)
    print(f"  {len(head_map):,} active cells with head values")

    nodes_path = Path(args.nodes)
    if not nodes_path.exists():
        print(f"Error: node CSV not found: {nodes_path}", file=sys.stderr)
        sys.exit(1)

    node_df = pd.read_csv(nodes_path)
    print(f"  {len(node_df):,} node-mapping rows loaded")

    sigma_range: Optional[List[int]] = None
    if args.sigma_range:
        sigma_range = [int(s.strip()) for s in args.sigma_range.split(",")]

    grid = GridGeometry(
        ncol=args.ncol,
        nrow=args.nrow,
        nlay=args.nlay,
        dx=args.dx,
        dy=args.dy,
        rotation_deg=args.rotation,
        origin_x=args.origin_x,
        origin_y=args.origin_y,
    )
    print(f"  Grid: {grid.ncol} x {grid.nrow} x {grid.nlay}, "
          f"dx={grid.dx}, dy={grid.dy}, rotation={grid.rotation_deg} deg")

    # --- Run estimation --------------------------------------------------
    result = estimate_smoothing_scale(head_map, grid, node_df, sigma_range)

    # --- Print formatted table -------------------------------------------
    _print_results(result)


def _print_results(result: Dict[str, Any]) -> None:
    """Pretty-print the sigma estimation results as a formatted table."""
    sigma_range = result["diagnostics"]["sigma_range"]

    # Gather per-method data
    cv = result["method_1_circular_variance"]
    ms = result["method_2_magnitude_stability"]
    iqr = result["method_3_iqr_crossover"]

    variances = cv.get("variances", [])
    magnitudes = ms.get("magnitudes", [])
    iqr_vals = iqr.get("iqr_values", [])
    n_projs = iqr.get("n_projections", [])

    print()
    print("=" * 90)
    print("  SIGMA ESTIMATION RESULTS")
    print("=" * 90)
    print()

    # Summary header
    header = (
        f"{'sigma':>6s}  "
        f"{'circ_var':>10s}  "
        f"{'median_mag':>12s}  "
        f"{'iqr':>12s}  "
        f"{'n_proj':>8s}"
    )
    print(header)
    print("-" * len(header))

    for i, sigma in enumerate(sigma_range):
        cv_str = f"{variances[i]:.6f}" if i < len(variances) and not math.isnan(variances[i]) else "      N/A"
        ms_str = f"{magnitudes[i]:.6e}" if i < len(magnitudes) and not math.isnan(magnitudes[i]) else "       N/A"
        iqr_str = f"{iqr_vals[i]:.6e}" if i < len(iqr_vals) and not math.isnan(iqr_vals[i]) else "        N/A"
        np_str = f"{n_projs[i]:d}" if i < len(n_projs) else "    N/A"

        print(f"{sigma:6d}  {cv_str:>10s}  {ms_str:>12s}  {iqr_str:>12s}  {np_str:>8s}")

    print()
    print("-" * 90)
    print(f"  Method 1 (Circular Variance Elbow) sigma  = {cv.get('elbow_sigma', 'N/A')}")
    print(f"  Method 2 (Magnitude Stability) sigma       = {ms.get('stability_sigma', 'N/A')}")
    print(f"  Method 3 (IQR Crossover) sigma             = {iqr.get('crossover_sigma', 'N/A')}")
    print(f"  ---")
    print(f"  RECOMMENDED sigma = {result['recommended_sigma']}")
    print(f"  Binding method(s): {', '.join(result['binding_method'])}")
    print("=" * 90)


if __name__ == "__main__":
    _cli()
