"""Hydraulic Hodge diagnostics for head-constrained topology inference.

This layer is intentionally separate from the geochemical affine sheaf.
It treats hydraulic head as a scalar potential on nodes and provides
graph-level cost terms based on physically meaningful diagnostics:

- **Local head-plane residuals**: penalises nodes whose head values don't
  fit a local planar trend, flagging possible data errors or mis-placed
  topology.

The original scalar-potential "obstruction energy" (d² = 0 for an exact
1-form) is retained behind an ``_experimental_`` prefix for future
multi-field sheaf work where α ≠ 1 and offsets are genuinely non-conservative.

The projected-gradient alignment term that existed in earlier versions has
been **removed** — directional prior responsibility lives entirely in
:mod:`~hydrosheaf.graph.flow_direction`.
"""

from __future__ import annotations

import math
from statistics import median
from typing import Any, Dict, Iterable, Mapping, Optional, Sequence, Tuple

import numpy as np

from ..config import Config
from ..graph.types import Edge
from ..log import get_logger
from .cohomology import (
    build_coboundary_matrix,
    build_rhs_vector,
    compute_cohomology,
    compute_cycle_obstruction,
    compute_edge_leverage,
)
from .directed_section import DirectedEdgeMap

logger = get_logger("sheaf.hydraulic_hodge")


def _safe_float(value: object) -> Optional[float]:
    if value is None:
        return None
    try:
        result = float(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return None
    if not math.isfinite(result):
        return None
    return result


# TODO: consolidate with flow_direction._sigmoid into shared utility
def _sigmoid(x: float) -> float:
    """Numerically stable logistic sigmoid."""
    if x >= 0:
        return 1.0 / (1.0 + math.exp(-x))
    ex = math.exp(x)
    return ex / (1.0 + ex)


def _sample_map(samples: object) -> Dict[str, Mapping[str, object]]:
    if isinstance(samples, Mapping):
        return {str(k): v for k, v in samples.items()}
    if isinstance(samples, Sequence):
        mapping: Dict[str, Mapping[str, object]] = {}
        for row in samples:
            if not isinstance(row, Mapping):
                continue
            site_id = row.get("site_id") or row.get("node_id") or row.get("sample_id")
            if site_id is None:
                continue
            mapping[str(site_id)] = row
        return mapping
    raise TypeError("Unsupported samples input type.")


def _resolve_head_value(sample: Mapping[str, object], config: Config) -> Optional[float]:
    primary = getattr(config, "hydraulic_hodge_head_key", "hydraulic_head")
    for key in (
        primary,
        "hydraulic_head",
        "head",
        "head_meas",
    ):
        value = _safe_float(sample.get(key))
        if value is not None:
            return value
    if getattr(config, "hydraulic_hodge_fallback_to_elevation", True):
        return _safe_float(sample.get("elevation"))
    return None


def extract_node_heads(
    samples: object,
    config: Config,
) -> Dict[str, float]:
    sample_map = _sample_map(samples)
    heads: Dict[str, float] = {}
    for node_id, sample in sample_map.items():
        head = _resolve_head_value(sample, config)
        if head is not None:
            heads[node_id] = head
    return heads


def _edge_distance_km(
    edge: Edge,
    sample_map: Mapping[str, Mapping[str, object]],
) -> Optional[float]:
    attrs = edge.attrs or {}
    distance = _safe_float(attrs.get("distance_km"))
    if distance is not None and distance > 0:
        return distance

    sample_u = sample_map.get(edge.u)
    sample_v = sample_map.get(edge.v)
    if sample_u is None or sample_v is None:
        return None

    x_u = _safe_float(sample_u.get("lon"))
    y_u = _safe_float(sample_u.get("lat"))
    x_v = _safe_float(sample_v.get("lon"))
    y_v = _safe_float(sample_v.get("lat"))
    if None in (x_u, y_u, x_v, y_v):
        return None

    assert x_u is not None and y_u is not None and x_v is not None and y_v is not None
    # Haversine — convert lon/lat degree-difference to km
    phi1 = math.radians(y_u)
    phi2 = math.radians(y_v)
    dphi = math.radians(y_v - y_u)
    dlam = math.radians(x_v - x_u)
    a = math.sin(dphi / 2) ** 2 + math.cos(phi1) * math.cos(phi2) * math.sin(dlam / 2) ** 2
    return 2.0 * 6371.0 * math.asin(math.sqrt(a))


def infer_reference_distance_km(
    edges: Iterable[Edge],
    sample_map: Mapping[str, Mapping[str, object]],
    config: Config,
) -> float:
    fixed = _safe_float(getattr(config, "hydraulic_hodge_reference_distance_km", 0.0))
    if fixed is not None and fixed > 0:
        return fixed

    distances = []
    for edge in edges:
        distance = _edge_distance_km(edge, sample_map)
        if distance is not None and distance > 0:
            distances.append(distance)

    if distances:
        return float(median(distances))
    return 1.0


# ---------------------------------------------------------------------------
# Experimental: scalar-potential obstruction energy (d² = 0 by construction)
# ---------------------------------------------------------------------------


def _experimental_obstruction_energy(edge_maps: Sequence[DirectedEdgeMap]) -> float:
    """Experimental obstruction energy — NOT used in the production cost path.

    For a scalar potential (head) with α = 1.0 and offset = -(Δh/d)·d_ref,
    the coboundary equation Dx - b = 0 is solved exactly by x_i = h_i.
    The residual is zero by construction for any graph on any head values.

    This function is preserved for future multi-field sheaf work where
    α ≠ 1 and offsets encode non-conservative transport/reaction relations.
    """
    if not edge_maps:
        return 0.0
    D = build_coboundary_matrix(edge_maps, dim=1)
    b = build_rhs_vector(edge_maps, dim=1)
    if D.shape[0] == 0:
        return 0.0
    try:
        x_sol, _, _, _ = np.linalg.lstsq(D.toarray(), b, rcond=None)
        residual = D @ x_sol - b
        return float(np.sum(residual ** 2))
    except Exception:
        return float(np.sum(b ** 2))


# ---------------------------------------------------------------------------
# Local head-plane residuals
# ---------------------------------------------------------------------------


def _node_coordinates(
    node_id: str,
    sample_map: Mapping[str, Mapping[str, object]],
) -> Optional[Tuple[float, float]]:
    """Extract (x, y) in consistent units from a sample dict."""
    sample = sample_map.get(node_id)
    if sample is None:
        return None
    for xk in ("lon", "easting", "x"):
        for yk in ("lat", "northing", "y"):
            rx = sample.get(xk)
            ry = sample.get(yk)
            if rx is not None and ry is not None:
                try:
                    return float(rx), float(ry)  # type: ignore[arg-type]
                except (TypeError, ValueError):
                    pass
    return None


def compute_local_head_plane_residuals(
    head_map: Dict[str, float],
    sample_map: Mapping[str, Mapping[str, object]],
    n_neighbors: int = 8,
) -> Dict[str, float]:
    """Compute local head-plane residuals for every node.

    For each node, fit a plane :math:`h(x,y) = a x + b y + c` to the *k*
    nearest neighbours' head values via least squares.  The residual is the
    absolute difference between the node's actual head and the planar
    prediction, normalised by the local head standard deviation.

    A large residual means the node's head doesn't fit the local spatial
    trend — possible data error, local sink/source, or incorrect topology
    placement.

    This is an **external reference**: the plane is fit from *neighbouring*
    data, not from the node itself, so it avoids the d² = 0 degeneracy that
    plagues self-consistent scalar-potential diagnostics.

    Parameters
    ----------
    head_map:
        ``{node_id: head_value}`` for all nodes.
    sample_map:
        ``{node_id: sample_dict}`` providing coordinate data.
    n_neighbors:
        Number of nearest spatial neighbours to use for the plane fit
        (default 8).  Must be ≥ 3.

    Returns
    -------
    ``{node_id: residual}`` where *residual* is in units of local head
    standard deviation.  Nodes without enough neighbours get NaN.
    """
    if n_neighbors < 3:
        n_neighbors = 3

    # Build coordinate + head arrays
    node_ids: list[str] = []
    xy: list[Tuple[float, float]] = []
    h_vals: list[float] = []
    for nid, h in head_map.items():
        c = _node_coordinates(nid, sample_map)
        if c is not None and math.isfinite(h):
            node_ids.append(nid)
            xy.append(c)
            h_vals.append(h)

    if len(node_ids) < n_neighbors + 1:
        return {}

    xy_arr = np.array(xy, dtype=np.float64)
    h_arr = np.array(h_vals, dtype=np.float64)

    # For each node, find k nearest spatial neighbours
    residuals: Dict[str, float] = {}
    for i, nid in enumerate(node_ids):
        # Euclidean distances to all other nodes
        dists = np.sqrt(np.sum((xy_arr - xy_arr[i]) ** 2, axis=1))
        # Exclude self
        order = np.argsort(dists)
        # Skip self (distance 0)
        nn_idx = [j for j in order if j != i][:n_neighbors]
        if len(nn_idx) < 3:
            continue

        # Fit plane: h = a*x + b*y + c
        X_nn = np.column_stack([
            xy_arr[nn_idx, 0],
            xy_arr[nn_idx, 1],
            np.ones(len(nn_idx)),
        ])
        h_nn = h_arr[nn_idx]

        try:
            coeffs, _, _, _ = np.linalg.lstsq(X_nn, h_nn, rcond=None)
        except np.linalg.LinAlgError:
            continue

        # Predicted head at the node
        h_pred = float(coeffs[0] * xy_arr[i, 0] + coeffs[1] * xy_arr[i, 1] + coeffs[2])
        h_actual = h_arr[i]

        # Normalise by local standard deviation
        local_std = float(np.std(h_nn))
        if local_std < 1e-8:
            local_std = 1.0

        residuals[nid] = abs(h_actual - h_pred) / local_std

    return residuals


# ---------------------------------------------------------------------------
# Graph cost function
# ---------------------------------------------------------------------------


def head_hodge_graph_cost(
    edges: Sequence[Edge],
    sample_map: Mapping[str, Mapping[str, object]],
    config: Config,
    reference_distance_km: Optional[float] = None,
    gradient_map: Optional[Mapping[str, tuple[float, float, float]]] = None,
    local_residuals: Optional[Dict[str, float]] = None,
) -> float:
    """Graph cost based on hydraulic-head consistency.

    Replaces the degenerate scalar-potential obstruction energy with
    physically meaningful terms:

    1. **Local head-plane residuals** (when *local_residuals*
       provided): sum of incident-node residuals, penalising nodes whose
       head values don't fit the local spatial trend.

    2. **Fallback**: zero-information uniform cost when no diagnostic
       is available.

    The uphill penalty that existed in earlier versions is **removed** —
    directional prior responsibility lives entirely in
    :mod:`~hydrosheaf.graph.flow_direction`.

    .. note::

       *gradient_map* is accepted for signature compatibility but is
       **not consumed** — the directional prior responsibility lives in
       :mod:`~hydrosheaf.graph.flow_direction`, which encodes edge
       direction as ``edge_confidence``.

    .. note::

       *local_residuals* are static node-level diagnostics computed once
       before MCMC (not per-graph). They should be interpreted as a
       **data-quality flag** rather than a dynamic edge-rejection signal.
       See :func:`compute_local_head_plane_residuals` for computational
       details — the plane is fit from neighbouring data, not the node
       itself, but the residual does not change when edges are rewired
       during MCMC.

    Parameters
    ----------
    edges:
        Current graph edges.
    sample_map:
        ``{node_id: sample_dict}`` for coordinate lookups.
    config:
        Hydrosheaf Config.
    reference_distance_km:
        Not used in the new cost path (kept for signature compat).
    gradient_map:
        ``{node_id: (gx, gy, gz)}`` — **not consumed** (kept for
        signature compatibility; see note above).
    local_residuals:
        ``{node_id: residual}`` from :func:`compute_local_head_plane_residuals`.

    Returns
    -------
    scalar cost (non-negative).
    """
    if not edges:
        return 0.0

    residual_weight = float(getattr(config, "head_plane_residual_weight", 1.0) or 1.0)
    n_edges = max(1, len(edges))

    cost = 0.0

    # --- Local head-plane residuals ---
    if local_residuals is not None and residual_weight > 0:
        residual_cost = 0.0
        for edge in edges:
            r_u = local_residuals.get(str(edge.u), 0.0)
            r_v = local_residuals.get(str(edge.v), 0.0)
            residual_cost += (r_u + r_v) / 2.0
        cost += residual_weight * residual_cost / n_edges

    return cost


# ---------------------------------------------------------------------------
# Diagnostics (preserve structure, update content)
# ---------------------------------------------------------------------------


def _build_head_maps(
    edges: Sequence[Edge],
    sample_map: Mapping[str, Mapping[str, object]],
    config: Config,
    reference_distance_km: Optional[float] = None,
) -> Dict[str, DirectedEdgeMap]:
    """Build 1-D head edge maps for diagnostic purposes.

    Kept for backward compatibility with :func:`compute_head_hodge_diagnostics`
    and any callers that rely on the head-map attributes being attached to edges.
    The *min_drop* clipping is removed — edge maps encode raw head drops.
    """
    heads = extract_node_heads(sample_map, config)
    ref_distance = reference_distance_km
    if ref_distance is None or ref_distance <= 0:
        ref_distance = infer_reference_distance_km(edges, sample_map, config)

    maps: Dict[str, DirectedEdgeMap] = {}
    for edge in edges:
        h_u = heads.get(edge.u)
        h_v = heads.get(edge.v)
        if h_u is None or h_v is None:
            continue

        distance = _edge_distance_km(edge, sample_map)
        if distance is None or distance <= 0:
            continue

        delta_h = h_u - h_v
        # Raw head drop — no min_drop clipping to avoid artificial inconsistency
        standardized_drop = (delta_h / distance) * ref_distance

        attrs = edge.attrs or {}
        weight = _safe_float(attrs.get("edge_confidence"))
        if weight is None:
            weight = _safe_float(attrs.get("p_uv"))
        if weight is None:
            weight = 1.0
        weight = min(1.0, max(1e-6, weight))

        maps[edge.edge_id] = DirectedEdgeMap(
            edge=edge,
            alpha=1.0,
            offset=[-standardized_drop],
            weight=weight,
            objective=0.0,
            transport_model="head_hodge",
            endmember_id=None,
            residual_norm=0.0,
        )

        attrs = dict(attrs)
        attrs["head_hodge_delta_h"] = delta_h
        attrs["head_hodge_distance_km"] = distance
        attrs["head_hodge_standardized_drop"] = standardized_drop
        attrs["head_hodge_reference_distance_km"] = ref_distance
        attrs["head_hodge_direction_valid"] = delta_h > 0
        edge.attrs = attrs

    return maps


def compute_head_hodge_diagnostics(
    edges: Sequence[Edge],
    sample_map: Mapping[str, Mapping[str, object]],
    config: Config,
    reference_distance_km: Optional[float] = None,
    compute_leverage: bool = True,
) -> Dict[str, object]:
    """Compute hydraulic-head diagnostics for a graph.

    The obstruction-related fields are labelled ``experimental_*`` to
    signal that they are degenerate for scalar potentials (d² = 0) and
    should not be used as production physics terms.

    Returns
    -------
    dict with keys:
        reference_distance_km : float
        direction_violation_count : int
            Number of edges with Δh ≤ 0 (uphill or flat).
        local_plane_residual_mean : float
            Mean local head-plane residual, if computed.
        local_plane_residual_max : float
        experimental_obstruction_energy : float
            Scalar-potential obstruction energy (degenerate; d² = 0).
        experimental_cycle_obstruction_max : float
        experimental_cycle_count : int
        experimental_obstruction_leverage : dict
        normalized_cost : float
            Cost from :func:`head_hodge_graph_cost`.
    """
    empty = {
        "reference_distance_km": 0.0,
        "direction_violation_count": 0,
        "local_plane_residual_mean": 0.0,
        "local_plane_residual_max": 0.0,
        "experimental_obstruction_energy": 0.0,
        "experimental_cycle_obstruction_max": 0.0,
        "experimental_cycle_count": 0,
        "experimental_obstruction_leverage": {},
        "normalized_cost": 0.0,
    }

    if not edges:
        return empty

    edge_maps = _build_head_maps(edges, sample_map, config, reference_distance_km)
    maps_list = list(edge_maps.values())
    if not maps_list:
        return empty

    ref_distance = _safe_float(
        (maps_list[0].edge.attrs or {}).get("head_hodge_reference_distance_km")
    ) or 0.0

    # --- Experimental: scalar-potential obstruction (d² = 0) ---
    coh = compute_cohomology(maps_list, dim=1)
    cycle_info = compute_cycle_obstruction(maps_list, dim=1)
    leverage: Dict[str, float] = {}
    if compute_leverage and len(maps_list) <= 200:
        leverage = compute_edge_leverage(maps_list, dim=1)
    elif compute_leverage:
        logger.info(
            "Skipping experimental Hodge leverage (too many edges: %d).",
            len(maps_list),
        )

    # --- Direction violations ---
    direction_violation_count = 0
    for edge in edges:
        delta_h = _safe_float((edge.attrs or {}).get("head_hodge_delta_h"))
        if delta_h is not None and delta_h <= 0:
            direction_violation_count += 1

    # --- Local plane residuals ---
    head_map = extract_node_heads(sample_map, config)
    residuals = compute_local_head_plane_residuals(
        head_map, sample_map,
        n_neighbors=int(getattr(config, "head_plane_residual_neighbors", 8) or 8),
    )
    res_vals = list(residuals.values())
    local_mean = float(np.mean(res_vals)) if res_vals else 0.0
    local_max = float(np.max(res_vals)) if res_vals else 0.0

    normalized_cost = head_hodge_graph_cost(
        edges, sample_map, config, reference_distance_km=ref_distance,
    )

    return {
        "reference_distance_km": ref_distance,
        "direction_violation_count": direction_violation_count,
        "local_plane_residual_mean": local_mean,
        "local_plane_residual_max": local_max,
        "experimental_obstruction_energy": float(coh["obstruction_energy"]),
        "experimental_cycle_obstruction_max": float(cycle_info["cycle_obstruction_max"]),
        "experimental_cycle_count": int(cycle_info["cycle_count"]),
        "experimental_obstruction_leverage": leverage,
        "normalized_cost": normalized_cost,
    }


def compute_head_hodge_edge_penalties(
    edges: Sequence[Edge],
    sample_map: Mapping[str, Mapping[str, object]],
    config: Config,
    reference_distance_km: Optional[float] = None,
    local_residuals: Optional[Dict[str, float]] = None,
) -> Dict[str, float]:
    """Compute per-edge penalties from local head-plane residuals.

    Each edge inherits the average residual of its incident nodes, scaled
    by *hydraulic_hodge_leverage_weight*.  When *local_residuals* is not
    provided it is computed on the fly.

    Returns
    -------
    ``{edge_id: penalty}`` (non-negative).
    """
    if local_residuals is None:
        head_map = extract_node_heads(sample_map, config)
        n_neighbors = int(getattr(config, "head_plane_residual_neighbors", 8) or 8)
        local_residuals = compute_local_head_plane_residuals(
            head_map, sample_map, n_neighbors=n_neighbors,
        )

    weight = float(getattr(config, "hydraulic_hodge_leverage_weight", 0.5) or 0.5)
    penalties: Dict[str, float] = {}
    for edge in edges:
        r_u = local_residuals.get(str(edge.u), 0.0) if local_residuals else 0.0
        r_v = local_residuals.get(str(edge.v), 0.0) if local_residuals else 0.0
        penalty = max(0.0, (r_u + r_v) / 2.0) * weight
        if penalty > 0:
            penalties[edge.edge_id] = penalty
    return penalties


def attach_head_hodge_attrs(
    selected_edges: Sequence[Edge],
    sample_map: Mapping[str, Mapping[str, object]],
    config: Config,
    reference_distance_km: Optional[float] = None,
    compute_leverage: bool = True,
) -> None:
    """Compute head diagnostics and attach per-edge attributes to selected edges.

    Attributes written (in-place via ``edge.attrs``):
        head_hodge_direction_violation_count : int
        head_hodge_normalized_cost : float
        head_hodge_reference_distance_km : float
        head_hodge_local_plane_residual_mean : float
        head_hodge_local_plane_residual_max : float
        head_hodge_experimental_obstruction_energy : float
        head_hodge_experimental_obstruction_leverage : float or None
        head_hodge_experimental_cycle_obstruction_max : float
        head_hodge_experimental_cycle_count : int
    """
    diagnostics = compute_head_hodge_diagnostics(
        list(selected_edges),
        sample_map,
        config,
        reference_distance_km=reference_distance_km,
        compute_leverage=compute_leverage,
    )
    leverage = diagnostics.get("experimental_obstruction_leverage", {})
    if not isinstance(leverage, Mapping):
        leverage = {}

    for edge in selected_edges:
        attrs = dict(edge.attrs or {})
        attrs["head_hodge_direction_violation_count"] = diagnostics["direction_violation_count"]
        attrs["head_hodge_normalized_cost"] = diagnostics["normalized_cost"]
        attrs["head_hodge_reference_distance_km"] = diagnostics["reference_distance_km"]
        attrs["head_hodge_local_plane_residual_mean"] = diagnostics["local_plane_residual_mean"]
        attrs["head_hodge_local_plane_residual_max"] = diagnostics["local_plane_residual_max"]
        attrs["head_hodge_experimental_obstruction_energy"] = diagnostics["experimental_obstruction_energy"]
        attrs["head_hodge_experimental_obstruction_leverage"] = leverage.get(edge.edge_id)
        attrs["head_hodge_experimental_cycle_obstruction_max"] = diagnostics["experimental_cycle_obstruction_max"]
        attrs["head_hodge_experimental_cycle_count"] = diagnostics["experimental_cycle_count"]
        edge.attrs = attrs
