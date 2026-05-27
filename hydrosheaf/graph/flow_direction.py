"""Steepest-descent flow-direction priors for topology inference.

Assigns each candidate directed edge a Bernoulli prior probability based on
how well the edge direction aligns with the steepest head-descent direction at
the source node.  The resulting scores are written to ``edge.attrs["edge_confidence"]``
so they feed directly into the Metropolis-Hastings topology posterior.

Physics rationale
-----------------
MODPATH advective particles follow the hydraulic head gradient exactly.
A false-positive edge u→v is one that points laterally or nearly laterally —
it passes the simple "downhill" test but is not on the primary flow path.
The steepest-descent score penalises such edges by comparing their direction
against the maximum-gradient direction at the source node.
"""

from __future__ import annotations

import math
from typing import Any, Dict, Mapping, Optional, Sequence

from .types import Edge


_PRIOR_FLOOR = 0.05   # minimum prior for any edge (keeps MCMC ergodic)
_PRIOR_CEIL = 0.95    # maximum prior (avoids log(1) = 0 trap)
_UNINFORMATIVE = 0.5  # returned when head data are absent for a node


def _get_config_float(config: Any, name: str, default: float) -> float:
    val = getattr(config, name, default)
    if val is None:
        return default
    if isinstance(val, dict):
        for key in ("value", "weight", "val", "default"):
            if key in val:
                try:
                    return float(val[key])
                except (ValueError, TypeError):
                    pass
        for v in val.values():
            try:
                return float(v)
            except (ValueError, TypeError):
                pass
        return default
    try:
        return float(val)
    except (ValueError, TypeError):
        return default


def _resolve_head(sample: Mapping[str, object], config: Any) -> Optional[float]:
    """Extract a head value from a sample dict, with elevation fallback."""
    primary = getattr(config, "hydraulic_hodge_head_key", "hydraulic_head")
    for key in (primary, "hydraulic_head", "head", "head_meas"):
        raw = sample.get(key)
        if raw is not None:
            try:
                v = float(raw)  # type: ignore[arg-type]
                if math.isfinite(v):
                    return v
            except (TypeError, ValueError):
                pass
    if getattr(config, "hydraulic_hodge_fallback_to_elevation", True):
        raw = sample.get("elevation")
        if raw is not None:
            try:
                v = float(raw)  # type: ignore[arg-type]
                if math.isfinite(v):
                    return v
            except (TypeError, ValueError):
                pass
    return None


def _resolve_coords(sample: Mapping[str, object]) -> Optional[tuple[float, float, float]]:
    """Extract (x, y, z) from a sample dict."""
    for xk in ("x", "easting", "lon"):
        for yk in ("y", "northing", "lat"):
            rx, ry = sample.get(xk), sample.get(yk)
            if rx is None or ry is None:
                continue
            try:
                xv, yv = float(rx), float(ry)  # type: ignore[arg-type]
            except (TypeError, ValueError):
                continue
            zv = 0.0
            for zk in ("elevation", "z", "depth"):
                rz = sample.get(zk)
                if rz is not None:
                    try:
                        zv = float(rz)  # type: ignore[arg-type]
                        break
                    except (TypeError, ValueError):
                        pass
            return xv, yv, zv
    return None


def _dot3(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def compute_steepest_descent_priors(
    candidate_edges: Sequence[Edge],
    sample_map: Mapping[str, Mapping[str, object]],
    config: Any,
) -> Dict[str, float]:
    """Compute steepest-descent alignment scores for every candidate edge.

    Parameters
    ----------
    candidate_edges:
        Universe of directed Edge objects.
    sample_map:
        ``{node_id: sample_dict}`` providing head and coordinate data.
    config:
        Hydrosheaf Config (reads ``steepest_descent_angular_weight``,
        ``steepest_descent_head_weight``, ``hydraulic_hodge_head_key``,
        ``hydraulic_hodge_fallback_to_elevation``).

    Returns
    -------
    dict mapping ``edge_id -> float`` score in [0.05, 0.95].
    High score = well-aligned with steepest head-descent at source node u.
    """
    angular_w = _get_config_float(config, "steepest_descent_angular_weight", 1.0)
    head_w = _get_config_float(config, "steepest_descent_head_weight", 1.0)

    # --- Resolve head and coordinates for every node in the sample map ---
    heads: Dict[str, float] = {}
    coords: Dict[str, tuple[float, float, float]] = {}
    for node_id, sample in sample_map.items():
        h = _resolve_head(sample, config)
        if h is not None:
            heads[node_id] = h
        c = _resolve_coords(sample)
        if c is not None:
            coords[node_id] = c

    # If no head data at all, return uninformative priors for everything
    if not heads:
        return {
            (e.edge_id if hasattr(e, "edge_id") else f"{e.u}->{e.v}"): _UNINFORMATIVE
            for e in candidate_edges
        }

    # --- Group edges by source node to find the steepest direction per node ---
    # steepest_dir[u] = unit 3-vector toward the neighbor with max head-drop-per-km
    # max_grad[u]     = that maximum normalised gradient (head drop / distance)
    steepest_dir: Dict[str, tuple[float, float, float]] = {}
    max_grad: Dict[str, float] = {}

    edges_from: Dict[str, list[Edge]] = {}
    for edge in candidate_edges:
        edges_from.setdefault(str(edge.u), []).append(edge)

    for u, out_edges in edges_from.items():
        h_u = heads.get(u)
        c_u = coords.get(u)
        if h_u is None or c_u is None:
            continue

        best_grad = -math.inf
        best_uvec: Optional[tuple[float, float, float]] = None

        for edge in out_edges:
            v = str(edge.v)
            h_v = heads.get(v)
            c_v = coords.get(v)
            if h_v is None or c_v is None:
                continue
            drop = h_u - h_v
            if drop <= 0:
                continue
            dx = c_v[0] - c_u[0]
            dy = c_v[1] - c_u[1]
            dz = c_v[2] - c_u[2]
            dist = math.sqrt(dx * dx + dy * dy + dz * dz)
            if dist < 1e-9:
                continue
            grad = drop / dist
            if grad > best_grad:
                best_grad = grad
                best_uvec = (dx / dist, dy / dist, dz / dist)

        if best_uvec is not None:
            steepest_dir[u] = best_uvec
            max_grad[u] = best_grad

    # --- Score every candidate edge ---
    scores: Dict[str, float] = {}
    for edge in candidate_edges:
        eid = edge.edge_id if hasattr(edge, "edge_id") else f"{edge.u}->{edge.v}"
        u = str(edge.u)
        v = str(edge.v)
        h_u = heads.get(u)
        h_v = heads.get(v)

        if h_u is None or h_v is None or u not in coords or v not in coords:
            scores[eid] = _UNINFORMATIVE
            continue

        delta_h = h_u - h_v
        if delta_h <= 0.0:
            # Uphill or flat — strongly discourage
            scores[eid] = _PRIOR_FLOOR
            continue

        if u not in steepest_dir:
            # Node has exactly one downhill neighbour — it IS the steepest path
            scores[eid] = 0.80
            continue

        # Angular alignment with steepest direction at u
        c_u = coords[u]
        c_v = coords[v]
        dx = c_v[0] - c_u[0]
        dy = c_v[1] - c_u[1]
        dz = c_v[2] - c_u[2]
        dist = math.sqrt(dx * dx + dy * dy + dz * dz)
        if dist < 1e-9:
            scores[eid] = _UNINFORMATIVE
            continue

        uv_unit = (dx / dist, dy / dist, dz / dist)
        cos_theta = max(0.0, min(1.0, _dot3(uv_unit, steepest_dir[u])))

        # Relative head-gradient: how large is this edge's drop-per-unit-length
        # compared to the steepest direction at this node
        rel_grad = min(1.0, (delta_h / dist) / max_grad[u])

        # Combined score: angular alignment × relative gradient, both in [0,1]
        raw = (cos_theta ** angular_w) * (rel_grad ** head_w)
        scores[eid] = _PRIOR_FLOOR + (_PRIOR_CEIL - _PRIOR_FLOOR) * raw

    return scores


def apply_steepest_descent_priors(
    candidate_edges: Sequence[Edge],
    sample_map: Mapping[str, Mapping[str, object]],
    config: Any,
) -> None:
    """Mutate ``edge.attrs["edge_confidence"]`` on every candidate edge in-place.

    This is the injection point called from ``infer_topology_map_edges()`` when
    ``config.steepest_descent_enabled`` is True.  After this call the MCMC
    Bernoulli prior automatically uses the steepest-descent scores.

    For backward compatibility, this function always runs steepest-descent
    regardless of the config flag (the flag gates whether it is *called*,
    not what it does internally).
    """
    scores = compute_steepest_descent_priors(candidate_edges, sample_map, config)
    threshold = _get_config_float(config, "steepest_descent_threshold", 0.0)

    for edge in candidate_edges:
        if not hasattr(edge, "attrs"):
            continue
        eid = edge.edge_id if hasattr(edge, "edge_id") else f"{edge.u}->{edge.v}"
        score = scores.get(eid, _UNINFORMATIVE)
        attrs = dict(edge.attrs or {})
        attrs["edge_confidence"] = max(threshold, score) if threshold > 0 else score
        attrs["steepest_descent_score"] = score
        edge.attrs = attrs


# ---------------------------------------------------------------------------
# Projected-gradient priors (continuous head field)
# ---------------------------------------------------------------------------


def _sigmoid(x: float) -> float:
    """Numerically stable logistic sigmoid."""
    if x >= 0:
        return 1.0 / (1.0 + math.exp(-x))
    ex = math.exp(x)
    return ex / (1.0 + ex)


def compute_projected_gradient_priors(
    candidate_edges: Sequence[Edge],
    sample_map: Mapping[str, Mapping[str, object]],
    config: Any,
    gradient_map: Mapping[str, tuple[float, float, float]],
) -> Dict[str, float]:
    """Score each candidate edge by projecting the continuous head gradient
    at the source node onto the edge direction.

    For each edge u→v:

    1. Get the flow-direction vector at u: ``g = -grad_h(u)``.
       g points in the direction water actually flows (downhill).
    2. Compute the unit vector d_uv from u to v.
    3. Project ``proj = dot(g, d_uv)`` — this is the component of the
       gradient along the edge direction.
       - proj > 0: edge aligns with flow (high prior)
       - proj ≈ 0: edge is perpendicular to flow (low prior)
       - proj < 0: edge points against flow (floor prior)
    4. Map proj → [0.05, 0.95] via a sigmoid centred on the median
       projection across all edges.

    This replaces the discretized ``cos_θ^w1 × rel_grad^w2`` scoring.
    The projection naturally encodes both directional alignment AND
    gradient magnitude without double-counting.

    Parameters
    ----------
    candidate_edges:
        Universe of candidate directed Edge objects.
    sample_map:
        ``{node_id: sample_dict}`` providing coordinate data.
    config:
        Hydrosheaf Config (reads ``projected_gradient_sharpness``).
    gradient_map:
        ``{node_id: (gx, gy, gz)}`` from :func:`~hydrosheaf.physics.modflow_head.map_gradient_to_nodes`.

    Returns
    -------
    dict mapping ``edge_id → float`` score in [0.05, 0.95].
    """
    sharpness = _get_config_float(config, "projected_gradient_sharpness", 10.0)

    # Resolve coordinates
    coords: Dict[str, tuple[float, float, float]] = {}
    for node_id, sample in sample_map.items():
        c = _resolve_coords(sample)
        if c is not None:
            coords[node_id] = c

    # Compute raw projections for all edges that have gradient + coordinate data
    projections: list[float] = []
    edge_projections: Dict[str, float] = {}
    missing_nodes: set[str] = set()

    for edge in candidate_edges:
        eid = edge.edge_id if hasattr(edge, "edge_id") else f"{edge.u}->{edge.v}"
        u = str(edge.u)
        v = str(edge.v)

        grad = gradient_map.get(u)
        c_u = coords.get(u)
        c_v = coords.get(v)

        if grad is None:
            missing_nodes.add(u)
            continue
        if c_u is None or c_v is None:
            continue

        # Flow direction vector (downhill, water flows this way)
        gx, gy, gz = grad
        g = (-gx, -gy, -gz)

        # Edge unit vector
        dx = c_v[0] - c_u[0]
        dy = c_v[1] - c_u[1]
        dz = c_v[2] - c_u[2]
        dist = math.sqrt(dx * dx + dy * dy + dz * dz)
        if dist < 1e-9:
            continue

        d_uv = (dx / dist, dy / dist, dz / dist)
        proj = _dot3(g, d_uv)
        edge_projections[eid] = proj
        projections.append(proj)

    # Adaptive centering and scaling.
    # When the edge set has sufficient diversity (IQR > threshold), use
    # median-centring and IQR-scaling so the sigmoid distinguishes
    # better/worse-aligned edges within the candidate set.
    # When diversity is low (few edges or near-identical projections),
    # the data do not support directional discrimination — return
    # uninformative priors for all edges.
    if not projections:
        # No usable gradient data for any edge → uninformative for all
        return {
            (e.edge_id if hasattr(e, "edge_id") else f"{e.u}->{e.v}"): _UNINFORMATIVE
            for e in candidate_edges
        }

    sorted_proj = sorted(projections)
    n = len(sorted_proj)
    median_proj = sorted_proj[n // 2]
    q1 = sorted_proj[n // 4] if n >= 4 else sorted_proj[0]
    q3 = sorted_proj[(3 * n) // 4] if n >= 4 else sorted_proj[-1]
    iqr = q3 - q1

    if iqr > 1e-8 and n >= 4:
        # High diversity: rank edges relative to the set
        scale = iqr
        center = median_proj
    else:
        # Low diversity: all edges get uninformative prior
        # (flat terrain — data don't support directional discrimination)
        return {
            (e.edge_id if hasattr(e, "edge_id") else f"{e.u}->{e.v}"): _UNINFORMATIVE
            for e in candidate_edges
        }

    # Map projections to scores
    scores: Dict[str, float] = {}
    for edge in candidate_edges:
        eid = edge.edge_id if hasattr(edge, "edge_id") else f"{edge.u}->{edge.v}"

        if eid in edge_projections:
            proj = edge_projections[eid]
            if proj <= 0.0:
                # Anti-aligned with flow — floor prior
                scores[eid] = _PRIOR_FLOOR
            else:
                # Sigmoid mapping: dimensionless via (proj - center) / scale
                z = sharpness * (proj - center) / scale
                raw = _sigmoid(z)
                scores[eid] = _PRIOR_FLOOR + (_PRIOR_CEIL - _PRIOR_FLOOR) * raw
        elif str(edge.u) in missing_nodes or str(edge.v) in missing_nodes:
            scores[eid] = _UNINFORMATIVE
        else:
            # No gradient data but coordinates available → uninformative
            scores[eid] = _UNINFORMATIVE

    return scores


def apply_flow_direction_priors(
    candidate_edges: Sequence[Edge],
    sample_map: Mapping[str, Mapping[str, object]],
    config: Any,
    gradient_map: Optional[Mapping[str, tuple[float, float, float]]] = None,
) -> None:
    """Unified dispatch for flow-direction edge priors.

    Mutates ``edge.attrs["edge_confidence"]`` on every candidate edge in-place.

    Dispatch logic:
    - If ``config.projected_gradient_enabled`` is True and *gradient_map* is
      provided → :func:`compute_projected_gradient_priors` (continuous field).
    - Elif ``config.steepest_descent_enabled`` is True →
      :func:`compute_steepest_descent_priors` (discretised, backward-compat).
    - Otherwise → no-op (leaves existing ``edge_confidence`` unchanged).

    This is the canonical injection point for topology-posterior callers.
    """
    if getattr(config, "projected_gradient_enabled", False) and gradient_map is not None:
        scores = compute_projected_gradient_priors(
            candidate_edges, sample_map, config, gradient_map
        )
        attr_name = "projected_gradient_score"
    elif getattr(config, "steepest_descent_enabled", False):
        scores = compute_steepest_descent_priors(candidate_edges, sample_map, config)
        attr_name = "steepest_descent_score"
    else:
        return

    threshold = _get_config_float(
        config,
        "steepest_descent_threshold"
        if attr_name == "steepest_descent_score"
        else "projected_gradient_threshold",
        0.0,
    )

    for edge in candidate_edges:
        if not hasattr(edge, "attrs"):
            continue
        eid = edge.edge_id if hasattr(edge, "edge_id") else f"{edge.u}->{edge.v}"
        score = scores.get(eid, _UNINFORMATIVE)
        attrs = dict(edge.attrs or {})
        attrs["edge_confidence"] = max(threshold, score) if threshold > 0 else score
        attrs[attr_name] = score
        edge.attrs = attrs
