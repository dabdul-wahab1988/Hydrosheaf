"""FloPy-based readers for MODPATH particle tracking outputs (optional)."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import math

from ..graph.types import Edge
from .priors import PhysicsPrior


@dataclass(frozen=True)
class NodeCoord:
    node_id: str
    x: float
    y: float
    z: Optional[float] = None


def _safe_float(value: object) -> Optional[float]:
    try:
        if value in (None, ""):
            return None
        return float(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return None


def _euclid(a: NodeCoord, x: float, y: float, z: Optional[float] = None) -> float:
    dx = a.x - x
    dy = a.y - y
    dz = 0.0
    if z is not None and a.z is not None:
        dz = a.z - z
    return math.sqrt(dx * dx + dy * dy + dz * dz)


def _nearest_node(
    nodes: Sequence[NodeCoord], x: float, y: float, z: Optional[float], max_dist: float
) -> Optional[str]:
    best_id = None
    best_d = float("inf")
    for node in nodes:
        d = _euclid(node, x, y, z)
        if d < best_d:
            best_d = d
            best_id = node.node_id
    if best_id is None or best_d > max_dist:
        return None
    return best_id


def node_coords_from_samples(
    samples: Iterable[Mapping[str, object]],
    *,
    node_id_key: str = "site_id",
    x_key: str = "x",
    y_key: str = "y",
    z_key: Optional[str] = None,
) -> List[NodeCoord]:
    coords: List[NodeCoord] = []
    for row in samples:
        node_id = row.get(node_id_key)
        if node_id in (None, ""):
            continue
        x = _safe_float(row.get(x_key))
        y = _safe_float(row.get(y_key))
        if x is None or y is None:
            continue
        z = _safe_float(row.get(z_key)) if z_key else None
        coords.append(NodeCoord(node_id=str(node_id), x=float(x), y=float(y), z=z))
    return coords


def priors_from_modpath_endpoints(
    endpoints_path: str,
    nodes: Sequence[NodeCoord],
    *,
    max_match_distance: float = 500.0,
    time_unit_days: float = 1.0,
    source: str = "modpath_endpoints",
) -> List[PhysicsPrior]:
    """Build PhysicsPrior edges from a MODPATH endpoints file.

    Assumes endpoint coordinates are in the same coordinate system as `nodes`.
    """
    try:
        import flopy  # type: ignore
    except Exception as exc:  # pragma: no cover
        raise RuntimeError(
            "FloPy is required to read MODPATH endpoints. Install flopy to use this feature."
        ) from exc

    ep = flopy.utils.EndpointFile(endpoints_path)
    data = ep.get_alldata()

    # Aggregations: (u,v) -> [times]
    times_by_uv: Dict[Tuple[str, str], List[float]] = {}
    counts_from_u: Dict[str, int] = {}

    # Common MODPATH fields: x0,y0,z0 (start), x,y,z (end), time
    for rec in data:
        x0 = float(rec["x0"]) if "x0" in rec.dtype.names else None
        y0 = float(rec["y0"]) if "y0" in rec.dtype.names else None
        z0 = float(rec["z0"]) if "z0" in rec.dtype.names else None
        x1 = float(rec["x"]) if "x" in rec.dtype.names else None
        y1 = float(rec["y"]) if "y" in rec.dtype.names else None
        z1 = float(rec["z"]) if "z" in rec.dtype.names else None
        t = float(rec["time"]) if "time" in rec.dtype.names else None
        if x0 is None or y0 is None or x1 is None or y1 is None or t is None:
            continue
        u = _nearest_node(nodes, x0, y0, z0, max_match_distance)
        v = _nearest_node(nodes, x1, y1, z1, max_match_distance)
        if u is None or v is None or u == v:
            continue
        tau_days = float(t) * float(time_unit_days)
        times_by_uv.setdefault((u, v), []).append(tau_days)
        counts_from_u[u] = counts_from_u.get(u, 0) + 1

    priors: List[PhysicsPrior] = []
    for (u, v), times in times_by_uv.items():
        times_sorted = sorted(times)
        n = len(times_sorted)
        if n == 0:
            continue
        mean = sum(times_sorted) / n
        p10 = times_sorted[int(0.1 * (n - 1))]
        p90 = times_sorted[int(0.9 * (n - 1))]
        # crude std
        var = sum((tt - mean) ** 2 for tt in times_sorted) / max(1, n - 1)
        std = math.sqrt(var)
        p_uv = float(n) / float(max(1, counts_from_u.get(u, n)))
        priors.append(
            PhysicsPrior(
                u=u,
                v=v,
                p_uv=p_uv,
                tt_mean_days=float(mean),
                tt_p10_days=float(p10),
                tt_p90_days=float(p90),
                tt_std_days=float(std),
                source=source,
            )
        )
    return priors


def edges_from_priors(priors: Iterable[PhysicsPrior]) -> List[Edge]:
    return [Edge(edge_id=p.edge_id(), u=p.u, v=p.v, attrs=p.attrs()) for p in priors]
