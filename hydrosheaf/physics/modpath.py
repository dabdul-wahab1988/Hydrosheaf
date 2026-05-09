"""Readers for MODPATH particle tracking outputs."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import flopy

from ..graph.types import Edge

from .priors import PhysicsPrior


@dataclass(frozen=True)
class NodeCoord:
    node_id: str
    x: float
    y: float
    z: Optional[float] = None


@dataclass(frozen=True)
class ModpathEndpoint:
    particle_id: int
    x0: float
    y0: float
    z0: Optional[float]
    x: float
    y: float
    z: Optional[float]
    time: float
    initial_cell: Optional[int] = None
    final_cell: Optional[int] = None
    status: Optional[int] = None


@dataclass(frozen=True)
class ModpathPathlinePoint:
    particle_id: int
    x: float
    y: float
    z: Optional[float]
    time: float
    cell: Optional[int] = None
    sequence: int = 0


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


def _field(rec: object, names: Sequence[str]) -> Optional[float]:
    dtype_names = getattr(getattr(rec, "dtype", None), "names", None) or ()
    for name in names:
        if name in dtype_names:
            return float(rec[name])
    return None


def _load_flopy_endpoint_records(endpoints_path: str) -> List[ModpathEndpoint]:
    ep = flopy.utils.EndpointFile(endpoints_path)
    data = ep.get_alldata()

    records: List[ModpathEndpoint] = []
    for idx, rec in enumerate(data, start=1):
        x0 = _field(rec, ("x0",))
        y0 = _field(rec, ("y0",))
        z0 = _field(rec, ("z0", "zloc0"))
        x = _field(rec, ("x",))
        y = _field(rec, ("y",))
        z = _field(rec, ("z", "zloc"))
        time = _field(rec, ("time",))
        if x0 is None or y0 is None or x is None or y is None or time is None:
            continue
        particle_id = int(_field(rec, ("particleid",)) or idx)
        initial_cell = _field(rec, ("node0",))
        final_cell = _field(rec, ("node",))
        status = _field(rec, ("status", "ipcode"))
        records.append(
            ModpathEndpoint(
                particle_id=particle_id,
                x0=float(x0),
                y0=float(y0),
                z0=z0,
                x=float(x),
                y=float(y),
                z=z,
                time=float(time),
                initial_cell=int(initial_cell) if initial_cell is not None else None,
                final_cell=int(final_cell) if final_cell is not None else None,
                status=int(status) if status is not None else None,
            )
        )
    return records


def _load_compact_modpath5_endpoint_records(endpoints_path: str) -> List[ModpathEndpoint]:
    """Load the 14-column compact MODPATH 5 endpoint layout.

    MODPATH 5 compact endpoint files used by some USGS archives omit row/column
    fields and store start/final model-cell numbers instead. The columns needed
    for Hydrosheaf priors are stable in those files:
    final coordinates/time followed by starting coordinates/cell metadata.
    """
    records: List[ModpathEndpoint] = []
    with open(endpoints_path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            text = line.strip()
            if not text or text.startswith("@") or text.startswith("#"):
                continue
            parts = text.split()
            if len(parts) != 14:
                continue
            try:
                values = [float(part) for part in parts]
            except ValueError:
                continue
            records.append(
                ModpathEndpoint(
                    particle_id=len(records) + 1,
                    final_cell=int(values[1]),
                    x=float(values[2]),
                    y=float(values[3]),
                    z=float(values[4]),
                    time=float(values[5]),
                    x0=float(values[6]),
                    y0=float(values[7]),
                    z0=float(values[8]),
                    initial_cell=int(values[9]),
                    status=int(values[12]),
                )
            )
    if not records:
        raise ValueError(
            "Could not parse a supported MODPATH endpoint layout. Expected a "
            "FloPy-readable endpoint file or 14-column compact MODPATH 5 file."
        )
    return records


def load_modpath_endpoint_records(endpoints_path: str) -> List[ModpathEndpoint]:
    """Load MODPATH endpoint records with a compact MODPATH 5 fallback."""
    try:
        return _load_flopy_endpoint_records(endpoints_path)
    except Exception:
        return _load_compact_modpath5_endpoint_records(endpoints_path)


def _load_flopy_pathline_points(pathline_path: str) -> List[ModpathPathlinePoint]:
    pathline = flopy.utils.PathlineFile(pathline_path)
    data = pathline.get_alldata()
    points: List[ModpathPathlinePoint] = []
    for rec in data:
        rows = rec if hasattr(rec, "__iter__") and not isinstance(rec, tuple) else [rec]
        for row in rows:
            particle_id = _field(row, ("particleid", "particleidloc", "particlegroup"))
            x = _field(row, ("x",))
            y = _field(row, ("y",))
            z = _field(row, ("z", "zloc"))
            time = _field(row, ("time",))
            cell = _field(row, ("node", "cell"))
            if particle_id is None or x is None or y is None or time is None:
                continue
            points.append(
                ModpathPathlinePoint(
                    particle_id=int(particle_id),
                    x=float(x),
                    y=float(y),
                    z=z,
                    time=float(time),
                    cell=int(cell) if cell is not None else None,
                    sequence=len(points),
                )
            )
    if not points:
        raise ValueError("No pathline points were parsed by FloPy.")
    return points


def _load_compact_modpath5_pathline_points(pathline_path: str) -> List[ModpathPathlinePoint]:
    """Load the 8-column compact MODPATH 5 pathline layout."""
    points: List[ModpathPathlinePoint] = []
    sequence_by_particle: Dict[int, int] = {}
    with open(pathline_path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            text = line.strip()
            if not text or text.startswith("@") or text.startswith("#"):
                continue
            parts = text.split()
            if len(parts) != 8:
                continue
            try:
                values = [float(part) for part in parts]
            except ValueError:
                continue
            particle_id = int(values[0])
            sequence = sequence_by_particle.get(particle_id, 0)
            sequence_by_particle[particle_id] = sequence + 1
            points.append(
                ModpathPathlinePoint(
                    particle_id=particle_id,
                    x=float(values[1]),
                    y=float(values[2]),
                    z=float(values[3]),
                    time=float(values[4]),
                    cell=int(values[6]),
                    sequence=sequence,
                )
            )
    if not points:
        raise ValueError(
            "Could not parse a supported MODPATH pathline layout. Expected a "
            "FloPy-readable pathline file or 8-column compact MODPATH 5 file."
        )
    return points


def load_modpath_pathline_points(pathline_path: str) -> List[ModpathPathlinePoint]:
    """Load MODPATH pathline points with a compact MODPATH 5 fallback."""
    try:
        return _load_flopy_pathline_points(pathline_path)
    except Exception:
        return _load_compact_modpath5_pathline_points(pathline_path)


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
    data = load_modpath_endpoint_records(endpoints_path)

    # Aggregations: (u,v) -> [times]
    times_by_uv: Dict[Tuple[str, str], List[float]] = {}
    counts_from_u: Dict[str, int] = {}

    # Common MODPATH fields: x0,y0,z0 (start), x,y,z (end), time
    for rec in data:
        x0 = rec.x0
        y0 = rec.y0
        z0 = rec.z0
        x1 = rec.x
        y1 = rec.y
        z1 = rec.z
        t = rec.time
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
