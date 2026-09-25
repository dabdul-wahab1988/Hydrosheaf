"""USGS-facing, truth-blind inputs for the redesigned M4 benchmark.

This module is deliberately narrower than the general HydroSheaf inference
pipeline.  It converts the public Savage archive inputs into the observation
schema consumed by :mod:`hydrosheaf.validation.topology_v2` and keeps the
MODPATH reference graph out of that conversion path.

The archive stores the Savage endpoint coordinates in a projected model frame
whose native unit is feet.  Pairwise distances are therefore Euclidean in that
frame (after a feet-to-metres conversion), never haversine distances.  The
origin/rotation metadata is retained in the frame object for provenance; a
rigid translation/rotation cannot change pairwise distances, but making the
transform explicit prevents projected coordinates from being mistaken for
latitude/longitude.
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
import math
from pathlib import Path
from typing import Any


# Exact US survey foot conversion used by State Plane New Hampshire metadata.
FEET_TO_METRES = 1200.0 / 3937.0
M4_A = "M4-A_sparse_observation"
M4_B = "M4-B_archive_informed"


@dataclass(frozen=True)
class SavageProjectedFrame:
    """Declared coordinate frame for the public Savage model archive."""

    coordinate_system: str = "NAD83 / State Plane New Hampshire (US Feet)"
    native_units: str = "US survey/model feet"
    output_units: str = "metres"
    rotation_deg: float = -12.0
    origin_x_ft: float = 961030.4
    origin_y_ft: float = 112955.0
    transform_source: str = (
        "M4 phase2_savage_pipeline/phase2b grid metadata; ANGROT and NH State Plane origin"
    )

    def project(self, x_native: float, y_native: float) -> tuple[float, float]:
        """Transform model-local projected coordinates to projected metres."""

        theta = math.radians(float(self.rotation_deg))
        cosine = math.cos(theta)
        sine = math.sin(theta)
        world_x_ft = self.origin_x_ft + cosine * float(x_native) - sine * float(y_native)
        world_y_ft = self.origin_y_ft + sine * float(x_native) + cosine * float(y_native)
        return world_x_ft * FEET_TO_METRES, world_y_ft * FEET_TO_METRES

    def to_dict(self) -> dict[str, object]:
        return {
            "coordinate_system": self.coordinate_system,
            "native_units": self.native_units,
            "output_units": self.output_units,
            "rotation_deg": self.rotation_deg,
            "origin_x_ft": self.origin_x_ft,
            "origin_y_ft": self.origin_y_ft,
            "transform_source": self.transform_source,
            "distance_method": "Euclidean in projected coordinates",
            "latitude_longitude_interpretation": False,
        }


def _finite(value: object) -> float | None:
    if value is None:
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _records(nodes: Any) -> list[dict[str, object]]:
    if hasattr(nodes, "to_dict"):
        rows = nodes.to_dict("records")
    else:
        rows = [dict(row) for row in nodes]
    result: list[dict[str, object]] = []
    for row in rows:
        node_id = str(row.get("node_id", row.get("site_id", ""))).strip()
        if not node_id:
            raise ValueError("Savage node records require node_id.")
        x = _finite(row.get("x"))
        y = _finite(row.get("y"))
        if x is None or y is None:
            raise ValueError(f"Savage node {node_id} is missing projected x/y coordinates.")
        copied = dict(row)
        copied["node_id"] = node_id
        copied["x"] = x
        copied["y"] = y
        result.append(copied)
    result.sort(key=lambda row: str(row["node_id"]))
    if len({str(row["node_id"]) for row in result}) != len(result):
        raise ValueError("Savage node records contain duplicate node IDs.")
    return result


def build_savage_observations(
    nodes: Any,
    *,
    mode: str,
    heads: Mapping[str, object] | None = None,
    budget_context: Mapping[str, Mapping[str, object]] | None = None,
    frame: SavageProjectedFrame | None = None,
    head_scale_to_m: float = FEET_TO_METRES,
    head_sigma_m: float = 0.05,
) -> tuple[dict[str, object], ...]:
    """Build truth-blind M4 observations from nodes and public model fields.

    ``heads`` and ``budget_context`` are used only for M4-B.  M4-B fails
    closed if either channel is absent or incomplete, because silently using
    M4-A evidence would make the archive-informed label misleading.
    """

    resolved_frame = frame or SavageProjectedFrame()
    rows = _records(nodes)
    resolved_heads = {str(key): _finite(value) for key, value in (heads or {}).items()}
    if mode not in {M4_A, M4_B}:
        raise ValueError(f"Unknown M4 mode: {mode}")
    if mode == M4_B:
        missing_heads = [
            str(row["node_id"])
            for row in rows
            if resolved_heads.get(str(row["node_id"])) is None
        ]
        if missing_heads:
            raise ValueError(
                "M4-B requires an actual FHD head for every benchmark node; "
                f"missing {len(missing_heads)} nodes."
            )
        if budget_context is None:
            raise ValueError("M4-B requires aggregated public CBC source/sink context.")

    observations: list[dict[str, object]] = []
    for row in rows:
        node_id = str(row["node_id"])
        x_m, y_m = resolved_frame.project(float(row["x"]), float(row["y"]))
        observation: dict[str, object] = {
            "site_id": node_id,
            "x_m": x_m,
            "y_m": y_m,
            "coordinate_frame": resolved_frame.coordinate_system,
        }
        elevation = _finite(row.get("z"))
        if elevation is not None:
            # Endpoint z is retained as a sparse vertical/elevation proxy.  It
            # is intentionally not represented as hydraulic head.
            observation["elevation"] = elevation
            observation["elevation_source"] = "MODPATH endpoint z proxy"

        if mode == M4_B:
            head_ft = resolved_heads[node_id]
            assert head_ft is not None
            observation["hydraulic_head"] = float(head_ft) * float(head_scale_to_m)
            observation["head_sigma_m"] = float(head_sigma_m)
            context = dict(budget_context.get(node_id, {}))
            # The CBC aggregation covers all nodes.  A zero means that the
            # node had no record in the public budget channel, not that the
            # channel was unavailable.
            for key in (
                "well_rate",
                "river_leakage",
                "recharge",
                "head_boundary_flux",
            ):
                value = _finite(context.get(key))
                observation[key] = 0.0 if value is None else value
            observation["source_sink_context_source"] = "MODFLOW CBC public budget records"
        observations.append(observation)
    return tuple(observations)


def _sum_node_records(
    budget_file: Any,
    record_name: str,
    *,
    counts: dict[str, int],
) -> dict[int, float]:
    """Sum a structured CBC node/q record by one-based MODFLOW node."""

    values: defaultdict[int, float] = defaultdict(float)
    try:
        blocks = budget_file.get_data(text=record_name, full3D=False)
    except Exception:
        counts[record_name] = 0
        return {}
    count = 0
    for block in blocks or ():
        names = getattr(getattr(block, "dtype", None), "names", None) or ()
        if "node" not in names or "q" not in names:
            continue
        for record in block:
            node = int(record["node"])
            values[node] += float(record["q"])
            count += 1
    counts[record_name] = count
    return dict(values)


def _node_index(node_id: str, *, nrow: int, ncol: int) -> tuple[int, int, int] | None:
    text = str(node_id)
    if "_" not in text:
        return None
    try:
        flat = int(text.rsplit("_", 1)[1]) - 1
    except ValueError:
        return None
    if flat < 0:
        return None
    layer_size = int(nrow) * int(ncol)
    layer, remainder = divmod(flat, layer_size)
    row, col = divmod(remainder, int(ncol))
    return layer, row, col


def aggregate_savage_cbc_context(
    cbc_path: str | Path,
    node_ids: Iterable[str],
    *,
    nrow: int = 202,
    ncol: int = 183,
) -> tuple[dict[str, dict[str, float]], dict[str, object]]:
    """Aggregate public MODFLOW CBC source/sink records by node.

    The MODFLOW record signs are retained.  ``well_rate`` is negative for
    extraction at the Savage receptor wells; the topology-v2 feature builder
    converts signed rates into graded source/sink strengths.  No MODPATH
    endpoint or edge data are read here.
    """

    from flopy.utils import CellBudgetFile

    path = Path(cbc_path)
    if not path.exists():
        raise FileNotFoundError(path)
    budget_file = CellBudgetFile(str(path))
    counts: dict[str, int] = {}
    wells = _sum_node_records(budget_file, "WELLS", counts=counts)
    rivers = _sum_node_records(budget_file, "RIVER LEAKAGE", counts=counts)
    boundaries = _sum_node_records(budget_file, "HEAD DEP BOUNDS", counts=counts)

    recharge_by_node: dict[int, float] = {}
    try:
        recharge_blocks = budget_file.get_data(text="RECHARGE", full3D=False)
        arrays = [
            block
            for block in (recharge_blocks or ())
            if getattr(block, "ndim", 0) == 2
        ]
        if arrays:
            recharge_array = arrays[-1]
            for node_id in node_ids:
                index = _node_index(str(node_id), nrow=nrow, ncol=ncol)
                if index is None:
                    continue
                layer, row, col = index
                # MODFLOW recharge enters the upper active layer.  The CBC
                # array is represented on the row/column grid, so only layer 0
                # receives the direct recharge feature.
                if layer == 0 and row < recharge_array.shape[0] and col < recharge_array.shape[1]:
                    recharge_by_node[int(str(node_id).rsplit("_", 1)[1])] = float(
                        recharge_array[row, col]
                    )
    except Exception:
        recharge_by_node = {}
    counts["RECHARGE"] = len(recharge_by_node)

    context: dict[str, dict[str, float]] = {}
    for node_id in sorted({str(node) for node in node_ids}):
        try:
            node = int(node_id.rsplit("_", 1)[1])
        except (IndexError, ValueError):
            node = -1
        context[node_id] = {
            "well_rate": float(wells.get(node, 0.0)),
            "river_leakage": float(rivers.get(node, 0.0)),
            "recharge": float(recharge_by_node.get(node, 0.0)),
            "head_boundary_flux": float(boundaries.get(node, 0.0)),
        }
    metadata = {
        "cbc_path": str(path),
        "record_names_used": ["WELLS", "RIVER LEAKAGE", "HEAD DEP BOUNDS", "RECHARGE"],
        "record_counts": dict(sorted(counts.items())),
        "node_count": len(context),
        "truth_blind": True,
    }
    return context, metadata


__all__ = [
    "FEET_TO_METRES",
    "M4_A",
    "M4_B",
    "SavageProjectedFrame",
    "aggregate_savage_cbc_context",
    "build_savage_observations",
]
