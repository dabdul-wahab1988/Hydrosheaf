"""File-based physics priors for edges (p_uv, travel time, etc.)."""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Tuple

from ..graph.types import Edge


@dataclass(frozen=True)
class PhysicsPrior:
    u: str
    v: str
    p_uv: Optional[float] = None
    tt_mean_days: Optional[float] = None
    tt_p10_days: Optional[float] = None
    tt_p90_days: Optional[float] = None
    tt_std_days: Optional[float] = None
    source: str = "physics"

    def edge_id(self) -> str:
        return f"{self.u}->{self.v}"

    def attrs(self) -> Dict[str, object]:
        attrs: Dict[str, object] = {"physics_source": self.source}
        if self.p_uv is not None:
            attrs["p_uv"] = float(self.p_uv)
            attrs["edge_confidence"] = float(self.p_uv)
        if self.tt_mean_days is not None:
            attrs["physics_tau_mean_days"] = float(self.tt_mean_days)
            # Also provide generic keys some pipelines may check.
            attrs["edge_residence_time_days"] = float(self.tt_mean_days)
        if self.tt_std_days is not None:
            attrs["physics_tau_std_days"] = float(self.tt_std_days)
        if self.tt_p10_days is not None:
            attrs["physics_tau_p10_days"] = float(self.tt_p10_days)
        if self.tt_p90_days is not None:
            attrs["physics_tau_p90_days"] = float(self.tt_p90_days)
        return attrs


def _safe_float(value: object) -> Optional[float]:
    try:
        if value in (None, ""):
            return None
        return float(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return None


def load_physics_priors(path: str) -> List[PhysicsPrior]:
    """Load priors from CSV or JSON.

    CSV columns (recommended):
      - u, v (required)
      - p_uv (optional)
      - tt_mean_days (optional)
      - tt_std_days (optional) or tt_p10_days/tt_p90_days (optional)

    JSON formats supported:
      - list of objects with the same keys as the CSV
    """
    ext = Path(path).suffix.lower()
    if ext in {".json"}:
        with open(path, "r", encoding="utf-8") as handle:
            data = json.load(handle)
        if not isinstance(data, list):
            raise ValueError("physics priors JSON must be a list of objects.")
        priors: List[PhysicsPrior] = []
        for row in data:
            if not isinstance(row, Mapping):
                continue
            u = str(row.get("u") or "")
            v = str(row.get("v") or "")
            if not u or not v:
                continue
            priors.append(
                PhysicsPrior(
                    u=u,
                    v=v,
                    p_uv=_safe_float(row.get("p_uv")),
                    tt_mean_days=_safe_float(
                        row.get("tt_mean_days") or row.get("tau_mean_days")
                    ),
                    tt_p10_days=_safe_float(row.get("tt_p10_days")),
                    tt_p90_days=_safe_float(row.get("tt_p90_days")),
                    tt_std_days=_safe_float(
                        row.get("tt_std_days") or row.get("tau_std_days")
                    ),
                    source=str(row.get("source") or "physics"),
                )
            )
        return priors

    # CSV
    with open(path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        priors = []
        for row in reader:
            u = str(row.get("u") or "")
            v = str(row.get("v") or "")
            if not u or not v:
                continue
            priors.append(
                PhysicsPrior(
                    u=u,
                    v=v,
                    p_uv=_safe_float(row.get("p_uv") or row.get("edge_confidence")),
                    tt_mean_days=_safe_float(
                        row.get("tt_mean_days")
                        or row.get("tau_mean_days")
                        or row.get("edge_residence_time_days")
                    ),
                    tt_p10_days=_safe_float(row.get("tt_p10_days")),
                    tt_p90_days=_safe_float(row.get("tt_p90_days")),
                    tt_std_days=_safe_float(
                        row.get("tt_std_days") or row.get("tau_std_days")
                    ),
                    source=str(row.get("source") or "physics"),
                )
            )
        return priors


def apply_physics_priors(
    edges: List[Edge],
    priors: Iterable[PhysicsPrior],
    *,
    mode: str = "override",
) -> List[Edge]:
    """Apply physics priors to edges.

    Modes
    -----
    - override: keep existing edge set, overwrite/attach attrs where priors exist
    - merge: union of existing edges and physics edges
    - only: ignore existing edges; use physics edges only
    """
    if mode not in {"override", "merge", "only"}:
        raise ValueError("mode must be one of: override, merge, only")

    prior_by_uv: Dict[Tuple[str, str], PhysicsPrior] = {(p.u, p.v): p for p in priors}
    if mode == "only":
        merged: List[Edge] = []
        for p in prior_by_uv.values():
            merged.append(Edge(edge_id=p.edge_id(), u=p.u, v=p.v, attrs=p.attrs()))
        return merged

    merged_by_id: Dict[str, Edge] = {e.edge_id: e for e in edges}
    for (u, v), p in prior_by_uv.items():
        edge_id = f"{u}->{v}"
        if edge_id not in merged_by_id:
            if mode == "merge":
                merged_by_id[edge_id] = Edge(edge_id=edge_id, u=u, v=v, attrs=p.attrs())
            continue
        edge = merged_by_id[edge_id]
        attrs = dict(edge.attrs or {})
        attrs.update(p.attrs())
        edge.attrs = attrs

    return list(merged_by_id.values())
