from __future__ import annotations

import csv
import json
from dataclasses import asdict
from datetime import datetime
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

from .contracts import VadoseForcingSample, VadoseLayer, VadoseLinksRow, VadoseProfile

try:
    import yaml  # type: ignore
except ModuleNotFoundError:  # pragma: no cover
    yaml = None


def load_forcing_csv(
    path: str,
    *,
    time_column: str = "timestamp",
    time_format: str = "%Y-%m-%d",
    precipitation_column: str = "P_mm_day",
    et0_column: str = "ET0_mm_day",
    irrigation_column: str = "I_mm_day",
) -> List[VadoseForcingSample]:
    forcing: List[VadoseForcingSample] = []
    with open(path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            ts_raw = (row.get(time_column) or "").strip()
            if not ts_raw:
                continue
            try:
                ts = datetime.strptime(ts_raw, time_format)
            except ValueError:
                # Try common alternatives.
                for fmt in ["%Y-%m-%d %H:%M:%S", "%Y/%m/%d", "%m/%d/%Y"]:
                    try:
                        ts = datetime.strptime(ts_raw, fmt)
                        break
                    except ValueError:
                        continue
                else:
                    continue

            def _f(key: str) -> float:
                raw = row.get(key)
                if raw in (None, ""):
                    return 0.0
                try:
                    return float(raw)
                except ValueError:
                    return 0.0

            forcing.append(
                VadoseForcingSample(
                    timestamp=ts,
                    precipitation_mm_day=_f(precipitation_column),
                    et0_mm_day=_f(et0_column),
                    irrigation_mm_day=_f(irrigation_column),
                )
            )
    forcing.sort(key=lambda s: s.timestamp)
    return forcing


def load_links_csv(path: str) -> List[VadoseLinksRow]:
    links: List[VadoseLinksRow] = []
    with open(path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            u = str(row.get("u") or "").strip()
            v = str(row.get("v") or "").strip()
            if not u or not v:
                continue
            depth = row.get("u_depth_m") or row.get("depth_m")
            if depth in (None, ""):
                continue
            try:
                u_depth_m = float(depth)
            except ValueError:
                continue
            p_uv = 1.0
            if row.get("p_uv") not in (None, ""):
                try:
                    p_uv = float(row["p_uv"])
                except ValueError:
                    p_uv = 1.0
            links.append(VadoseLinksRow(u=u, v=v, u_depth_m=u_depth_m, p_uv=p_uv))
    return links


def load_water_table_csv(
    path: str,
    *,
    time_column: str = "timestamp",
    time_format: str = "%Y-%m-%d",
    wt_depth_column: str = "water_table_depth_m",
) -> List[Tuple[datetime, float]]:
    series: List[Tuple[datetime, float]] = []
    with open(path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            ts_raw = (row.get(time_column) or "").strip()
            if not ts_raw:
                continue
            try:
                ts = datetime.strptime(ts_raw, time_format)
            except ValueError:
                continue
            raw = row.get(wt_depth_column)
            if raw in (None, ""):
                continue
            try:
                depth = float(raw)
            except ValueError:
                continue
            series.append((ts, depth))
    series.sort(key=lambda p: p[0])
    return series


def load_profile(path: str) -> VadoseProfile:
    """Load vadose profile from JSON or YAML (if PyYAML is available).

    JSON/YAML schema:
      {
        "profile_id": "siteA",
        "depth_m": 5.0,
        "root_depth_m": 1.0,  # optional
        "layers": [
          {"thickness_m": 0.3, "texture": "sandy_loam"},
          {"thickness_m": 0.7, "theta_r":..., "theta_s":..., "alpha_1_m":..., "n":..., "ks_m_day":...}
        ]
      }
    """
    ext = Path(path).suffix.lower()
    if ext in {".yaml", ".yml"}:
        if yaml is None:
            raise RuntimeError("PyYAML not installed; provide profile as JSON or install pyyaml.")
        with open(path, "r", encoding="utf-8") as handle:
            data = yaml.safe_load(handle)
    else:
        with open(path, "r", encoding="utf-8") as handle:
            data = json.load(handle)

    if not isinstance(data, dict):
        raise ValueError("profile must be an object/dict")
    profile_id = str(data.get("profile_id") or "vadose_profile")
    depth_m = float(data.get("depth_m"))
    root_depth_m = data.get("root_depth_m")
    layers_raw = data.get("layers")
    if not isinstance(layers_raw, list) or not layers_raw:
        raise ValueError("profile.layers must be a non-empty list")
    layers: List[VadoseLayer] = []
    for entry in layers_raw:
        if not isinstance(entry, dict):
            continue
        layers.append(
            VadoseLayer(
                thickness_m=float(entry.get("thickness_m")),
                texture=entry.get("texture"),
                theta_r=entry.get("theta_r"),
                theta_s=entry.get("theta_s"),
                alpha_1_m=entry.get("alpha_1_m"),
                n=entry.get("n"),
                ks_m_day=entry.get("ks_m_day"),
                l=float(entry.get("l", 0.5)),
            )
        )
    return VadoseProfile(
        profile_id=profile_id,
        depth_m=float(depth_m),
        layers=layers,
        root_depth_m=float(root_depth_m) if root_depth_m is not None else None,
    )


