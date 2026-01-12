"""Build reaction bounds from PHREEQC outputs."""

from typing import Dict, Iterable, List, Mapping, Tuple

from ..config import Config
from ..graph.build import build_edges

SI_KEYS = {
    "calcite": "si_calcite",
    "dolomite": "si_dolomite",
    "gypsum": "si_gypsum",
    "halite": "si_halite",
    "fluorite": "si_fluorite",
    "siderite": "si_siderite",
    "apatite": "si_apatite",
    "goethite": "si_goethite",
    "rhodochrosite": "si_rhodochrosite",
    "sylvite": "si_sylvite",
}


def _si_value(sample: Mapping[str, object], label: str) -> float:
    key = SI_KEYS.get(label)
    if key is None:
        return float("nan")
    value = sample.get(key)
    if value is None and label == "gypsum":
        value = sample.get("si_anhydrite")
    return float(value) if value is not None else float("nan")


def build_edge_bounds(
    phreeqc_results: Mapping[str, Mapping[str, object]],
    edges: Iterable[object],
    labels: List[str],
    mineral_mask: List[bool],
    config: Config,
) -> Dict[str, Dict[str, object]]:
    bounds: Dict[str, Dict[str, object]] = {}
    tau = config.si_threshold_tau
    inf = float("inf")

    for edge in build_edges(edges):
        entry: Dict[str, object] = {
            "lb": [-inf] * len(labels),
            "ub": [inf] * len(labels),
            "constraints_active": {},
            "si_u": {},
            "si_v": {},
            "phreeqc_ok": False,
            "charge_error": None,
            "skipped_reason": None,
        }
        sample_u = phreeqc_results.get(edge.u)
        sample_v = phreeqc_results.get(edge.v)

        if not sample_u or not sample_v:
            entry["skipped_reason"] = "phreeqc_missing"
        elif not sample_u.get("phreeqc_ok") or not sample_v.get("phreeqc_ok"):
            entry["skipped_reason"] = "phreeqc_unavailable"
        else:
            entry["phreeqc_ok"] = True
            entry["charge_error"] = sample_v.get("charge_balance")

        for idx, label in enumerate(labels):
            if mineral_mask[idx]:
                si_u = _si_value(sample_u or {}, label)
                si_v = _si_value(sample_v or {}, label)
                entry["si_u"][label] = si_u
                entry["si_v"][label] = si_v
                if not entry["phreeqc_ok"]:
                    if label in config.signed_reaction_labels:
                        entry["lb"][idx] = -inf
                        entry["ub"][idx] = inf
                        entry["constraints_active"][label] = "free"
                    else:
                        entry["lb"][idx] = 0.0
                        entry["ub"][idx] = inf
                        entry["constraints_active"][label] = "dissolution_only"
                elif si_u < -tau and si_v < -tau:
                    entry["lb"][idx] = 0.0
                    entry["ub"][idx] = inf
                    entry["constraints_active"][label] = "dissolution_only"
                elif entry["phreeqc_ok"] and si_v > tau:
                    entry["lb"][idx] = -inf
                    entry["ub"][idx] = 0.0
                    entry["constraints_active"][label] = "precipitation_only"
                else:
                    entry["lb"][idx] = -inf
                    entry["ub"][idx] = inf
                    entry["constraints_active"][label] = "free"
            else:
                if label in config.signed_reaction_labels:
                    entry["lb"][idx] = -inf
                    entry["ub"][idx] = inf
                    entry["constraints_active"][label] = "free"
                else:
                    entry["lb"][idx] = 0.0
                    entry["ub"][idx] = inf
                    entry["constraints_active"][label] = "dissolution_only"

        if not config.constraints_hard:
            entry["lb"] = None
            entry["ub"] = None

        bounds[edge.edge_id] = entry

    return bounds
