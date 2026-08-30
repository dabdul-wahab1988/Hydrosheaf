"""Run the Ghana field demonstration through the pipeline the manuscript describes.

The field results reported at submission were built by
``scripts/analysis/run_m2_field_benchmarks.py``, which connects each site to its
two nearest neighbours with a 2-D Euclidean KD-tree on (lon, lat) and applies no
directional test at all: no p_ij, no edge_p_min, no edge_gradient_min, no
edge_radius_km, and no sheaf refinement (the pipeline stage report records
``sheaf_refinement: not_requested`` and returns as many fitted edges as it was
given). Section 2.3 of the manuscript describes a different construction.

This script runs both so they can be compared on the same data:

  A. as-run      -- the KD-tree kNN construction that produced the reported
                    numbers, reproduced here for reference
  B. documented  -- infer_edges(method="probabilistic") with elevation as the
                    head proxy, followed by refine_edges_with_sheaf with
                    sheaf_cohomology_enabled=True, i.e. the construction the
                    manuscript actually describes

It reports, for each: candidate and retained edge counts, how many retained
edges run downhill / uphill / flat on the elevation proxy, how many are
reciprocal (both A->B and B->A present), the chemistry-fit distribution, and --
for B -- the sheaf cohomology diagnostics (H0, H1, rank D, affine obstruction
energy, per-edge leverage, cycle obstruction).

Outputs
    M2/m2_benchmark/results/revision/field_pipeline_comparison.csv    per-edge
    M2/m2_benchmark/results/revision/field_pipeline_comparison_summary.csv
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from hydrosheaf import Config  # noqa: E402
from hydrosheaf.api import fit_network_pipeline  # noqa: E402
from hydrosheaf.data.units import mgL_to_mmolL  # noqa: E402
from hydrosheaf.graph.types import Edge  # noqa: E402
from hydrosheaf.inference.network_fit import infer_edges  # noqa: E402

OUT_DIR = ROOT / "M2" / "m2_benchmark" / "results" / "revision"


def _manu_samples() -> list[dict[str, Any]]:
    df = pd.read_csv(ROOT / "data" / "FieldData" / "LowerAnayari" / "manu.csv")
    ions = ["Ca", "Mg", "Na", "K", "HCO3", "Cl", "SO4", "NO3", "F"]
    out = []
    for _, row in df.iterrows():
        s: dict[str, Any] = {
            "site_id": str(row["Sample ID"]),
            "x": row["X coordinate"], "y": row["Y coordinate"], "z": row["Elevation"],
            # keys the documented graph builder reads
            "lon": float(row["X coordinate"]), "lat": float(row["Y coordinate"]),
            "elevation": float(row["Elevation"]),
            "pH": row["pH"], "temp_C": row["Temp"],
            "18O": row["d18O"], "2H": row["d2H"], "type": "sample",
        }
        for ion in ions:
            val = row[ion]
            if isinstance(val, str) and "<" in val:
                val = 0.0005
            s[ion] = mgL_to_mmolL(float(val), ion)
        s["Fe"] = mgL_to_mmolL(float(row["Fe"]), "Fe")
        out.append(s)
    return out


def _talensi_samples() -> list[dict[str, Any]]:
    df = pd.read_csv(ROOT / "data" / "FieldData" / "Talensi_MiningArea" / "talensi.csv")
    ions = ["Ca", "Mg", "Na", "K", "HCO3", "Cl", "SO4", "NO3"]
    out = []
    for _, row in df.iterrows():
        s: dict[str, Any] = {
            "site_id": str(row["Code"]),
            "x": row["Longitude"], "y": row["Latitude"], "z": row["Elevation"],
            "lon": float(row["Longitude"]), "lat": float(row["Latitude"]),
            "elevation": float(row["Elevation"]),
            "pH": row["pH"], "temp_C": row["Temp"],
            "18O": row["d18O"], "2H": row["d2H"], "type": "sample",
        }
        for ion in ions:
            s[ion] = mgL_to_mmolL(float(row[ion]), ion)
        s["Fe"] = mgL_to_mmolL(float(row["Fe"]), "Fe")
        s["F"] = 0.0
        out.append(s)
    return out


def _config(site: str, samples: list[dict[str, Any]], *, documented: bool) -> Config:
    minerals = (["calcite", "dolomite", "gypsum", "albite", "halite", "fluorite"]
                if site == "Manu"
                else ["calcite", "dolomite", "pyrite_oxidation_aerobic", "albite", "halite"])
    cfg = Config(
        ion_order=["Ca", "Mg", "Na", "K", "HCO3", "Cl", "SO4", "NO3", "F", "Fe"],
        weights=[1.0] * 10,
        conservative_weights=[0.01] * 10,
        lmwl_a=7.22, lmwl_b=8.66,
        isotope_enabled=True,
        active_minerals=minerals,
        exchange_enabled=True,
        honest_modeling=True,
        geologic_bias="crystalline",
    )
    # the field sites carry no nuclear tracers; leaving the age stage on makes
    # the pipeline try to treat Cl as a residence-time tracer and fail
    cfg.sheaf_age_enabled = False
    if documented:
        cfg.sheaf_cohomology_enabled = True
    _local_dilute_endmember(cfg, samples, site.lower())
    return cfg


def _local_dilute_endmember(config: Config, samples: list[dict[str, Any]], label: str) -> None:
    candidates = []
    for sample in samples:
        vec, ok = [], True
        for ion in config.ion_order:
            v = sample.get(ion)
            if v is None or not np.isfinite(float(v)):
                ok = False
                break
            vec.append(float(v))
        if ok:
            candidates.append((sum(max(v, 0.0) for v in vec), str(sample["site_id"]), vec, sample))
    if not candidates:
        return
    _, site_id, vec, sample = min(candidates, key=lambda i: i[0])
    name = f"{label}_local_dilute_{site_id}"
    config.mixing_endmembers = {name: vec}
    if "18O" in sample and "2H" in sample:
        config.mixing_endmembers_isotopes = {name: (float(sample["18O"]), float(sample["2H"]))}


def _knn_edges(samples: list[dict[str, Any]], site: str) -> list[Edge]:
    """Exactly the construction used for the submitted field results."""
    from scipy.spatial import cKDTree

    coords = np.array([[s["x"], s["y"]] for s in samples])
    tree = cKDTree(coords)
    edges = []
    for i, s in enumerate(samples):
        _d, idx = tree.query([s["x"], s["y"]], k=3)
        for j in idx:
            if i != int(j):
                edges.append(Edge(u=samples[i]["site_id"], v=samples[int(j)]["site_id"],
                                  edge_id=f"{site}_{i}_{j}"))
    return edges


def _direction_stats(rows: list[dict[str, Any]], elev: dict[str, float]) -> dict[str, Any]:
    down = up = flat = 0
    pairs = {(r["u"], r["v"]) for r in rows}
    for r in rows:
        zu, zv = elev.get(str(r["u"])), elev.get(str(r["v"]))
        if zu is None or zv is None:
            continue
        if zu > zv:
            down += 1
        elif zu < zv:
            up += 1
        else:
            flat += 1
    reciprocal = sum(1 for (u, v) in pairs if (v, u) in pairs)
    return {"downhill": down, "uphill": up, "flat": flat, "reciprocal_edges": reciprocal}


def run_site(site: str, samples: list[dict[str, Any]]) -> tuple[list[dict], list[dict]]:
    elev = {str(s["site_id"]): float(s["elevation"]) for s in samples}
    per_edge: list[dict] = []
    summary: list[dict] = []

    for mode in ("as_run", "documented"):
        cfg = _config(site, samples, documented=mode == "documented")
        if mode == "as_run":
            cand = _knn_edges(samples, site)
            n_cand = len(cand)
            results, report = fit_network_pipeline(samples, cand, cfg)
        else:
            cand = infer_edges(samples, method="probabilistic", config=cfg)
            n_cand = len(cand)
            if not cand:
                summary.append({"site": site, "mode": mode, "n_candidates": 0,
                                "n_retained": 0, "note": "probabilistic builder returned no edges"})
                continue
            results, report = fit_network_pipeline(
                samples, cand, cfg, sheaf_refinement_enabled=True)

        stages = {k: (v.get("status") if isinstance(v, dict) else v)
                  for k, v in (report.get("stage_status", {}) or {}).items()}
        # EdgeResult carries the cohomology diagnostics as fields; p_uv and
        # edge_confidence live on the candidate Edge's attrs (set by the
        # probabilistic builder only).
        cand_attrs = {e.edge_id: (e.attrs or {}) for e in cand}
        rows = []
        for r in results:
            attrs = cand_attrs.get(r.edge_id, {})
            row = {"site": site, "mode": mode, "edge_id": r.edge_id,
                   "u": str(r.u), "v": str(r.v),
                   "chemistry_r2": r.chemistry_r2,
                   "edge_confidence": attrs.get("edge_confidence"),
                   "p_uv": attrs.get("p_uv"),
                   "head_gradient": attrs.get("head_gradient")}
            for k in ("sheaf_h0_dim", "sheaf_h1_dim", "sheaf_obstruction_energy",
                      "sheaf_obstruction_leverage", "sheaf_cycle_obstruction_max",
                      "sheaf_cycle_count", "edge_sheaf_score_global"):
                row[k] = getattr(r, k, None)
            rows.append(row)
        per_edge.extend(rows)

        chem_vals = [r["chemistry_r2"] for r in rows
                     if r["chemistry_r2"] is not None and np.isfinite(float(r["chemistry_r2"]))]
        s = {"site": site, "mode": mode, "n_candidates": n_cand, "n_retained": len(rows),
             "rejected": n_cand - len(rows),
             "median_chemistry_r2": float(np.median(chem_vals)) if chem_vals else None,
             "sheaf_refinement_stage": stages.get("sheaf_refinement"),
             **_direction_stats(rows, elev)}
        for k in ("sheaf_h0_dim", "sheaf_h1_dim", "sheaf_obstruction_energy",
                  "sheaf_cycle_count", "sheaf_cycle_obstruction_max"):
            vals = [r[k] for r in rows if r.get(k) is not None]
            s[k] = vals[0] if vals else None
        lev = [r["sheaf_obstruction_leverage"] for r in rows
               if r.get("sheaf_obstruction_leverage") is not None]
        s["leverage_max"] = float(np.max(lev)) if lev else None
        s["leverage_median"] = float(np.median(lev)) if lev else None
        summary.append(s)
        print(f"  {site:<8} {mode:<11} candidates={n_cand:>4} retained={len(rows):>4} "
              f"rejected={n_cand - len(rows):>4} stage={stages.get('sheaf_refinement')}")
    return per_edge, summary


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    all_edges: list[dict] = []
    all_summary: list[dict] = []
    for site, samples in (("Manu", _manu_samples()), ("Talensi", _talensi_samples())):
        e, s = run_site(site, samples)
        all_edges.extend(e)
        all_summary.extend(s)

    pd.DataFrame(all_edges).to_csv(OUT_DIR / "field_pipeline_comparison.csv", index=False)
    pd.DataFrame(all_summary).to_csv(
        OUT_DIR / "field_pipeline_comparison_summary.csv", index=False)
    print("\n" + json.dumps(all_summary, indent=2, default=str))


if __name__ == "__main__":
    main()
