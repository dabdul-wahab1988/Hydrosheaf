"""M3 prototype: graph topology as a bomb-peak BRANCH SELECTOR on real USGS data.

Motivation
==========
The graph prior evaluated in the main M3 benchmark regularises the age *field*:
for a directed edge u -> v it pulls the downstream age toward the upstream age
(Section 2.7, eq. 5).  That is a smoother.  A smoother cannot add information; it
trades variance for bias, which is why the corrected-forcing benchmark finds it
neutral at best (weak parameter smoothing -0.006 log10) and why it does not
improve withheld-tritium prediction (+12.3%).

The tritium input function is *double-valued* across the bomb peak: a single
measured concentration aliases to a young and an old solution (e.g. ~8 TU maps to
roughly 3 yr or roughly 50 yr).  This is a discrete ambiguity, not a continuous
smoothness problem.  Flow ordering can resolve it: if an upstream neighbour is
confidently old, the downstream node's young alias is infeasible.

`run_m3_network_dating_demo.py` demonstrates this on a synthetic twin with known
ages (within-factor-2 0.63 -> 0.84 on ambiguous nodes).  This script is the
real-data test: it applies the same selection rule to the USGS public-supply
benchmark, against the published USGS/TracerLPM reference ages, with a randomised
graph as the negative control.

Design
======
Per site:
  1. Enumerate the alias set of piston-flow ages consistent with the measured 3H,
     using that site's own WISER input history (continued to the present).
  2. `age_single`  = the alias with the best tritium match (single-node baseline).
  3. `age_graph`   = topological root->leaf selection over the depth-constrained
     graph, restricting each node's aliases to those no younger than its selected
     upstream neighbour, then taking the best tritium match among the feasible set.
  4. `age_random`  = the same rule over a randomised within-study-unit graph
     (negative control, seed 42), matched in edge count.

Reported on the AMBIGUOUS subset (>= 2 aliases separated by > 15 yr), because that
is the only subset where the graph can act.  Sites with a unique alias are
identical by construction under all three rules.

This is a capability test on real data, not a claim of field validation: the
reference ages remain USGS LPM products, not observed ages.
"""
from __future__ import annotations

import argparse
import json
import sys
import warnings
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[3]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
sys.path.insert(0, str(Path(__file__).resolve().parent))

warnings.filterwarnings("ignore")

from hydrosheaf.nuclear.lpm import convolve_input
from hydrosheaf.nuclear.tracer_inputs import SiteInputContext, build_site_tracer_histories

RESULTS = Path(__file__).resolve().parents[1] / "results"
RESULTS.mkdir(parents=True, exist_ok=True)

SEED = 42
LAMBDA_Y = np.log(2.0) / 12.32
AGE_GRID = np.linspace(0.5, 90.0, 460)
AMBIGUITY_SEPARATION_YR = 15.0
ORDERING_TOLERANCE_YR = 2.0
MAX_EDGE_KM = 25.0
MIN_DEPTH_DROP_M, MAX_DEPTH_DROP_M = 5.0, 200.0


def _haversine_km(lat1, lon1, lat2, lon2):
    r = 6371.0
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dp, dl = p2 - p1, np.radians(lon2 - lon1)
    a = np.sin(dp / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2) ** 2
    return float(2 * r * np.arcsin(np.sqrt(a)))


def _tritium_curve(history, sample_year):
    """Forward-model 3H for every age on the grid, for one site."""
    years = np.asarray(history.years, dtype=float)
    values = np.asarray(history.values, dtype=float)
    return np.array(
        [convolve_input(sample_year, float(t), years, values, LAMBDA_Y, model_type="PFM")
         for t in AGE_GRID]
    )


def candidate_ages(obs, curve):
    """Distinct piston-flow ages consistent with a measured tritium value."""
    tol = max(0.4, 0.12 * obs)
    hits = AGE_GRID[np.abs(curve - obs) <= tol]
    if hits.size == 0:
        return [float(AGE_GRID[int(np.argmin(np.abs(curve - obs)))])]
    clusters, cur = [], [hits[0]]
    for h in hits[1:]:
        if h - cur[-1] <= 3.0:
            cur.append(h)
        else:
            clusters.append(cur)
            cur = [h]
    clusters.append(cur)
    reps = []
    for cl in clusters:
        cl = np.asarray(cl)
        pred = np.interp(cl, AGE_GRID, curve)
        reps.append(float(cl[int(np.argmin(np.abs(pred - obs)))]))
    return sorted({round(r, 1) for r in reps})


def build_depth_graph(frame):
    """Depth-constrained edges: same study unit, downstream deeper, short distance.

    Mirrors the graph-family rules used by run_m3_real_usgs_graph_benchmark.py so
    that this prototype is comparable with the main graph benchmark.
    """
    edges = []
    for _, grp in frame.groupby("StudyUnit"):
        if len(grp) < 5:
            continue
        rows = grp.sort_values("depth_m").to_dict("records")
        for i, up in enumerate(rows):
            best = None
            for down in rows[i + 1:]:
                drop = down["depth_m"] - up["depth_m"]
                if not (MIN_DEPTH_DROP_M <= drop <= MAX_DEPTH_DROP_M):
                    continue
                d = _haversine_km(up["lat"], up["lon"], down["lat"], down["lon"])
                if d <= MAX_EDGE_KM and (best is None or d < best[1]):
                    best = (down["site_id"], d)
            if best is not None:
                edges.append((up["site_id"], best[0]))
    return edges


def build_random_graph(frame, n_edges, rng):
    """Within-study-unit random edges, matched in count (negative control)."""
    pool = []
    for _, grp in frame.groupby("StudyUnit"):
        ids = grp["site_id"].tolist()
        if len(ids) < 5:
            continue
        for _ in range(len(ids)):
            u, v = rng.choice(ids, size=2, replace=False)
            pool.append((str(u), str(v)))
    rng.shuffle(pool)
    return pool[:n_edges]


def select_along_graph(nodes, edges):
    """Root->leaf alias selection enforcing downstream >= upstream age."""
    children, parents = {}, {}
    for u, v in edges:
        if u in nodes and v in nodes and v not in parents:
            children.setdefault(u, []).append(v)
            parents[v] = u
    selected = {}
    roots = [n for n in nodes if n not in parents]
    stack = list(roots)
    seen = set()
    while stack:
        n = stack.pop()
        if n in seen:
            continue
        seen.add(n)
        p = parents.get(n)
        floor = selected[p] - ORDERING_TOLERANCE_YR if p in selected else -np.inf
        cands = nodes[n]["candidates"]
        feasible = [c for c in cands if c >= floor]
        pool = feasible if feasible else cands
        selected[n] = min(pool, key=lambda c: abs(np.interp(c, AGE_GRID, nodes[n]["curve"]) - nodes[n]["obs"]))
        stack.extend(children.get(n, []))
    for n in nodes:  # nodes in cycles or unreachable keep their single-node choice
        selected.setdefault(n, nodes[n]["age_single"])
    return selected


MRVA_DIR = (REPO / "M2" / "m2_benchmark" / "external" / "usgs_age" / "input"
            / "MRVA_GroundwaterAge_2018_20")


def load_mrva_dataset():
    """Mississippi River Valley alluvial aquifer release (Gratzer et al., 2025a).

    A single contiguous alluvial aquifer, so well spacing is far tighter than the
    national public-supply network (median nearest-neighbour ~15 km against
    40-95 km), and the release includes nested piezometers at common locations.
    That is the regime in which a flow-ordering prior can actually act.

    Well attributes come from the Table1_Wells shapefile's DBF, which is attached
    to the ScienceBase item as a shapefile facet rather than as a listed file.
    """
    wells = pd.read_csv(MRVA_DIR / "Table1_Wells_attributes.csv")
    ages = pd.read_csv(MRVA_DIR / "Table2_MeanAgeSummary.csv", encoding="utf-8-sig")
    trit = pd.read_csv(MRVA_DIR / "Table3_tritium.csv", encoding="utf-8-sig")

    for col in ("ALT_VA", "WELL_DEPTH", "DEC_LAT_VA", "DEC_LONG_V"):
        wells[col] = pd.to_numeric(wells[col], errors="coerce")
    # -9999 is the release's NoData flag.
    wells = wells[(wells.WELL_DEPTH > 0) & (wells.ALT_VA > -9000)]

    trit["Trit_TU"] = pd.to_numeric(trit["Trit_TU"], errors="coerce")
    trit["year"] = pd.to_datetime(trit["date"], errors="coerce").dt.year
    trit = (trit.dropna(subset=["Trit_TU"])
                .groupby("siteag", as_index=False)
                .agg(tritium_TU=("Trit_TU", "median"), sample_year=("year", "median")))

    ages["avg_mean_age"] = pd.to_numeric(ages["avg_mean_age"], errors="coerce")

    df = (wells.merge(ages[["siteag", "avg_mean_age"]], on="siteag", how="inner")
               .merge(trit, on="siteag", how="inner"))
    df = df.rename(columns={
        "siteag": "site_id", "DEC_LAT_VA": "lat", "DEC_LONG_V": "lon",
        "WELL_DEPTH": "depth_m", "avg_mean_age": "reference_age_years",
    })
    # One contiguous aquifer: a single grouping unit, unlike the national release.
    df["StudyUnit"] = "Mississippi River Valley alluvial aquifer (MRVA)"
    df["sample_year"] = df["sample_year"].fillna(2019.0)
    return df


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--max-sites", type=int, default=0, help="0 = all eligible sites")
    ap.add_argument("--dataset", choices=("national", "mrva"), default="national")
    args = ap.parse_args(argv)

    if args.dataset == "mrva":
        df = load_mrva_dataset()
    else:
        import run_m3_usgs_benchmark as bench
        df = bench.load_benchmark_dataset()
    need = ["site_id", "lat", "lon", "sample_year", "depth_m", "StudyUnit",
            "tritium_TU", "reference_age_years"]
    frame = df[[c for c in need if c in df.columns]].copy()
    # Several USGS columns arrive as object dtype with sentinel strings; coerce
    # before filtering so comparisons do not raise on mixed types.
    for col in ("lat", "lon", "sample_year", "depth_m", "tritium_TU", "reference_age_years"):
        frame[col] = pd.to_numeric(frame[col], errors="coerce")
    frame = frame.dropna(subset=need)
    frame = frame[(frame.tritium_TU > 0.05) & (frame.reference_age_years > 0)]
    if args.max_sites:
        frame = frame.head(args.max_sites)
    print(f"eligible sites with 3H + depth + reference age: {len(frame)}", flush=True)

    nodes = {}
    for rec in frame.to_dict("records"):
        ctx = SiteInputContext(str(rec["site_id"]), float(rec["sample_year"]),
                               float(rec["lat"]), float(rec["lon"]))
        hist = build_site_tracer_histories(ctx)["3H"]
        curve = _tritium_curve(hist, float(rec["sample_year"]))
        obs = float(rec["tritium_TU"])
        cands = candidate_ages(obs, curve)
        single = min(cands, key=lambda c: abs(np.interp(c, AGE_GRID, curve) - obs))
        nodes[str(rec["site_id"])] = {
            "curve": curve, "obs": obs, "candidates": cands, "age_single": single,
            "ref": float(rec["reference_age_years"]),
            "ambiguous": int(len(cands) >= 2 and (max(cands) - min(cands) > AMBIGUITY_SEPARATION_YR)),
        }

    depth_edges = build_depth_graph(frame)
    rng = np.random.default_rng(SEED)
    rand_edges = build_random_graph(frame, len(depth_edges), rng)
    print(f"depth-constrained edges: {len(depth_edges)} | randomised control edges: {len(rand_edges)}",
          flush=True)

    sel_graph = select_along_graph(nodes, depth_edges)
    sel_rand = select_along_graph(nodes, rand_edges)

    rows = []
    for sid, n in nodes.items():
        rows.append({
            "site_id": sid, "tritium_TU": n["obs"], "n_alias": len(n["candidates"]),
            "ambiguous": n["ambiguous"], "reference_age_years": n["ref"],
            "age_single": n["age_single"],
            "age_graph": sel_graph[sid], "age_random": sel_rand[sid],
        })
    out = pd.DataFrame(rows)
    tag = "" if args.dataset == "national" else f"_{args.dataset}"
    out.to_csv(RESULTS / f"m3_branch_selection_benchmark{tag}.csv", index=False)

    def metrics(sub, col):
        est = sub[col].clip(lower=0.01).to_numpy()
        ref = sub["reference_age_years"].to_numpy()
        e = np.log10(est) - np.log10(ref)
        return {
            "n": int(len(sub)),
            "median_abs_log10_error": float(np.median(np.abs(e))),
            "log10_rmse": float(np.sqrt(np.mean(e ** 2))),
            "within_factor_2": float(np.mean(np.abs(e) <= np.log10(2))),
        }

    summary = []
    for label, sub in (("all_sites", out), ("ambiguous_only", out[out.ambiguous == 1])):
        for mode, col in (("single_node", "age_single"),
                          ("graph_branch_selection", "age_graph"),
                          ("randomised_control", "age_random")):
            summary.append({"subset": label, "mode": mode, **metrics(sub, col)})
    summ = pd.DataFrame(summary)
    summ.to_csv(RESULTS / f"m3_branch_selection_summary{tag}.csv", index=False)

    print("\n" + summ.to_string(index=False), flush=True)

    (RESULTS / f"m3_branch_selection_manifest{tag}.json").write_text(json.dumps({
        "dataset": args.dataset,
        "run_utc": datetime.now(timezone.utc).isoformat(),
        "n_sites": int(len(out)),
        "n_ambiguous": int(out.ambiguous.sum()),
        "n_depth_edges": len(depth_edges),
        "n_random_edges": len(rand_edges),
        "seed": SEED,
        "ambiguity_separation_yr": AMBIGUITY_SEPARATION_YR,
        "ordering_tolerance_yr": ORDERING_TOLERANCE_YR,
    }, indent=2), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
