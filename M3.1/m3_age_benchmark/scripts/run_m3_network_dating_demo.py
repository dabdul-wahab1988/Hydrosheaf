"""M3 network-enhanced dating — controlled ambiguity demonstration.

Real-data finding (m3_real_usgs_graph_benchmark): graph priors do NOT improve age
RMSE on well-constrained USGS tracers; they only enforce age-ordering consistency
(severe violations 306 -> ~0). This benchmark shows the CONDITION under which
network-enhanced dating genuinely improves *accuracy*: the tritium bomb-peak
AMBIGUITY regime, where a single measured value aliases to two very different ages
(e.g. 8.4 TU -> ~3 yr OR ~50 yr) and single-node inversion is degenerate. Flow-
ordering along the graph (downstream must be at least as old as upstream) selects
the correct alias.

Controlled twin with KNOWN ages; drives the real Hydrosheaf tritium model
(nuclear.lpm.convolve_input, nuclear.input_history). Deterministic (seed 1234).
This is a capability demonstration, not field validation.
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))
from hydrosheaf.nuclear.lpm import convolve_input
from hydrosheaf.nuclear.input_history import build_default_tritium_input

RESULTS = Path(__file__).resolve().parents[1] / "results"
RESULTS.mkdir(parents=True, exist_ok=True)
SEED = 1234
SAMPLE_DATE = 2023.0
HIST = build_default_tritium_input("global")
LAMBDA_Y = np.log(2) / 12.32
GRID = np.linspace(0.5, 90.0, 460)
GRID_TU = np.array([convolve_input(SAMPLE_DATE, t, HIST.years, HIST.values,
                                   LAMBDA_Y, model_type="PFM") for t in GRID])


def tritium(age):
    return float(convolve_input(SAMPLE_DATE, max(age, 0.01), HIST.years, HIST.values,
                                LAMBDA_Y, model_type="PFM"))


def candidate_ages(obs):
    """Distinct tritium-consistent ages (alias solutions) for a measured TU."""
    tol = max(0.4, 0.12 * obs)
    hits = GRID[np.abs(GRID_TU - obs) <= tol]
    if len(hits) == 0:
        return [float(GRID[int(np.argmin(np.abs(GRID_TU - obs)))])]
    clusters, cur = [], [hits[0]]
    for h in hits[1:]:
        if h - cur[-1] <= 3.0:
            cur.append(h)
        else:
            clusters.append(cur); cur = [h]
    clusters.append(cur)
    reps = []
    for cl in clusters:
        cl = np.array(cl)
        reps.append(float(cl[int(np.argmin(np.abs([tritium(c) for c in cl] - np.array(obs))))]))
    return sorted(set(round(r, 1) for r in reps))


def generate(rng, n_chains=9):
    """Flow chains bracketing ambiguous-zone nodes between clearer neighbours."""
    nodes, edges = {}, []
    nid = 0
    for c in range(n_chains):
        # young root -> ambiguous mid (bomb-peak alias) -> older tail
        root_age = float(rng.uniform(4, 14))
        mid_age = float(rng.uniform(45, 56))     # aliases to a young (~3-9 yr) solution
        tail_age = float(rng.uniform(62, 80))    # clearly old (near/after bomb peak)
        chain_ages = [root_age, mid_age, tail_age]
        if rng.random() < 0.5:  # some chains add a 4th, unambiguous young/old node
            chain_ages.insert(1, float(rng.uniform(15, 28)))
        prev = None
        for a in chain_ages:
            name = f"W{nid:02d}"; nid += 1
            nodes[name] = {"true_age": a, "3H_obs": tritium(a) * (1 + rng.normal(0, 0.05))}
            if prev is not None:
                edges.append((prev, name))
            prev = name
    return nodes, edges


def main():
    rng = np.random.default_rng(SEED)
    nodes, edges = generate(rng)
    ids = list(nodes)

    # candidate alias sets + ambiguity flag
    for n in ids:
        cands = candidate_ages(nodes[n]["3H_obs"])
        nodes[n]["candidates"] = cands
        nodes[n]["ambiguous"] = int(len(cands) >= 2 and (max(cands) - min(cands) > 15))
        # single-node estimate = global best tritium match
        nodes[n]["age_single"] = min(cands, key=lambda c: abs(tritium(c) - nodes[n]["3H_obs"]))

    # network estimate: topological (root->leaf) selection enforcing downstream >= upstream
    children = {n: [] for n in ids}
    parents = {n: None for n in ids}
    for (u, v) in edges:
        children[u].append(v); parents[v] = u
    roots = [n for n in ids if parents[n] is None]
    order = []
    stack = list(roots)
    while stack:
        x = stack.pop(); order.append(x); stack.extend(children[x])
    for n in order:
        p = parents[n]
        floor = nodes[p]["age_network"] - 2.0 if p is not None else -np.inf
        feas = [c for c in nodes[n]["candidates"] if c >= floor]
        pool = feas if feas else nodes[n]["candidates"]
        nodes[n]["age_network"] = min(pool, key=lambda c: abs(tritium(c) - nodes[n]["3H_obs"]))

    rows = [{"node": n, "true_age": nodes[n]["true_age"], "tritium_TU": nodes[n]["3H_obs"],
             "n_alias": len(nodes[n]["candidates"]), "ambiguous": nodes[n]["ambiguous"],
             "age_single": nodes[n]["age_single"], "age_network": nodes[n]["age_network"]}
            for n in ids]
    df = pd.DataFrame(rows)
    df.to_csv(RESULTS / "m3_network_dating_demo.csv", index=False)

    def wf2(est, tru):
        return float(np.mean(np.abs(np.log10(est / tru)) < np.log10(2)))
    def rmse_log(est, tru):
        return float(np.sqrt(np.mean((np.log10(est) - np.log10(tru)) ** 2)))

    def summ(sub, label):
        t = sub.true_age.to_numpy()
        return {"subset": label, "n": len(sub),
                "wf2_single": wf2(sub.age_single.to_numpy(), t),
                "wf2_network": wf2(sub.age_network.to_numpy(), t),
                "rmse_log_single": rmse_log(sub.age_single.to_numpy(), t),
                "rmse_log_network": rmse_log(sub.age_network.to_numpy(), t)}
    out = pd.DataFrame([summ(df, "all"),
                        summ(df[df.ambiguous == 1], "ambiguous"),
                        summ(df[df.ambiguous == 0], "unambiguous")])
    out.to_csv(RESULTS / "m3_network_dating_demo_summary.csv", index=False)

    # age-ordering violations (consistency) single vs network
    def violations(col):
        return int(sum(1 for (u, v) in edges if nodes[v][col] < nodes[u][col] - 2.0))
    cons = pd.DataFrame([{"method": "single_node", "ordering_violations": violations("age_single")},
                         {"method": "network", "ordering_violations": violations("age_network")}])
    cons.to_csv(RESULTS / "m3_network_dating_demo_consistency.csv", index=False)

    print("=== ACCURACY (single vs network) ===")
    print(out.round(3).to_string(index=False))
    print("\n=== ORDERING VIOLATIONS ===")
    print(cons.to_string(index=False))
    print(f"\nnodes {len(df)}, ambiguous {int(df.ambiguous.sum())}, edges {len(edges)}")
    print("M3 network-dating demo ->", RESULTS)


if __name__ == "__main__":
    main()
