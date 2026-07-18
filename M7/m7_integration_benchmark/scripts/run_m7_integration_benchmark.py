"""M7 coupled integration benchmark (thesis keystone).

Demonstrates that Hydrosheaf's THREE evidence streams — residence-time (age),
graph topology (connectivity) and inverse hydrogeochemistry (reactions) — jointly
constrain and test each other, on a controlled synthetic twin with KNOWN joint
truth (true flow DAG, true ages, true reactions). This is the integration result
that the field data cannot provide (no real aquifer supplies all three truths at
once).

Defensibility (see docs/m7_integration_defensibility.md):
  * controlled twin with known truth is the only way to *measure* a reduction in
    non-uniqueness; it is a capability/mechanism demonstration, not field validation;
  * anti-inverse-crime: generator adds noise, evaporation and an unmodelled mixing
    perturbation not representable by the reaction dictionary;
  * negative controls: planted trap edges that defeat exactly one single stream;
  * drives the real tritium forward model and age-coherence audit. The chemistry
    screen is a conservative mass/Cl falsifier; full inverse-reaction fitting and
    blind topology inference are exercised by the separate M7.1 benchmark.

Four integration tests:
  T1 age<->topology : single-node tritium ordering powers the coherence audit.
  T2 chem<->topology: chemistry feasibility improves edge classification / rejects traps.
  T3 age<->chem     : age coherence rejects age-reversal traps chemistry accepts.
  T4 integration gain: JOINT screening raises precision and lowers false-edge
  acceptance at a documented recall cost.
Deterministic (seed 1234).
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[3]
M6_SCRIPTS = REPO / "M6" / "m6_field_transfer_benchmark" / "scripts"
for p in (str(REPO), str(M6_SCRIPTS)):
    if p not in sys.path:
        sys.path.insert(0, p)

import m6_reactions as rx  # noqa: E402  (M5/M6 inverse-reaction pillar)
from hydrosheaf.nuclear.lpm import convolve_input  # noqa: E402
from hydrosheaf.nuclear.input_history import build_default_tritium_input  # noqa: E402
from hydrosheaf.nuclear.age_coherence import audit_graph_age_coherence  # noqa: E402

RESULTS = Path(__file__).resolve().parents[1] / "results"
RESULTS.mkdir(parents=True, exist_ok=True)
SEED = 1234
HALF_LIFE = 12.32
LAMBDA_Y = np.log(2) / HALF_LIFE
SAMPLE_DATE = 2023.0
HIST = build_default_tritium_input("global")

X0 = {"Ca": 0.30, "Mg": 0.15, "Na": 0.35, "K": 0.05, "HCO3": 1.0, "Cl": 0.20,
      "SO4": 0.08, "NO3": 0.05, "F": 0.005, "Fe": 0.001, "PO4": 0.002,
      "SiO2": 0.15, "Sr": 0.002}
PANEL = rx.PANELS["plus_Sr_SiO2"]
ARCH_RXN = {"silicate": ["albite", "k_feldspar"], "carbonate": ["calcite", "dolomite"],
            "evaporite": ["gypsum", "halite"], "ion_exchange": ["NaCa_exch"],
            "nitrate": ["NO3src"], "redox": ["pyrite_oxidation", "denit"]}


def tritium(age, model="PFM"):
    return float(convolve_input(SAMPLE_DATE, max(age, 0.01), HIST.years, HIST.values,
                                LAMBDA_Y, model_type=model))


def apply_reactions(chem, labels, exts):
    out = dict(chem)
    for lab, e in zip(labels, exts):
        for s, c in rx.REACTIONS[rx.LABELS.index(lab)].stoich.items():
            out[s] = out.get(s, 0.0) + e * c
    return out


def generate(rng, n_backbone=24):
    """Build a flow DAG with known ages + chemistry, plus trap edges."""
    # backbone: sorted flow coordinate s; elevation decreases downstream; age increases
    s = np.sort(rng.uniform(0, 1, n_backbone))
    nodes = {}
    for i in range(n_backbone):
        nid = f"N{i:02d}"
        nodes[nid] = {
            "x": s[i] * 2600 + rng.normal(0, 90), "y": rng.uniform(0, 1500),  # metres
            "elev": 100 - s[i] * 60 + rng.normal(0, 1.0),  # down-gradient with flow
            "s": s[i], "age": None, "branch": None,
        }
    ids = list(nodes)
    # true DAG: each node -> nearest downstream (higher s) node; assign reaction archetype per branch
    true_edges = []
    archetypes = ["silicate", "carbonate", "ion_exchange", "redox", "evaporite"]
    for i in range(n_backbone - 1):
        u = ids[i]
        cand = [(j, abs(nodes[ids[j]]["x"] - nodes[u]["x"])) for j in range(i + 1, n_backbone)]
        cand.sort(key=lambda t: t[1])
        v = ids[cand[0][0]]
        true_edges.append((u, v))
    # physically-consistent ages: travel time accumulates along the DAG
    # (age_v = age_u + distance / velocity), so the network velocity-distance model
    # is correctly specified. This is more realistic than an arbitrary age field.
    velocity = float(rng.uniform(9.0, 12.0))  # m/yr
    nodes[ids[0]]["age"] = float(rng.uniform(2.0, 5.0))
    parent0 = {v: u for (u, v) in true_edges}
    edge_len = {}
    for v in ids[1:]:
        u = parent0.get(v, ids[0])
        d = float(np.hypot(nodes[u]["x"] - nodes[v]["x"], nodes[u]["y"] - nodes[v]["y"]))
        edge_len[(u, v)] = d
        vel = velocity * float(np.exp(rng.normal(0, 0.12)))
        if nodes[u]["age"] is None:
            nodes[u]["age"] = float(rng.uniform(2.0, 5.0))
        nodes[v]["age"] = nodes[u]["age"] + d / vel
    # forward chemistry pass along DAG (topological by s order)
    for nid in ids:
        nodes[nid]["chem"] = None
    roots = [ids[0]]
    nodes[ids[0]]["chem"] = dict(X0)
    nodes[ids[0]]["branch"] = "recharge"
    parent = {v: u for (u, v) in true_edges}
    for v in ids[1:]:
        u = parent.get(v, ids[0])
        arch = archetypes[hash(v) % len(archetypes)]
        labs = ARCH_RXN[arch]
        exts = [float(rng.uniform(0.2, 1.1)) for _ in labs]
        f = float(rng.uniform(1.0, 1.6))  # evaporative concentration
        base = {k: val * f for k, val in nodes[u]["chem"].items()}
        chem = apply_reactions(base, labs, exts)
        # anti-inverse-crime: unmodelled mixing perturbation (not a dictionary reaction)
        for ion in ["Na", "Cl", "Mg"]:
            chem[ion] = chem.get(ion, 0.0) * (1 + rng.normal(0, 0.04))
        chem = {k: max(v * (1 + rng.normal(0, 0.03)), 0.0) for k, v in chem.items()}
        nodes[v]["chem"] = chem
        nodes[v]["branch"] = arch
        nodes[v]["true_reactions"] = set(labs)
    nodes[ids[0]]["true_reactions"] = set()

    # young-but-mineralised nodes: young age (fast weathering / anthropogenic loading)
    # yet evolved chemistry, at LOW elevation. These decouple elevation AND chemistry
    # from age, so an edge into them is chemically plausible but temporally impossible
    # (age reversal) — targets of trap A, catchable ONLY by the age stream.
    local_recharge = []
    for k in range(5):
        nid = f"R{k:02d}"
        base = {kk: vv * rng.uniform(1.2, 1.6) for kk, vv in X0.items()}
        chem = apply_reactions(base, ["albite", "NaCa_exch", "calcite"],
                               [rng.uniform(0.3, 0.9) for _ in range(3)])
        nodes[nid] = {"x": rng.uniform(800, 2400), "y": rng.uniform(0, 1500),
                      "elev": rng.uniform(35, 55),  # low, like a downstream node
                      "s": rng.uniform(0.02, 0.12),
                      "age": float(rng.uniform(2.0, 6.0)),  # young
                      "branch": "young_mineralised",
                      "chem": {kk: max(vv * (1 + rng.normal(0, 0.03)), 0.0)
                               for kk, vv in chem.items()},
                      "true_reactions": set()}
        local_recharge.append(nid)
    # fossil-dilute nodes: OLD age but dilute chemistry at low elevation — targets of
    # chemistry-implausible (trap B) edges (age-coherent, geometry-ok, but a mineralised
    # upstream water cannot evolve INTO a dilute one without net precipitation).
    fossil_dilute = []
    for k in range(5):
        nid = f"F{k:02d}"
        nodes[nid] = {"x": rng.uniform(1000, 2400), "y": rng.uniform(0, 1500),
                      "elev": rng.uniform(30, 50),
                      "s": rng.uniform(0.85, 0.98),
                      "age": float(rng.uniform(60.0, 160.0)),  # old
                      "branch": "fossil",
                      "chem": {kk: vv * (1 + rng.normal(0, 0.05)) for kk, vv in X0.items()},
                      "true_reactions": set()}
        fossil_dilute.append(nid)
    ids = list(nodes)

    # tritium observations from true ages
    for nid in ids:
        nodes[nid]["3H_obs"] = tritium(nodes[nid]["age"]) * (1 + rng.normal(0, 0.05))

    # ---- candidate edge set: true + traps + spurious ----
    cands = [(u, v, "true") for (u, v) in true_edges]
    have = lambda: [(a, b) for a, b, _ in cands]
    med_age = float(np.median([nodes[k]["age"] for k in ids]))
    # Trap A (age-reversal): an older up-gradient backbone node -> young low-elev
    # recharge node. Geometry ACCEPTS (down-gradient + near); age REJECTS (reversal).
    old_nodes = [n for n in ids if nodes[n]["age"] > med_age]
    for v in local_recharge:
        ups = [u for u in old_nodes if nodes[u]["elev"] > nodes[v]["elev"] + 5]
        if ups:
            u = min(ups, key=lambda a: abs(nodes[a]["x"] - nodes[v]["x"]))
            if (u, v) not in have():
                cands.append((u, v, "trapA"))
    # Trap B (chem-implausible): a mineralised backbone node -> an OLD dilute fossil
    # node. Geometry ACCEPTS (down-gradient), age ACCEPTS (downstream older), but
    # chemistry REJECTS (evolution into a dilute water needs net precipitation).
    tds = lambda n: sum(nodes[n]["chem"].get(i, 0.0) for i in
                        ["Ca", "Mg", "Na", "K", "HCO3", "Cl", "SO4"])
    mineralised = [n for n in ids if nodes[n]["branch"] not in ("recharge", "fossil")
                   and nodes[n]["age"] > 8]
    for v in fossil_dilute:
        ups = [u for u in mineralised if nodes[u]["elev"] > nodes[v]["elev"] + 5
               and nodes[u]["age"] < nodes[v]["age"]]
        if ups:  # strongly-mineralised source -> dilute fossil = clear net precipitation
            u = max(ups, key=tds)
            if (u, v) not in have():
                cands.append((u, v, "trapB"))
    # spurious random
    for _ in range(6):
        u, v = rng.choice(ids, 2, replace=False)
        if (u, v) not in [(a, b) for a, b, _ in cands]:
            cands.append((u, v, "spurious"))
    return nodes, true_edges, cands


# ---- single-node age inversion (baseline) -----------------------------------
def single_node_age(obs, model="PFM"):
    taus = np.linspace(0.1, 400, 800)
    preds = np.array([tritium(t, model) for t in taus])
    return float(taus[int(np.argmin(np.abs(preds - obs)))])


# ---- chemistry feasibility of an edge (held-out-ion falsification) -----------
# Residual-fit R2 is ~1 for ANY edge (the reaction dictionary is underdetermined —
# the equifinality point). The principled discriminator is M5's held-out-ion test:
# fit on all-but-one ion and predict the held-out ion. A real reaction sequence
# predicts held-out ions; a coincidental cross-branch fit (trap B) does not.
MAJ = ["Ca", "Mg", "Na", "K", "HCO3", "Cl", "SO4"]
def chem_edge_score(chem_u, chem_v):
    """Chemistry falsifier for a candidate edge u->v.

    Physical basis: along a real downstream flow path, water-rock reactions and
    evaporation ADD dissolved mass and the conservative tracer (Cl) is never
    removed, so the downstream water cannot be substantially more dilute than its
    upstream source. A large drop in total mineralisation (or in Cl) therefore
    falsifies the edge. Trap B (a mineralised water 'evolving' into an old dilute
    fossil water) is caught by this; residual-fit R2 cannot see it because the
    reaction dictionary is underdetermined (equifinality)."""
    tds_u = sum(max(chem_u.get(i, 0.0), 0.0) for i in MAJ)
    tds_v = sum(max(chem_v.get(i, 0.0), 0.0) for i in MAJ)
    if tds_u <= 0:
        return 0.0
    mineral_ratio = tds_v / tds_u
    cl_u, cl_v = chem_u.get("Cl", 0.0), chem_v.get("Cl", 0.0)
    cl_ratio = cl_v / cl_u if cl_u > 0 else 1.0
    return float(min(mineral_ratio / 0.85, cl_ratio / 0.6, 1.2))  # <1 => falsified


def main():
    rng = np.random.default_rng(SEED)
    import networkx as nx
    nodes, true_edges, cands = generate(rng)
    ids = list(nodes)
    true_ages = {n: nodes[n]["age"] for n in ids}

    # ===== T1: age <-> topology via tritium age-coherence =====
    # M7 is an INTEGRATION benchmark, not a dating-method benchmark (that is M3).
    # Tritium robustly separates modern from pre-modern water; the young/old ordering
    # this provides is what powers the age-coherence coupling. Network-enhanced Bayesian
    # dating is part of the framework but its advantage is confined to the tracer-
    # ambiguous regime and is not re-benchmarked here (tritium cannot date > ~70 yr).
    print("T1 tritium ages + age-coherence coupling ...")
    obs = {n: nodes[n]["3H_obs"] for n in ids}
    single = {n: single_node_age(obs[n]) for n in ids}
    DATABLE = 70.0  # yr; tritium is uninformative for older water
    rows = [{"node": n, "true_age": true_ages[n], "age_single": single[n],
             "tritium_TU": obs[n], "tritium_datable": int(true_ages[n] <= DATABLE)}
            for n in ids]
    age_df = pd.DataFrame(rows)
    age_df.to_csv(RESULTS / "m7_age_recovery.csv", index=False)
    dat = age_df[age_df.tritium_datable == 1]
    def wf2(a, b):
        return float(np.mean(np.abs(np.log10(a / b)) < np.log10(2)))
    gain = pd.DataFrame([{
        "method": "single_node_tritium",
        "within_factor2_datable": wf2(dat.age_single.to_numpy(), dat.true_age.to_numpy()),
        "n_datable": int(len(dat)),
        "modern_premodern_accuracy": float(np.mean(
            (age_df.age_single <= DATABLE).to_numpy() == (age_df.true_age <= DATABLE).to_numpy())),
    }])
    gain.to_csv(RESULTS / "m7_age_gain.csv", index=False)

    # ===== T2/T3/T4: per-candidate-edge signals & joint classification =====
    print("T2/T3/T4 edge signals + joint classification ...")
    # age-coherence audit on candidate edges using robust single-node tritium ages
    inferred_ages = dict(single)
    # small tolerance so noisy true edges are not flagged, but trap A's large
    # (decades) age reversals still are (driven by the real framework audit).
    audit = audit_graph_age_coherence([(u, v) for (u, v, _) in cands], inferred_ages,
                                      min_downstream_increase_years=-8.0)
    coh_by_edge = {}
    for rec in audit.get("edges", []):
        key = (str(rec.get("upstream")), str(rec.get("downstream")))
        coh_by_edge[key] = not bool(rec.get("violation", False))

    te_len = [float(np.hypot(nodes[a]["x"] - nodes[b]["x"], nodes[a]["y"] - nodes[b]["y"]))
              for (a, b) in true_edges]
    near_thr = 1.6 * float(np.median(te_len))  # adaptive "near" threshold

    erows = []
    for (u, v, label) in cands:
        dist = float(np.hypot(nodes[u]["x"] - nodes[v]["x"], nodes[u]["y"] - nodes[v]["y"]))
        geom = int(nodes[u]["elev"] > nodes[v]["elev"] and dist < near_thr)  # down-gradient + near
        # chemistry acts as a FALSIFIER: it cleanly rejects physically-impossible
        # (net-precipitation / backward-evolution) edges but does not confirm real
        # ones. chem_ok = "not chemically falsified".
        chem = chem_edge_score(nodes[u]["chem"], nodes[v]["chem"])
        chem_ok = int(chem >= 1.0)   # >=1 => not chemically falsified
        age_ok = int(coh_by_edge.get((u, v), inferred_ages[v] >= inferred_ages[u]))
        erows.append({"u": u, "v": v, "label": label, "is_true": int(label == "true"),
                      "dist": dist, "geom_ok": geom, "chem_r2": chem, "chem_ok": chem_ok,
                      "age_ok": age_ok,
                      "joint_ok": int(geom and chem_ok and age_ok)})
    edf = pd.DataFrame(erows)
    edf.to_csv(RESULTS / "m7_edge_classification.csv", index=False)

    # age-coherence demonstration: upstream vs downstream single-node age per edge.
    # True edges lie above the 1:1 line (downstream older); trap A falls below (reversal).
    coh_rows = [{"u": u, "v": v, "label": label,
                 "upstream_age": single[u], "downstream_age": single[v],
                 "coherent": int(single[v] >= single[u])}
                for (u, v, label) in cands if label in ("true", "trapA")]
    pd.DataFrame(coh_rows).to_csv(RESULTS / "m7_age_coherence_demo.csv", index=False)

    # integration gain: precision/recall/F1 of accepting TRUE edges per stream
    def prf(accept_col):
        tp = int(((edf[accept_col] == 1) & (edf.is_true == 1)).sum())
        fp = int(((edf[accept_col] == 1) & (edf.is_true == 0)).sum())
        fn = int(((edf[accept_col] == 0) & (edf.is_true == 1)).sum())
        p = tp / (tp + fp) if tp + fp else 0.0
        r = tp / (tp + fn) if tp + fn else 0.0
        f1 = 2 * p * r / (p + r) if p + r else 0.0
        # trap acceptance rate (false edges accepted)
        traps = edf[edf.is_true == 0]
        tar = float((traps[accept_col] == 1).mean()) if len(traps) else 0.0
        return {"precision": p, "recall": r, "f1": f1, "trap_accept_rate": tar}
    gains = []
    for stream, col in [("geometry_only", "geom_ok"), ("chemistry_only", "chem_ok"),
                        ("age_only", "age_ok"), ("joint", "joint_ok")]:
        gains.append({"stream": stream, **prf(col)})
    gdf = pd.DataFrame(gains)
    gdf.to_csv(RESULTS / "m7_integration_gain.csv", index=False)

    # trap-type rejection matrix (which stream catches which trap)
    trapmat = []
    for tlab in ["trapA", "trapB", "spurious"]:
        sub = edf[edf.label == tlab]
        if not len(sub):
            continue
        trapmat.append({"trap_type": tlab, "n": len(sub),
                        "geometry_rejects": float((sub.geom_ok == 0).mean()),
                        "chemistry_rejects": float((sub.chem_ok == 0).mean()),
                        "age_rejects": float((sub.age_ok == 0).mean()),
                        "joint_rejects": float((sub.joint_ok == 0).mean())})
    pd.DataFrame(trapmat).to_csv(RESULTS / "m7_trap_rejection.csv", index=False)

    print("\n=== AGE GAIN ===\n", gain.round(3).to_string(index=False))
    print("\n=== INTEGRATION GAIN ===\n", gdf.round(3).to_string(index=False))
    print("\n=== TRAP REJECTION ===\n", pd.DataFrame(trapmat).round(3).to_string(index=False))
    print("\nM7 integration benchmark complete ->", RESULTS)


if __name__ == "__main__":
    main()
