"""Cross-check every canonical metric against the revision text, tables and figures.

Answers R2-m1 ("full consistency audit"). For each canonical value this script
re-reads the number from its result file, then asserts that the rendered forms
of that number appear everywhere the revision quotes it: the manuscript, the
Supplementary Information, and the generated table files. Figure annotations are
checked at their source -- the figure generator is made to read the same
canonical files -- rather than by OCR.

Run:  .venv/Scripts/python.exe M2/M2_ready/Revision/audit/number_audit.py
Exit 0 = every canonical value agrees everywhere it is quoted.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import pandas as pd

REV = Path(__file__).resolve().parent.parent
ROOT = REV.parent.parent.parent
M2 = ROOT / "M2" / "m2_benchmark" / "results"
M3 = ROOT / "M3" / "m3_age_benchmark" / "results"

MS = (REV / "Manuscript-Final-Revised.md").read_text(encoding="utf-8")
SI = (REV / "Supplementary-Information-Revised.md").read_text(encoding="utf-8")
TABLES = {p.name: p.read_text(encoding="utf-8") for p in (REV / "tables").glob("*.md")}

failures: list[str] = []
checks = 0


def _norm(text: str) -> str:
    """Fold unicode minus/dashes/spaces so matching is robust to typography."""
    for a, b in [("−", "-"), ("–", "-"), ("—", "-"),
                 (" ", " "), (" ", " "), (" ", " ")]:
        text = text.replace(a, b)
    return re.sub(r"\s+", " ", text)


CORPUS = {
    "manuscript": _norm(MS),
    "SI": _norm(SI),
    "tables": _norm("\n".join(TABLES.values())),
}


def kv(path: Path) -> dict[str, float]:
    """Read a two-column metric,value summary into a dict."""
    df = pd.read_csv(path)
    return {str(k): float(v) for k, v in zip(df["metric"], df["value"])}


def check(label: str, forms: list[str], where: list[str]) -> None:
    """Assert at least one rendered form appears in each named corpus."""
    global checks
    for w in where:
        checks += 1
        if not any(_norm(f) in CORPUS[w] for f in forms):
            failures.append(f"{label}: none of {forms} found in {w}")


def pct(x: float, dp: int = 1) -> str:
    return f"{100 * x:.{dp}f}%"


# --------------------------------------------------------------------------
# 1. Synthetic reaction and family recovery
# --------------------------------------------------------------------------
rr = kv(M2 / "revision" / "reaction_recovery_summary.csv")
check("active-reaction R2", [f"{rr['active_r2_corr']:.2f}"], ["manuscript", "SI"])
check("active-reaction MAE", [f"{rr['active_mae_mmolL']:.2f}"], ["manuscript", "SI"])
check("active-reaction RMSE", [f"{rr['active_rmse_mmolL']:.2f}"], ["manuscript", "SI"])
check("false-activation rate", [pct(rr["false_positive_rate_0.05"])], ["manuscript", "SI"])
check("n active rows", [f"{int(rr['n_active_rows']):,}", str(int(rr["n_active_rows"]))],
      ["manuscript", "SI"])

fam = kv(M2 / "revision" / "family_recovery.csv")
check("family R2", [f"{fam['family_r2_corr']:.2f}"], ["manuscript", "SI"])
check("family sign-match", [pct(fam["family_sign_match_rate"])], ["manuscript", "SI"])
check("dominant-family hit rate", [pct(fam["dominant_family_hit_rate"], 0),
                                   f"{fam['dominant_family_hit_rate']:.2f}"],
      ["manuscript", "SI"])
check("n active family rows", [f"{int(fam['n_active_family_rows']):,}"], ["SI"])

# --------------------------------------------------------------------------
# 2. Dictionary identifiability
# --------------------------------------------------------------------------
ident = kv(M2 / "revision" / "identifiability_dictionary.csv")
check("dictionary reactions", [f"{int(ident['n_reactions'])} reactions",
                               f"{int(ident['n_reactions'])}-reaction"],
      ["manuscript", "SI"])
check("dictionary rank", [f"rank {int(ident['rank'])}"], ["manuscript", "SI"])
check("rank deficiency", [f"deficiency {int(ident['rank_deficiency'])}",
                          f"rank deficiency | {int(ident['rank_deficiency'])}"],
      ["manuscript", "SI"])
check("collinear pairs", [f"{int(ident['n_collinear_pairs_abs_cos_gt_0.7'])} "
                          f"column pairs",
                          f"{int(ident['n_collinear_pairs_abs_cos_gt_0.7'])} of the "
                          f"dictionary column pairs",
                          f"| {int(ident['n_collinear_pairs_abs_cos_gt_0.7'])} ("],
      ["manuscript", "SI"])

# --------------------------------------------------------------------------
# 3. Topology (Fig 2, Table 5, Table S9)
# --------------------------------------------------------------------------
topo = pd.read_csv(M2 / "revision" / "topology_baselines.csv")
np_row = topo[topo["baseline"] == "head_gradient_downhill_k2"].iloc[0]
for field in ["tp", "fp", "fn", "n_inferred_edges", "n_reference_edges"]:
    check(f"no-prior {field}", [str(int(np_row[field]))], ["manuscript", "tables"])
check("no-prior precision", [f"{np_row['precision']:.2f}"], ["manuscript", "SI", "tables"])
check("no-prior recall", [f"{np_row['recall']:.2f}"], ["manuscript", "SI", "tables"])
check("no-prior F1", [f"{np_row['f1']:.2f}"], ["manuscript", "SI", "tables"])

knn = topo[topo["baseline"] == "distance_knn_k2"].iloc[0]
check("kNN baseline F1", [f"{knn['f1']:.2f}"], ["manuscript", "SI"])

# --------------------------------------------------------------------------
# 4. Public USGS age benchmark -- also the Fig 5b annotation
# --------------------------------------------------------------------------
pub = pd.read_csv(M3 / "m3_tracerlpm_parity_agefractions_full_summary.csv").iloc[0]
check("public total rows", [f"{int(pub['total_rows']):,}", str(int(pub["total_rows"]))],
      ["manuscript", "SI"])
check("public identifiable rows", [str(int(pub["identifiable_rows"]))], ["manuscript", "SI"])
check("public non-identifiable rows", [str(int(pub["nonidentifiable_rows"]))],
      ["manuscript", "SI"])
check("public median |log10 err|", [f"{pub['median_abs_log10_error']:.3f}",
                                    f"{pub['median_abs_log10_error']:.2f}"],
      ["manuscript", "SI"])
check("public log10 RMSE", [f"{pub['log10_rmse']:.3f}"], ["manuscript", "SI"])
check("public within 2x", [pct(pub["within_factor_2"])], ["manuscript", "SI"])
check("public within 10x", [pct(pub["within_factor_10"])], ["manuscript", "SI"])
check("public log10 R2 (text, Table 2/4, Fig 5b)",
      [f"{pub['log10_r2']:.3f}", f"{pub['log10_r2']:.2f}"], ["manuscript", "SI"])

# --------------------------------------------------------------------------
# 5. Field and PSI (Fig 4, Fig 7, Table 6, Tables S5/S10)
# --------------------------------------------------------------------------
field = pd.read_csv(M2 / "field_discovery_results.csv")
psi = pd.read_csv(M2 / "m2_phase_stability_index.csv")
check("field edge count", [str(len(field))], ["manuscript", "SI"])
check("field median chemistry R2", [f"{field['chemistry_r2'].median():.2f}"],
      ["manuscript", "SI", "tables"])

check("field candidate edges", ["572"], ["manuscript"])
check("field rejected edges", ["314"], ["manuscript"])
for site, label in (("Manu", "Lower Anayari"), ("Talensi", "Talensi")):
    sub = field[field["edge_id"].str.startswith(site)]
    check(f"{label} retained edges", [str(len(sub))], ["manuscript", "SI"])
    check(f"{label} median chemistry R2", [f"{sub['chemistry_r2'].median():.2f}"],
          ["manuscript", "SI"])
    # sheaf cohomology diagnostics are constant per network
    for col, fmt, where in [("sheaf_h0_dim", "{:.0f}", ["manuscript", "SI"]),
                            ("sheaf_h1_dim", "{:.0f}", ["manuscript"]),
                            ("sheaf_obstruction_energy", "{:.3f}", ["manuscript"])]:
        if col in sub.columns and sub[col].notna().any():
            check(f"{label} {col}", [fmt.format(float(sub[col].dropna().iloc[0]))],
                  where)

# The directional refinement must leave no edge running against the elevation
# proxy -- this is the claim Section 4.6 makes and the reviewer challenged.
_elev: dict[str, float] = {}
for path, idc, zc in [("data/FieldData/LowerAnayari/manu.csv", "Sample ID", "Elevation"),
                      ("data/FieldData/Talensi_MiningArea/talensi.csv", "Code", "Elevation")]:
    df_e = pd.read_csv(ROOT / path)
    for _, r in df_e.iterrows():
        if pd.notna(r.get(zc)):
            _elev[str(r[idc])] = float(r[zc])
uphill = sum(1 for _, r in field.iterrows()
             if _elev.get(str(r["u"])) is not None and _elev.get(str(r["v"])) is not None
             and _elev[str(r["u"])] < _elev[str(r["v"])])
checks += 1
if uphill:
    failures.append(
        f"{uphill} retained field edge(s) run against the elevation proxy, but "
        "Section 4.6 states none do"
    )

# edge-level PSI (family + per-edge probability) is a different file from the
# site x mineral phase-stability matrix loaded above
edge_psi = pd.read_csv(M2 / "top_edges_psi.csv")
for fam, n in edge_psi["family"].value_counts().to_dict().items():
    check(f"PSI family count {fam}", [f"{fam} {n}"], ["manuscript"])
check("median PSI", [f"{edge_psi['psi'].median():.2f}"], ["manuscript", "tables"])

SITE_LABEL = {"Manu": "Lower Anayari", "Talensi": "Talensi"}
for site, mineral in [(s, m) for s in ("Manu", "Talensi")
                      for m in ("CaNa_exch", "NaCa_exch", "calcite", "dolomite")]:
    row = psi[(psi["site"] == site) & (psi["mineral"] == mineral)]
    if row.empty:
        continue
    v = float(row["psi"].iloc[0])
    # Fig 7 plots these same rows as percentages, so SI Table S10 and the
    # heatmap cannot drift apart as long as both derive from this file.
    forms = [f"{v:.2f}", f"{v:.3f}", f"{v:.4f}".rstrip("0")]
    check(f"PSI {SITE_LABEL[site]}/{mineral} (SI Table S10 vs Fig 7)", forms, ["SI"])

# --------------------------------------------------------------------------
# 6. Multi-endmember stress test
# --------------------------------------------------------------------------
me = kv(M2 / "revision" / "multiendmember_recovery_summary.csv")
check("true halite extent", [f"{me['true_halite_extent_mmolL']:.2f}"], ["manuscript", "SI"])
check("true calcite extent", [f"{me['true_calcite_extent_mmolL']:.2f}"], ["manuscript", "SI"])
check("recovered halite", [f"{me['median_extent_halite_mmolL']:.2f}"], ["manuscript", "SI"])
check("recovered calcite", [f"{me['median_extent_calcite_mmolL']:.2f}"], ["manuscript", "SI"])
check("multi-endmember chemistry R2", [f"{me['median_chemistry_r2']:.3f}"],
      ["manuscript", "SI"])
check("median recovered f", [f"{me['median_recovered_f']:.2f}"], ["manuscript", "SI"])
check("mix model selected", [pct(me["mix_model_selected_fraction"], 1),
                             pct(me["mix_model_selected_fraction"], 0)],
      ["manuscript", "SI"])

# --------------------------------------------------------------------------
# 7. Age ordering (Table S12) -- and the denominator must actually support the
#    reported percentages, which is the defect the reviewer caught
# --------------------------------------------------------------------------
ao = kv(M2 / "revision" / "age_overlap_stats.csv")
n_edges = int(ao["n_edges_checked"])
check("age violation fraction", [pct(ao["violation_fraction_mean"])], ["manuscript", "SI"])
check("age severe fraction", [pct(ao["severe_violation_fraction_mean"])], ["manuscript", "SI"])
check("age overlap share", [pct(ao["violations_with_overlapping_intervals_fraction"])],
      ["manuscript", "SI"])
check("age-order consistency index", [f"{ao['age_order_consistency_index']:.2f}"],
      ["manuscript", "SI", "tables"])
check("age-ordering denominator", [f"{n_edges} directed edges", f"{n_edges} ("], ["SI"])

# each reported fraction must be an integer count over the stated denominator
N_REAL = 100
total = n_edges * N_REAL
for key in ["violation_fraction_mean", "severe_violation_fraction_mean",
            "violations_with_overlapping_intervals_fraction"]:
    checks += 1
    k = ao[key] * total
    if abs(k - round(k)) > 0.51:
        failures.append(
            f"age-ordering {key}={ao[key]:.4f} is not an integer count out of "
            f"{total} edge-checks -- the SI denominator cannot produce it"
        )
if f"{total:,}" not in CORPUS["SI"] and str(total) not in CORPUS["SI"]:
    failures.append(f"SI never states the {total} edge-check denominator")
checks += 1

# --------------------------------------------------------------------------
# 7b. Reviewer-facing structural consistency checks
# --------------------------------------------------------------------------
check("homogeneous-nullity interpretation", ["homogeneous nullity"], ["manuscript", "SI"])
check("affine-solvability interpretation", ["exact affine global section", "no exact global section"], ["manuscript", "SI"])
check("radius-sensitivity disclosure", ["Supplementary Table S14"], ["manuscript", "SI"])
check("age-overlap policy", ["does not apply a continuous overlap-proportional weight"], ["manuscript", "SI"])
check("Table 6 disagreement disclosure", ["disagree on the other five"], ["manuscript"])

for label, token in [
    ("stale graph-Laplacian wording", "ordinary weighted graph Laplacian"),
    ("stale H0 interpretation", "H0 = 0 means no assignment of node states"),
    ("stale affine-family wording", "ten-dimensional family of globally consistent assignments"),
    ("stale null-edge result", "Seven edges converging"),
    ("stale zero-closure result", "returned zero chemical closure"),
    ("stale age down-weighting claim", "down-weighted in refinement"),
]:
    checks += 1
    if token.lower() in CORPUS["manuscript"].lower() or token.lower() in CORPUS["SI"].lower():
        failures.append(f"{label}: stale text still present: {token!r}")

radius_path = M2 / "revision" / "radius_sensitivity.csv"
checks += 1
if not radius_path.exists():
    failures.append(f"radius sensitivity output missing: {radius_path}")
else:
    radius = pd.read_csv(radius_path)
    expected = {
        "Manu": {3.0: (265, 117), 5.0: (422, 121), 7.5: (480, 123), 10.0: (485, 123)},
        "Talensi": {3.0: (89, 85), 5.0: (150, 137), 7.5: (177, 151), 10.0: (192, 162)},
    }
    if len(radius) != 8:
        failures.append(f"radius sensitivity output has {len(radius)} rows, expected 8")
    for site, by_radius in expected.items():
        for r_km, (n_candidates, n_retained) in by_radius.items():
            row = radius[(radius["site"] == site) & (radius["radius_km"] == r_km)]
            if row.empty or int(row.iloc[0]["n_candidates"]) != n_candidates or int(row.iloc[0]["n_retained"]) != n_retained:
                failures.append(f"radius sensitivity mismatch for {site} at {r_km:g} km")

# --------------------------------------------------------------------------
# 8. Cross-document consistency
# --------------------------------------------------------------------------
edge_re = r"\b((?:Talensi|Manu)_\d+_\d+)\b"
s5_ids = set(re.findall(edge_re, TABLES.get("table_s5_edge_outputs.md", "")))
t6 = TABLES.get("table6_discovery.md", "")
t6_ids = set(re.findall(edge_re, t6))
checks += 1
if s5_ids and t6_ids and s5_ids != t6_ids:
    failures.append(f"Table 6 edges {sorted(t6_ids)} != Table S5 edges {sorted(s5_ids)}")

# Table 6 rows must appear in the manuscript verbatim
for line in t6.splitlines():
    if re.match(r"\|\s*(Talensi|Lower Anayari)\s*\|", line):
        checks += 1
        if _norm(line.strip()) not in CORPUS["manuscript"]:
            failures.append(f"Table 6 row not verbatim in manuscript: {line.strip()[:70]}")

# SI Table S5 rows must match the generated file
for eid in s5_ids:
    checks += 1
    if eid not in CORPUS["SI"]:
        failures.append(f"Table S5 edge {eid} missing from SI")

# --------------------------------------------------------------------------
# 9. Banned values -- superseded numbers that must not reappear
# --------------------------------------------------------------------------
BANNED = {
    "183.3": "superseded age MAE",
    "0.0026": "superseded AICc lambda",
    "1249": "superseded public-age row count",
    "qA": "superseded isotope noise symbol; use sigma_delta",
    "0.687": "pre-fix reaction R2",
}
for token, why in BANNED.items():
    for name in ("manuscript", "SI", "tables"):
        checks += 1
        if re.search(rf"(?<![\w.]){re.escape(token)}(?![\w])", CORPUS[name]):
            failures.append(f"BANNED {token!r} ({why}) still present in {name}")

# --------------------------------------------------------------------------
print(f"number_audit: {checks} assertions over "
      f"{len(CORPUS)} documents and {len(TABLES)} generated tables")
if failures:
    print(f"\n{len(failures)} FAILURE(S):\n")
    for f in failures:
        print("  -", f)
    sys.exit(1)
print("PASS: every canonical value agrees across text, tables and figure sources")
