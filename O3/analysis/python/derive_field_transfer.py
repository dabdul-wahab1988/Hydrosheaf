"""FIG-6: field- and archive-transfer scope across the three components.

M5 field rows are the per-edge median evidence-lifted resolution index (ELRI),
grouped by dataset from the already-published edge-level export; the group
medians reproduce M5's own docs/m5_results_summary.md values exactly
(NorthernGhana 0.072, Talensi 0.147, LowerAnayari 0.072), confirmed before
this script was written. M4 rows are the three public MODFLOW/MODPATH
archives. M3 has no field-transfer row: its only external reference is the
public USGS national release, disclosed explicitly rather than left as a
silent gap in the figure.

Run:  .venv/Scripts/python.exe O3/analysis/python/derive_field_transfer.py
"""

from __future__ import annotations

import pandas as pd

from _common import M4, M5, write


def m5_field() -> pd.DataFrame:
    elri = pd.read_csv(M5 / "results/external_field_evidence_lifted_resolution.csv")
    ghana = pd.read_csv(M5 / "results/ghana_field_hydrosheaf_core_evidence.csv")
    med = elri.groupby("dataset").agg(
        n_edges=("edge_id", "nunique"),
        median_elri=("evidence_lifted_resolution_index", "median"),
    ).reset_index()
    med["dataset"] = med["dataset"].replace({"NorthernGhana.xlsx": "Northern Ghana"})
    med["component"] = "Reaction"
    med["site_type"] = "field chemistry (wet/dry candidate edges)"
    med["score_name"] = "median evidence-lifted resolution index (ELRI)"
    med["score_value"] = med["median_elri"]
    med["claim_scope"] = ("Field plausibility / candidate-class audit, not "
                          "reaction-truth or flow-path validation")
    return med[["component", "dataset", "site_type", "n_edges", "score_name",
                "score_value", "claim_scope"]]


def m4_field() -> pd.DataFrame:
    rows = []
    names = {"tier_1_savage": "Savage (NH)", "tier_2_great_miami": "Great Miami (OH)",
             "tier_3_long_island": "Long Island (NY)"}
    for tier, label in names.items():
        p = M4 / f"results/{tier}_archive_summary.csv"
        if not p.exists():
            continue
        r = pd.read_csv(p).iloc[0]
        active = pd.notna(r.get("n_particles")) and r.get("n_particles", 0) > 0
        rows.append(dict(
            component="Topology", dataset=label,
            site_type="public MODFLOW/MODPATH archive",
            n_edges=int(r["n_endpoint_edges"]) if active and pd.notna(r["n_endpoint_edges"]) else 0,
            score_name="edge F1 (endpoint/pathline self-consistency)",
            score_value=float(r["edge_f1"]) if active and pd.notna(r.get("edge_f1")) else float("nan"),
            claim_scope=("Diagnostic/prior-assisted self-consistency, not "
                        "independent field-truth validation" if active else
                        "Documented fallback stub; no active validation rows")))
    return pd.DataFrame(rows)


def m3_field() -> pd.DataFrame:
    return pd.DataFrame([dict(
        component="Age / residence time", dataset="(none)",
        site_type="no field-transfer benchmark",
        n_edges=0, score_name="n/a", score_value=float("nan"),
        claim_scope="Only external reference is the public USGS national "
                   "groundwater-age release; no independent field-site "
                   "transfer benchmark exists for the age layer")])


def main() -> None:
    df = pd.concat([m5_field(), m4_field(), m3_field()], ignore_index=True)
    write(df, "field_transfer.csv")


if __name__ == "__main__":
    main()
