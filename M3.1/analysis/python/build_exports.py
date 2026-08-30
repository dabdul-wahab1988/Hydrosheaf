"""Build tidy, read-only CSV exports for the M3.1 R figure/map pipeline.

Python owns computation; this script performs no new scientific inference. It
copies already-validated `tables/Manuscript_Ready/*.csv` outputs of the
regenerated `M3.1/m3_age_benchmark` pipeline, joins site metadata for the new
study-area map, aggregates the tracer-withholding cross-validation results
into the summary the manuscript reports, and reshapes the identified-TTD
infeasibility-audit JSON into tidy tables. R (`M3.1/analysis/r/`) consumes
only the CSVs written here and recomputes no reported statistic.

Run from the repository root:
    .venv\\Scripts\\python.exe M3.1\\analysis\\python\\build_exports.py
"""
from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT))

BENCH = REPO_ROOT / "M3.1" / "m3_age_benchmark"
RESULTS = BENCH / "results"
TABLES = BENCH / "tables" / "Manuscript_Ready"
DATA_OUT = REPO_ROOT / "M3.1" / "manuscript" / "artifacts" / "data"
DATA_OUT.mkdir(parents=True, exist_ok=True)

SITES_TXT = (
    REPO_ROOT / "M2" / "m2_benchmark" / "external" / "usgs_age" / "input"
    / "DataForNationalGroundwaterAge_1_1" / "Table_1_Sites.txt"
)

CV_TRACERS = ("3H", "SF6", "14C", "CFC11", "CFC12")
CV_UNITS = {"3H": "TU", "SF6": "pptv", "14C": "pmC", "CFC11": "pptv", "CFC12": "pptv"}


def _write(df: pd.DataFrame, name: str) -> None:
    path = DATA_OUT / name
    df.to_csv(path, index=False)
    print(f"wrote {path} ({len(df)} rows)")


def copy_manuscript_ready_tables() -> None:
    """Carry the already-validated, regenerated Manuscript_Ready tables into
    the R-facing evidence interface unchanged (byte-for-byte copy, no
    recomputation)."""
    for src in sorted(TABLES.glob("*.csv")):
        shutil.copyfile(src, DATA_OUT / src.name)
        print(f"copied {src.name}")


def build_site_locations() -> None:
    sites = pd.read_csv(SITES_TXT, sep="\t", low_memory=False, encoding="latin-1")
    sites = sites.rename(columns={
        "SampleID": "site_id", "LatDD83": "lat", "LongDD83": "lon",
    })
    sites = sites[["site_id", "lat", "lon", "State", "StudyUnit", "AqGroup"]].drop_duplicates("site_id")

    design = pd.read_csv(RESULTS / "m3_design_matrix_results.csv", low_memory=False)
    strict = design[design["scenario_id"] == "tracerlpm_strict_parity"]
    reportable_ids = set(strict.loc[strict["est_age_multi"].notna(), "site_id"].astype(str))
    used_ids = set(design["site_id"].astype(str).unique())

    sites["site_id"] = sites["site_id"].astype(str)
    sites = sites[sites["site_id"].isin(used_ids)].copy()
    sites["reportable_strict_parity"] = sites["site_id"].isin(reportable_ids)
    sites = sites.dropna(subset=["lat", "lon"])
    _write(sites, "site_locations.csv")


def build_atmospheric_input_histories() -> None:
    from hydrosheaf.nuclear.input_history import build_default_tritium_input
    from hydrosheaf.nuclear.multi_tracer import build_atmospheric_tracer_input

    rows = []
    trit = build_default_tritium_input()
    for year, value in zip(trit.years, trit.values):
        rows.append({"tracer": "3H", "recharge_year": float(year), "value": float(value), "unit": "TU"})
    for tracer, unit in (("SF6", "pptv"), ("CFC12", "pptv")):
        hist = build_atmospheric_tracer_input(tracer)
        for year, value in zip(hist.years, hist.values):
            rows.append({"tracer": tracer, "recharge_year": float(year), "value": float(value), "unit": unit})
    _write(pd.DataFrame(rows), "atmospheric_input_histories.csv")


def build_strict_parity_scatter() -> None:
    design = pd.read_csv(RESULTS / "m3_design_matrix_results.csv", low_memory=False)
    strict = design[design["scenario_id"] == "tracerlpm_strict_parity"].copy()
    strict = strict[strict["est_age_multi"].notna() & strict["ref_age"].notna()]
    out = strict[["site_id", "StudyUnit", "AqGroup", "ref_age", "est_age_multi", "age_class",
                  "within_factor_2", "within_factor_10"]].copy()
    out["log10_ref_age"] = np.log10(out["ref_age"])
    out["log10_est_age"] = np.log10(out["est_age_multi"])
    _write(out, "strict_parity_scatter.csv")


def build_edge_geometry() -> None:
    src = RESULTS / "m3_real_usgs_graph_edges.csv"
    df = pd.read_csv(src)
    _write(df, "edge_geometry.csv")


def build_tracer_withholding_summary() -> None:
    rows = []
    for tracer in CV_TRACERS:
        path = RESULTS / f"m3_cv_benchmark_{tracer}.csv"
        if not path.exists():
            continue
        df = pd.read_csv(path, low_memory=False)
        eligible_n = len(df)
        finite = df.dropna(subset=["pred_single"])
        reportable_n = len(finite)

        def _rmse(col: str) -> float:
            e = finite["true_val"] - finite[col]
            return float(np.sqrt(np.mean(np.square(e)))) if len(e) else float("nan")

        rows.append({
            "tracer": tracer,
            "unit": CV_UNITS[tracer],
            "eligible_rows": eligible_n,
            "reportable_rows": reportable_n,
            "rmse_single_baseline": _rmse("pred_single"),
            "rmse_hydraulic_graph": _rmse("pred_phys_hyd"),
            "rmse_depth_graph": _rmse("pred_phys_dep"),
            "rmse_randomised_control": _rmse("pred_rand_graph"),
        })
    _write(pd.DataFrame(rows), "tracer_withholding_summary.csv")


def build_infeasibility_exports() -> None:
    path = RESULTS / "m3_infeasibility_diagnostics.json"
    if not path.exists():
        print(f"skip: {path} not found (identified-TTD diagnostic not rerun)")
        return
    diag = json.loads(path.read_text(encoding="utf-8"))
    # Finding 4 of docs/m3_identified_ttd_infeasibility_audit_20260731.md ran
    # the pairwise/minimal-infeasible-subset test at k=6.0 specifically ("for
    # every site with at least three constraints... tested for feasibility at
    # k=6.0"), not at the k=1.96 baseline used for the headline 27.85% rate.
    # Values below are cross-checked against that document's tables (e.g.
    # CFC11+CFC12 32.7%, 14C+3H 0.8%, MIS sizes 19/188/17/2).
    by_k = diag["subset"]["by_k"]["k=6"]

    pairwise = by_k["pairwise"]
    rows = []
    for pair_key, stats in pairwise.items():
        if stats.get("seen", 0) < by_k.get("pairwise_reporting_minimum", 0):
            continue
        a, b = pair_key.split("|") if "|" in pair_key else pair_key.split("+")
        rows.append({
            "tracer_a": a.strip(), "tracer_b": b.strip(),
            "seen": stats["seen"], "infeasible": stats["infeasible"],
            "rate": stats["infeasible"] / stats["seen"] if stats["seen"] else float("nan"),
        })
    _write(pd.DataFrame(rows), "infeasibility_pairwise.csv")

    mis = by_k["minimal_infeasible_subset_size_counts"]
    mis_rows = [{"mis_size": int(k), "count": v} for k, v in mis.items()]
    _write(pd.DataFrame(mis_rows), "infeasibility_mis_size.csv")

    # Conditioning-tracer-count vs. infeasibility rate. Transcribed from the
    # reviewed finding in docs/m3_identified_ttd_infeasibility_audit_20260731.md
    # Finding 1 (denominators cross-checked against
    # m3_identified_ttd_development_full_20260731.csv n_conditioning_tracers
    # value_counts: 594, 1626, 768, 303, 138 - exact match).
    conditioning = pd.DataFrame([
        {"n_conditioning_tracers": 1, "infeasible": 9, "total": 594},
        {"n_conditioning_tracers": 2, "infeasible": 327, "total": 1626},
        {"n_conditioning_tracers": 3, "infeasible": 347, "total": 768},
        {"n_conditioning_tracers": 4, "infeasible": 188, "total": 303},
        {"n_conditioning_tracers": 5, "infeasible": 104, "total": 138},
    ])
    conditioning["rate"] = conditioning["infeasible"] / conditioning["total"]
    _write(conditioning, "infeasibility_by_conditioning_size.csv")


def main() -> int:
    copy_manuscript_ready_tables()
    build_site_locations()
    build_atmospheric_input_histories()
    build_strict_parity_scatter()
    build_edge_geometry()
    build_tracer_withholding_summary()
    build_infeasibility_exports()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
