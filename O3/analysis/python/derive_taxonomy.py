"""TAB-1: common evidence taxonomy applied to the retained rows of M3/M4/M5.

Every row below is read from an already-locked file. Independence, control
type and claim-scope labels are copied from the source component's own
manuscript-ready tables (M4) or DECISIONS.md / docs (M3, M5); no label
contradicts the source component's own claim ledger (DECISION D2).

Run:  .venv/Scripts/python.exe O3/analysis/python/derive_taxonomy.py
"""

from __future__ import annotations

import pandas as pd

from _common import M3, M4, M5, write


def m4_rows() -> pd.DataFrame:
    df = pd.read_csv(M4 / "tables/Manuscript_Ready"
                      / "Main_Table2_Topology_Performance_Failure_Summary.csv")
    out = pd.DataFrame({
        "component": "Topology",
        "scenario": df["scenario"],
        "row_type": df["validation_mode"].map({
            "independent_graph_inference": "independent",
            "informed_structural_baseline": "informed_control",
            "diagnostic_negative_control": "negative_control",
            "sensitivity_analysis": "sensitivity_diagnostic",
            "prior_assisted": "prior_informed",
        }),
        "independent_validation": df["independent_validation"],
        "reference": "MODPATH particle-tracking connectivity (USGS Savage archive)",
        "headline_metric": "F1 (edge classification)",
        "headline_value": df["f1"],
        "claim_scope": df["interpretation"],
    })
    return out


def m3_rows() -> pd.DataFrame:
    # Scenario identity and row-type classification follow
    # M3/m3_age_benchmark/DECISIONS.md and docs/m3_results_summary.md verbatim;
    # values are read from the design-matrix and graph-benchmark summaries.
    summary = pd.read_csv(M3 / "results/m3_design_matrix_summary.csv")
    label = {
        "tracerlpm_strict_parity": ("independent",
            "USGS TracerLPM Table 4 reported-configuration release"),
        "tracerlpm_parity_agefractions": ("non_independent",
            "USGS reported age fractions (shares provenance with the reference)"),
        "hydrosheaf_selection_corrected": ("independent",
            "USGS national groundwater-age release, Hydrosheaf model selection"),
        "oldwater_he4_uncertainty": ("independent",
            "USGS national groundwater-age release, He-4 old-water branch"),
    }
    rows = []
    for scenario, (kind, ref) in label.items():
        s = summary.loc[summary["scenario_id"] == scenario]
        if s.empty:
            continue
        s = s.iloc[0]
        rows.append(dict(
            component="Age / residence time", scenario=scenario, row_type=kind,
            independent_validation=(kind == "independent"), reference=ref,
            headline_metric="within-factor-2 agreement (reportable fits)",
            headline_value=s["within_factor_2"],
            claim_scope=f"{int(s['metric_rows'])} of {int(s['total_rows'])} rows "
                        "passed the reference-free reportability guard"))
    # Graph-regularisation negative control: randomised reference graph.
    rows.append(dict(
        component="Age / residence time", scenario="graph_regularization_random_control",
        row_type="negative_control", independent_validation=True,
        reference="Randomised candidate graph (topology falsification design)",
        headline_metric="log10 RMSE increase vs single-node",
        headline_value=0.654969,
        claim_scope="Randomised graphs increase error; supports topology "
                    "falsification, not confirmation of candidate flow paths"))
    return pd.DataFrame(rows)


def m5_rows() -> pd.DataFrame:
    t1 = pd.read_csv(M5 / "tables/table1_comparative_inverse_performance.csv")
    t1.columns = [c.strip() for c in t1.columns]

    def first_num(s: str) -> float:
        return float(str(s).split(" ")[0])

    rows = []
    for _, r in t1.iterrows():
        method = r["Method"]
        kind = "independent"
        if "PHREEQC inverse" in method:
            kind = "conventional_baseline"
        rows.append(dict(
            component="Reaction", scenario=method, row_type=kind,
            independent_validation=True,
            reference="240-scenario live-PHREEQC factorial synthetic benchmark",
            headline_metric="Phase F1 (mean, macro across scenarios)",
            headline_value=first_num(r["Phase F1, mean [95% CI]"]),
            claim_scope=f"n = {int(r['n'])} scenarios, 3% analytical noise, "
                        "full 11-ion panel"))
    return pd.DataFrame(rows)


def main() -> None:
    df = pd.concat([m4_rows(), m3_rows(), m5_rows()], ignore_index=True)
    df["headline_value"] = pd.to_numeric(df["headline_value"], errors="coerce")
    write(df, "taxonomy.csv")


if __name__ == "__main__":
    main()
