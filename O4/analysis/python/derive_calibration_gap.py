"""TAB-2/FIG-5 (M8 rows): nominal 95% interval coverage (internal) vs
realised parameter-recovery error against known truth (external), under the
matched calibration model and under an independently generated numerical
truth, plus the kinetic rank-one structural result.

Reads M8's own registered manuscript-artifact CSVs
(`M8/m8_calibration_benchmark/manuscript/artifacts/`) and one locked
provenance result (`RUN-M8-INDEPENDENT-20260728-01/strategy_summary.csv`)
with no recomputation of the underlying calibration.

Run:  .venv/Scripts/python.exe O4/analysis/python/derive_calibration_gap.py
"""

from __future__ import annotations

import pandas as pd

from _common import M8, write


def matched_model_extremes() -> pd.DataFrame:
    """Fixed 16-design sweep: success rate is 1.0 and mean_coverage_95pc
    stays within a narrow band in every design (the internal signal), while
    the realised recovery error spans more than an order of magnitude (the
    external signal), at constant success."""
    df = pd.read_csv(M8 / "manuscript" / "artifacts" / "m8_transport_parameter_summary.csv")
    rows = []
    for parameter in ["dispersivity", "decay"]:
        g = df.loc[df["parameter"] == parameter]
        best = g.loc[g["median_abs_log10_error"].idxmin()]
        worst = g.loc[g["median_abs_log10_error"].idxmax()]
        for tag, row in [("lowest_error_design", best), ("highest_error_design", worst)]:
            rows.append(
                dict(
                    component="M8 calibration",
                    axis_x=f"{parameter}: {tag} ({row['design']})",
                    internal_signal_name="success_rate",
                    internal_signal_value=float(row["success_rate"]),
                    internal_signal_2_name="mean_coverage_95pc",
                    internal_signal_2_value=float(row["mean_coverage_95pc"]),
                    external_signal_name="median_abs_log10_error",
                    external_signal_value=float(row["median_abs_log10_error"]),
                    external_ci_low=float(row["ci95_low"]),
                    external_ci_high=float(row["ci95_high"]),
                    condition_number=float(row["condition"]),
                    n=int(row["n"]),
                    condition="matched_model_fixed_design_sweep",
                    source="manuscript/artifacts/m8_transport_parameter_summary.csv",
                )
            )
    return pd.DataFrame(rows)


def model_form_shift() -> pd.DataFrame:
    """The same 50 d observation under the matched analytical model versus
    an independently generated 240-cell numerical truth: point-error and
    coverage move in different directions for the same parameter."""
    matched = pd.read_csv(M8 / "manuscript" / "artifacts" / "m8_transport_oed_summary_reviewed.csv")
    matched_local = matched.loc[(matched["start_regime"] == "local") & (matched["strategy"] == "A_optimal")]
    indep = pd.read_csv(M8 / "results" / "RUN-M8-INDEPENDENT-20260728-01" / "strategy_summary.csv")
    indep_eopt = indep.loc[indep["strategy"] == "E_optimal"].iloc[0]
    indep_none = indep.loc[indep["strategy"] == "no_new_measurement"].iloc[0]

    rows = []
    for parameter in ["dispersivity", "decay"]:
        m = matched_local.loc[matched_local["parameter"] == parameter].iloc[0]
        rows.append(
            dict(
                component="M8 calibration", axis_x=f"{parameter}: matched model, 50d vs none",
                model_form="matched_analytical", parameter=parameter,
                internal_signal_name="mean_coverage_95pc", internal_signal_value=float(m["mean_coverage_95pc"]),
                external_signal_name="median_abs_log10_error", external_signal_value=float(m["median_abs_log10_error"]),
                n=int(m["n"]), source="manuscript/artifacts/m8_transport_oed_summary_reviewed.csv",
            )
        )
        err_col = f"median_abs_log10_error_{parameter}"
        cov_col = f"coverage_95pc_{parameter}"
        rows.append(
            dict(
                component="M8 calibration", axis_x=f"{parameter}: independent truth, 50d",
                model_form="independent_numerical", parameter=parameter,
                internal_signal_name="coverage_95pc", internal_signal_value=float(indep_eopt[cov_col]),
                external_signal_name="median_abs_log10_error", external_signal_value=float(indep_eopt[err_col]),
                n=int(indep_eopt["n"]), source="results/RUN-M8-INDEPENDENT-20260728-01/strategy_summary.csv",
            )
        )
        rows.append(
            dict(
                component="M8 calibration", axis_x=f"{parameter}: independent truth, no new measurement",
                model_form="independent_numerical_baseline", parameter=parameter,
                internal_signal_name="coverage_95pc", internal_signal_value=float(indep_none[cov_col]),
                external_signal_name="median_abs_log10_error", external_signal_value=float(indep_none[err_col]),
                n=int(indep_none["n"]), source="results/RUN-M8-INDEPENDENT-20260728-01/strategy_summary.csv",
            )
        )
    return pd.DataFrame(rows)


def kinetic_structural_confound() -> pd.DataFrame:
    df = pd.read_csv(M8 / "manuscript" / "artifacts" / "m8_kinetics_structural_summary.csv")
    df["component"] = "M8 calibration"
    df["condition"] = "kinetic_rate_law_structural_confound"
    df["source"] = "manuscript/artifacts/m8_kinetics_structural_summary.csv"
    return df


def frontier_active_learning() -> pd.DataFrame:
    path = M8 / "provenance" / "runs" / "RUN-M8-FRONTIER-AL-20260728-01" / "strategy_summary.csv"
    df = pd.read_csv(path)
    df["component"] = "M8 calibration"
    df["condition"] = "frontier_active_learning_qualification"
    df["source"] = "provenance/runs/RUN-M8-FRONTIER-AL-20260728-01/strategy_summary.csv"
    return df


def main() -> None:
    write(matched_model_extremes(), "calibration_matched_extremes.csv")
    write(model_form_shift(), "calibration_model_form_shift.csv")
    write(kinetic_structural_confound(), "calibration_kinetic_structural.csv")
    write(frontier_active_learning(), "calibration_frontier_active_learning.csv")


if __name__ == "__main__":
    main()
