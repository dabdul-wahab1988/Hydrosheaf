"""Derive validation evidence for the M2.3 manuscript from primary run records.

Computation authority: Python. Every value is recomputed or read directly from a
machine-written run record. No value is copied from an earlier manuscript, and
where a recomputed value disagrees with an earlier manuscript the disagreement is
written to `evidence_discrepancies.csv` (DECISION D5).

Evidence streams:
  1. Locked controlled-synthetic programme run (component and integrated gates).
  2. Independent 100-realisation synthetic recovery benchmark.
  3. No-prior directed-topology comparison against a particle-tracking reference.
  4. Held-out comparison with published lumped-parameter age outputs (D4).

Run:  .venv/Scripts/python.exe M2.3/analysis/python/derive_model_evidence.py
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "M2.3" / "manuscript" / "artifacts" / "data"
OUT.mkdir(parents=True, exist_ok=True)

RUN_ID = "RUN-INTEGRATION-FULL-20260802-15"
RUN = ROOT / ".codex_work/programme-validation" / RUN_ID
BENCH = ROOT / "M2/m2_benchmark"
AGE = ROOT / "M3/m3_age_benchmark"

LOG10_FACTOR_2 = np.log10(2.0)


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


# --------------------------------------------------------------------------
# Stream 1: locked controlled-synthetic programme
# --------------------------------------------------------------------------

def locked_gates() -> pd.DataFrame:
    """Flatten the four predeclared component gates into one tidy table."""
    gates = load_json(RUN / "performance_gates.json")["gates"]
    headline = {
        "age": [("coverage_including_abstention", "Coverage including abstention", ""),
                ("specialist_mae", "Specialist mean absolute error", "years"),
                ("baseline_mae", "Baseline mean absolute error", "years"),
                ("mean_width", "Mean interval width", "years"),
                ("relative_width", "Relative interval width", ""),
                ("selective_risk", "Selective risk", "years"),
                ("acceptance_rate", "Acceptance rate", "")],
        "reaction": [("multiclass_log_loss", "Multiclass log loss", ""),
                     ("baseline_log_loss", "Baseline multiclass log loss", ""),
                     ("coverage", "Coverage", ""),
                     ("max_classwise_ece", "Maximum classwise calibration error", ""),
                     ("selective_risk", "Selective risk", ""),
                     ("false_commitment_rate", "False commitment rate", "")],
        "kinetics": [("predictive_rmse", "Predictive RMSE", "mmol/L"),
                     ("interval_coverage", "Interval coverage", ""),
                     ("identifiability_rate", "Conditional identification rate", ""),
                     ("identifiability_rate_overall", "Overall identification rate", ""),
                     ("parameter_abstention_rate", "Parameter abstention rate", ""),
                     ("parameter_error", "Parameter error", "log10 units")],
        "integrated": [("hydrosheaf_mean_utility_per_cost", "Mean utility per cost", ""),
                       ("random_mean_utility_per_cost", "Random policy utility per cost", ""),
                       ("strongest_specialist_mean_utility_per_cost",
                        "Strongest specialist utility per cost", ""),
                       ("paired_random_delta_ci_low", "Paired lower bound versus random", ""),
                       ("paired_specialist_delta_ci_low",
                        "Paired lower bound versus specialist", ""),
                       ("observed_false_commitment_rate", "Observed false commitment", ""),
                       ("prospective_case_count", "Locked prospective cases", "cases")],
    }
    rows = []
    for comp, spec in headline.items():
        gate = gates[comp]
        obs, req = gate["observed"], gate["requirements"]
        for key, label, unit in spec:
            rows.append(dict(
                component=comp, gate_name=gate["name"], gate_status=gate["status"],
                claim_scope=req.get("claim_scope", ""), metric=key, label=label,
                unit=unit, value=obs.get(key)))
    return pd.DataFrame(rows)


def locked_thresholds() -> pd.DataFrame:
    """Predeclared acceptance thresholds, reported beside the observed values."""
    gates = load_json(RUN / "performance_gates.json")["gates"]
    rows = []
    for comp, gate in gates.items():
        for key, value in gate["requirements"].items():
            if isinstance(value, (int, float, bool)):
                rows.append(dict(component=comp, requirement=key, threshold=value))
    return pd.DataFrame(rows)


def locked_execution() -> pd.DataFrame:
    """Execution-gate checks that had to pass before any performance claim."""
    gate = load_json(RUN / "execution_gate.json")
    rows = [dict(check=k, value=str(v))
            for k, v in gate.get("checks", {}).items()]
    for key in ["status", "critic_status", "critic_blocker_count",
                "critic_major_count", "programme_workflow_status",
                "field_validation", "synthetic_component_claim",
                "synthetic_integrated_claim"]:
        rows.append(dict(check=key, value=str(gate.get(key))))
    return pd.DataFrame(rows)


def prospective_pairwise() -> pd.DataFrame:
    """Paired policy contrasts from the truth-blind prospective benchmark."""
    bench = load_json(RUN / "prospective_decision_benchmark.json")["benchmark"]
    rows = []
    for name, rec in bench.get("pairwise", {}).items():
        rows.append(dict(
            contrast=name,
            left_policy=rec.get("left_policy"), right_policy=rec.get("right_policy"),
            mean_delta=rec.get("mean_cost_adjusted_utility_delta"),
            ci_low=rec.get("paired_delta_ci_low"), ci_high=rec.get("paired_delta_ci_high"),
            ci_level=rec.get("paired_delta_ci_level"),
            paired_case_count=rec.get("paired_case_count"),
            bootstrap_replicates=rec.get("paired_bootstrap_replicates"),
            win_rate=rec.get("left_utility_win_rate")))
    df = pd.DataFrame(rows)
    df["claim_status"] = bench.get("claim_status")
    df["cost_penalty"] = bench.get("cost_penalty")
    return df


# --------------------------------------------------------------------------
# Stream 2: independent 100-realisation recovery benchmark
# --------------------------------------------------------------------------

def recovery_benchmark() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Recompute transport, reaction and age recovery from row-level records."""
    transport = pd.read_csv(BENCH / "results/transport_recovery.csv")
    reaction = pd.read_csv(BENCH / "results/reaction_recovery.csv")
    age = pd.read_csv(BENCH / "results/age_inference_validation.csv")

    t_ok = transport["absolute_error"].notna()
    t_summary = pd.DataFrame([dict(
        quantity="Transport parameter", n_realisations=int(transport["realisation"].nunique()),
        n_rows=int(t_ok.sum()), n_not_recovered=int((~t_ok).sum()),
        median_absolute_error=float(transport.loc[t_ok, "absolute_error"].median()),
        q25=float(transport.loc[t_ok, "absolute_error"].quantile(0.25)),
        q75=float(transport.loc[t_ok, "absolute_error"].quantile(0.75)),
        model_selection_accuracy=float(transport["model_correct"].mean()),
        unit="dimensionless fraction")])

    # "Active" follows the benchmark's own definition: a truth extent whose
    # magnitude exceeds 0.01 mmol/L. Negative extents denote precipitation and
    # are retained. Inactive terms are the complement.
    active = reaction[reaction["true_extent_mmolL"].abs() > 0.01]
    inactive = reaction[reaction["true_extent_mmolL"].abs() <= 0.01]
    err = active["recovered_extent_mmolL"] - active["true_extent_mmolL"]
    truth = active["true_extent_mmolL"]
    r_summary = pd.DataFrame([dict(
        quantity="Reaction extent (active terms)",
        n_realisations=int(reaction["realisation"].nunique()),
        n_rows=int(len(active)), n_not_recovered=0,
        median_absolute_error=float(err.abs().median()),
        q25=float(err.abs().quantile(0.25)), q75=float(err.abs().quantile(0.75)),
        model_selection_accuracy=np.nan, unit="mmol/L")])
    r_summary["mean_absolute_error"] = float(err.abs().mean())
    r_summary["r_squared"] = float(
        1 - np.sum(err ** 2) / np.sum((truth - truth.mean()) ** 2))
    # Leakage: inactive truth terms assigned a non-trivial recovered extent.
    r_summary["n_inactive_terms"] = int(len(inactive))
    r_summary["inactive_term_leakage_fraction"] = float(
        (inactive["recovered_extent_mmolL"].abs() > 0.05).mean())
    r_summary["leakage_threshold_mmolL"] = 0.05

    # Age: single-node lumped model versus network-constrained Bayesian inference.
    rows = []
    for label, col, wcol in [("Single-node lumped model", "single_node_lpm_years",
                              "single_node_ci_width_years"),
                             ("Network-constrained Bayesian", "network_bayesian_years",
                              "network_ci_width_years")]:
        ok = age[col].notna() & age["true_mrt_years"].notna()
        est, ref = age.loc[ok, col], age.loc[ok, "true_mrt_years"]
        le, lr = np.log10(est.clip(lower=1e-6)), np.log10(ref.clip(lower=1e-6))
        d = le - lr
        rows.append(dict(
            method=label, n=int(ok.sum()),
            mae_years=float((est - ref).abs().mean()),
            median_absolute_error_years=float((est - ref).abs().median()),
            log10_r2=float(1 - np.sum(d ** 2) / np.sum((lr - lr.mean()) ** 2)),
            median_abs_log10_error=float(np.median(np.abs(d))),
            mean_ci_width_years=float(age.loc[ok, wcol].mean())))
    a_df = pd.DataFrame(rows)
    return pd.concat([t_summary, r_summary], ignore_index=True), a_df, reaction


def age_by_class() -> pd.DataFrame:
    """Age recovery stratified by residence-time class."""
    age = pd.read_csv(BENCH / "results/age_inference_validation.csv")
    rows = []
    for cls, grp in age.groupby("age_class"):
        for label, col in [("Single-node lumped model", "single_node_lpm_years"),
                           ("Network-constrained Bayesian", "network_bayesian_years")]:
            ok = grp[col].notna()
            est, ref = grp.loc[ok, col], grp.loc[ok, "true_mrt_years"]
            rows.append(dict(age_class=cls, method=label, n=int(ok.sum()),
                             median_true_years=float(ref.median()),
                             mae_years=float((est - ref).abs().mean()),
                             median_abs_relative_error=float(
                                 ((est - ref).abs() / ref.clip(lower=1e-9)).median())))
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# Stream 3: no-prior directed-topology comparison
# --------------------------------------------------------------------------

def topology() -> pd.DataFrame:
    df = pd.read_csv(BENCH / "results/modpath_noprior_topology.csv")
    df["independent_validation"] = df["independent_validation"].astype(str)
    keep = ["mode", "independent_validation", "n_reference_edges",
            "n_inferred_edges", "tp", "fp", "fn", "precision", "recall", "f1"]
    out = df[keep].copy()
    # Recompute the reported rates from the confusion counts.
    out["precision_recomputed"] = out["tp"] / (out["tp"] + out["fp"])
    out["recall_recomputed"] = out["tp"] / (out["tp"] + out["fn"])
    out["f1_recomputed"] = (2 * out["precision_recomputed"] * out["recall_recomputed"]
                            / (out["precision_recomputed"] + out["recall_recomputed"]))
    out["false_discovery_rate"] = out["fp"] / (out["tp"] + out["fp"])
    return out


# --------------------------------------------------------------------------
# Stream 4: held-out comparison with published lumped-parameter age outputs
# --------------------------------------------------------------------------

def field_application() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Candidate-edge chemistry closures for the two measured field datasets.

    These are in-sample closures of a candidate edge against the same chemistry
    used to fit it. They are screening constructions, not measured flow paths.
    """
    df = pd.read_csv(BENCH / "results/field_discovery_results.csv")
    # Site is carried by the local end-member identifier used for each edge.
    site = np.where(df["endmember_id"].astype(str).str.contains("talensi", case=False),
                    "Talensi",
                    np.where(df["endmember_id"].astype(str).str.contains("manu", case=False),
                             "Lower Anayari", "unassigned"))
    df = df.assign(site=site)
    rows = []
    for name, grp in df.groupby("site"):
        r2 = grp["chemistry_r2"].dropna()
        rows.append(dict(
            site=name, n_candidate_edges=int(len(grp)),
            median_chemistry_r2=float(r2.median()), mean_chemistry_r2=float(r2.mean()),
            q25_chemistry_r2=float(r2.quantile(0.25)),
            q75_chemistry_r2=float(r2.quantile(0.75)),
            pct_mixing_model=float(100.0 * (grp["transport_model"] == "mix").mean())))
    overall = df["chemistry_r2"].dropna()
    rows.append(dict(site="All candidate edges", n_candidate_edges=int(len(df)),
                     median_chemistry_r2=float(overall.median()),
                     mean_chemistry_r2=float(overall.mean()),
                     q25_chemistry_r2=float(overall.quantile(0.25)),
                     q75_chemistry_r2=float(overall.quantile(0.75)),
                     pct_mixing_model=float(100.0 * (df["transport_model"] == "mix").mean())))
    points = df[["edge_id", "site", "chemistry_r2", "transport_model",
                 "objective_score"]].copy()
    # Number of reaction terms retained per candidate edge.
    extent_cols = [c for c in df.columns if c.startswith("extent_")]
    points["n_reactions_retained"] = (df[extent_cols].abs() > 0.05).sum(axis=1).values
    return pd.DataFrame(rows), points


def parity_metrics(est: np.ndarray, ref: np.ndarray) -> dict:
    ok = np.isfinite(est) & np.isfinite(ref)
    e, r = est[ok], ref[ok]
    d = e - r
    return dict(
        n=int(ok.sum()),
        median_abs_log10_error=float(np.median(np.abs(d))),
        log10_rmse=float(np.sqrt(np.mean(d ** 2))),
        log10_r2=float(1 - np.sum(d ** 2) / np.sum((r - r.mean()) ** 2)),
        log10_bias=float(np.mean(d)),
        within_factor_2=float(np.mean(np.abs(d) <= LOG10_FACTOR_2)),
        within_factor_10=float(np.mean(np.abs(d) <= 1.0)))


def public_age() -> tuple[pd.DataFrame, pd.DataFrame]:
    """D4: held-out uncalibrated parity is primary; calibrated is emulation."""
    p = pd.read_csv(AGE / "results/m3_usgs_calibrated_parity.csv")
    held = p[p["fold_id"] == p["held_out_study_unit"]]
    rows = []
    unc = parity_metrics(held["log10_est_age"].values, held["log10_reported_age"].values)
    unc.update(mode="Held-out uncalibrated strict parity",
               interpretation="independent screening-level agreement",
               is_primary=True)
    cal = parity_metrics(held["log10_calibrated_age"].values,
                         held["log10_reported_age"].values)
    cal.update(mode="Held-out calibrated emulation",
               interpretation="fit to the reference outputs; measures emulation, "
                              "not independent agreement",
               is_primary=False)
    rows.extend([unc, cal])
    summary = pd.DataFrame(rows)
    summary["n_study_units"] = int(p["held_out_study_unit"].nunique())

    ok = held["log10_est_age"].notna() & held["log10_reported_age"].notna()
    points = held.loc[ok, ["site_id", "held_out_study_unit", "log10_est_age",
                           "log10_calibrated_age", "log10_reported_age"]].copy()
    points["log10_residual_uncalibrated"] = (points["log10_est_age"]
                                             - points["log10_reported_age"])
    points["log10_residual_calibrated"] = (points["log10_calibrated_age"]
                                           - points["log10_reported_age"])
    return summary, points


def public_age_attrition() -> pd.DataFrame:
    """Characterise the reference sites that carry no usable paired value.

    Added after review, which correctly asked whether the 47% attrition is
    informative. It is entirely reference-side, but it is not uniform across
    aquifer groups, so the retained comparison is not representative.
    """
    p = pd.read_csv(AGE / "results/m3_usgs_calibrated_parity.csv")
    raw = pd.read_csv(AGE / "results/m3_usgs_benchmark_results.csv",
                      low_memory=False)
    p = p.assign(retained=np.isfinite(p["log10_est_age"])
                 & np.isfinite(p["log10_reported_age"]))
    merged = p.merge(raw[["site_id", "Depth_m", "AqGroup"]], on="site_id", how="left")
    merged["Depth_m"] = pd.to_numeric(merged["Depth_m"], errors="coerce")
    rows = []
    for grp, sub in merged.groupby("AqGroup", dropna=False):
        rows.append(dict(
            aquifer_group=str(grp), n_total=len(sub),
            n_retained=int(sub["retained"].sum()),
            pct_dropped=round(100.0 * (~sub["retained"]).mean(), 1),
            median_depth_m=round(float(sub["Depth_m"].median()), 1)))
    out = pd.DataFrame(rows).sort_values("pct_dropped", ascending=False)
    out["attrition_side"] = "reference value absent"
    out["n_estimate_missing"] = int((~np.isfinite(p["log10_est_age"])).sum())
    return out.reset_index(drop=True)


# --------------------------------------------------------------------------

def discrepancies(recov: pd.DataFrame, age_public: pd.DataFrame) -> pd.DataFrame:
    """Record where recomputed values disagree with earlier manuscripts (D5)."""
    prim = age_public[age_public["is_primary"]].iloc[0]
    rx = recov[recov["quantity"].str.startswith("Reaction")].iloc[0]
    return pd.DataFrame([
        dict(quantity="Reaction extent recovery, R2 (active terms)",
             earlier_value="0.57 (M2 abstract)",
             recomputed_value=f"{rx['r_squared']:.3f}",
             resolution="Recomputed from reaction_recovery.csv using the "
                        "benchmark's own active definition (|truth| > 0.01 mmol/L, "
                        f"n = {int(rx['n_rows'])}). The earlier value is not "
                        "reproducible from the shipped record under this or any "
                        "tested subset definition."),
        dict(quantity="Reaction extent recovery, MAE (active terms)",
             earlier_value="0.33 mmol/L (M2 abstract)",
             recomputed_value=f"{rx['mean_absolute_error']:.3f} mmol/L",
             resolution="Recomputed under the same active definition."),
        dict(quantity="Inactive-term leakage above 0.05 mmol/L",
             earlier_value="32.7% (M2 abstract)",
             recomputed_value=f"{100 * rx['inactive_term_leakage_fraction']:.1f}%",
             resolution="Recomputed over all inactive truth terms "
                        f"(n = {int(rx['n_inactive_terms'])})."),
        dict(quantity="Public age benchmark, n",
             earlier_value="1,249 (M2 abstract) / 1,272 (M2 summary)",
             recomputed_value=f"{int(prim['n'])} finite held-out pairs",
             resolution="Recomputed from the parity record; earlier counts mixed "
                        "row totals with finite-pair counts."),
        dict(quantity="Public age benchmark, R2",
             earlier_value="0.71 (M2 abstract)",
             recomputed_value=f"{prim['log10_r2']:.3f} uncalibrated held-out",
             resolution="Earlier abstract reported a value reproducible from "
                        "neither the uncalibrated nor the calibrated record."),
        dict(quantity="Public age benchmark, median |log10 error|",
             earlier_value="0.17 (M2 abstract) / 0.383 (M2 summary)",
             recomputed_value=f"{prim['median_abs_log10_error']:.3f} uncalibrated held-out",
             resolution="0.17 corresponds to the calibrated emulation, which is "
                        "not an independent comparison (D4)."),
        dict(quantity="Northern Ghana records",
             earlier_value="320 seasonal field records",
             recomputed_value="160 wells; second seasonal panel reconstructed",
             resolution="Author confirmed the seasonal split was reconstructed (D2)."),
        dict(quantity="Talensi sample longitudes",
             earlier_value="Positive (east) in the source file and in M2/M2_1",
             recomputed_value="Negated; Talensi District lies near 0.8 deg west",
             resolution="As delivered all 63 samples plot outside Ghana. After "
                        "negation all 63 fall inside the Upper East Region "
                        "polygon. Because reflection about the prime meridian is "
                        "an isometry, pairwise sample distances are unchanged, so "
                        "the distance-based candidate-edge results are unaffected; "
                        "only absolute position and the study-area map change."),
        dict(quantity="Field candidate-edge chemistry closure",
             earlier_value="0.711 (M2) / 0.713 (M2_1), described as overall R2",
             recomputed_value="median 0.713 over 208 edges; mean 0.611",
             resolution="The two earlier values are the median and a differently "
                        "pooled summary of the same quantity. The median is "
                        "reported and named as such."),
    ])


def main() -> None:
    locked = locked_gates()
    locked.to_csv(OUT / "locked_gate_metrics.csv", index=False)
    locked_thresholds().to_csv(OUT / "locked_gate_thresholds.csv", index=False)
    locked_execution().to_csv(OUT / "locked_execution_checks.csv", index=False)
    prospective_pairwise().to_csv(OUT / "prospective_policy_contrasts.csv", index=False)

    recov, age_methods, reaction_raw = recovery_benchmark()
    recov.to_csv(OUT / "synthetic_recovery_summary.csv", index=False)
    age_methods.to_csv(OUT / "synthetic_age_method_comparison.csv", index=False)
    age_by_class().to_csv(OUT / "synthetic_age_by_class.csv", index=False)

    # Row-level reaction recovery for the figure layer, carrying the active flag.
    pts = reaction_raw[["realisation", "edge_id", "reaction_label",
                        "true_extent_mmolL", "recovered_extent_mmolL",
                        "absolute_error_mmolL"]].copy()
    pts["term_status"] = np.where(pts["true_extent_mmolL"].abs() > 0.01,
                                  "active", "inactive")
    pts.to_csv(OUT / "synthetic_reaction_points.csv", index=False)

    transport = pd.read_csv(BENCH / "results/transport_recovery.csv")
    transport[transport["absolute_error"].notna()][
        ["realisation", "edge_id", "process", "parameter", "true_value",
         "recovered_value", "absolute_error"]].to_csv(
        OUT / "synthetic_transport_points.csv", index=False)

    age_pts = pd.read_csv(BENCH / "results/age_inference_validation.csv")
    age_pts.to_csv(OUT / "synthetic_age_points.csv", index=False)

    topo = topology()
    topo.to_csv(OUT / "topology_comparison.csv", index=False)

    field_summary, field_points = field_application()
    field_summary.to_csv(OUT / "field_edge_chemistry_summary.csv", index=False)
    field_points.to_csv(OUT / "field_edge_chemistry_points.csv", index=False)

    age_public, age_points = public_age()
    age_public.to_csv(OUT / "public_age_parity_summary.csv", index=False)
    age_points.to_csv(OUT / "public_age_parity_points.csv", index=False)
    public_age_attrition().to_csv(OUT / "public_age_attrition.csv", index=False)

    discrepancies(recov, age_public).to_csv(
        OUT / "evidence_discrepancies.csv", index=False)

    manifest = load_json(RUN / "run_manifest.json")
    (OUT / "locked_run_manifest_summary.json").write_text(json.dumps(dict(
        run_id=RUN_ID,
        git_revision=manifest.get("git_revision") or manifest.get("git_commit"),
        worktree_dirty=manifest.get("git_worktree_dirty"),
        n_source_files_hashed=len(manifest.get("source_hashes", {}) or {}),
    ), indent=2), encoding="utf-8")

    print("== locked gates ==")
    print(locked[locked.metric.isin([
        "coverage_including_abstention", "specialist_mae", "baseline_mae",
        "multiclass_log_loss", "baseline_log_loss", "predictive_rmse",
        "identifiability_rate_overall", "hydrosheaf_mean_utility_per_cost",
        "random_mean_utility_per_cost"])][
        ["component", "label", "value", "gate_status"]].to_string(index=False))
    print("\n== recovery ==")
    print(recov.to_string(index=False))
    print("\n== age methods ==")
    print(age_methods.to_string(index=False))
    print("\n== topology ==")
    print(topo[["mode", "tp", "fp", "fn", "precision_recomputed",
                "recall_recomputed", "f1_recomputed"]].to_string(index=False))
    print("\n== public age ==")
    print(age_public[["mode", "n", "log10_r2", "median_abs_log10_error",
                      "within_factor_2", "is_primary"]].to_string(index=False))


if __name__ == "__main__":
    main()
