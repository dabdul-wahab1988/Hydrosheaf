"""Reviewer-requested, deterministic post-hoc checks for M5--M7.

The script never rewrites locked primary outputs.  It reads their canonical
CSV files and writes clearly labelled post-hoc review-sensitivity tables.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def wilson(successes: int, total: int, z: float = 1.959963984540054) -> tuple[float, float]:
    p = successes / total
    denominator = 1.0 + z * z / total
    centre = (p + z * z / (2.0 * total)) / denominator
    half = z * math.sqrt(p * (1.0 - p) / total + z * z / (4.0 * total * total)) / denominator
    return centre - half, centre + half


def run_m5() -> None:
    results = ROOT / "M5" / "m5_inverse_reaction_benchmark" / "results"
    output = results / "m5_review_sensitivity.csv"
    fits = pd.read_csv(results / "benchmark_fits.csv")
    rows: list[dict[str, object]] = []

    bounded = fits[fits["method"] == "bounded_ls"]
    for label, subset in (
        ("all", bounded),
        ("converged", bounded[bounded["converged"]]),
        ("non_converged", bounded[~bounded["converged"]]),
    ):
        rows.extend(
            [
                {"analysis": "bounded_ls", "subset": label, "metric": "phase_f1", "value": subset["phase_f1"].mean(), "n": len(subset)},
                {"analysis": "bounded_ls", "subset": label, "metric": "false_discovery_rate", "value": subset["false_discovery_rate"].mean(), "n": len(subset)},
                {"analysis": "bounded_ls", "subset": label, "metric": "class_f1", "value": subset["class_f1"].mean(), "n": len(subset)},
            ]
        )

    constrained = fits[
        (fits["method"] == "thermo_elastic_net")
        & np.isclose(fits["noise_level"], 0.03)
        & (fits["panel"] == "full_11")
    ]
    low = int((constrained["class_f1"] < 0.80).sum())
    low_ci, high_ci = wilson(low, len(constrained))
    rows.extend(
        [
            {"analysis": "thermodynamic_screening", "subset": "3pct_noise_full_11", "metric": "fraction_class_f1_below_0.80", "value": low / len(constrained), "n": len(constrained)},
            {"analysis": "thermodynamic_screening", "subset": "3pct_noise_full_11", "metric": "wilson_ci95_low", "value": low_ci, "n": len(constrained)},
            {"analysis": "thermodynamic_screening", "subset": "3pct_noise_full_11", "metric": "wilson_ci95_high", "value": high_ci, "n": len(constrained)},
        ]
    )

    scores = pd.read_csv(results / "mechanism_resolution_scores.csv")
    mixed = scores[scores["archetype"] == "mixed"]
    counts = mixed["true_resolution_class"].value_counts()
    rows.extend(
        [
            {"analysis": "mrs_transfer", "subset": "mixed_heldout", "metric": "accuracy", "value": mixed["resolution_class_correct"].mean(), "n": len(mixed)},
            {"analysis": "mrs_transfer", "subset": "mixed_heldout", "metric": "uniform_random_baseline", "value": 0.25, "n": len(mixed)},
            {"analysis": "mrs_transfer", "subset": "mixed_heldout", "metric": "majority_class_baseline", "value": counts.max() / len(mixed), "n": len(mixed)},
        ]
    )

    measurement = pd.read_csv(results / "next_best_measurement.csv")
    top_indexes = measurement.groupby("scenario_id")["measurement_value_score"].idxmax()
    top = measurement.loc[top_indexes]
    alternatives = measurement.drop(top_indexes)
    for label, subset in (("top_ranked", top), ("other_candidates", alternatives)):
        rows.extend(
            [
                {"analysis": "next_best_measurement", "subset": label, "metric": "mean_realised_class_f1_gain", "value": subset["realised_class_f1_gain"].mean(), "n": len(subset)},
                {"analysis": "next_best_measurement", "subset": label, "metric": "mean_support_change", "value": subset["support_change"].mean(), "n": len(subset)},
            ]
        )

    grid = pd.read_csv(results / "hyperparameter_selection.csv")
    for column in ("lambda_l1", "lambda_l2", "support_threshold_mmolL"):
        values = np.sort(grid[column].unique())
        rows.extend(
            [
                {"analysis": "regularisation_grid", "subset": column, "metric": "minimum", "value": values[0], "n": len(values)},
                {"analysis": "regularisation_grid", "subset": column, "metric": "maximum", "value": values[-1], "n": len(values)},
                {"analysis": "regularisation_grid", "subset": column, "metric": "n_values", "value": len(values), "n": len(values)},
            ]
        )

    pd.DataFrame(rows).to_csv(output, index=False)
    print(f"wrote {output}")


def _m6_import():
    scripts = ROOT / "M6" / "m6_field_transfer_benchmark" / "scripts"
    sys.path.insert(0, str(scripts))
    import m6_common  # type: ignore
    return m6_common


def _custom_corroborated(family: str, row: pd.Series, lifters: list[str], profile: dict[str, float]) -> bool:
    available = set(lifters)
    if family == "silicate":
        return "sr_sio2" in available and np.isfinite(row.get("SiO2_mgL", np.nan)) and row["SiO2_mgL"] >= profile["sio2"]
    if family == "carbonate":
        si_ok = "si" in available and np.isfinite(row.get("Calcite_SI", np.nan))
        sr_ok = "sr_sio2" in available and np.isfinite(row.get("Sr_mgL", np.nan)) and row["Sr_mgL"] >= profile["sr"]
        return bool(si_ok or sr_ok)
    if family == "evaporite":
        if "iso" not in available:
            return False
        d18, d2 = row.get("d18O", np.nan), row.get("d2H", np.nan)
        return bool(np.isfinite(d18) and np.isfinite(d2) and (d2 - 8.0 * d18) < profile["d_excess"])
    if family == "trace_mineral":
        return bool(np.isfinite(row.get("F", np.nan)) and row["F"] > profile["f"])
    return True


def run_m6() -> None:
    m6 = _m6_import()
    module = ROOT / "M6" / "m6_field_transfer_benchmark"
    results = module / "results"
    output = results / "m6_review_sensitivity.csv"
    rows: list[dict[str, object]] = []
    ng = m6.load_northern_ghana()
    pairs = m6.seasonal_well_pairs(ng)
    dry_lookup = {well: dry for well, (_wet, dry) in pairs.items()}

    gate = pd.read_csv(results / "m6_field_gate_perwell.csv")
    profiles = {
        "low": {"sio2": 8.0, "sr": 0.08, "f": 0.04, "d_excess": 10.0},
        "predeclared": {"sio2": 10.0, "sr": 0.10, "f": 0.05, "d_excess": 8.0},
        "high": {"sio2": 12.0, "sr": 0.12, "f": 0.06, "d_excess": 6.0},
    }
    gated_families = set(m6.TRACER_GATED)
    for profile_name, profile in profiles.items():
        for tier, subset in gate.groupby("tier", sort=False):
            lifters = list(m6.TIERS[tier]["lifters"])
            revised: list[str] = []
            for record in subset.itertuples(index=False):
                base = record.ungated_class
                family = record.dominant_family
                corroborated = _custom_corroborated(family, dry_lookup[record.well], lifters, profile)
                if family in gated_families and not corroborated:
                    index = max(0, m6.RESOLUTION_ORDER.index(base) - 1)
                    revised.append(m6.RESOLUTION_ORDER[index])
                else:
                    revised.append(base)
            rows.append(
                {
                    "analysis": "combined_gate_thresholds",
                    "setting": profile_name,
                    "scope": tier,
                    "metric": "fraction_non_identifiable",
                    "value": float(np.mean(np.asarray(revised) == "non_identifiable")),
                    "n": len(revised),
                }
            )

    selected = list(sorted(pairs.items()))[:40]
    for tier in ("tier2_fluoride", "tier4_full_metadata"):
        cfg = m6.TIERS[tier]
        for well, (wet, dry) in selected:
            residual = m6.residual_vector(wet, dry, cfg["ions"])
            method = "thermo" if "si" in cfg["lifters"] else "nonneg"
            supports: dict[float, set[str]] = {}
            extents: dict[float, np.ndarray] = {}
            for high in (10.0, 20.0, 50.0):
                fit = m6.m5.fit_inverse(
                    residual,
                    cfg["ions"],
                    method=method,
                    lambda_l1=0.01,
                    lambda_l2=0.001,
                    upstream_si=m6._upstream_si(dry) if method == "thermo" else {},
                    penalty_scales=m6._penalty_for_panel(cfg["ions"], high=high),
                )
                extents[high] = np.asarray(fit["extents"], float)
                supports[high] = m6.m5.support_from_extents(extents[high])
            for high in (10.0, 50.0):
                rows.extend(
                    [
                        {
                            "analysis": "unmeasured_ion_penalty",
                            "setting": f"{high:g}_vs_20",
                            "scope": tier,
                            "metric": "support_jaccard",
                            "value": m6.m5.jaccard(supports[20.0], supports[high]),
                            "n": 1,
                        },
                        {
                            "analysis": "unmeasured_ion_penalty",
                            "setting": f"{high:g}_vs_20",
                            "scope": tier,
                            "metric": "max_absolute_extent_change_mmolL",
                            "value": float(np.max(np.abs(extents[20.0] - extents[high]))),
                            "n": 1,
                        },
                    ]
                )

    rng = np.random.default_rng(9191)
    cfg = m6.TIERS["tier4_full_metadata"]
    for well, (wet, dry) in selected[:20]:
        residual = m6.residual_vector(wet, dry, cfg["ions"])
        penalty = m6._penalty_for_panel(cfg["ions"])
        base_fit = m6.m5.fit_inverse(
            residual,
            cfg["ions"],
            method="thermo",
            lambda_l1=0.01,
            lambda_l2=0.001,
            upstream_si=m6._upstream_si(dry),
            penalty_scales=penalty,
        )
        base_support = m6.m5.support_from_extents(np.asarray(base_fit["extents"], float))
        jaccards: list[float] = []
        for _ in range(200):
            perturbed = residual + rng.normal(0, 0.05 * np.maximum(np.abs(residual), 1.0e-6))
            fit = m6.m5.fit_inverse(
                perturbed,
                cfg["ions"],
                method="thermo",
                lambda_l1=0.01,
                lambda_l2=0.001,
                upstream_si=m6._upstream_si(dry),
                penalty_scales=penalty,
            )
            support = m6.m5.support_from_extents(np.asarray(fit["extents"], float))
            jaccards.append(m6.m5.jaccard(base_support, support))
        standard_deviation = float(np.std(jaccards, ddof=1))
        for count in (8, 48, 200):
            rows.append(
                {
                    "analysis": "bootstrap_stability_precision",
                    "setting": f"B={count}",
                    "scope": well,
                    "metric": "estimated_standard_error",
                    "value": standard_deviation / math.sqrt(count),
                    "n": count,
                }
            )

    frame = pd.DataFrame(rows)
    # Keep unit-level evidence and append compact aggregate rows for manuscript use.
    aggregates = (
        frame[frame["analysis"].isin(["unmeasured_ion_penalty", "bootstrap_stability_precision"])]
        .groupby(["analysis", "setting", "metric"], as_index=False)
        .agg(value=("value", "mean"), n=("n", "sum"))
    )
    frame = pd.concat([frame, aggregates.assign(scope="aggregate")], ignore_index=True)
    frame.to_csv(output, index=False)
    print(f"wrote {output}")


def run_m7() -> None:
    scripts = ROOT / "M7" / "m7_nonuniqueness_benchmark" / "scripts"
    sys.path.insert(0, str(scripts))
    import m7_3_analysis as analysis  # type: ignore
    from sklearn.linear_model import LogisticRegression

    module = ROOT / "M7" / "m7_nonuniqueness_benchmark"
    results = module / "results" / "m7_3_locked"
    output = results / "review_sensitivity.csv"
    development = pd.read_csv(results / "development_edge_features.csv")
    test = pd.read_csv(results / "locked_test_edge_features.csv")
    rows: list[dict[str, object]] = []

    topology = pd.read_csv(results / "topology_age_sensitivity.csv")
    for threshold in (200.0, 400.0, 1000.0):
        for (regime, condition), subset in topology.groupby(["tracer_regime", "graph_condition"]):
            rows.append(
                {
                    "analysis": "ess_threshold",
                    "setting": f"ESS>={threshold:g}",
                    "scope": f"{regime}:{condition}",
                    "metric": "stable_case_fraction",
                    "value": float((subset["importance_ess"] >= threshold).mean()),
                    "n": len(subset),
                }
            )

    def fit_models(c_value: float) -> dict[str, dict[str, object]]:
        truth = development["is_true_edge"].to_numpy(int)
        fitted: dict[str, dict[str, object]] = {}
        for panel, names in analysis.EVIDENCE_FEATURES.items():
            x = development[list(names)].to_numpy(float)
            means = x.mean(axis=0)
            scales = x.std(axis=0)
            scales[scales < 1.0e-8] = 1.0
            estimator = LogisticRegression(C=c_value, max_iter=2000, random_state=7300 + len(panel), solver="lbfgs")
            estimator.fit((x - means) / scales, truth)
            fitted[panel] = {
                "panel": panel,
                "feature_names": list(names),
                "means": means.tolist(),
                "scales": scales.tolist(),
                "coefficients": estimator.coef_[0].tolist(),
                "intercept": float(estimator.intercept_[0]),
                "regularization_C": c_value,
            }
        return fitted

    base_models: dict[str, dict[str, object]] | None = None
    for c_value in (0.1, 0.3, 1.0, 3.0, 10.0):
        models = fit_models(c_value)
        if c_value == 1.0:
            base_models = models
        summary, _case, _conflict = analysis.evaluate_evidence_conditions(test, models)
        native = summary[(summary["condition"] == "native") & (summary["panel"] == "HAC")].iloc[0]
        for metric in ("pr_auc", "brier", "log_loss", "mean_edge_entropy"):
            rows.append(
                {
                    "analysis": "l2_sensitivity",
                    "setting": f"C={c_value:g}",
                    "scope": "native:HAC",
                    "metric": metric,
                    "value": float(native[metric]),
                    "n": len(test),
                }
            )
        rows.append(
            {
                "analysis": "l2_sensitivity",
                "setting": f"C={c_value:g}",
                "scope": "native:HAC",
                "metric": "coefficients_json",
                "value": json.dumps(models["HAC"]["coefficients"]),
                "n": len(test),
            }
        )

    assert base_models is not None
    for condition in analysis.CONDITIONS:
        conditioned = analysis.apply_evidence_condition(test, condition)
        probabilities = {
            panel: analysis.predict_evidence_model(conditioned, model)
            for panel, model in base_models.items()
        }
        univariate = np.column_stack([probabilities["H"], probabilities["A"], probabilities["C"]])
        span = univariate.max(axis=1) - univariate.min(axis=1)
        integrated = probabilities["HAC"]
        wrong = (integrated >= 0.5) != conditioned["is_true_edge"].to_numpy(int)
        for threshold in (0.10, 0.20, 0.30, 0.50):
            conflict = span >= threshold
            rows.extend(
                [
                    {
                        "analysis": "conflict_threshold",
                        "setting": f"span>={threshold:.2f}",
                        "scope": condition,
                        "metric": "conflict_fraction",
                        "value": float(conflict.mean()),
                        "n": len(test),
                    },
                    {
                        "analysis": "conflict_threshold",
                        "setting": f"span>={threshold:.2f}",
                        "scope": condition,
                        "metric": "error_rate_within_conflict" if conflict.any() else "error_rate_within_conflict_unavailable",
                        "value": float(wrong[conflict].mean()) if conflict.any() else np.nan,
                        "n": int(conflict.sum()),
                    },
                ]
            )

    reaction = pd.read_csv(results / "reaction_edge_nonuniqueness.csv")
    by_case = (
        reaction.groupby(["seed", "tier"], as_index=False)
        .agg(accuracy=("modal_family_correct", "mean"))
    )
    tier_names = set(by_case["tier"].astype(str))
    core_name = "core_8" if "core_8" in tier_names else "core"
    enhanced_name = "enhanced_12" if "enhanced_12" in tier_names else "enhanced"
    core = by_case[by_case["tier"] == core_name].set_index("seed")["accuracy"]
    enhanced = by_case[by_case["tier"] == enhanced_name].set_index("seed")["accuracy"]
    common = sorted(set(core.index) & set(enhanced.index))
    paired = enhanced.loc[common].to_numpy(float) - core.loc[common].to_numpy(float)
    rng = np.random.default_rng(7427)
    sampled = rng.choice(paired, size=(10_000, len(paired)), replace=True).mean(axis=1)
    rows.extend(
        [
            {"analysis": "posthoc_modal_accuracy", "setting": "enhanced_minus_core", "scope": "case_block", "metric": "mean_difference", "value": float(paired.mean()), "n": len(paired)},
            {"analysis": "posthoc_modal_accuracy", "setting": "enhanced_minus_core", "scope": "case_block", "metric": "ci95_low", "value": float(np.quantile(sampled, 0.025)), "n": len(paired)},
            {"analysis": "posthoc_modal_accuracy", "setting": "enhanced_minus_core", "scope": "case_block", "metric": "ci95_high", "value": float(np.quantile(sampled, 0.975)), "n": len(paired)},
        ]
    )

    pd.DataFrame(rows).to_csv(output, index=False)
    print(f"wrote {output}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("module", choices=("m5", "m6", "m7", "all"))
    args = parser.parse_args()
    if args.module in {"m5", "all"}:
        run_m5()
    if args.module in {"m6", "all"}:
        run_m6()
    if args.module in {"m7", "all"}:
        run_m7()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
