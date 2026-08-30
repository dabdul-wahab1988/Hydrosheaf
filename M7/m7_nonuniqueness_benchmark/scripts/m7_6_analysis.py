"""Analysis primitives for the M7.6 controlled-synthetic M3 mechanism audit.

The module is deliberately independent of HydroSheaf's production inference
code.  It supplies a small, auditable fusion layer for the four synthetic
evidence streams and a local transit-time-distribution feasibility diagnostic.
The latter is the claim-bearing part of M7.6: it asks whether declared shared
tracer nuisance or CFC degradation can produce an M3-like, pairwise family
conflict.  It never treats synthetic output as field validation.
"""

from __future__ import annotations

from itertools import combinations
import math
from typing import Dict, Iterable, Mapping, Sequence

import numpy as np
import pandas as pd
from scipy.optimize import linprog
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import (
    average_precision_score,
    brier_score_loss,
    log_loss,
    mean_absolute_error,
    roc_auc_score,
)

from independent_modflow_generator_v2 import (
    CARBON14_HALF_LIFE_YEARS,
    C14_ANALYTICAL_SIGMA_PMC,
    D18O_ANALYTICAL_SIGMA_PERMIL,
    DEUTERIUM_EXCESS_SIGMA_PERMIL,
    HALOCARBON_ANALYTICAL_CV,
    HALOCARBON_GASES,
    HE4_ACCUMULATION_CCPG_PER_YEAR,
    HE4_AIR_SATURATED_CCPG,
    TRITIUM_HALF_LIFE_YEARS,
    atmospheric_mixing_ratio,
    halocarbon_panel,
    helium_panel,
    radiocarbon_pmc,
)

STREAMS = ("H", "C", "N", "E")
TARGETS = ("T1", "T2", "T3")
SUBSETS = tuple(
    "".join(stream for stream in STREAMS if mask & (1 << index))
    for mask in range(1, 1 << len(STREAMS))
    for index in [0]
    if True
)
# The compact expression above is intentionally replaced with lexical order;
# it makes the persisted panel order stable across Python versions.
SUBSETS = tuple(
    "".join(stream for index, stream in enumerate(STREAMS) if mask & (1 << index))
    for mask in range(1, 1 << len(STREAMS))
)

ION_FEATURES = (
    "Ca",
    "Mg",
    "Na",
    "K",
    "HCO3",
    "Cl",
    "SO4",
    "NO3",
    "F",
    "Fe",
    "PO4",
    "SiO2",
)
NUCLEAR_TRACERS = (
    "c14_pmc",
    "cfc11_pptv",
    "cfc12_pptv",
    "cfc113_pptv",
    "sf6_pptv",
    "he4_ccpg",
    "h3_he3_TU",
    "tritium_TU",
    "argon39_pmc",
)
E_TRACERS = ("d18O_permil", "d2H_permil", "sr87_sr86")
REDOX_PROCESSES = frozenset(("sulfate_reduction", "iron_reduction"))
M3_CFC_TRACERS = ("cfc11_pptv", "cfc12_pptv", "cfc113_pptv")
M3_TRITIUM_PAIR = ("tritium_TU", "h3_he3_TU")


def _clip_probability(value: np.ndarray | float) -> np.ndarray:
    return np.clip(np.asarray(value, dtype=float), 1.0e-7, 1.0 - 1.0e-7)


def _site_particle_milestone(site_id: str) -> tuple[str, int]:
    text = str(site_id)
    particle, milestone = text.rsplit("_M", 1)
    return particle, int(milestone)


def _particle_key(site_id: str) -> str:
    return _site_particle_milestone(site_id)[0]


def _normalise(values: Sequence[float]) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    scale = float(np.nanstd(values))
    if not np.isfinite(scale) or scale < 1.0e-9:
        scale = 1.0
    return (values - float(np.nanmean(values))) / scale


def _safe_log(value: float) -> float:
    return float(np.log(max(1.0e-8, float(value))))


def nuclear_age_proxy(row: Mapping[str, object]) -> float:
    """Robust age proxy used only as an evidence feature, never as truth."""

    candidates: list[float] = []
    tritium = max(0.02, float(row["tritium_TU"]))
    candidates.append(max(0.0, -17.0 * _safe_log(tritium / 6.2)))
    argon = max(0.05, float(row["argon39_pmc"]))
    candidates.append(max(0.0, -269.0 / math.log(2.0) * _safe_log(argon / 100.0)))
    c14 = max(0.05, float(row["c14_pmc"]))
    candidates.append(
        max(0.0, -CARBON14_HALF_LIFE_YEARS / math.log(2.0) * _safe_log(c14 / 100.0))
    )
    he4 = float(row["he4_ccpg"])
    candidates.append(
        max(0.0, (he4 - HE4_AIR_SATURATED_CCPG) / HE4_ACCUMULATION_CCPG_PER_YEAR)
    )
    # The 3H/3He partner is an ingrowth clock and saturates at 70 years in the
    # generator; using the closed-form inverse keeps this feature transparent.
    h3 = max(1.0e-8, float(row["tritium_TU"]))
    he3 = max(0.0, float(row["h3_he3_TU"]))
    candidates.append(
        max(0.0, TRITIUM_HALF_LIFE_YEARS / math.log(2.0) * _safe_log(1.0 + he3 / (0.90 * h3)))
    )
    # Atmospheric tracers are converted by a nearest-history search.  Their
    # common-mode temperature nuisance is therefore visible to the diagnostic.
    sample_date = float(row.get("sample_date", 2025.5))
    for gas in HALOCARBON_GASES:
        observed = max(0.0, float(row[gas]))
        ages = np.linspace(0.25, 100.0, 250)
        expected = np.asarray(
            [atmospheric_mixing_ratio(gas, sample_date - age) for age in ages]
        )
        candidates.append(float(ages[int(np.argmin(abs(expected - observed)))]))
    finite = np.asarray([value for value in candidates if np.isfinite(value)])
    return float(np.clip(np.median(finite), 0.0, 200.0))


def _chemistry_score(left: Mapping[str, object], right: Mapping[str, object]) -> float:
    distances = []
    for ion in ION_FEATURES:
        left_value = max(1.0e-8, float(left[ion]))
        right_value = max(1.0e-8, float(right[ion]))
        distances.append(abs(math.log(left_value) - math.log(right_value)))
    distances.extend(
        [
            abs(float(left["pH"]) - float(right["pH"])),
            0.1 * abs(float(left["temp_c"]) - float(right["temp_c"])),
        ]
    )
    return float(-np.mean(distances))


def _isotope_score(left: Mapping[str, object], right: Mapping[str, object]) -> float:
    d18 = (float(left["d18O_permil"]) - float(right["d18O_permil"])) / 0.20
    d2h = (float(left["d2H_permil"]) - float(right["d2H_permil"])) / 2.0
    sr = (float(left["sr87_sr86"]) - float(right["sr87_sr86"])) / 4.0e-4
    return float(-math.sqrt(d18 * d18 + d2h * d2h + sr * sr))


def _particle_records(observations: Sequence[Mapping[str, object]]) -> dict[str, list[Mapping[str, object]]]:
    grouped: dict[str, list[Mapping[str, object]]] = {}
    for row in observations:
        grouped.setdefault(_particle_key(str(row["site_id"])), []).append(row)
    for rows in grouped.values():
        rows.sort(key=lambda row: _site_particle_milestone(str(row["site_id"]))[1])
    return grouped


def build_edge_features(
    observations: Sequence[Mapping[str, object]],
    *,
    true_edges: Sequence[tuple[str, str]] | None = None,
) -> pd.DataFrame:
    """Build one deterministic candidate-edge row per directed within-particle pair."""

    rows = {str(row["site_id"]): row for row in observations}
    true = {f"{u}->{v}" for u, v in (true_edges or ())}
    case_seed = int(str(next(iter(rows))).split("MF", 1)[1].split("_", 1)[0])
    grouped = _particle_records(observations)
    records: list[dict[str, object]] = []
    age_proxy = {site: nuclear_age_proxy(row) for site, row in rows.items()}
    for particle_rows in grouped.values():
        for left in particle_rows:
            for right in particle_rows:
                left_id = str(left["site_id"])
                right_id = str(right["site_id"])
                if left_id == right_id:
                    continue
                left_m = _site_particle_milestone(left_id)[1]
                right_m = _site_particle_milestone(right_id)[1]
                head_drop = float(left["hydraulic_head"]) - float(right["hydraulic_head"])
                distance = max(1.0, math.hypot(float(left["x_m"]) - float(right["x_m"]), float(left["y_m"]) - float(right["y_m"])))
                edge_id = f"{left_id}->{right_id}"
                records.append(
                    {
                        "seed": case_seed,
                        "edge_id": edge_id,
                        "u": left_id,
                        "v": right_id,
                        "milestone_delta": right_m - left_m,
                        "H_score": head_drop / distance,
                        "C_score": _chemistry_score(left, right),
                        "N_score": -abs(age_proxy[left_id] - age_proxy[right_id]) / 20.0,
                        "E_score": _isotope_score(left, right),
                        "is_true_edge": int(edge_id in true),
                    }
                )
    return pd.DataFrame.from_records(records)


def _node_feature_row(row: Mapping[str, object], *, mixing_fraction: float | None = None) -> dict[str, object]:
    age_proxy = nuclear_age_proxy(row)
    ion_total = float(sum(float(row[ion]) for ion in ION_FEATURES))
    return {
        "seed": int(str(row["site_id"]).split("MF", 1)[1].split("_", 1)[0]),
        "site_id": str(row["site_id"]),
        "particle": _particle_key(str(row["site_id"])),
        "milestone": _site_particle_milestone(str(row["site_id"]))[1],
        "H_head": float(row["hydraulic_head"]),
        "H_elevation": float(row["elevation"]),
        "C_log_ion_total": math.log1p(max(0.0, ion_total)),
        "C_pH": float(row["pH"]),
        "C_temp": float(row["temp_c"]),
        "N_age_proxy": age_proxy,
        "N_log_c14": _safe_log(float(row["c14_pmc"])),
        "N_log_cfc11": _safe_log(float(row["cfc11_pptv"])),
        "N_log_cfc12": _safe_log(float(row["cfc12_pptv"])),
        "N_log_cfc113": _safe_log(float(row["cfc113_pptv"])),
        "N_log_sf6": _safe_log(float(row["sf6_pptv"])),
        "N_log_he4": _safe_log(float(row["he4_ccpg"])),
        "N_log_h3he3": _safe_log(float(row["h3_he3_TU"])),
        "N_log_tritium": _safe_log(float(row["tritium_TU"])),
        "N_log_argon39": _safe_log(float(row["argon39_pmc"])),
        "E_d18O": float(row["d18O_permil"]),
        "E_d2H": float(row["d2H_permil"]),
        "E_sr": float(row["sr87_sr86"]),
        "true_mixing_fraction": mixing_fraction,
    }


def build_node_features(
    observations: Sequence[Mapping[str, object]],
    *,
    true_mixing_fraction: Mapping[str, float] | None = None,
    true_ages: Mapping[str, float] | None = None,
) -> pd.DataFrame:
    records = []
    for row in observations:
        site = str(row["site_id"])
        record = _node_feature_row(
            row,
            mixing_fraction=(true_mixing_fraction or {}).get(site),
        )
        if true_ages is not None:
            record["true_age_years"] = float(true_ages[site])
        records.append(record)
    return pd.DataFrame.from_records(records)


TARGET_FEATURES: dict[str, dict[str, tuple[str, ...]]] = {
    "T2": {
        "H": ("H_head", "H_elevation"),
        "C": ("C_log_ion_total", "C_pH", "C_temp"),
        "N": (
            "N_age_proxy",
            "N_log_c14",
            "N_log_cfc11",
            "N_log_cfc12",
            "N_log_cfc113",
            "N_log_sf6",
            "N_log_he4",
            "N_log_h3he3",
            "N_log_tritium",
            "N_log_argon39",
        ),
        # Environmental isotopes are conservative recharge/source tracers.  A
        # target-conditional age mask prevents them from becoming an accidental
        # age predictor while retaining them in the panel-combination matrix.
        "E": (),
    },
    "T3": {
        "H": ("H_head", "H_elevation"),
        "C": ("C_log_ion_total", "C_pH", "C_temp"),
        "N": ("N_log_c14", "N_log_cfc11", "N_log_cfc12", "N_log_cfc113", "N_log_sf6"),
        "E": ("E_d18O", "E_d2H", "E_sr"),
    },
}


def _feature_names(target: str, panel: str) -> list[str]:
    if target == "T1":
        return [f"{stream}_score" for stream in panel]
    return [feature for stream in panel for feature in TARGET_FEATURES[target][stream]]


def fit_models(
    development_nodes: pd.DataFrame,
    development_edges: pd.DataFrame,
) -> dict[str, dict[str, dict[str, object]]]:
    """Fit all development-only panel models for T1--T3."""

    models: dict[str, dict[str, dict[str, object]]] = {target: {} for target in TARGETS}
    for panel in SUBSETS:
        names = _feature_names("T1", panel)
        x = development_edges[names].to_numpy(float)
        y = development_edges["is_true_edge"].to_numpy(int)
        estimator = make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000, random_state=7601))
        estimator.fit(x, y)
        models["T1"][panel] = {"kind": "logistic", "features": names, "estimator": estimator}
        for target, outcome in (("T2", "true_age_years"), ("T3", "true_mixing_fraction")):
            names = _feature_names(target, panel)
            if names:
                estimator = make_pipeline(StandardScaler(), Ridge(alpha=1.0))
                estimator.fit(development_nodes[names].to_numpy(float), development_nodes[outcome].to_numpy(float))
                residual = development_nodes[outcome].to_numpy(float) - estimator.predict(development_nodes[names].to_numpy(float))
                residual_sd = float(np.std(residual, ddof=1)) if len(residual) > 1 else 1.0
                models[target][panel] = {"kind": "ridge", "features": names, "estimator": estimator, "residual_sd": max(1.0e-6, residual_sd)}
            else:
                mean = float(development_nodes[outcome].mean())
                residual_sd = float(development_nodes[outcome].std(ddof=1)) if len(development_nodes) > 1 else 1.0
                models[target][panel] = {"kind": "mean", "features": [], "mean": mean, "residual_sd": max(1.0e-6, residual_sd)}
    return models


def model_summary(models: Mapping[str, Mapping[str, Mapping[str, object]]]) -> dict[str, dict[str, dict[str, object]]]:
    """Return JSON-safe model metadata without serialising sklearn objects."""

    output: dict[str, dict[str, dict[str, object]]] = {}
    for target, panels in models.items():
        output[target] = {}
        for panel, model in panels.items():
            item = {key: value for key, value in model.items() if key != "estimator"}
            estimator = model.get("estimator")
            if estimator is not None:
                regressor = estimator[-1]
                item["coefficients"] = np.asarray(regressor.coef_, dtype=float).reshape(-1).tolist()
                item["intercept"] = float(np.asarray(regressor.intercept_, dtype=float).reshape(-1)[0])
            output[target][panel] = item
    return output


def _permute_stream_columns(frame: pd.DataFrame, panel: str, *, seed: int, salt: int) -> pd.DataFrame:
    output = frame.copy()
    for stream in panel:
        columns = [column for column in output.columns if column.startswith(f"{stream}_")]
        if not columns:
            continue
        rng = np.random.default_rng(int(seed) * 1009 + int(salt) + ord(stream))
        order = rng.permutation(len(output))
        for column in columns:
            output[column] = output[column].to_numpy()[order]
    return output


def predict_panel(
    frame: pd.DataFrame,
    model: Mapping[str, object],
    *,
    target: str,
) -> np.ndarray:
    if model["kind"] == "mean":
        return np.full(len(frame), float(model["mean"]), dtype=float)
    names = list(model["features"])
    return np.asarray(model["estimator"].predict(frame[names].to_numpy(float)), dtype=float)


def edge_metrics(frame: pd.DataFrame, probability: np.ndarray) -> dict[str, float]:
    truth = frame["is_true_edge"].to_numpy(int)
    probability = _clip_probability(probability)
    return {
        "pr_auc": float(average_precision_score(truth, probability)),
        "roc_auc": float(roc_auc_score(truth, probability)) if len(np.unique(truth)) == 2 else float("nan"),
        "brier": float(brier_score_loss(truth, probability)),
        "log_loss": float(log_loss(truth, probability, labels=[0, 1])),
        "mean_edge_entropy": float(np.mean(-(probability * np.log(probability) + (1.0 - probability) * np.log(1.0 - probability)) / math.log(2.0))),
        "overconfident_error_fraction": float(np.mean(((probability <= 0.1) | (probability >= 0.9)) & (truth != (probability >= 0.5)))),
    }


def predict_case_panels(
    nodes: pd.DataFrame,
    edges: pd.DataFrame,
    models: Mapping[str, Mapping[str, Mapping[str, object]]],
    *,
    condition: str = "native",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Predict one case without reading target truth."""

    if condition not in {"native", "permuted", "identity_relation"}:
        raise ValueError(f"Unknown evidence condition {condition!r}")
    edge_rows: list[dict[str, object]] = []
    node_rows: list[dict[str, object]] = []
    seed = int(edges["seed"].iloc[0])
    for panel in SUBSETS:
        edge_input = edges if condition == "native" else _permute_stream_columns(edges, panel, seed=seed, salt=91)
        if condition == "identity_relation":
            probability = np.where(edge_input["milestone_delta"].to_numpy(int) == 1, 0.90, 0.10)
        else:
            probability = predict_panel(edge_input, models["T1"][panel], target="T1")
        for edge_id, value in zip(edges["edge_id"], probability):
            edge_rows.append({"seed": seed, "target": "T1", "condition": condition, "panel": panel, "edge_id": edge_id, "probability": float(value)})
        for target in ("T2", "T3"):
            node_input = nodes if condition == "native" else _permute_stream_columns(nodes, panel, seed=seed, salt=191 + ord(target[1]))
            if condition == "identity_relation":
                if target == "T2":
                    prediction = 1.5 + 15.0 * node_input["milestone"].to_numpy(float)
                else:
                    prediction = np.full(len(node_input), 0.5, dtype=float)
                residual_sd = 1.0
            else:
                model = models[target][panel]
                prediction = predict_panel(node_input, model, target=target)
                residual_sd = float(model["residual_sd"])
            lower = prediction - 1.96 * residual_sd
            upper = prediction + 1.96 * residual_sd
            if target == "T2":
                lower = np.maximum(0.0, lower)
            else:
                lower = np.clip(lower, 0.0, 1.0)
                upper = np.clip(upper, 0.0, 1.0)
            for site_id, pred, low, high in zip(nodes["site_id"], prediction, lower, upper):
                node_rows.append({
                    "seed": seed,
                    "target": target,
                    "condition": condition,
                    "panel": panel,
                    "site_id": site_id,
                    "prediction": float(pred),
                    "lower": float(low),
                    "upper": float(high),
                })
    return pd.DataFrame(edge_rows), pd.DataFrame(node_rows)


def score_case_predictions(
    node_truth: pd.DataFrame,
    edge_truth: pd.DataFrame,
    edge_predictions: pd.DataFrame,
    node_predictions: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Join truth after prediction and calculate case-level metrics."""

    edge_rows: list[dict[str, object]] = []
    node_rows: list[dict[str, object]] = []
    for keys, group in edge_predictions.groupby(["seed", "condition", "panel"], sort=True):
        seed, condition, panel = keys
        truth = edge_truth[edge_truth["seed"] == seed].set_index("edge_id").loc[group["edge_id"], "is_true_edge"].to_numpy(int)
        metric_frame = pd.DataFrame({"is_true_edge": truth})
        edge_rows.append({"seed": int(seed), "target": "T1", "condition": condition, "panel": panel, **edge_metrics(metric_frame, group["probability"].to_numpy(float))})
    for keys, group in node_predictions.groupby(["seed", "target", "condition", "panel"], sort=True):
        seed, target, condition, panel = keys
        truth_column = "true_age_years" if target == "T2" else "true_mixing_fraction"
        truth_map = node_truth[node_truth["seed"] == seed].set_index("site_id")[truth_column]
        truth = truth_map.loc[group["site_id"]].to_numpy(float)
        prediction = group["prediction"].to_numpy(float)
        lower = group["lower"].to_numpy(float)
        upper = group["upper"].to_numpy(float)
        node_rows.append({
            "seed": int(seed),
            "target": target,
            "condition": condition,
            "panel": panel,
            "mae": float(mean_absolute_error(truth, prediction)),
            "mean_interval_width": float(np.mean(upper - lower)),
            "coverage": float(np.mean((truth >= lower) & (truth <= upper))),
            "mean_prediction": float(np.mean(prediction)),
        })
    return pd.DataFrame(edge_rows), pd.DataFrame(node_rows)


def evaluate_case_panels(
    nodes: pd.DataFrame,
    edges: pd.DataFrame,
    models: Mapping[str, Mapping[str, Mapping[str, object]]],
    *,
    condition: str = "native",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Convenience wrapper for development-only evaluation."""

    edge_predictions, node_predictions = predict_case_panels(
        nodes.drop(columns=["true_age_years", "true_mixing_fraction"], errors="ignore"),
        edges.drop(columns=["is_true_edge"], errors="ignore"),
        models,
        condition=condition,
    )
    return score_case_predictions(nodes, edges, edge_predictions, node_predictions)


def bootstrap_case_contrasts(
    metrics: pd.DataFrame,
    *,
    metric_columns: Sequence[str],
    comparisons: Sequence[tuple[str, str, str]],
    n_bootstrap: int = 10_000,
    random_seed: int = 7607,
) -> pd.DataFrame:
    """Paired case-block bootstrap for a complete predeclared contrast family."""

    if metrics.empty:
        return pd.DataFrame()
    rng = np.random.default_rng(int(random_seed))
    rows: list[dict[str, object]] = []
    group_columns = [column for column in ("target", "condition", "nuisance_level") if column in metrics.columns]
    for group_values, group in metrics.groupby(group_columns, sort=True):
        if not isinstance(group_values, tuple):
            group_values = (group_values,)
        selectors = dict(zip(group_columns, group_values))
        for left_panel, right_panel, label in comparisons:
            left = group[group["panel"] == left_panel].set_index("seed")
            right = group[group["panel"] == right_panel].set_index("seed")
            common = sorted(set(left.index) & set(right.index))
            if not common:
                continue
            for metric in metric_columns:
                paired = left.loc[common, metric].to_numpy(float) - right.loc[common, metric].to_numpy(float)
                sampled = rng.choice(paired, size=(int(n_bootstrap), len(paired)), replace=True).mean(axis=1)
                rows.append({
                    **selectors,
                    "contrast": label,
                    "left_panel": left_panel,
                    "right_panel": right_panel,
                    "metric": metric,
                    "mean_difference": float(np.mean(paired)),
                    "ci95_low": float(np.quantile(sampled, 0.025)),
                    "ci95_high": float(np.quantile(sampled, 0.975)),
                    "n_cases": len(common),
                    "n_bootstrap": int(n_bootstrap),
                    "resampling_unit": "independent_MODFLOW_case",
                })
    return pd.DataFrame(rows)


def _expected_nuclear_matrix(
    row: Mapping[str, object],
    ages: np.ndarray,
) -> np.ndarray:
    """Reference tracer responses used by the local feasibility LP."""

    sample_date = float(row.get("sample_date", 2025.5))
    tritium = 6.2 * np.exp(-ages / 17.0)
    argon = 100.0 * np.exp(-np.log(2.0) * ages / 269.0) * (1.0 - 0.05 * (1.0 - np.exp(-ages / 55.0)))
    d13c = float(row.get("d13C_permil", -23.0))
    c14 = np.asarray([radiocarbon_pmc(float(age), d13c)["c14_pmc"] for age in ages])
    cfc = {gas: np.asarray([atmospheric_mixing_ratio(gas, sample_date - age) for age in ages]) for gas in HALOCARBON_GASES}
    he4 = HE4_AIR_SATURATED_CCPG + HE4_ACCUMULATION_CCPG_PER_YEAR * ages
    he3 = np.asarray([helium_panel(float(age), float(t))["h3_he3_TU"] for age, t in zip(ages, tritium)])
    columns = {
        "c14_pmc": c14,
        "cfc11_pptv": cfc["cfc11_pptv"],
        "cfc12_pptv": cfc["cfc12_pptv"],
        "cfc113_pptv": cfc["cfc113_pptv"],
        "sf6_pptv": cfc["sf6_pptv"],
        "he4_ccpg": he4,
        "h3_he3_TU": he3,
        "tritium_TU": tritium,
        "argon39_pmc": argon,
    }
    return np.vstack([columns[tracer] for tracer in NUCLEAR_TRACERS])


def _tracer_tolerance(tracer: str, observed: float, k_sigma: float) -> float:
    base = {
        "c14_pmc": C14_ANALYTICAL_SIGMA_PMC,
        "tritium_TU": 0.10,
        "argon39_pmc": 1.20,
        "h3_he3_TU": max(1.0e-4, HALOCARBON_ANALYTICAL_CV * abs(observed)),
        "he4_ccpg": max(1.0e-10, HALOCARBON_ANALYTICAL_CV * abs(observed)),
        "cfc11_pptv": max(1.0e-3, HALOCARBON_ANALYTICAL_CV * abs(observed)),
        "cfc12_pptv": max(1.0e-3, HALOCARBON_ANALYTICAL_CV * abs(observed)),
        "cfc113_pptv": max(1.0e-3, HALOCARBON_ANALYTICAL_CV * abs(observed)),
        "sf6_pptv": max(1.0e-3, HALOCARBON_ANALYTICAL_CV * abs(observed)),
    }
    return float(k_sigma * base[tracer])


def _feasible_convex_ttd(
    expected: np.ndarray,
    observed: np.ndarray,
    tolerances: np.ndarray,
    indexes: Sequence[int],
) -> bool:
    selected = np.asarray(list(indexes), dtype=int)
    matrix = expected[selected]
    obs = observed[selected]
    tol = tolerances[selected]
    n_ages = matrix.shape[1]
    upper = matrix @ np.ones(n_ages)  # shape check without creating a second matrix
    del upper
    a_ub = np.vstack((matrix, -matrix))
    b_ub = np.concatenate((obs + tol, -(obs - tol)))
    result = linprog(
        np.zeros(n_ages, dtype=float),
        A_ub=a_ub,
        b_ub=b_ub,
        A_eq=np.ones((1, n_ages), dtype=float),
        b_eq=np.ones(1, dtype=float),
        bounds=[(0.0, 1.0)] * n_ages,
        method="highs",
    )
    return bool(result.success)


def ttd_feasibility_diagnostics(
    observations: Sequence[Mapping[str, object]],
    *,
    redox_by_node: Mapping[str, str],
    k_sigma: float = 1.96,
    age_grid: np.ndarray | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return node-level and pair-level feasibility results.

    The pairwise audit is intentionally limited to the nine synthetic nuclear
    tracers so that its geometry is comparable to the M3 local-TTD diagnostic.
    Truth labels are not used in the feasibility LP.
    """

    ages = np.asarray(age_grid if age_grid is not None else np.linspace(0.25, 100.0, 240), dtype=float)
    node_rows: list[dict[str, object]] = []
    pair_rows: list[dict[str, object]] = []
    tracer_index = {tracer: index for index, tracer in enumerate(NUCLEAR_TRACERS)}
    all_indexes = tuple(range(len(NUCLEAR_TRACERS)))
    pairs = list(combinations(NUCLEAR_TRACERS, 2))
    for row in observations:
        site = str(row["site_id"])
        expected = _expected_nuclear_matrix(row, ages)
        observed = np.asarray([float(row[tracer]) for tracer in NUCLEAR_TRACERS], dtype=float)
        tolerances = np.asarray([_tracer_tolerance(tracer, float(row[tracer]), float(k_sigma)) for tracer in NUCLEAR_TRACERS], dtype=float)
        full_feasible = _feasible_convex_ttd(expected, observed, tolerances, all_indexes)
        pair_count = 0
        cfc_pair_count = 0
        tritium_pair_count = 0
        for left, right in pairs:
            indexes = (tracer_index[left], tracer_index[right])
            feasible = _feasible_convex_ttd(expected, observed, tolerances, indexes)
            if not feasible:
                pair_count += 1
                if left in M3_CFC_TRACERS and right in M3_CFC_TRACERS:
                    cfc_pair_count += 1
                if {left, right} == set(M3_TRITIUM_PAIR):
                    tritium_pair_count += 1
            pair_rows.append({
                "site_id": site,
                "redox_class": redox_by_node.get(site, "unknown"),
                "tracer_left": left,
                "tracer_right": right,
                "pair_family": "cfc_cfc" if left in M3_CFC_TRACERS and right in M3_CFC_TRACERS else "tritium_pair" if {left, right} == set(M3_TRITIUM_PAIR) else "other",
                "cfc11_in_pair": int(left == "cfc11_pptv" or right == "cfc11_pptv"),
                "cfc12_in_pair": int(left == "cfc12_pptv" or right == "cfc12_pptv"),
                "pair_infeasible": int(not feasible),
                "k_sigma": float(k_sigma),
            })
        node_rows.append({
            "site_id": site,
            "redox_class": redox_by_node.get(site, "unknown"),
            "full_panel_infeasible": int(not full_feasible),
            "infeasible_pair_count": int(pair_count),
            "cfc_cfc_infeasible_pair_count": int(cfc_pair_count),
            "tritium_pair_infeasible": int(tritium_pair_count > 0),
            "nuclear_tracer_count": len(NUCLEAR_TRACERS),
            "k_sigma": float(k_sigma),
            "age_grid_min_years": float(ages.min()),
            "age_grid_max_years": float(ages.max()),
            "age_grid_n": int(len(ages)),
        })
    return pd.DataFrame(node_rows), pd.DataFrame(pair_rows)


def summarise_m3_mechanism(
    node_diagnostics: pd.DataFrame,
    pair_diagnostics: pd.DataFrame,
    *,
    nuisance_level: str,
) -> pd.DataFrame:
    """Aggregate M3-like signatures without using hidden truth as a predictor."""

    rows: list[dict[str, object]] = []
    for redox_class, group in node_diagnostics.groupby("redox_class", sort=True):
        rows.append({
            "nuisance_level": nuisance_level,
            "redox_class": redox_class,
            "n_nodes": int(len(group)),
            "full_infeasibility_rate": float(group["full_panel_infeasible"].mean()),
            "mean_infeasible_pair_count": float(group["infeasible_pair_count"].mean()),
            "cfc_cfc_infeasibility_rate": float(group["cfc_cfc_infeasible_pair_count"].gt(0).mean()),
            "tritium_pair_infeasibility_rate": float(group["tritium_pair_infeasible"].mean()),
        })
    for family, group in pair_diagnostics.groupby("pair_family", sort=True):
        if family == "other":
            continue
        rows.append({
            "nuisance_level": nuisance_level,
            "redox_class": "ALL",
            "pair_family": family,
            "n_pairs": int(len(group)),
            "pair_infeasibility_rate": float(group["pair_infeasible"].mean()),
        })
    return pd.DataFrame(rows)
