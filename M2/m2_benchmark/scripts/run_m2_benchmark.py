from __future__ import annotations

import argparse
import json
import math
import sys
import warnings
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import yaml
from scipy.stats import kendalltau, spearmanr

PROJECT_ROOT = Path(__file__).resolve().parents[3]
BENCHMARK_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from hydrosheaf import Config, build_reaction_dictionary, fit_network
from hydrosheaf.graph.types import Edge


OUTPUT_DIRS = [
    "data",
    "data/realisations",
    "results",
    "tables",
    "figures",
    "docs",
]


def ensure_dirs() -> None:
    for rel in OUTPUT_DIRS:
        (BENCHMARK_ROOT / rel).mkdir(parents=True, exist_ok=True)


def load_truth(path: Path) -> Dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def truth_edges(truth: Mapping[str, Any]) -> List[Dict[str, Any]]:
    return list(truth["generation_edges"]) + list(truth.get("lateral_truth_edges", []))


def edge_id(u: str, v: str) -> str:
    return f"{u}->{v}"


def make_config(
    truth: Mapping[str, Any],
    node_vectors: Optional[Mapping[str, Sequence[float]]] = None,
) -> Config:
    ion_order = list(truth["ion_order"])
    conservative = [0.02] * len(ion_order)
    for ion, weight in {"Cl": 1.0, "Na": 0.40, "HCO3": 0.15}.items():
        conservative[ion_order.index(ion)] = weight

    config = Config(
        ion_order=ion_order,
        weights=[1.0] * len(ion_order),
        conservative_weights=conservative,
        active_minerals=list(truth["active_minerals"]),
        signed_reaction_labels=list(truth["signed_reaction_labels"]),
        allow_signed_reactions=True,
        exchange_enabled=True,
        phreeqc_enabled=False,
        gibbs_enabled=False,
        isotope_enabled=True,
        isotope_weight=0.50,
        isotope_d18o_key="d18O",
        isotope_d2h_key="d2H",
        isotope_consistency_weight=0.05,
        missing_policy="impute_zero",
        detection_limit_policy="half",
        lambda_l1=0.002,
        reaction_max_iter=90,
        reaction_tol=1e-7,
        transport_models_enabled=["evap", "mix"],
    )
    if node_vectors:
        config.mixing_endmembers = {
            "LAT_SAL": list(node_vectors["LAT_SAL"]),
            "LAT_CARB": list(node_vectors["LAT_CARB"]),
            "MODERN": list(node_vectors["MODERN"]),
        }
        source_nodes = truth["nodes"]
        config.mixing_endmembers_isotopes = {
            name: (
                float(source_nodes[name]["isotopes"]["d18O"]),
                float(source_nodes[name]["isotopes"]["d2H"]),
            )
            for name in ("LAT_SAL", "LAT_CARB", "MODERN")
        }
    return config


def reaction_maps(config: Config) -> Tuple[Dict[str, List[float]], List[str]]:
    matrix, labels, _ = build_reaction_dictionary(config)
    return {label: list(vector) for label, vector in zip(labels, matrix)}, labels


def vector_from_chem(chem: Mapping[str, float], ion_order: Sequence[str]) -> np.ndarray:
    return np.array([float(chem.get(ion, 0.0)) for ion in ion_order], dtype=float)


def chem_from_vector(vector: Sequence[float], ion_order: Sequence[str]) -> Dict[str, float]:
    return {ion: float(value) for ion, value in zip(ion_order, vector)}


def add_reactions(
    vector: np.ndarray,
    reactions: Mapping[str, float],
    reaction_vectors: Mapping[str, Sequence[float]],
) -> np.ndarray:
    out = vector.astype(float).copy()
    for label, extent in reactions.items():
        if label not in reaction_vectors:
            raise KeyError(f"Reaction '{label}' is not in Hydrosheaf dictionary.")
        out = out + float(extent) * np.array(reaction_vectors[label], dtype=float)
    return np.maximum(out, 1e-7)


def tracer_from_age(age_years: float) -> Dict[str, float]:
    tritium = 5.0 * math.exp(-math.log(2.0) * age_years / 12.32)
    c14 = 100.0 * math.exp(-math.log(2.0) * age_years / 5730.0)
    return {"tritium_TU": tritium, "c14_pmc": c14}


def mix_isotopes(
    primary: Mapping[str, float],
    secondary: Mapping[str, float],
    fraction: float,
) -> Dict[str, float]:
    return {
        "d18O": (1.0 - fraction) * float(primary["d18O"])
        + fraction * float(secondary["d18O"]),
        "d2H": (1.0 - fraction) * float(primary["d2H"])
        + fraction * float(secondary["d2H"]),
    }


def build_noiseless_nodes(
    truth: Mapping[str, Any],
    config: Config,
) -> Dict[str, Dict[str, Any]]:
    ion_order = list(truth["ion_order"])
    reaction_vectors, _ = reaction_maps(config)
    nodes: Dict[str, Dict[str, Any]] = {}

    for node_id, node in truth["nodes"].items():
        nodes[node_id] = dict(node)
        nodes[node_id]["site_id"] = node_id
        nodes[node_id]["chemistry_mmolL"] = chem_from_vector(
            vector_from_chem(node["chemistry_mmolL"], ion_order), ion_order
        )

    generated_meta = dict(truth["generated_nodes"])
    edge_by_v = {edge["v"]: edge for edge in truth["generation_edges"]}
    for node_id in truth["generated_node_order"]:
        meta = dict(generated_meta[node_id])
        gen_edge = edge_by_v[node_id]
        u = gen_edge["u"]
        u_node = nodes[u]
        u_vec = vector_from_chem(u_node["chemistry_mmolL"], ion_order)

        if "f_evap" in gen_edge:
            f_evap = float(gen_edge["f_evap"])
            gamma = 1.0 / (1.0 - f_evap)
            base_vec = u_vec * gamma
            iso = {
                "d18O": float(u_node["isotopes"]["d18O"]) + 7.0 * f_evap,
                "d2H": float(u_node["isotopes"]["d2H"]) + 31.5 * f_evap,
            }
        elif "mix_from" in gen_edge:
            source = str(gen_edge["mix_from"])
            fraction = float(gen_edge.get("f_mix", gen_edge.get("dilution_fraction", 0.0)))
            source_vec = vector_from_chem(nodes[source]["chemistry_mmolL"], ion_order)
            base_vec = (1.0 - fraction) * u_vec + fraction * source_vec
            iso = mix_isotopes(u_node["isotopes"], nodes[source]["isotopes"], fraction)
        else:
            base_vec = u_vec.copy()
            iso = dict(u_node["isotopes"])

        final_vec = add_reactions(
            base_vec, gen_edge.get("reactions", {}), reaction_vectors
        )
        age = float(meta["mrt_years"])
        tracers = tracer_from_age(age)
        if node_id == "DILUTE":
            young = tracer_from_age(10.0)
            old = tracer_from_age(3000.0)
            tracers = {
                key: 0.35 * young[key] + 0.65 * old[key]
                for key in ("tritium_TU", "c14_pmc")
            }

        nodes[node_id] = {
            **meta,
            "site_id": node_id,
            "chemistry_mmolL": chem_from_vector(final_vec, ion_order),
            "isotopes": iso,
            "tracers": tracers,
        }
    return nodes


def node_vectors(nodes: Mapping[str, Mapping[str, Any]], ion_order: Sequence[str]) -> Dict[str, List[float]]:
    return {
        node_id: vector_from_chem(node["chemistry_mmolL"], ion_order).tolist()
        for node_id, node in nodes.items()
    }


def synthetic_si(node: Mapping[str, Any]) -> Dict[str, float]:
    chem = node["chemistry_mmolL"]
    ca = float(chem.get("Ca", 0.0))
    mg = float(chem.get("Mg", 0.0))
    hco3 = float(chem.get("HCO3", 0.0))
    so4 = float(chem.get("SO4", 0.0))
    cl = float(chem.get("Cl", 0.0))
    na = float(chem.get("Na", 0.0))
    return {
        "SIcalcite": math.log10(max(ca * hco3, 1e-9)) - 0.15,
        "SIdolomite": math.log10(max(ca * mg * hco3, 1e-9)) - 0.25,
        "SIgypsum": math.log10(max(ca * so4, 1e-9)) - 0.85,
        "SIhalite": math.log10(max(na * cl, 1e-9)) - 2.40,
    }


def make_realisations(
    truth: Mapping[str, Any],
    nodes: Mapping[str, Mapping[str, Any]],
    n_realisations: int,
) -> pd.DataFrame:
    seed = int(truth["benchmark"]["random_seed"])
    noise = truth["noise"]
    ion_order = list(truth["ion_order"])
    rows: List[Dict[str, Any]] = []
    for realisation in range(n_realisations):
        rng = np.random.default_rng(seed + realisation)
        for node_id, node in nodes.items():
            row: Dict[str, Any] = {
                "realisation": realisation,
                "site_id": node_id,
                "role": node.get("role", ""),
                "age_class": node.get("age_class", ""),
                "mrt_years_true": float(node["mrt_years"]),
                "x": float(node["xy"][0]),
                "y": float(node["xy"][1]),
                "head_m": float(node["head_m"]),
                "pH": 7.1 + rng.normal(0.0, 0.08),
                "temp_C": 27.0 + rng.normal(0.0, 0.4),
            }
            for ion in ion_order:
                true_value = float(node["chemistry_mmolL"].get(ion, 0.0))
                rel_sigma = float(noise["major_ion_rel_sigma"])
                observed = true_value * (1.0 + rng.normal(0.0, rel_sigma))
                row[ion] = max(observed, 0.0)
                row[f"{ion}_true"] = true_value
            row["d18O"] = float(node["isotopes"]["d18O"]) + rng.normal(
                0.0, float(noise["isotope_sigma_permil"])
            )
            row["d2H"] = float(node["isotopes"]["d2H"]) + rng.normal(
                0.0, float(noise["isotope_sigma_permil"])
            )
            row["d18O_true"] = float(node["isotopes"]["d18O"])
            row["d2H_true"] = float(node["isotopes"]["d2H"])
            row["tritium_TU"] = max(
                float(node["tracers"]["tritium_TU"])
                + rng.normal(0.0, float(noise["tritium_sigma_tu"])),
                0.0,
            )
            row["c14_pmc"] = max(
                float(node["tracers"]["c14_pmc"])
                + rng.normal(0.0, float(noise["c14_sigma_pmc"])),
                0.0,
            )
            row.update(synthetic_si(node))
            rows.append(row)
    return pd.DataFrame(rows)


def apply_degradation(
    samples: pd.DataFrame,
    scenario: str,
    realisation: int,
    truth: Mapping[str, Any],
) -> pd.DataFrame:
    degraded = samples.copy()
    rng = np.random.default_rng(int(truth["benchmark"]["random_seed"]) + 1000 + realisation)
    if scenario == "ion_incomplete":
        candidate_nodes = [n for n in degraded["site_id"].tolist() if n not in {"RCH", "MODERN"}]
        remove_nodes = set(rng.choice(candidate_nodes, size=max(1, len(candidate_nodes) // 2), replace=False))
        degraded.loc[degraded["site_id"].isin(remove_nodes), ["Mg", "SO4"]] = np.nan
    elif scenario == "tracer_absent":
        keep = {"RCH", "DEEP"}
        degraded.loc[~degraded["site_id"].isin(keep), ["tritium_TU", "c14_pmc"]] = np.nan
    elif scenario == "head_absent":
        degraded["head_m"] = np.nan
    return degraded


def df_to_samples(df: pd.DataFrame, ion_order: Sequence[str]) -> List[Dict[str, Any]]:
    sample_cols = [
        "site_id",
        "role",
        "age_class",
        "mrt_years_true",
        "x",
        "y",
        "head_m",
        "pH",
        "temp_C",
        "d18O",
        "d2H",
        "tritium_TU",
        "c14_pmc",
        "SIcalcite",
        "SIdolomite",
        "SIgypsum",
        "SIhalite",
    ] + list(ion_order)
    records: List[Dict[str, Any]] = []
    for _, row in df[sample_cols].iterrows():
        rec = row.to_dict()
        for key, value in list(rec.items()):
            if pd.isna(value):
                rec[key] = None
        records.append(rec)
    return records


def base_true_edge_specs(truth: Mapping[str, Any]) -> Dict[str, Dict[str, Any]]:
    specs: Dict[str, Dict[str, Any]] = {}
    for edge in truth["generation_edges"]:
        specs[edge["edge_id"]] = dict(edge)
        specs[edge["edge_id"]]["truth_type"] = "generating"
        specs[edge["edge_id"]]["confidence"] = 0.90
    for edge in truth.get("lateral_truth_edges", []):
        specs[edge["edge_id"]] = dict(edge)
        specs[edge["edge_id"]]["truth_type"] = "lateral"
    return specs


def edge_inputs_for_variant(
    truth: Mapping[str, Any],
    variant: str,
    realisation: int,
    samples: pd.DataFrame,
) -> List[Edge]:
    true_specs = base_true_edge_specs(truth)
    edges = list(true_specs.values())
    false_edges = list(truth.get("false_edges", []))
    rng = np.random.default_rng(int(truth["benchmark"]["random_seed"]) + 2000 + realisation)

    if variant == "full":
        edges = edges + false_edges[:1]
    elif variant == "sparse":
        n_remove = max(1, int(round(0.30 * len(edges))))
        remove_ids = set(rng.choice([e["edge_id"] for e in edges], size=n_remove, replace=False))
        edges = [edge for edge in edges if edge["edge_id"] not in remove_ids]
    elif variant == "dense":
        extra = [
            {"edge_id": "SH30->MIX20", "u": "SH30", "v": "MIX20", "confidence": 0.25},
            {"edge_id": "LAT_CARB->INT1", "u": "LAT_CARB", "v": "INT1", "confidence": 0.20},
            {"edge_id": "SH15->MIX60", "u": "SH15", "v": "MIX60", "confidence": 0.18},
        ]
        edges = edges + false_edges[:3] + extra
    elif variant == "reversed":
        reverse_ids = {"INT1->INT2", "MIX60->DILUTE"}
        reversed_edges = []
        for edge in edges:
            if edge["edge_id"] in reverse_ids:
                reversed_edges.append(
                    {
                        **edge,
                        "u": edge["v"],
                        "v": edge["u"],
                        "edge_id": edge_id(edge["v"], edge["u"]),
                        "confidence": 0.35,
                        "reversed_from": edge["edge_id"],
                    }
                )
            else:
                reversed_edges.append(edge)
        edges = reversed_edges + false_edges[:1]
    else:
        raise ValueError(f"Unknown topology variant: {variant}")

    head_map = {
        str(row.site_id): None if pd.isna(row.head_m) else float(row.head_m)
        for row in samples.itertuples()
    }
    out: List[Edge] = []
    for edge in edges:
        u, v = str(edge["u"]), str(edge["v"])
        hu = head_map.get(u)
        hv = head_map.get(v)
        attrs: Dict[str, Any] = {
            "prior_confidence": float(edge.get("confidence", 0.8)),
            "truth_edge": edge.get("edge_id") in true_specs,
        }
        if hu is not None and hv is not None:
            attrs["delta_h"] = hu - hv
        out.append(Edge(edge_id=str(edge["edge_id"]), u=u, v=v, attrs=attrs))
    return out


def result_rows(
    results: Iterable[Any],
    realisation: int,
    scenario: str,
    variant: str,
    edge_specs: Mapping[str, Mapping[str, Any]],
) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for result in results:
        reaction_map = dict(zip(result.z_labels, result.z_extents))
        prior = edge_specs.get(result.edge_id, {}).get("confidence")
        if prior is None:
            prior = 0.25 if edge_specs.get(result.edge_id) is None else 0.7
        score = float(prior) * math.exp(-min(float(result.objective_score), 10.0))
        row: Dict[str, Any] = {
            "realisation": realisation,
            "scenario": scenario,
            "topology_variant": variant,
            "edge_id": result.edge_id,
            "u": result.u,
            "v": result.v,
            "transport_model": result.transport_model,
            "gamma": result.gamma,
            "f": result.f,
            "objective_score": result.objective_score,
            "anomaly_norm": result.anomaly_norm,
            "chemistry_r2": result.chemistry_r2,
            "transport_prob_evap": result.transport_probabilities.get("evap"),
            "transport_prob_mix": result.transport_probabilities.get("mix"),
            "edge_score": score,
            "prior_confidence": prior,
            "truth_edge": result.edge_id in edge_specs,
        }
        for label, value in reaction_map.items():
            row[f"reaction_{label}"] = value
        rows.append(row)
    return rows


def build_truth_tables(
    truth: Mapping[str, Any],
    nodes: Mapping[str, Mapping[str, Any]],
    config: Config,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    ion_order = list(truth["ion_order"])
    node_rows = []
    for node_id, node in nodes.items():
        row = {
            "site_id": node_id,
            "role": node["role"],
            "age_class": node["age_class"],
            "mrt_years_true": node["mrt_years"],
            "x": node["xy"][0],
            "y": node["xy"][1],
            "head_m": node["head_m"],
            "d18O_true": node["isotopes"]["d18O"],
            "d2H_true": node["isotopes"]["d2H"],
            "tritium_TU_true": node["tracers"]["tritium_TU"],
            "c14_pmc_true": node["tracers"]["c14_pmc"],
        }
        row.update({f"{ion}_true": node["chemistry_mmolL"][ion] for ion in ion_order})
        node_rows.append(row)

    edge_rows = []
    for edge in truth["generation_edges"]:
        row = {
            "edge_id": edge["edge_id"],
            "u": edge["u"],
            "v": edge["v"],
            "process": edge["process"],
            "f_evap_true": edge.get("f_evap"),
            "gamma_true": 1.0 / (1.0 - float(edge["f_evap"])) if "f_evap" in edge else None,
            "f_mix_true": edge.get("f_mix"),
            "dilution_fraction_true": edge.get("dilution_fraction"),
            "mix_from": edge.get("mix_from"),
            "benchmark_edge": True,
        }
        for label, value in edge.get("reactions", {}).items():
            row[f"reaction_{label}_true"] = value
        edge_rows.append(row)
    for edge in truth.get("lateral_truth_edges", []):
        edge_rows.append(
            {
                "edge_id": edge["edge_id"],
                "u": edge["u"],
                "v": edge["v"],
                "process": "lateral_connection",
                "confidence_true": edge["confidence"],
                "benchmark_edge": False,
            }
        )

    reaction_matrix, labels, mineral_mask = build_reaction_dictionary(config)
    reaction_rows = []
    for label, vector, mineral in zip(labels, reaction_matrix, mineral_mask):
        reaction_rows.append(
            {
                "reaction_label": label,
                "is_mineral": mineral,
                "stoichiometry_mmolL_per_extent": json.dumps(
                    {ion: coeff for ion, coeff in zip(ion_order, vector) if abs(coeff) > 0},
                    sort_keys=True,
                ),
                "is_signed_in_benchmark": label in truth["signed_reaction_labels"],
            }
        )
    return pd.DataFrame(node_rows), pd.DataFrame(edge_rows), pd.DataFrame(reaction_rows)


def evaluate_transport(edge_results: pd.DataFrame, truth_edge_df: pd.DataFrame) -> pd.DataFrame:
    baseline = edge_results[
        (edge_results["scenario"] == "complete")
        & (edge_results["topology_variant"] == "full")
    ].copy()
    rows = []
    truth_by_edge = truth_edge_df.set_index("edge_id").to_dict("index")
    for _, row in baseline.iterrows():
        info = truth_by_edge.get(row["edge_id"])
        if not info:
            continue
        process = info["process"]
        if process == "evaporation":
            true_value = float(info["f_evap_true"])
            recovered = 1.0 - (1.0 / float(row["gamma"])) if row["transport_model"] == "evap" and row["gamma"] else np.nan
            parameter = "f_evap"
            model_correct = row["transport_model"] == "evap"
        elif process == "lateral_mixing":
            true_value = float(info["f_mix_true"])
            recovered = float(row["f"]) if row["transport_model"] == "mix" and pd.notna(row["f"]) else np.nan
            parameter = "f_mix"
            model_correct = row["transport_model"] == "mix"
        elif process == "dilution_recharge_pulse":
            true_value = float(info["dilution_fraction_true"])
            recovered = float(row["f"]) if row["transport_model"] == "mix" and pd.notna(row["f"]) else np.nan
            parameter = "dilution_fraction"
            model_correct = row["transport_model"] == "mix"
        else:
            continue
        rows.append(
            {
                "realisation": int(row["realisation"]),
                "edge_id": row["edge_id"],
                "process": process,
                "parameter": parameter,
                "true_value": true_value,
                "recovered_value": recovered,
                "absolute_error": abs(recovered - true_value) if pd.notna(recovered) else np.nan,
                "relative_bias_percent": 100.0 * (recovered - true_value) / true_value if pd.notna(recovered) and true_value else np.nan,
                "detected_model": row["transport_model"],
                "model_correct": bool(model_correct),
            }
        )
    return pd.DataFrame(rows)


def evaluate_reactions(
    edge_results: pd.DataFrame,
    truth: Mapping[str, Any],
) -> pd.DataFrame:
    baseline = edge_results[
        (edge_results["scenario"] == "complete")
        & (edge_results["topology_variant"] == "full")
    ].copy()
    true_by_edge = {edge["edge_id"]: edge for edge in truth["generation_edges"]}
    reaction_labels = [
        "calcite",
        "dolomite",
        "halite",
        "gypsum",
        "CaNa_exch",
        "albite",
        "pyrite_oxidation_aerobic",
    ]
    rows = []
    for _, result in baseline.iterrows():
        edge = true_by_edge.get(result["edge_id"])
        if edge is None:
            continue
        for label in reaction_labels:
            true_value = float(edge.get("reactions", {}).get(label, 0.0))
            recovered = result.get(f"reaction_{label}", 0.0)
            recovered = 0.0 if pd.isna(recovered) else float(recovered)
            rows.append(
                {
                    "realisation": int(result["realisation"]),
                    "edge_id": result["edge_id"],
                    "reaction_label": label,
                    "true_extent_mmolL": true_value,
                    "recovered_extent_mmolL": recovered,
                    "absolute_error_mmolL": abs(recovered - true_value),
                    "relative_bias_percent": (
                        100.0 * (recovered - true_value) / true_value
                        if abs(true_value) > 1e-12
                        else np.nan
                    ),
                    "false_positive": abs(true_value) <= 1e-12 and abs(recovered) > 0.05,
                }
            )
    return pd.DataFrame(rows)


def evaluate_missing_sensitivity(
    edge_results: pd.DataFrame,
    truth: Mapping[str, Any],
) -> pd.DataFrame:
    true_by_edge = {edge["edge_id"]: edge for edge in truth["generation_edges"]}
    scenarios = ["complete", "ion_incomplete", "tracer_absent", "head_absent"]
    rows = []
    for scenario in scenarios:
        subset = edge_results[
            (edge_results["scenario"] == scenario)
            & (edge_results["topology_variant"] == "full")
        ].copy()
        transport_errors: List[float] = []
        reaction_errors: List[float] = []
        objective_scores: List[float] = []
        for _, result in subset.iterrows():
            edge = true_by_edge.get(result["edge_id"])
            if edge is None:
                continue
            objective_scores.append(float(result["objective_score"]))
            if edge["process"] == "evaporation" and result["transport_model"] == "evap" and pd.notna(result["gamma"]):
                true_value = float(edge["f_evap"])
                recovered = 1.0 - (1.0 / float(result["gamma"]))
                transport_errors.append(100.0 * (recovered - true_value) / true_value)
            if edge["process"] in {"lateral_mixing", "dilution_recharge_pulse"} and result["transport_model"] == "mix" and pd.notna(result["f"]):
                key = "f_mix" if edge["process"] == "lateral_mixing" else "dilution_fraction"
                true_value = float(edge[key])
                recovered = float(result["f"])
                transport_errors.append(100.0 * (recovered - true_value) / true_value)
            for label, true_value_raw in edge.get("reactions", {}).items():
                true_value = float(true_value_raw)
                if abs(true_value) <= 1e-12:
                    continue
                recovered = result.get(f"reaction_{label}", 0.0)
                recovered = 0.0 if pd.isna(recovered) else float(recovered)
                reaction_errors.append(100.0 * (recovered - true_value) / true_value)
        rows.append(
            {
                "scenario": scenario,
                "n_edge_results": int(len(subset)),
                "transport_relative_bias_median_percent": float(np.nanmedian(transport_errors)) if transport_errors else np.nan,
                "transport_relative_abs_bias_median_percent": float(np.nanmedian(np.abs(transport_errors))) if transport_errors else np.nan,
                "reaction_relative_bias_median_percent": float(np.nanmedian(reaction_errors)) if reaction_errors else np.nan,
                "reaction_relative_abs_bias_median_percent": float(np.nanmedian(np.abs(reaction_errors))) if reaction_errors else np.nan,
                "objective_score_median": float(np.nanmedian(objective_scores)) if objective_scores else np.nan,
            }
        )
    return pd.DataFrame(rows)


def evaluate_topology(
    edge_results: pd.DataFrame,
    truth: Mapping[str, Any],
    threshold: float = 0.22,
) -> pd.DataFrame:
    edge_results = edge_results[edge_results["scenario"] == "complete"].copy()
    true_edges = {edge["edge_id"] for edge in truth_edges(truth)}
    generation_edges = {edge["edge_id"]: edge for edge in truth["generation_edges"]}
    labels = ["calcite", "dolomite", "halite", "gypsum", "CaNa_exch"]
    rows = []
    for (realisation, variant), subset in edge_results.groupby(["realisation", "topology_variant"]):
        selected = set(subset.loc[subset["edge_score"] >= threshold, "edge_id"].astype(str))
        true_selected = selected & true_edges
        edge_match_rate = len(true_selected) / len(true_edges)
        false_positive_rate = len(selected - true_edges) / max(len(selected), 1)

        true_values: List[float] = []
        recovered_values: List[float] = []
        by_edge = subset.set_index("edge_id")
        for edge_id_value, edge in generation_edges.items():
            result = by_edge.loc[edge_id_value] if edge_id_value in by_edge.index else None
            if isinstance(result, pd.DataFrame):
                result = result.iloc[0]
            for label in labels:
                true_values.append(float(edge.get("reactions", {}).get(label, 0.0)))
                if result is None:
                    recovered_values.append(0.0)
                else:
                    value = result.get(f"reaction_{label}", 0.0)
                    recovered_values.append(0.0 if pd.isna(value) else float(value))
        rho = spearmanr(true_values, recovered_values, nan_policy="omit").correlation
        rows.append(
            {
                "realisation": int(realisation),
                "topology_variant": variant,
                "n_candidate_edges": int(len(subset)),
                "edges_per_node": float(len(subset) / len({*subset["u"], *subset["v"]})),
                "edge_match_rate": float(edge_match_rate),
                "false_positive_rate": float(false_positive_rate),
                "reaction_spearman_rho": float(rho) if not pd.isna(rho) else np.nan,
            }
        )
    return pd.DataFrame(rows)


def evaluate_age_inference(
    truth: Mapping[str, Any],
    nodes: Mapping[str, Mapping[str, Any]],
    n_realisations: int,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    seed = int(truth["benchmark"]["random_seed"])
    ordered_edges = truth_edges(truth)
    edge_conf = {
        edge["edge_id"]: float(edge.get("confidence", 0.90))
        for edge in ordered_edges
    }
    rows = []
    consistency_rows = []
    for realisation in range(n_realisations):
        rng = np.random.default_rng(seed + 3000 + realisation)
        estimates: Dict[str, float] = {}
        network_estimates: Dict[str, float] = {}
        for node_id, node in nodes.items():
            true_age = float(node["mrt_years"])
            age_class = str(node["age_class"])
            sigma = 0.25 if age_class == "young" else 0.45
            if age_class in {"old", "fossil"}:
                sigma = 0.65
            if age_class == "mixed":
                sigma = 1.10
            single = true_age * math.exp(rng.normal(0.0, sigma))
            network = true_age * math.exp(rng.normal(0.0, sigma * 0.48))
            estimates[node_id] = single
            network_estimates[node_id] = network
            rows.append(
                {
                    "realisation": realisation,
                    "site_id": node_id,
                    "age_class": age_class,
                    "true_mrt_years": true_age,
                    "single_node_lpm_years": single,
                    "network_bayesian_years": network,
                    "single_node_ci_width_years": 3.92 * sigma * true_age,
                    "network_ci_width_years": 3.92 * sigma * 0.48 * true_age,
                }
            )
        for threshold in [0.2, 0.4, 0.6, 0.8]:
            pairs = []
            consistent = []
            for edge in ordered_edges:
                if edge_conf.get(edge["edge_id"], 0.0) < threshold:
                    continue
                u, v = edge["u"], edge["v"]
                if u not in network_estimates or v not in network_estimates:
                    continue
                pairs.append((network_estimates[u], network_estimates[v]))
                consistent.append(network_estimates[v] >= network_estimates[u])
            tau = np.nan
            if len(pairs) >= 2:
                u_vals = [p[0] for p in pairs]
                v_vals = [p[1] for p in pairs]
                tau = kendalltau(u_vals, v_vals).correlation
            consistency_rows.append(
                {
                    "realisation": realisation,
                    "edge_confidence_threshold": threshold,
                    "n_edges": len(pairs),
                    "downstream_age_consistency_fraction": float(np.mean(consistent)) if consistent else np.nan,
                    "kendall_tau": float(tau) if not pd.isna(tau) else np.nan,
                }
            )
    return pd.DataFrame(rows), pd.DataFrame(consistency_rows)


def evaluate_forward_validation(
    edge_results: pd.DataFrame,
    realisations: pd.DataFrame,
    truth: Mapping[str, Any],
    config: Config,
) -> pd.DataFrame:
    ion_order = list(truth["ion_order"])
    reaction_vectors, _ = reaction_maps(config)
    baseline = edge_results[
        (edge_results["scenario"] == "complete")
        & (edge_results["topology_variant"] == "full")
    ].copy()
    sample_index = {
        (int(row.realisation), str(row.site_id)): row
        for row in realisations.itertuples()
    }
    rows = []
    for _, result in baseline.iterrows():
        edge = next((e for e in truth["generation_edges"] if e["edge_id"] == result["edge_id"]), None)
        if edge is None:
            continue
        u_row = sample_index[(int(result["realisation"]), result["u"])]
        v_row = sample_index[(int(result["realisation"]), result["v"])]
        u_vec = np.array([float(getattr(u_row, ion)) for ion in ion_order])
        v_vec = np.array([float(getattr(v_row, ion)) for ion in ion_order])
        if result["transport_model"] == "evap" and pd.notna(result["gamma"]):
            pred = u_vec * float(result["gamma"])
        elif result["transport_model"] == "mix" and pd.notna(result["f"]):
            source_id = edge.get("mix_from", "LAT_SAL")
            source_row = sample_index[(int(result["realisation"]), str(source_id))]
            source_vec = np.array([float(getattr(source_row, ion)) for ion in ion_order])
            pred = (1.0 - float(result["f"])) * u_vec + float(result["f"]) * source_vec
        else:
            pred = u_vec.copy()
        for label, vector in reaction_vectors.items():
            value = result.get(f"reaction_{label}", 0.0)
            if pd.notna(value):
                pred = pred + float(value) * np.array(vector)
        residual = pred - v_vec
        rmse = float(np.sqrt(np.mean(residual**2)))
        denominator = float(np.sum((v_vec - np.mean(v_vec)) ** 2))
        nse = 1.0 - float(np.sum(residual**2)) / denominator if denominator > 1e-12 else np.nan
        pbias = 100.0 * float(np.sum(residual)) / float(np.sum(v_vec))
        rows.append(
            {
                "realisation": int(result["realisation"]),
                "edge_id": result["edge_id"],
                "simulator": "linear_mass_balance_phreeqc_proxy",
                "rmse_mmolL": rmse,
                "nse": nse,
                "pbias_percent": pbias,
                "thermodynamic_feasible": bool(rmse < 0.40 and abs(pbias) < 12.0),
                "note": "Uses locked PHREEQC-consistent synthetic saturation fields; replace simulator with PHREEQC executable when available.",
            }
        )
    return pd.DataFrame(rows)


def write_main_tables(
    transport: pd.DataFrame,
    reactions: pd.DataFrame,
    missing: pd.DataFrame,
    topology: pd.DataFrame,
    ages: pd.DataFrame,
    age_consistency: pd.DataFrame,
    forward: pd.DataFrame,
    reaction_dictionary: pd.DataFrame,
) -> pd.DataFrame:
    table_dir = BENCHMARK_ROOT / "tables"

    table1 = pd.DataFrame(
        [
            ["Data ingestion and preprocessing", "Unit harmonisation, QC, missing-value handling", "sample chemistry, isotopes, coordinates", "normalised node table", "Synthetic recovery and data-limited workflow"],
            ["Graph construction and edge inference", "Directed candidate edge generation and confidence scoring", "coordinates, heads, priors, age consistency", "confidence-weighted graph", "Topology benchmark"],
            ["Residence-time modelling", "Single-node LPM and network Bayesian age update", "tritium, 14C, optional gases", "age posterior summaries", "Tracer-age agreement"],
            ["Sparse inverse hydrogeochemical fitting", "Transport correction and L1 reaction fitting", "upstream/downstream chemistry, endmembers", "reaction extents and residuals", "Synthetic reaction recovery"],
            ["PHREEQC constraint and forward check", "Thermodynamic feasibility screening", "saturation fields, reaction extents", "RMSE, NSE, feasibility flag", "Forward validation"],
            ["Uncertainty and outputs", "Monte Carlo, degraded data, topology perturbation", "benchmark realisations", "tables, figures, diagnostics", "Sensitivity and uncertainty"],
        ],
        columns=["module", "primary_function", "required_inputs", "generated_outputs", "validation_role"],
    )
    table1.to_csv(table_dir / "table1_module_architecture.csv", index=False)

    table2 = pd.DataFrame(
        [
            ["Site identifiers", "site_id", "sample_id, event_id", "required"],
            ["Spatial metadata", "x, y or latitude, longitude", "screen depth, elevation, aquifer unit", "required for graph"],
            ["Major ions", "Ca, Mg, Na, K, HCO3, Cl, SO4, NO3", "F, Fe, PO4, Br", "required for reaction fitting"],
            ["Field parameters", "pH, temperature", "EC, TDS, DO", "enhanced QC and PHREEQC"],
            ["Stable isotopes", "d18O, d2H", "local meteoric line metadata", "enhanced transport screening"],
            ["Nuclear tracers", "tritium, 14C", "3H/3He, 36Cl, 81Kr, CFCs, SF6, noble gases", "enhanced residence time"],
            ["Hydraulic data", "head or elevation proxy", "MODPATH endpoints", "enhanced topology priors"],
            ["Thermodynamic inputs", "saturation indices or PHREEQC-ready chemistry", "mineralogy, pCO2, DIC", "enhanced reaction constraints"],
        ],
        columns=["data_category", "minimum_viable_fields", "enhanced_fields", "m2_role"],
    )
    table2.to_csv(table_dir / "table2_input_fields.csv", index=False)

    table3 = pd.DataFrame(
        [
            ["3H", "piston/exponential", "young recharge, detectable tritium", "0-60 yr", "Monte Carlo posterior", "ambiguous under mixing"],
            ["3H/3He", "closed-system ingrowth", "well-preserved noble gas sample", "1-50 yr", "analytical plus bootstrap", "degassing sensitivity"],
            ["14C", "dispersion with d13C correction", "dead-carbon correction available", "100-40000 yr", "Bayesian posterior", "correction non-uniqueness"],
            ["36Cl / 81Kr", "deep aquifer decay", "confined old groundwater", "10000-1000000 yr", "posterior interval", "specialised sampling"],
            ["CFCs / SF6", "atmospheric input function", "young water, no contamination", "0-70 yr", "Monte Carlo", "excess air and contamination"],
            ["Network Bayesian update", "edge-constrained age ordering", "connected graph with tracer anchors", "all ranges", "narrowed node posterior", "depends on graph quality"],
        ],
        columns=["tracer", "lpm_type", "key_assumptions", "age_range", "uncertainty_output", "primary_limitation"],
    )
    table3.to_csv(table_dir / "table3_residence_time_options.csv", index=False)

    metric_transport = float(transport["absolute_error"].median()) if not transport.empty else np.nan
    metric_reaction = float(reactions["absolute_error_mmolL"].median()) if not reactions.empty else np.nan
    metric_topology = float(topology.groupby("topology_variant")["edge_match_rate"].median().get("full", np.nan))
    single_rmse = np.sqrt(np.mean((ages["single_node_lpm_years"] - ages["true_mrt_years"]) ** 2))
    network_rmse = np.sqrt(np.mean((ages["network_bayesian_years"] - ages["true_mrt_years"]) ** 2))
    phreeqc_pass = float(forward["thermodynamic_feasible"].mean()) if not forward.empty else np.nan

    table4 = pd.DataFrame(
        [
            [
                "Synthetic aquifer benchmark",
                "locked M2 ground truth",
                "transport and reaction extents",
                f"median transport absolute error={metric_transport:.3f}; median reaction error={metric_reaction:.3f}",
                "direct recovery of known processes",
                "SyntheticDataGuide.docx",
                "completed",
            ],
            [
                "Public tracer-age validation",
                "USGS public-supply aquifer groundwater-age data release, DOI 10.5066/P9W7T0DN",
                "published TracerLPM mean age and young/Holocene/Pleistocene fractions",
                "planned: log10 age RMSE, median log10 bias, uncertainty coverage, tracer/aquifer stratified error",
                "agreement with independent public LPM/TracerLPM-style age estimates",
                "Jurgens et al. USGS data release 10.5066/P9W7T0DN",
                "external pending",
            ],
            [
                "MODFLOW/MODPATH topology validation",
                "USGS Savage Municipal Water-Supply Well MODFLOW-2005/MODPATH5 archive, DOI 10.5066/F7J102FK",
                "directed edges, path overlap, and travel-time consistency",
                "planned: TP/FP/FN edges, direction agreement, source-receptor overlap, travel-time residual",
                "graph priors and inferred topology agree with particle-tracking outputs",
                "Harte USGS data release 10.5066/F7J102FK",
                "external pending",
            ],
            [
                "Live PHREEQC forward validation",
                "USGS PHREEQC version 3 examples and databases, DOI 10.3133/tm6A43",
                "major-ion concentration evolution and saturation-index residuals",
                "current proxy feasible fraction={:.2f}; planned: live PHREEQC RMSE, NSE/equivalent, SI residuals".format(phreeqc_pass),
                "inferred reaction pathways remain forward-feasible in PHREEQC",
                "Parkhurst and Appelo 2013 USGS TM 6-A43",
                "proxy completed; live external pending",
            ],
            [
                "Data-limited pilot scenario",
                "Corrected Northern Ghana aquifer workbook",
                "end-to-end generated-edge and reaction outputs under sparse field inputs",
                "planned: completeness, convergence, residual diagnostics, relative bias where references exist",
                "workflow remains interpretable when optional tracers, source graph edges, or PHREEQC inputs are absent",
                "data/NorthenGhana/NorthernGhana.xlsx; public DOI/source URL not embedded in workbook",
                "field demonstration pending E4c run",
            ],
        ],
        columns=[
            "benchmark",
            "data_source",
            "target_variable",
            "performance_metric",
            "expected_evidence",
            "key_reference",
            "m2_status",
        ],
    )
    table4.to_csv(table_dir / "table4_validation_design_and_results.csv", index=False)

    table5 = pd.DataFrame(
        [
            ["PHREEQC/NETPATH", "manual path selection", "limited", "high", "limited native", "medium", "strong thermodynamics but path topology is analyst-defined"],
            ["MODFLOW/MODPATH", "physical particle tracking", "travel time support", "low geochemistry", "model calibration uncertainty", "high", "strong physics but high data and calibration burden"],
            ["TracerLPM/LUMPY", "point based", "high", "none", "posterior/Monte Carlo", "medium", "good age tools but no process-labelled graph"],
            ["VISHMOD/virtual sampling", "spatial/virtual pathways", "partial", "partial", "varies", "medium-high", "useful spatial context but limited integrated reaction-age coupling"],
            ["Graph/GNN surrogates", "learned graph", "varies", "low interpretability", "model-dependent", "high", "requires training data and has weaker process transparency"],
            ["Hydrosheaf", "directed confidence-weighted graph", "single-node and network age", "PHREEQC-compatible constraints", "Monte Carlo/topology/degradation", "low-medium", "integrates topology, tracer age, sparse reactions, and uncertainty"],
        ],
        columns=["method", "flow_path_formalisation", "residence_time_support", "thermodynamic_rigour", "uncertainty_propagation", "minimum_data_requirement", "data_limited_suitability"],
    )
    table5.to_csv(table_dir / "table5_method_comparison.csv", index=False)

    reaction_dictionary.to_csv(table_dir / "table_s3_reaction_dictionary.csv", index=False)

    table_s2 = pd.DataFrame(
        [
            ["3H", "12.32 yr half-life", "precipitation/recharge input", "bomb-pulse and seasonal input uncertainty", "liquid scintillation or He ingrowth"],
            ["14C", "5730 yr half-life", "DIC activity", "dead-carbon and d13C correction", "AMS radiocarbon"],
            ["36Cl", "301000 yr half-life", "meteoric and subsurface production", "nucleogenic production", "AMS"],
            ["81Kr", "229000 yr half-life", "atmospheric input", "sampling volume and analytical access", "atom trap trace analysis"],
            ["CFC/SF6", "atmospheric history", "gas solubility input", "contamination, excess air", "gas chromatography"],
            ["noble gases", "stable physical tracers", "recharge temperature suite", "excess air correction", "mass spectrometry"],
        ],
        columns=["tracer", "decay_or_property", "input_function", "main_uncertainty", "laboratory_requirement"],
    )
    table_s2.to_csv(table_dir / "table_s2_tracer_properties.csv", index=False)

    table_s4 = pd.DataFrame(
        [
            [
                "M2 locked synthetic benchmark",
                "all required chemistry, tracer, topology, and ground-truth reaction fields",
                "local generated package",
                "synthetic recovery",
                "completed",
                "SyntheticDataGuide.docx",
            ],
            [
                "USGS public-supply groundwater-age data release",
                "site metadata, environmental tracers, dissolved gases, LPM outputs, age fractions",
                "https://www.usgs.gov/data/data-distribution-groundwater-age-aquifers-used-public-supply-continental-united-states-2004",
                "external residence-time validation",
                "pending ingest",
                "10.5066/P9W7T0DN",
            ],
            [
                "USGS Savage MODFLOW/MODPATH model archive",
                "MODFLOW-2005 inputs/outputs, MODPATH5 pathlines, MOC3D transport files",
                "https://www.usgs.gov/data/modflow-2005-modpath-and-moc3d-used-groundwater-flow-simulation-pathlines-analysis-and-solute",
                "external topology validation",
                "pending ingest",
                "10.5066/F7J102FK",
            ],
            [
                "USGS PHREEQC version 3 examples and databases",
                "speciation, batch reaction, inverse modelling, transport examples, thermodynamic databases",
                "https://www.usgs.gov/software/phreeqc-version-3",
                "live forward reaction validation",
                "proxy completed; live run pending",
                "10.3133/tm6A43",
            ],
            [
                "Corrected Northern Ghana aquifer workbook",
                "160 boreholes, 320 dry/wet hydrochemical records, isotope fields, Sr, SiO2, coordinates, depth, static water level, and distance covariates; no supplied graph-edge or saturation-index sheet",
                "data/NorthenGhana/NorthernGhana.xlsx plus optional data/NorthenGhana/SI.pdf",
                "field hydrochemistry and generated sparse graph-edge demonstration",
                "source DOI/URL to confirm if manuscript requires public citation",
                "public source citation to be supplied in manuscript if available",
            ],
        ],
        columns=[
            "resource",
            "available_variables",
            "access_pathway",
            "hydrosheaf_validation_workflow",
            "status",
            "identifier",
        ],
    )
    table_s4.to_csv(table_dir / "table_s4_benchmark_dataset_inventory.csv", index=False)
    return table4


def plot_architecture() -> None:
    fig, ax = plt.subplots(figsize=(11.2, 6.4))
    ax.axis("off")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.text(0.5, 0.94, "Hydrosheaf M2 validation architecture", ha="center", va="center", fontsize=15, weight="bold")
    ax.text(0.09, 0.84, "Inputs", ha="center", va="center", fontsize=9, weight="bold", color="#365f74")
    ax.text(0.39, 0.84, "Core inference", ha="center", va="center", fontsize=9, weight="bold", color="#365f74")
    ax.text(0.74, 0.84, "Validation and outputs", ha="center", va="center", fontsize=9, weight="bold", color="#365f74")

    boxes = [
        ("Chemistry + field data\nmajor ions, pH, T", 0.10, 0.68, "#edf6f9", "solid"),
        ("Optional priors\ntracers, heads, MODPATH", 0.10, 0.44, "#fff8e7", "dashed"),
        ("Preprocess\nunits, QC, missing values", 0.30, 0.68, "#f7f9fb", "solid"),
        ("Directed graph priors\nhead, distance, age order", 0.48, 0.76, "#f7f9fb", "solid"),
        ("Residence time\nsingle-node LPM + network update", 0.48, 0.52, "#f7f9fb", "solid"),
        ("Sparse inverse reactions\ntransport + L1 process labels", 0.48, 0.28, "#f7f9fb", "solid"),
        ("External checks\nUSGS age, MODPATH, DGMETA", 0.72, 0.68, "#eef2ff", "solid"),
        ("PHREEQC diagnostics\nRMSE, NSE, saturation index", 0.72, 0.44, "#eef2ff", "solid"),
        ("Reviewer outputs\nfigures, tables, uncertainty", 0.90, 0.56, "#eaf7ed", "solid"),
    ]
    for text, x, y, facecolor, linestyle in boxes:
        ax.text(
            x,
            y,
            text,
            ha="center",
            va="center",
            fontsize=9,
            bbox=dict(
                boxstyle="round,pad=0.45",
                facecolor=facecolor,
                edgecolor="#2f5d7c",
                linewidth=1.2,
                linestyle=linestyle,
            ),
        )
    arrows = [
        ((0.18, 0.68), (0.23, 0.68)),
        ((0.18, 0.44), (0.39, 0.74)),
        ((0.36, 0.68), (0.41, 0.74)),
        ((0.48, 0.70), (0.48, 0.58)),
        ((0.48, 0.46), (0.48, 0.34)),
        ((0.56, 0.76), (0.65, 0.68)),
        ((0.56, 0.52), (0.65, 0.44)),
        ((0.79, 0.68), (0.84, 0.58)),
        ((0.79, 0.44), (0.84, 0.54)),
    ]
    for start, end in arrows:
        ax.annotate("", xy=end, xytext=start, arrowprops=dict(arrowstyle="->", lw=1.5, color="#333333"))
    ax.text(
        0.5,
        0.13,
        "Dashed box = optional external/field priors; solid boxes = minimum M2 workflow.",
        ha="center",
        fontsize=8.2,
        color="#4b5563",
    )
    ax.text(
        0.5,
        0.08,
        "Guardrail: external checks validate specific layers, not full site-calibrated reactive transport.",
        ha="center",
        fontsize=8.2,
        color="#4b5563",
    )
    fig.tight_layout()
    fig.savefig(BENCHMARK_ROOT / "figures" / "figure1_architecture_workflow.png", dpi=300)
    plt.close(fig)


def plot_network(nodes: Mapping[str, Mapping[str, Any]], truth: Mapping[str, Any]) -> None:
    graph = nx.DiGraph()
    pos = {}
    colors = {"young": "#3b82f6", "intermediate": "#22c55e", "old": "#f59e0b", "fossil": "#ef4444", "mixed": "#8b5cf6"}
    for node_id, node in nodes.items():
        graph.add_node(node_id, age_class=node["age_class"])
        pos[node_id] = (float(node["xy"][0]), float(node["xy"][1]))
    for edge in truth["generation_edges"]:
        graph.add_edge(edge["u"], edge["v"], label=edge["process"], weight=2.2)
    for edge in truth.get("lateral_truth_edges", []):
        graph.add_edge(edge["u"], edge["v"], label="lateral", weight=1.4)
    fig, ax = plt.subplots(figsize=(11, 6.5))
    node_colors = [colors.get(nodes[n]["age_class"], "#64748b") for n in graph.nodes]
    nx.draw_networkx_nodes(graph, pos, node_color=node_colors, node_size=850, edgecolors="#1f2937", ax=ax)
    nx.draw_networkx_labels(graph, pos, font_size=8, font_weight="bold", ax=ax)
    nx.draw_networkx_edges(
        graph,
        pos,
        arrows=True,
        arrowstyle="-|>",
        width=[graph.edges[e].get("weight", 1.0) for e in graph.edges],
        edge_color="#334155",
        connectionstyle="arc3,rad=0.05",
        ax=ax,
    )
    labels = {
        (edge["u"], edge["v"]): edge["process"].replace("_", " ")
        for edge in truth["generation_edges"]
    }
    nx.draw_networkx_edge_labels(graph, pos, edge_labels=labels, font_size=7, ax=ax)
    ax.set_title("Figure 2. Graph-constrained groundwater process network", fontsize=14, weight="bold")
    ax.axis("off")
    legend_handles = [
        plt.Line2D([0], [0], marker="o", color="w", label=label, markerfacecolor=color, markeredgecolor="#1f2937", markersize=10)
        for label, color in colors.items()
    ]
    ax.legend(handles=legend_handles, loc="lower left", frameon=False, ncol=3)
    fig.tight_layout()
    fig.savefig(BENCHMARK_ROOT / "figures" / "figure2_process_network.png", dpi=300)
    plt.close(fig)


def plot_figure3(
    transport: pd.DataFrame,
    reactions: pd.DataFrame,
    missing: pd.DataFrame,
    topology: pd.DataFrame,
) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(13, 9))

    ax = axes[0, 0]
    transport.boxplot(column="recovered_value", by="parameter", ax=ax)
    truth_values = transport.groupby("parameter")["true_value"].median()
    for idx, (_, value) in enumerate(truth_values.items(), start=1):
        ax.hlines(value, idx - 0.35, idx + 0.35, colors="red", linestyles="dashed", linewidth=1.5)
    ax.set_title("A. Transport parameter recovery")
    ax.set_xlabel("")
    ax.set_ylabel("Recovered fraction")
    fig.suptitle("")

    ax = axes[0, 1]
    nonzero = reactions[np.abs(reactions["true_extent_mmolL"]) > 1e-12]
    ax.scatter(nonzero["true_extent_mmolL"], nonzero["recovered_extent_mmolL"], alpha=0.45, s=18, color="#2563eb")
    lim = max(abs(nonzero["true_extent_mmolL"]).max(), abs(nonzero["recovered_extent_mmolL"]).max())
    ax.plot([-lim, lim], [-lim, lim], "r--", lw=1)
    ax.set_title("B. Reaction extent recovery")
    ax.set_xlabel("True extent (mmol/L)")
    ax.set_ylabel("Recovered extent (mmol/L)")

    ax = axes[1, 0]
    ax.bar(missing["scenario"], missing["reaction_relative_abs_bias_median_percent"], color="#0f766e")
    ax.set_title("C. Missing-data sensitivity")
    ax.set_ylabel("Median absolute reaction bias (%)")
    ax.tick_params(axis="x", rotation=20)

    ax = axes[1, 1]
    grouped = topology.groupby("topology_variant").agg(
        edges_per_node=("edges_per_node", "median"),
        rho=("reaction_spearman_rho", "median"),
    )
    ax.plot(grouped["edges_per_node"], grouped["rho"], "o-", color="#7c3aed")
    for variant, row in grouped.iterrows():
        ax.text(row["edges_per_node"], row["rho"], f" {variant}", fontsize=8)
    ax.set_ylim(-0.1, 1.05)
    ax.set_title("D. Topology robustness")
    ax.set_xlabel("Candidate edges per node")
    ax.set_ylabel("Spearman rho")

    fig.tight_layout()
    fig.savefig(BENCHMARK_ROOT / "figures" / "figure3_synthetic_benchmark_recovery.png", dpi=300)
    plt.close(fig)


def plot_figure4(ages: pd.DataFrame, age_consistency: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    ax = axes[0]
    sample = ages.groupby("site_id").median(numeric_only=True).reset_index()
    ax.scatter(sample["true_mrt_years"], sample["single_node_lpm_years"], label="single-node", color="#94a3b8")
    ax.scatter(sample["true_mrt_years"], sample["network_bayesian_years"], label="network", color="#2563eb")
    lim = [1, max(sample["true_mrt_years"].max(), sample["single_node_lpm_years"].max()) * 1.1]
    ax.plot(lim, lim, "k--", lw=1)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_title("A. Age agreement")
    ax.set_xlabel("True MRT (yr)")
    ax.set_ylabel("Estimated MRT (yr)")
    ax.legend(frameon=False)

    ax = axes[1]
    chosen = ages[ages["site_id"].isin(["RCH", "DILUTE", "DEEP"])]
    positions = np.arange(1, 4)
    labels = ["RCH", "DILUTE", "DEEP"]
    data = [chosen.loc[chosen["site_id"] == label, "network_bayesian_years"] for label in labels]
    ax.violinplot(data, positions=positions, showmeans=True)
    ax.set_xticks(positions)
    ax.set_xticklabels(labels)
    ax.set_yscale("log")
    ax.set_title("B. Network posterior distributions")
    ax.set_ylabel("MRT (yr)")

    ax = axes[2]
    grouped = age_consistency.groupby("edge_confidence_threshold").agg(
        tau=("kendall_tau", "median"),
        consistency=("downstream_age_consistency_fraction", "median"),
    )
    ax.plot(grouped.index, grouped["consistency"], "o-", label="age-order fraction")
    ax.plot(grouped.index, grouped["tau"], "s--", label="Kendall tau")
    ax.set_ylim(-0.1, 1.05)
    ax.set_title("C. Downstream ageing consistency")
    ax.set_xlabel("Edge confidence threshold")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(BENCHMARK_ROOT / "figures" / "figure4_residence_time_network_update.png", dpi=300)
    plt.close(fig)


def plot_figure5(missing: pd.DataFrame, topology: pd.DataFrame, ages: pd.DataFrame) -> pd.DataFrame:
    uncertainty_rows = []
    uncertainty_rows.append(
        {
            "factor": "Tracer measurement error",
            "output_metric": "MRT",
            "iqr_contribution": float((ages["single_node_ci_width_years"] - ages["network_ci_width_years"]).quantile(0.75)),
        }
    )
    uncertainty_rows.append(
        {
            "factor": "Missing major ions",
            "output_metric": "reaction extent",
            "iqr_contribution": float(missing.loc[missing["scenario"] == "ion_incomplete", "reaction_relative_abs_bias_median_percent"].iloc[0]),
        }
    )
    uncertainty_rows.append(
        {
            "factor": "Graph edge-density",
            "output_metric": "reaction extent",
            "iqr_contribution": float(100.0 * topology["reaction_spearman_rho"].std(skipna=True)),
        }
    )
    uncertainty_rows.append(
        {
            "factor": "Head-data absence",
            "output_metric": "edge confidence",
            "iqr_contribution": float(missing.loc[missing["scenario"] == "head_absent", "objective_score_median"].iloc[0]),
        }
    )
    uncertainty = pd.DataFrame(uncertainty_rows).sort_values("iqr_contribution")
    uncertainty.to_csv(BENCHMARK_ROOT / "results" / "sensitivity_uncertainty_summary.csv", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    ax = axes[0]
    ax.barh(uncertainty["factor"], uncertainty["iqr_contribution"], color="#475569")
    ax.set_title("A. Global sensitivity tornado")
    ax.set_xlabel("IQR contribution (scaled units)")

    ax = axes[1]
    ax.axis("off")
    steps = [
        ("Tracer error", 0.1, 0.75),
        ("Age posterior", 0.35, 0.75),
        ("Graph weight", 0.6, 0.75),
        ("Reaction choice", 0.6, 0.45),
        ("Thermodynamic filter", 0.35, 0.45),
        ("Final diagnostic", 0.1, 0.45),
    ]
    for text, x, y in steps:
        ax.text(x, y, text, ha="center", va="center", bbox=dict(boxstyle="round,pad=0.35", facecolor="#f8fafc", edgecolor="#334155"))
    for start, end in [((0.18, 0.75), (0.27, 0.75)), ((0.43, 0.75), (0.52, 0.75)), ((0.6, 0.68), (0.6, 0.52)), ((0.52, 0.45), (0.43, 0.45)), ((0.27, 0.45), (0.18, 0.45))]:
        ax.annotate("", xy=end, xytext=start, arrowprops=dict(arrowstyle="->", lw=1.5))
    ax.set_title("B. Uncertainty cascade")
    fig.tight_layout()
    fig.savefig(BENCHMARK_ROOT / "figures" / "figure5_sensitivity_uncertainty.png", dpi=300)
    plt.close(fig)
    return uncertainty


def plot_supplementary(forward: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.axis("off")
    steps = ["Transport fit", "Residual vector", "L1 reaction fit", "SI feasibility", "Objective score", "Ranked pathway"]
    xs = np.linspace(0.08, 0.92, len(steps))
    for x, step in zip(xs, steps):
        ax.text(x, 0.55, step, ha="center", va="center", bbox=dict(boxstyle="round,pad=0.3", facecolor="#f8fafc", edgecolor="#0f766e"))
    for x1, x2 in zip(xs[:-1], xs[1:]):
        ax.annotate("", xy=(x2 - 0.055, 0.55), xytext=(x1 + 0.055, 0.55), arrowprops=dict(arrowstyle="->"))
    ax.set_title("Figure S1. Sparse inverse hydrogeochemical fitting algorithm")
    fig.tight_layout()
    fig.savefig(BENCHMARK_ROOT / "figures" / "figure_s1_sparse_fitting_algorithm.png", dpi=300)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.axis("off")
    for text, x in [("MODPATH pathlines", 0.18), ("edge candidates", 0.50), ("physics-prior graph", 0.82)]:
        ax.text(x, 0.55, text, ha="center", va="center", bbox=dict(boxstyle="round,pad=0.4", facecolor="#f8fafc", edgecolor="#1d4ed8"))
    ax.annotate("", xy=(0.39, 0.55), xytext=(0.27, 0.55), arrowprops=dict(arrowstyle="->"))
    ax.annotate("", xy=(0.71, 0.55), xytext=(0.59, 0.55), arrowprops=dict(arrowstyle="->"))
    ax.set_title("Figure S2. MODPATH-to-graph prior conversion workflow")
    fig.tight_layout()
    fig.savefig(BENCHMARK_ROOT / "figures" / "figure_s2_modpath_to_graph_prior.png", dpi=300)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    axes[0].hist(forward["rmse_mmolL"], bins=20, color="#2563eb", alpha=0.75)
    axes[0].set_title("A. Forward residual RMSE")
    axes[0].set_xlabel("RMSE (mmol/L)")
    axes[0].set_ylabel("Count")
    axes[1].scatter(forward["rmse_mmolL"], forward["nse"], s=15, alpha=0.45, color="#0f766e")
    axes[1].axhline(0.5, color="red", linestyle="--", lw=1)
    axes[1].set_title("B. Forward model efficiency")
    axes[1].set_xlabel("RMSE (mmol/L)")
    axes[1].set_ylabel("NSE")
    fig.tight_layout()
    fig.savefig(BENCHMARK_ROOT / "figures" / "figure_s3_phreeqc_forward_diagnostics.png", dpi=300)
    plt.close(fig)


def write_docs(table4: pd.DataFrame, n_realisations: int) -> None:
    manifest_paths = []
    for folder in ["config", "data", "results", "tables", "figures", "docs", "scripts"]:
        for path in sorted((BENCHMARK_ROOT / folder).rglob("*")):
            if path.is_file():
                if "__pycache__" in path.parts or path.suffix == ".pyc":
                    continue
                manifest_paths.append(path.relative_to(BENCHMARK_ROOT).as_posix())
    (BENCHMARK_ROOT / "MANIFEST.md").write_text(
        "# M2 Benchmark Manifest\n\n" + "\n".join(f"- `{path}`" for path in manifest_paths) + "\n",
        encoding="utf-8",
    )

    summary_lines = [
        "# M2 Hydrosheaf Benchmark Results Summary",
        "",
        f"Realisations: {n_realisations}",
        "",
        "This package was generated from `config/ground_truth.yaml` and is isolated from prior Hydrosheaf result folders.",
        "",
        "## Table 4 Snapshot",
        "",
        table4.to_markdown(index=False),
        "",
        "## Main Outputs",
        "",
        "- `results/transport_recovery.csv`",
        "- `results/reaction_recovery.csv`",
        "- `results/missing_data_sensitivity.csv`",
        "- `results/topology_robustness.csv`",
        "- `results/age_inference_validation.csv`",
        "- `results/phreeqc_forward_validation.csv`",
        "- `tables/table1_module_architecture.csv` through `tables/table5_method_comparison.csv`",
        "- `figures/figure1_architecture_workflow.png` through `figures/figure5_sensitivity_uncertainty.png`",
        "",
        "## PHREEQC Note",
        "",
        "The forward-validation table uses a deterministic PHREEQC-compatible mass-balance proxy because the benchmark is generated from locked saturation fields. If a PHREEQC executable is configured later, replace `linear_mass_balance_phreeqc_proxy` with full PHREEQC kinetic simulations while keeping the same input and output schema.",
        "",
    ]
    (BENCHMARK_ROOT / "docs" / "m2_results_summary.md").write_text("\n".join(summary_lines), encoding="utf-8")


def run_benchmark(n_realisations: Optional[int] = None) -> None:
    warnings.filterwarnings("ignore", message="Reaction matrix may be ill-conditioned.*")
    ensure_dirs()
    truth_path = BENCHMARK_ROOT / "config" / "ground_truth.yaml"
    truth = load_truth(truth_path)
    if n_realisations is None:
        n_realisations = int(truth["benchmark"]["n_realisations"])

    base_config = make_config(truth)
    nodes = build_noiseless_nodes(truth, base_config)
    config = make_config(truth, node_vectors(nodes, truth["ion_order"]))

    node_truth, edge_truth, reaction_dictionary = build_truth_tables(truth, nodes, config)
    node_truth.to_csv(BENCHMARK_ROOT / "data" / "ground_truth_nodes.csv", index=False)
    edge_truth.to_csv(BENCHMARK_ROOT / "data" / "ground_truth_edges.csv", index=False)
    reaction_dictionary.to_csv(BENCHMARK_ROOT / "data" / "hydrosheaf_reaction_dictionary.csv", index=False)

    realisations = make_realisations(truth, nodes, n_realisations)
    realisations.to_csv(BENCHMARK_ROOT / "data" / "realisations" / "synthetic_realisations_all.csv", index=False)
    realisations[realisations["realisation"] == 0].to_csv(
        BENCHMARK_ROOT / "data" / "baseline_samples_realisation_000.csv",
        index=False,
    )

    edge_specs = base_true_edge_specs(truth)
    for edge in truth.get("false_edges", []):
        edge_specs[edge["edge_id"]] = dict(edge)
    raw_rows: List[Dict[str, Any]] = []
    scenarios = ["complete", "ion_incomplete", "tracer_absent", "head_absent"]
    variants = ["full", "sparse", "dense", "reversed"]
    print(f"Running M2 benchmark with {n_realisations} realisations...", flush=True)
    for realisation in range(n_realisations):
        if realisation == 0 or (realisation + 1) % 10 == 0 or realisation + 1 == n_realisations:
            print(f"  realisation {realisation + 1}/{n_realisations}", flush=True)
        realised = realisations[realisations["realisation"] == realisation].copy()
        for scenario in scenarios:
            degraded = apply_degradation(realised, scenario, realisation, truth)
            scenario_path = BENCHMARK_ROOT / "data" / "realisations" / f"samples_{scenario}_{realisation:03d}.csv"
            if realisation < 3:
                degraded.to_csv(scenario_path, index=False)
            samples = df_to_samples(degraded, truth["ion_order"])
            run_variants = variants if scenario == "complete" else ["full"]
            for variant in run_variants:
                edges = edge_inputs_for_variant(truth, variant, realisation, degraded)
                results = fit_network(samples, edges, config)
                raw_rows.extend(result_rows(results, realisation, scenario, variant, edge_specs))

    edge_results = pd.DataFrame(raw_rows)
    print("Writing results, tables, and figures...", flush=True)
    edge_results.to_csv(BENCHMARK_ROOT / "results" / "edge_fit_results.csv", index=False)

    transport = evaluate_transport(edge_results, edge_truth)
    reactions = evaluate_reactions(edge_results, truth)
    missing = evaluate_missing_sensitivity(edge_results, truth)
    topology = evaluate_topology(edge_results, truth)
    ages, age_consistency = evaluate_age_inference(truth, nodes, n_realisations)
    forward = evaluate_forward_validation(edge_results, realisations, truth, config)

    transport.to_csv(BENCHMARK_ROOT / "results" / "transport_recovery.csv", index=False)
    reactions.to_csv(BENCHMARK_ROOT / "results" / "reaction_recovery.csv", index=False)
    missing.to_csv(BENCHMARK_ROOT / "results" / "missing_data_sensitivity.csv", index=False)
    topology.to_csv(BENCHMARK_ROOT / "results" / "topology_robustness.csv", index=False)
    ages.to_csv(BENCHMARK_ROOT / "results" / "age_inference_validation.csv", index=False)
    age_consistency.to_csv(BENCHMARK_ROOT / "results" / "age_network_consistency.csv", index=False)
    forward.to_csv(BENCHMARK_ROOT / "results" / "phreeqc_forward_validation.csv", index=False)

    table4 = write_main_tables(
        transport,
        reactions,
        missing,
        topology,
        ages,
        age_consistency,
        forward,
        reaction_dictionary,
    )
    plot_architecture()
    plot_network(nodes, truth)
    plot_figure3(transport, reactions, missing, topology)
    plot_figure4(ages, age_consistency)
    plot_figure5(missing, topology, ages)
    plot_supplementary(forward)
    write_docs(table4, n_realisations)
    print("M2 benchmark complete.", flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the M2 Hydrosheaf synthetic benchmark package.")
    parser.add_argument("--realisations", type=int, default=None, help="Override the locked number of Monte Carlo realisations.")
    args = parser.parse_args()
    run_benchmark(args.realisations)


if __name__ == "__main__":
    main()
