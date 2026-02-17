"""Reaction hyperparameter tuning for interpretability."""

from dataclasses import dataclass, field, replace
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import numpy as np

from ..config import Config
from ..data.schema import vector_from_sample
from ..graph.build import EdgeInput, build_edges
from ..graph.types import Edge
from ..inference.edge_fit import EdgeResult, fit_edge
from ..inference.network_fit import estimate_edge_residence_time_days
from ..models.reactions import build_reaction_dictionary
from ..models.redox import get_redox_constraints
from ..phreeqc.constraints import build_edge_bounds
from ..phreeqc.runner import run_phreeqc


@dataclass
class EdgeFitData:
    edge: Edge
    x_u: List[float]
    x_v: List[float]
    obs_u: Mapping[str, float]
    obs_v: Mapping[str, float]
    bounds: Dict[str, object]
    residence_time_days: Optional[float]
    config_edge: Config


@dataclass
class EdgeStabilitySummary:
    edge_id: str
    u: str
    v: str
    mean_jaccard: float
    mean_selected_reactions: float
    median_r2: float
    selection_prob: Dict[str, float] = field(default_factory=dict)
    sign_prob: Dict[str, float] = field(default_factory=dict)
    stable_reactions: List[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {
            "edge_id": self.edge_id,
            "u": self.u,
            "v": self.v,
            "mean_jaccard": self.mean_jaccard,
            "mean_selected_reactions": self.mean_selected_reactions,
            "median_r2": self.median_r2,
            "selection_prob": self.selection_prob,
            "sign_prob": self.sign_prob,
            "stable_reactions": self.stable_reactions,
        }


@dataclass
class LambdaStabilitySummary:
    lambda_l1: float
    stability_mean: float
    stability_se: float
    mean_selected_reactions: float
    median_r2: float
    edges_considered: int
    selection_prob: Dict[str, float] = field(default_factory=dict)
    sign_prob: Dict[str, float] = field(default_factory=dict)
    edge_summaries: List[EdgeStabilitySummary] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {
            "lambda_l1": self.lambda_l1,
            "stability_mean": self.stability_mean,
            "stability_se": self.stability_se,
            "mean_selected_reactions": self.mean_selected_reactions,
            "median_r2": self.median_r2,
            "edges_considered": self.edges_considered,
            "selection_prob": self.selection_prob,
            "sign_prob": self.sign_prob,
            "edge_summaries": [summary.to_dict() for summary in self.edge_summaries],
        }


@dataclass
class TuningReport:
    recommended_lambda_l1: float
    selection_rule: str
    candidates: List[LambdaStabilitySummary]
    min_r2: Optional[float] = None
    max_reactions: Optional[float] = None

    def to_dict(self) -> dict:
        return {
            "recommended_lambda_l1": self.recommended_lambda_l1,
            "selection_rule": self.selection_rule,
            "min_r2": self.min_r2,
            "max_reactions": self.max_reactions,
            "candidates": [candidate.to_dict() for candidate in self.candidates],
        }


def _sample_map(samples: object) -> Dict[str, Mapping[str, object]]:
    if isinstance(samples, Mapping):
        return {str(k): v for k, v in samples.items()}
    if isinstance(samples, Sequence):
        mapping: Dict[str, Mapping[str, object]] = {}
        for row in samples:
            site_id = row.get("site_id")
            if site_id is None:
                raise ValueError("Each sample row must include site_id.")
            mapping[str(site_id)] = row
        return mapping
    raise TypeError("Unsupported samples input type.")


def _mean_jaccard(sets: List[set]) -> float:
    if len(sets) < 2:
        return 1.0 if sets else 0.0
    total = 0.0
    count = 0
    for i in range(len(sets)):
        for j in range(i + 1, len(sets)):
            union = sets[i] | sets[j]
            if not union:
                score = 1.0
            else:
                score = len(sets[i] & sets[j]) / len(union)
            total += score
            count += 1
    return total / max(1, count)


def _selection_summary(
    labels: List[str],
    extents: List[float],
    threshold: float,
) -> Tuple[set, Dict[str, int], Dict[str, int]]:
    selected = set()
    pos_counts: Dict[str, int] = {}
    neg_counts: Dict[str, int] = {}
    for label, extent in zip(labels, extents):
        if abs(extent) < threshold:
            continue
        selected.add(label)
        if extent > 0:
            pos_counts[label] = pos_counts.get(label, 0) + 1
        elif extent < 0:
            neg_counts[label] = neg_counts.get(label, 0) + 1
    return selected, pos_counts, neg_counts


def _run_fit(
    fit_data: EdgeFitData,
    config_edge: Config,
    x_u: List[float],
    x_v: List[float],
) -> Optional[EdgeResult]:
    try:
        return fit_edge(
            x_u,
            x_v,
            config_edge,
            edge_id=fit_data.edge.edge_id,
            u=fit_data.edge.u,
            v=fit_data.edge.v,
            obs_v=fit_data.obs_v,
            bounds=fit_data.bounds,
            obs_u=fit_data.obs_u,
            residence_time_days=fit_data.residence_time_days,
        )
    except Exception:
        return None


def _build_edge_data(
    samples: object,
    edges: Iterable[EdgeInput],
    config: Config,
    phreeqc_results: Optional[Mapping[str, Mapping[str, object]]] = None,
    residence_time_overrides: Optional[Mapping[str, float]] = None,
) -> List[EdgeFitData]:
    sample_map = _sample_map(samples)
    built_edges = build_edges(edges)

    if config.phreeqc_enabled and phreeqc_results is None:
        phreeqc_results = run_phreeqc(sample_map.values(), config)

    edge_data: List[EdgeFitData] = []
    for edge in built_edges:
        if edge.u not in sample_map or edge.v not in sample_map:
            continue
        sample_u = sample_map[edge.u]
        sample_v = sample_map[edge.v]
        x_u, sample_u_norm = vector_from_sample(
            sample_u,
            config.ion_order,
            config.missing_policy,
            config.detection_limit_policy,
        )
        if x_u is None:
            continue
        x_v, sample_v_norm = vector_from_sample(
            sample_v,
            config.ion_order,
            config.missing_policy,
            config.detection_limit_policy,
        )
        if x_v is None:
            continue

        layer_key = getattr(config, "layer_key", "aquifer_layer")
        layer_value = sample_v.get(layer_key)
        layer_idx = None
        try:
            if layer_value not in (None, ""):
                layer_idx = int(layer_value)
        except (TypeError, ValueError):
            layer_idx = None

        config_edge = config
        if layer_idx is not None:
            minerals_for_layer = getattr(config, "layer_mineral_map", {}).get(layer_idx)
            if minerals_for_layer:
                config_edge = replace(config, active_minerals=list(minerals_for_layer))

        reaction_matrix, reaction_labels, mineral_mask = build_reaction_dictionary(
            config_edge
        )

        edge_bounds: Dict[str, object] = {}
        if config_edge.phreeqc_enabled and phreeqc_results is not None:
            edge_bounds = build_edge_bounds(
                phreeqc_results,
                [edge],
                reaction_labels,
                mineral_mask,
                config_edge,
            ).get(edge.edge_id, {})

        if "lb" not in edge_bounds or not isinstance(edge_bounds.get("lb"), list):
            edge_bounds["lb"] = [None] * len(reaction_labels)
        if "ub" not in edge_bounds or not isinstance(edge_bounds.get("ub"), list):
            edge_bounds["ub"] = [None] * len(reaction_labels)
        if "constraints_active" not in edge_bounds or not isinstance(
            edge_bounds.get("constraints_active"), dict
        ):
            edge_bounds["constraints_active"] = {}

        redox_overrides = get_redox_constraints(sample_v_norm, reaction_labels)
        if redox_overrides:
            lb_list = edge_bounds["lb"]
            ub_list = edge_bounds["ub"]
            assert isinstance(lb_list, list)
            assert isinstance(ub_list, list)
            for i, label in enumerate(reaction_labels):
                if label in redox_overrides:
                    lb_val, ub_val = redox_overrides[label]
                    lb_list[i] = lb_val
                    ub_list[i] = ub_val
            ca_dict = edge_bounds["constraints_active"]
            assert isinstance(ca_dict, dict)
            ca_dict["redox"] = "active"

        tau_edge = None
        if residence_time_overrides is not None:
            override = residence_time_overrides.get(edge.edge_id)
            try:
                tau_edge = float(override) if override is not None else None
            except (TypeError, ValueError):
                tau_edge = None
        if tau_edge is None:
            tau_edge = estimate_edge_residence_time_days(edge.attrs or {}, config_edge)

        edge_data.append(
            EdgeFitData(
                edge=edge,
                x_u=x_u,
                x_v=x_v,
                obs_u=sample_u_norm,
                obs_v=sample_v_norm,
                bounds=edge_bounds,
                residence_time_days=tau_edge,
                config_edge=config_edge,
            )
        )
    return edge_data


def _evaluate_lambda(
    edge_data: List[EdgeFitData],
    lambda_l1: float,
    n_resamples: int,
    perturbation: str,
    input_uncertainty_pct: float,
    selection_threshold: float,
    rng: np.random.Generator,
    stable_prob_threshold: float,
) -> LambdaStabilitySummary:
    edge_jaccards: List[float] = []
    edge_sizes: List[float] = []
    edge_r2: List[float] = []
    selection_counts: Dict[str, int] = {}
    sign_counts: Dict[str, int] = {}
    total_samples = 0
    edge_summaries: List[EdgeStabilitySummary] = []

    use_methods = [perturbation]
    if perturbation == "both":
        use_methods = ["monte_carlo", "bootstrap"]

    for fit_data in edge_data:
        config_edge = replace(fit_data.config_edge, lambda_l1=float(lambda_l1))

        base_result = _run_fit(
            fit_data, config_edge, fit_data.x_u, fit_data.x_v
        )
        if base_result is None:
            continue

        base_sets: List[set] = []
        base_sizes: List[int] = []
        base_r2: List[float] = []
        label_counts: Dict[str, int] = {}
        label_pos: Dict[str, int] = {}
        label_neg: Dict[str, int] = {}

        residuals = np.array(base_result.residual_vector, dtype=float)
        fitted = np.array(fit_data.x_v, dtype=float) - residuals

        for method in use_methods:
            for _ in range(n_resamples):
                if method == "monte_carlo":
                    sigma_u = (input_uncertainty_pct / 100.0) * np.array(
                        fit_data.x_u, dtype=float
                    )
                    sigma_v = (input_uncertainty_pct / 100.0) * np.array(
                        fit_data.x_v, dtype=float
                    )
                    sigma_u = np.maximum(sigma_u, 1e-6)
                    sigma_v = np.maximum(sigma_v, 1e-6)
                    x_u_star = rng.normal(fit_data.x_u, sigma_u)
                    x_v_star = rng.normal(fit_data.x_v, sigma_v)
                    x_u_star = np.maximum(x_u_star, 0.0)
                    x_v_star = np.maximum(x_v_star, 0.0)
                else:
                    resample_idx = rng.integers(0, len(residuals), len(residuals))
                    resampled = residuals[resample_idx]
                    x_u_star = np.array(fit_data.x_u, dtype=float)
                    x_v_star = fitted + resampled

                result = _run_fit(
                    fit_data,
                    config_edge,
                    x_u_star.tolist(),
                    x_v_star.tolist(),
                )
                if result is None:
                    continue
                if not result.z_labels:
                    continue

                selected, pos_counts, neg_counts = _selection_summary(
                    result.z_labels, result.z_extents, selection_threshold
                )
                base_sets.append(selected)
                base_sizes.append(len(selected))
                base_r2.append(result.chemistry_r2 or 0.0)

                for label in result.z_labels:
                    if label in selected:
                        label_counts[label] = label_counts.get(label, 0) + 1
                    if label in pos_counts:
                        label_pos[label] = label_pos.get(label, 0) + pos_counts[label]
                    if label in neg_counts:
                        label_neg[label] = label_neg.get(label, 0) + neg_counts[label]
                total_samples += 1

        if not base_sets:
            continue

        edge_jaccard = _mean_jaccard(base_sets)
        edge_size = float(np.mean(base_sizes)) if base_sizes else 0.0
        edge_r2_val = float(np.median(base_r2)) if base_r2 else 0.0
        edge_jaccards.append(edge_jaccard)
        edge_sizes.append(edge_size)
        edge_r2.append(edge_r2_val)

        edge_selection_prob: Dict[str, float] = {}
        edge_sign_prob: Dict[str, float] = {}
        for label, count in label_counts.items():
            edge_selection_prob[label] = count / max(1, len(base_sets))
            pos = label_pos.get(label, 0)
            neg = label_neg.get(label, 0)
            total = pos + neg
            edge_sign_prob[label] = (pos / total) if total else 0.0

        stable_reactions = [
            label
            for label, prob in edge_selection_prob.items()
            if prob >= stable_prob_threshold
        ]

        edge_summaries.append(
            EdgeStabilitySummary(
                edge_id=fit_data.edge.edge_id,
                u=fit_data.edge.u,
                v=fit_data.edge.v,
                mean_jaccard=edge_jaccard,
                mean_selected_reactions=edge_size,
                median_r2=edge_r2_val,
                selection_prob=edge_selection_prob,
                sign_prob=edge_sign_prob,
                stable_reactions=stable_reactions,
            )
        )

        for label, count in label_counts.items():
            selection_counts[label] = selection_counts.get(label, 0) + count
        for label, count in label_pos.items():
            sign_counts[label] = sign_counts.get(label, 0) + count

    stability_mean = float(np.mean(edge_jaccards)) if edge_jaccards else 0.0
    stability_se = (
        float(np.std(edge_jaccards, ddof=1)) / np.sqrt(len(edge_jaccards))
        if len(edge_jaccards) > 1
        else 0.0
    )
    mean_selected_reactions = (
        float(np.mean(edge_sizes)) if edge_sizes else 0.0
    )
    median_r2 = float(np.median(edge_r2)) if edge_r2 else 0.0

    selection_prob: Dict[str, float] = {}
    sign_prob: Dict[str, float] = {}
    for label, count in selection_counts.items():
        selection_prob[label] = count / max(1, total_samples)
        sign_prob[label] = sign_counts.get(label, 0) / max(1, count)

    return LambdaStabilitySummary(
        lambda_l1=float(lambda_l1),
        stability_mean=stability_mean,
        stability_se=stability_se,
        mean_selected_reactions=mean_selected_reactions,
        median_r2=median_r2,
        edges_considered=len(edge_jaccards),
        selection_prob=selection_prob,
        sign_prob=sign_prob,
        edge_summaries=edge_summaries,
    )


def tune_reaction_hyperparameters(
    samples: object,
    edges: Iterable[EdgeInput],
    config: Config,
    lambda_grid: Sequence[float],
    n_resamples: int = 200,
    perturbation: str = "monte_carlo",
    input_uncertainty_pct: float = 5.0,
    selection_threshold: float = 1e-6,
    stable_prob_threshold: float = 0.8,
    min_r2: Optional[float] = None,
    max_reactions: Optional[float] = None,
    random_state: Optional[int] = None,
    phreeqc_results: Optional[Mapping[str, Mapping[str, object]]] = None,
    residence_time_overrides: Optional[Mapping[str, float]] = None,
) -> TuningReport:
    grid = [float(value) for value in lambda_grid if float(value) > 0]
    if not grid:
        raise ValueError("lambda_grid must include positive values.")

    rng = np.random.default_rng(random_state)
    edge_data = _build_edge_data(
        samples,
        edges,
        config,
        phreeqc_results,
        residence_time_overrides,
    )
    if not edge_data:
        raise ValueError("No edges with complete data available for tuning.")

    candidates: List[LambdaStabilitySummary] = []
    for lambda_l1 in sorted(grid):
        summary = _evaluate_lambda(
            edge_data,
            lambda_l1,
            n_resamples,
            perturbation,
            input_uncertainty_pct,
            selection_threshold,
            rng,
            stable_prob_threshold,
        )
        candidates.append(summary)

    eligible = candidates
    if min_r2 is not None:
        eligible = [cand for cand in eligible if cand.median_r2 >= min_r2]
    if max_reactions is not None:
        eligible = [
            cand
            for cand in eligible
            if cand.mean_selected_reactions <= max_reactions
        ]
    if not eligible:
        eligible = candidates

    best = max(eligible, key=lambda item: item.stability_mean)
    threshold = best.stability_mean - best.stability_se
    near_best = [
        cand for cand in eligible if cand.stability_mean >= threshold
    ]
    chosen = max(near_best, key=lambda item: item.lambda_l1)

    return TuningReport(
        recommended_lambda_l1=chosen.lambda_l1,
        selection_rule="1se_stability",
        candidates=candidates,
        min_r2=min_r2,
        max_reactions=max_reactions,
    )
