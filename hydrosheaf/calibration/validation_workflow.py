"""
Assumption-calibrated topology validation workflow.

Provides an end-to-end runner that calibrates topological assumptions on
one label set, applies the fitted settings to a separate held-out validation
set, and emits manuscript-safe metrics without data leakage.
"""

from __future__ import annotations

import json
import os
from dataclasses import replace
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional

import pandas as pd

from ..config import Config as HConfig
from ..graph.types import Edge
from ..log import get_logger
from .adapters import TopologyCalibrationAdapter, TopologyCalibrationObservation
from .factory import _map_float, _map_int, setup_topology_adapter
from .glm import PESTGLM
from .pestpp.runner import run_pestpp

logger = get_logger("calibration.validation_workflow")

# ── helpers ──────────────────────────────────────────────────────────

CONFIG_ONLY_THRESHOLDS = frozenset({
    "evidence_threshold_probable",
    "evidence_threshold_validated",
})


def _resolve_path(file_path: Any) -> str:
    """Normalise a file path to an absolute string for comparison."""
    if file_path is None:
        return ""
    return os.path.abspath(str(file_path))


def _load_topology_observations(
    obs_file: str,
) -> List[TopologyCalibrationObservation]:
    """Load topology observations from a CSV file."""
    df = pd.read_csv(obs_file)
    observations: List[TopologyCalibrationObservation] = []
    for _, row in df.iterrows():
        edge_id = str(row.get("edge_id"))
        if not edge_id or edge_id == "nan":
            continue
        value = row.get("observed_present", row.get("present", row.get("value")))
        if pd.isna(value):
            continue
        observations.append(
            TopologyCalibrationObservation(
                edge_id=edge_id,
                observed_present=float(value),
                weight=float(row.get("weight", 1.0)),
            )
        )
    return observations


def _extract_assumption_params_from_optimal(
    optimal_params: Dict[str, float],
    assumption_param_names: List[str],
) -> Dict[str, float]:
    """Extract only assumption parameters from the optimal parameters dict."""
    return {
        name: float(optimal_params[name])
        for name in assumption_param_names
        if name in optimal_params
    }


def _extract_config_only_thresholds(cfg: HConfig) -> Dict[str, float]:
    """Extract evidence thresholds that are config-only, not calibratable."""
    result: Dict[str, float] = {}
    for attr_name in CONFIG_ONLY_THRESHOLDS:
        if hasattr(cfg, attr_name):
            result[attr_name] = float(getattr(cfg, attr_name))
    return result


def _edge_id_set(observations: List[TopologyCalibrationObservation]) -> set:
    """Extract the set of edge IDs from a list of observations."""
    return {obs.edge_id for obs in observations}


def _independence_status(
    cal_obs_file: str,
    val_obs_file: str,
    settings: Mapping[str, Any],
    overlapping: set,
) -> tuple[bool, str]:
    """Determine whether validation is independent for the primary claim.

    Different filenames are not sufficient evidence of independence: the
    same edge ID in both files is a leaked labelled unit, even when the two
    labels disagree.  Overlapping IDs may only be accepted when both files
    carry explicit, disjoint group assignments under the following contract::

        {
            "kind": "grouped",
            "group_column": "split_group",
            "calibration_groups": ["calibration"],
            "validation_groups": ["validation"],
            "allow_overlapping_edge_ids": True,
        }

    The contract is checked against the source CSVs rather than trusted from
    filenames or from the observed labels themselves.
    """
    if not overlapping:
        return (
            True,
            "Calibration and validation edge IDs are disjoint; filenames "
            "are not used as the independence criterion.",
        )

    contract = settings.get("grouped_independence_contract")
    if not isinstance(contract, Mapping):
        return (
            False,
            "independent_validation=False: calibration and validation contain "
            f"{len(overlapping)} overlapping edge IDs. Different filenames "
            "(including reversed labels) do not establish independent "
            "validation; provide an explicit grouped_independence_contract "
            "with disjoint group assignments to authorize this design.",
        )

    if contract.get("kind") != "grouped":
        return (
            False,
            "independent_validation=False: overlapping edge IDs require a "
            "grouped_independence_contract with kind='grouped'.",
        )
    if contract.get("allow_overlapping_edge_ids") is not True:
        return (
            False,
            "independent_validation=False: overlapping edge IDs require "
            "grouped_independence_contract.allow_overlapping_edge_ids=True.",
        )

    group_column = contract.get("group_column")
    if not isinstance(group_column, str) or not group_column.strip():
        return (
            False,
            "independent_validation=False: grouped independence requires a "
            "non-empty group_column.",
        )

    def _group_set(value: Any) -> set[str]:
        if isinstance(value, str):
            values = [value]
        elif isinstance(value, (list, tuple, set, frozenset)):
            values = value
        else:
            values = []
        return {str(item).strip() for item in values if str(item).strip()}

    calibration_groups = _group_set(contract.get("calibration_groups"))
    validation_groups = _group_set(contract.get("validation_groups"))
    if not calibration_groups or not validation_groups:
        return (
            False,
            "independent_validation=False: grouped independence requires "
            "non-empty calibration_groups and validation_groups.",
        )
    if calibration_groups & validation_groups:
        return (
            False,
            "independent_validation=False: calibration_groups and "
            "validation_groups must be disjoint.",
        )

    try:
        cal_df = pd.read_csv(cal_obs_file)
        val_df = pd.read_csv(val_obs_file)
    except Exception as exc:
        return (
            False,
            "independent_validation=False: could not verify the grouped "
            f"independence contract from the label files ({exc}).",
        )

    if group_column not in cal_df.columns or group_column not in val_df.columns:
        return (
            False,
            "independent_validation=False: grouped independence requires "
            f"'{group_column}' in both label files.",
        )

    def _observed_groups(df: pd.DataFrame) -> set[str]:
        return {
            str(value).strip()
            for value in df[group_column].dropna().tolist()
            if str(value).strip()
        }

    observed_calibration_groups = _observed_groups(cal_df)
    observed_validation_groups = _observed_groups(val_df)
    if not observed_calibration_groups or not observed_validation_groups:
        return (
            False,
            "independent_validation=False: grouped independence requires a "
            f"non-empty '{group_column}' assignment in both label files.",
        )

    unknown_calibration_groups = (
        observed_calibration_groups - calibration_groups
    )
    unknown_validation_groups = (
        observed_validation_groups - validation_groups
    )
    if unknown_calibration_groups or unknown_validation_groups:
        return (
            False,
            "independent_validation=False: grouped independence contract does "
            "not cover all observed group assignments.",
        )

    return (
        True,
        "Calibration and validation share "
        f"{len(overlapping)} edge IDs, but independence is authorized by an "
        "explicit grouped_independence_contract with disjoint '"
        f"{group_column}' assignments.",
    )


def _make_edge(
    edge_id: str, u: str = "", v: str = ""
) -> Edge:
    """Build an Edge dataclass from an edge_id string."""
    if "->" in edge_id and (not u or not v):
        parts = edge_id.split("->", 1)
        u = u or parts[0]
        v = v or parts[1]
    return Edge(edge_id=edge_id, u=u or "", v=v or "")


def _build_selected_edge_set(
    adapter: TopologyCalibrationAdapter,
    optimal_params: Dict[str, float],
) -> set:
    """Get the set of selected edge IDs after calibration."""
    try:
        selected = adapter.selected_edges(optimal_params, threshold=0.5)
        return {str(edge.edge_id) for edge in selected}
    except Exception:
        # Fallback: use edge_probabilities with threshold
        try:
            probs = adapter.edge_probabilities(optimal_params)
            return {
                edge_id for edge_id, prob in probs.items()
                if prob >= 0.5
            }
        except Exception:
            return set()


# ── main validation workflow ─────────────────────────────────────────

def run_assumption_calibration_validation(
    config: Any,
    calibration_result: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Run assumption calibration on one label set then validate on another.

    Parameters
    ----------
    config : CalibrationConfig
        Must be a topology calibration config whose ``adapter_settings``
        includes ``validation_observations_file``.
    calibration_result : dict, optional
        Pre-computed calibration result from the CLI. When provided,
        calibration is skipped and this result's optimal parameters are
        used directly. This avoids a redundant second calibration run
        and guarantees that ``results.json`` and the validation report
        were produced from the same parameter set.

    Returns
    -------
    dict
        Validation results dictionary (also written to
        ``assumption_validation_results.json`` in ``output_dir``).

    Raises
    ------
    ValueError
        If calibration and validation label files are the same path.
    """
    settings = config.adapter_settings

    # ── 1. Load label files and enforce independence ──────────────────
    cal_obs_file = (
        settings.get("observations_file")
        or config.observations_file
    )
    val_obs_file = settings.get("validation_observations_file")
    if not cal_obs_file:
        raise ValueError(
            "Topology calibration requires observations_file "
            "(calibration labels)."
        )
    if not val_obs_file:
        raise ValueError(
            "Validation workflow requires "
            "model.validation_observations_file."
        )

    cal_path = _resolve_path(cal_obs_file)
    val_path = _resolve_path(val_obs_file)
    if cal_path == val_path:
        raise ValueError(
            "Calibration and validation label files must be different. "
            f"Both resolved to: {cal_path}"
        )

    cal_observations = _load_topology_observations(cal_obs_file)
    val_observations = _load_topology_observations(val_obs_file)
    cal_edge_ids = _edge_id_set(cal_observations)
    val_edge_ids = _edge_id_set(val_observations)
    overlapping = cal_edge_ids & val_edge_ids
    independent_validation, independence_reason = _independence_status(
        cal_obs_file,
        val_obs_file,
        settings,
        overlapping,
    )

    logger.info(
        "Calibration labels: %d, Validation labels: %d, Overlapping edge IDs: %d",
        len(cal_observations), len(val_observations), len(overlapping),
    )

    if not cal_observations:
        raise ValueError("No calibration observations loaded.")
    if not val_observations:
        raise ValueError("No validation observations loaded.")

    # ── 2. Build topology adapter with CALIBRATION labels only ────────
    problem = setup_topology_adapter(config)
    logger.info(
        "Topology adapter built with %d candidate edges, %d calibration obs.",
        len(problem.candidate_edges), len(problem.observations),
    )

    assumption_params = settings.get("assumption_params") or []

    # ── 3. Run calibration or reuse pre-computed result ───────────────
    if calibration_result is not None:
        logger.info(
            "Reusing calibration result from CLI (skipping redundant run)."
        )
        result = calibration_result
    elif config.engine.startswith("pestpp"):
        pestpp_opts = dict(config.pestpp_options)
        if config.engine == "pestpp-ies":
            pestpp_opts.update(config.ies_settings)
        elif config.engine == "pestpp-sen":
            pestpp_opts.update(config.sen_settings)
        elif config.engine == "pestpp-swp":
            pestpp_opts.update(config.swp_settings)
        elif config.engine in ("pestpp-mou", "pestpp-opt"):
            pestpp_opts.update(config.opt_settings)
        elif config.engine == "pestpp-da":
            pestpp_opts.update(config.da_settings)
        result = run_pestpp(
            problem=problem,
            engine=config.engine,
            work_dir=config.work_dir or "pest_workspace",
            case_name="assumption_calibration",
            max_nfev=config.max_nfev,
            n_workers=config.n_workers,
            pestpp_options=pestpp_opts,
            pestpp_version=config.pestpp_version,
        )
    else:
        from .cli import _resolve_internal_parameters
        pest_params = _resolve_internal_parameters(
            problem, config.parameters, logger,
        )
        pest = PESTGLM(
            parameters=pest_params,
            observations=problem.get_observations(),
            model_runner=problem.run_model,
            n_workers=config.n_workers,
            worker_type="thread",
            loss=config.loss,
        )
        result = pest.calibrate(max_nfev=config.max_nfev)

    if not result.get("success", False):
        logger.warning(
            "Calibration did not converge. Validation will still be "
            "attempted with the best-fit parameters found."
        )

    optimal_params: Dict[str, float] = result.get("optimal_parameters", {})
    calibrated_assumption_params = _extract_assumption_params_from_optimal(
        optimal_params, assumption_params,
    )

    # ── 4. Evaluate on validation labels ──────────────────────────────
    # Build a Config with calibrated assumption params applied
    from dataclasses import replace
    if problem.config is not None and isinstance(problem.config, HConfig):
        cal_cfg = replace(problem.config)
    else:
        cal_cfg = HConfig()

    for a_name in assumption_params:
        if a_name in optimal_params and hasattr(cal_cfg, a_name):
            setattr(cal_cfg, a_name, float(optimal_params[a_name]))

    if assumption_params:
        cal_cfg.assumption_calibration_enabled = True
        cal_cfg.null_model_enabled = True
        cal_cfg.evidence_ladder_enabled = True

    # Extract config-only evidence thresholds (not calibratable)
    config_only_thresholds = _extract_config_only_thresholds(cal_cfg)

    # Build the calibrated-edge set via sheaf refinement
    from ..sheaf.topology_refine import refine_edges_with_sheaf

    # Build samples list from the adapter
    samples = problem.samples

    # Apply edge probabilities from optimal params to candidate edges
    calibrated_candidates = []
    for edge in problem.candidate_edges:
        p_name = TopologyCalibrationAdapter._param_name(edge.edge_id)
        if p_name in optimal_params:
            prob = TopologyCalibrationAdapter._sigmoid(
                float(optimal_params[p_name])
            )
            new_attrs = dict(edge.attrs or {})
            new_attrs["edge_confidence"] = prob
            calibrated_candidates.append(
                Edge(
                    edge_id=edge.edge_id,
                    u=edge.u,
                    v=edge.v,
                    attrs=new_attrs,
                )
            )
        else:
            calibrated_candidates.append(edge)

    # Run sheaf refinement with calibrated config
    if samples is not None and assumption_params:
        selected_edges = refine_edges_with_sheaf(
            samples, calibrated_candidates, cal_cfg,
        )
    else:
        # Fallback: simple threshold on edge_probabilities
        probs = problem.edge_probabilities(optimal_params)
        selected_edges = [
            edge for edge in problem.candidate_edges
            if probs.get(str(edge.edge_id), 0.0) >= 0.5
        ]

    selected_ids = {str(e.edge_id) for e in selected_edges}

    # ── 5. Compute validation metrics on validation labels only ──────
    # Metrics are restricted to edges that appear in the
    # validation_observations_file. Model-selected edges that are not
    # explicitly labelled are excluded from the evaluation so they
    # cannot be mis-counted as false positives or true negatives.

    tp = 0
    fp = 0
    tn = 0
    fn = 0

    for obs in val_observations:
        is_selected = obs.edge_id in selected_ids
        if obs.observed_present >= 0.5:
            if is_selected:
                tp += 1
            else:
                fn += 1
        else:
            if is_selected:
                fp += 1
            else:
                tn += 1

    total = tp + fp + tn + fn
    precision = tp / (tp + fp) if tp + fp else 0.0
    recall = tp / (tp + fn) if tp + fn else 0.0
    f1 = (2.0 * precision * recall / (precision + recall)
          if precision + recall else 0.0)
    accuracy = (tp + tn) / total if total > 0 else float("nan")

    validation_confusion_matrix = {
        "true_positives": tp,
        "false_positives": fp,
        "true_negatives": tn,
        "false_negatives": fn,
    }

    validation_metrics = {
        "precision": precision,
        "recall": recall,
        "f1": f1,
        "accuracy": accuracy,
    }

    # ── 6. Emit manuscript-safe report JSON ───────────────────────────
    report = {
        "calibrated_assumption_parameters": calibrated_assumption_params,
        "fixed_config_only_thresholds": config_only_thresholds,
        "selected_edges_all": sorted(list(selected_ids)),
        "evaluated_validation_edge_ids": sorted(list(val_edge_ids)),
        "selected_edge_details": [
            {"edge_id": e.edge_id, "u": e.u, "v": e.v}
            for e in selected_edges
        ],
        "validation_confusion_matrix": validation_confusion_matrix,
        "validation_metrics": validation_metrics,
        "calibration_label_file": cal_obs_file,
        "validation_label_file": val_obs_file,
        "independent_validation": independent_validation,
        "independence_reason": independence_reason,
        "n_calibration_labels": len(cal_observations),
        "n_validation_labels": len(val_observations),
        "n_overlapping_edge_ids": len(overlapping),
        "calibration_phi": result.get("phi", float("nan")),
        "calibration_success": result.get("success", False),
    }

    # Write report
    os.makedirs(config.output_dir, exist_ok=True)
    report_path = Path(config.output_dir) / "assumption_validation_results.json"

    import numpy as np

    def _convert(obj: Any) -> Any:
        if isinstance(obj, (np.integer,)):
            return int(obj)
        if isinstance(obj, (np.floating,)):
            return float(obj)
        raise TypeError

    with open(report_path, "w") as f:
        json.dump(report, f, indent=2, default=_convert)

    logger.info(
        "Validation complete. Report saved to %s", report_path,
    )
    logger.info(
        "Validation metrics — precision=%.4f recall=%.4f f1=%.4f accuracy=%.4f",
        precision, recall, f1, accuracy,
    )

    return report
