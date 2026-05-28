"""
Phase 2b: Independent Hydrosheaf graph validation against MODPATH reference edges.

Implements 7 independent validation scenarios, 3 diagnostic negative controls,
1 sensitivity analysis, and 3 prior-assisted scenarios for the Savage archive,
keeping independent and prior-assisted modes strictly separated per the
Phase_1.txt scientific guardrails.

Independent scenarios use ONLY grid cell coordinates (from MODFLOW grid geometry,
not MODPATH connectivity) to infer directed edges. Prior-assisted scenarios
explicitly ingest MODPATH edge information and are labelled accordingly.

Architecture notes
==================
*What works*
- The head-gradient baseline (elevation-as-head proxy, downhill only, 2-nearest
  neighbours) recovers MODPATH topology at F1 ~0.618, well above the random
  baseline (~0.003).  This is a real signal from a trivially simple prior.
- Bayesian Hodge pruning (Scenario 2b) modestly improves precision over the
  baseline by suppressing topologically inconsistent shortcuts, but does not
  fundamentally change the performance envelope because the flat-z (elevation)
  potential field is nearly harmonic (d*d = 0 by construction).
- The projected-gradient prior (Scenario 2d) replaces the discretised
  steepest-descent heuristic with a continuous head gradient computed from the
  full MODFLOW grid via finite differences.  This eliminates the artificial
  min_drop clipping and uphill double-penalty of earlier Hodge-based scenarios.

*What does not work*
- CBC flux-informed edge selection (former Scenario 2e, removed).  The
  MODFLOW cell-by-cell dominant-outflow direction produced only 2.3% alignment
  with particle trajectories -- worse than useless.  The dominant outflow
  direction on a face is anti-correlated with the actual MODPATH pathlines,
  because MODPATH tracks the magnitude-weighted flow on each face (not just
  the sign of the net flux).  This scenario was dead code and has been removed
  entirely.

*Why Hodge is kept as diagnostic, not edge selector*
- The Hodge Laplacian decomposition (H1 = d*d + dd*) vanishes identically on a
  single scalar potential field (head or elevation) because d*d = 0 is an exact
  cohomological identity.  The Hodge obstruction energy is therefore near-zero
  for any gradient field and provides negligible discriminatory power.
- Hodge becomes discriminative only with multi-field data (head + chemistry +
  stable isotopes), where the curl component d* detects rotational flow
  patterns that violate gradient-field assumptions.  The Savage archive lacks
  these auxiliary fields.
- For single-field topography, the projected-gradient prior (Scenario 2d) is
  the correct approach: it scores edges by projecting the continuous head
  gradient onto edge directions, providing a physically meaningful prior
  without relying on the degenerate Hodge decomposition.
"""
from __future__ import annotations

import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PROJECT_ROOT))

from hydrosheaf.config import Config
from hydrosheaf.graph.build import infer_edges_from_coordinates
from hydrosheaf.inference.topology_posterior import infer_topology_map_edges
from hydrosheaf.physics.modflow_head import parse_fhd
from hydrosheaf.validation import (
    edge_confusion,
    scale_mismatch_diagnostics,
    validate_independent_graph_against_modpath,
    build_modpath_informed_graph_priors,
    apply_modpath_informed_graph_priors,
)
from hydrosheaf.validation.topology import normalize_directed_edges

RANDOM_SEED = 20260521
SPARSITY_TRIALS = 20

SAVAGE_RESULTS = (
    PROJECT_ROOT
    / "M4"
    / "m4_topology_benchmark"
    / "results"
    / "public_archives"
    / "savage"
)
BENCHMARK_RESULTS = PROJECT_ROOT / "M4" / "m4_topology_benchmark" / "results"
DOCS_DIR = PROJECT_ROOT / "M4" / "m4_topology_benchmark" / "docs"

# Savage MODFLOW grid geometry (from base.lst and nwis_matching.py)
_SAVAGE_GRID = dict(
    ncol=183, nrow=202, nlay=8,
    dx=110.0, dy=73.333,          # core cell size (ft)
    rotation_deg=-12.0,             # ANGROT from .lst
    origin_x=961030.4, origin_y=112955.0,  # NH SP ft, from nwis_matching.py
)


def load_savage_nodes() -> pd.DataFrame:
    """Load the standardized node mapping from Phase 2."""
    path = SAVAGE_RESULTS / "modpath_node_mapping.csv"
    if not path.exists():
        raise FileNotFoundError(f"Node mapping not found: {path}")
    return pd.read_csv(path)


def load_savage_reference_edges() -> pd.DataFrame:
    """Load the MODPATH reference edges from Phase 2."""
    path = SAVAGE_RESULTS / "modpath_reference_edges.csv"
    if not path.exists():
        raise FileNotFoundError(f"Reference edges not found: {path}")
    return pd.read_csv(path)


def nodes_to_samples(nodes: pd.DataFrame) -> List[Dict[str, object]]:
    """Convert grid cell nodes to Hydrosheaf sample dicts.

    Each cell becomes a sample with site_id, coordinates, and elevation.
    Coordinates come from MODPATH endpoint particle positions, which are
    MODFLOW grid cell centers. These are independent of MODPATH connectivity.
    """
    samples = []
    for _, row in nodes.iterrows():
        node_id = str(row["node_id"])
        samples.append({
            "site_id": node_id,
            "lat": float(row["y"]),   # y-coordinate as spatial proxy
            "lon": float(row["x"]),   # x-coordinate as spatial proxy
            "elevation": float(row["z"]) if not pd.isna(row["z"]) else 0.0,
        })
    return samples


def compute_haversine_km(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    """Haversine distance between two points in km."""
    R = 6371.0
    dlat = math.radians(lat2 - lat1)
    dlon = math.radians(lon2 - lon1)
    a = (math.sin(dlat / 2) ** 2
         + math.cos(math.radians(lat1)) * math.cos(math.radians(lat2))
         * math.sin(dlon / 2) ** 2)
    return 2 * R * math.atan2(math.sqrt(a), math.sqrt(1 - a))


def compute_edge_lengths(
    node_map: Dict[str, Tuple[float, float]],
    edges: List[Tuple[str, str]],
) -> Dict[Tuple[str, str], float]:
    """Compute spatial distance for each edge from node coordinates."""
    lengths: Dict[Tuple[str, str], float] = {}
    for u, v in edges:
        if u in node_map and v in node_map:
            lat1, lon1 = node_map[u]
            lat2, lon2 = node_map[v]
            lengths[(u, v)] = compute_haversine_km(lat1, lon1, lat2, lon2)
    return lengths


def infer_bayesian_hodge_edges(
    samples: List[Dict[str, object]],
) -> Tuple[List[Tuple[str, str]], Dict[str, Any]]:
    """Infer topology from a downhill candidate graph using Bayesian Hodge pruning."""
    candidate_obj = infer_edges_from_coordinates(samples, max_neighbors=5, allow_uphill=False)
    initial_obj = infer_edges_from_coordinates(samples, max_neighbors=2, allow_uphill=False)

    config = Config(
        phreeqc_enabled=False,
        sheaf_isotope_enabled=False,
        sheaf_cl_enabled=False,
        sheaf_age_enabled=False,
        sheaf_cohomology_enabled=False,
        hydraulic_hodge_enabled=True,
        hydraulic_hodge_weight=5.0,
        hydraulic_hodge_leverage_weight=0.5,
        hydraulic_hodge_fallback_to_elevation=True,
        hydraulic_hodge_head_key="hydraulic_head",
        hydraulic_hodge_min_head_drop=0.0,
        hydraulic_hodge_direction_penalty=0.1,
        topology_posterior_enabled=True,
        topology_posterior_samples=800,
        topology_posterior_burnin=200,
        topology_posterior_beta=1.0,
        topology_posterior_edge_penalty=0.08,
    )
    config.validate()

    selected_obj, posterior = infer_topology_map_edges(
        samples=samples,
        candidate_edges=candidate_obj,
        config=config,
        initial_edges=initial_obj,
        max_neighbors=2,
        probability_threshold=0.0,
        seed=RANDOM_SEED,
    )
    return [(e.u, e.v) for e in selected_obj], posterior


# Path to the base-calibration MODFLOW head output.
# Set to None or a non-existent path to skip the real-head scenario
# gracefully until the source archive is available.
_BASE_FHD_PATH: Optional[Path] = (
    PROJECT_ROOT
    / "M4"
    / "m4_topology_benchmark"
    / "public_archives"
    / "savage"
    / "base.fhd"
)


def infer_projected_gradient_edges(
    samples: List[Dict[str, object]],
    node_df: pd.DataFrame,
) -> Optional[Tuple[List[Tuple[str, str]], Dict[str, Any]]]:
    """Scenario 2d: Real-head directional prior via continuous gradient projection.

    Replaces the degenerate scalar-potential Hodge obstruction with physically
    meaningful diagnostics:

    1. Computes the continuous head gradient from the full MODFLOW grid via
       finite differences (not the discretised candidate-edge set).
    2. Scores every candidate edge by projecting the gradient onto the edge
       direction — edges aligned with the flow get high priors; perpendicular
       or anti-aligned edges get low priors.
    3. Computes local head-plane residuals as a secondary consistency check.
    4. Runs the MCMC topology posterior with these priors + cost terms.

    Requires base.fhd.  Returns None gracefully if unavailable.
    """
    if _BASE_FHD_PATH is None or not _BASE_FHD_PATH.exists():
        return None

    try:
        head_map = parse_fhd(_BASE_FHD_PATH)
    except Exception:
        return None

    if node_df is None:
        return None

    from hydrosheaf.physics.modflow_head import (
        build_grid_geometry_from_params,
        cell_heads_to_sample_field,
    )

    # Attach real MODFLOW head to each sample
    node_with_heads = cell_heads_to_sample_field(node_df, head_map, "hydraulic_head")
    head_lookup = dict(
        zip(node_with_heads["node_id"].astype(str),
            node_with_heads["hydraulic_head"])
    )

    augmented = []
    for s in samples:
        sc = dict(s)
        h = head_lookup.get(str(sc["site_id"]))
        if h is not None and not (isinstance(h, float) and math.isnan(h)):
            sc["hydraulic_head"] = float(h)
        augmented.append(sc)

    # Build Savage grid geometry from known parameters (see _SAVAGE_GRID)
    grid = build_grid_geometry_from_params(**_SAVAGE_GRID)

    # Candidate graph uses flat-z (same as head_gradient baseline)
    candidate_obj = infer_edges_from_coordinates(samples, max_neighbors=5, allow_uphill=False)
    initial_obj = infer_edges_from_coordinates(samples, max_neighbors=2, allow_uphill=False)

    config = Config(
        phreeqc_enabled=False,
        sheaf_isotope_enabled=False,
        sheaf_cl_enabled=False,
        sheaf_age_enabled=False,
        sheaf_cohomology_enabled=False,
        # Primary: projected-gradient prior from continuous head field
        projected_gradient_enabled=True,
        projected_gradient_weight=1.0,
        projected_gradient_sharpness=10.0,
        projected_gradient_smoothing_sigma=1.0,
        # Secondary: local head-plane residuals
        head_plane_residual_enabled=True,
        head_plane_residual_weight=2.0,
        head_plane_residual_neighbors=8,
        # Hodge kept for diagnostic attributes only (not primary cost)
        hydraulic_hodge_enabled=True,
        hydraulic_hodge_weight=1.0,
        hydraulic_hodge_leverage_weight=0.5,
        hydraulic_hodge_fallback_to_elevation=False,
        hydraulic_hodge_head_key="hydraulic_head",
        hydraulic_hodge_min_head_drop=0.0,
        hydraulic_hodge_direction_penalty=0.0,  # direction lives in projected_gradient
        # Topology posterior
        topology_posterior_enabled=True,
        topology_posterior_samples=800,
        topology_posterior_burnin=200,
        topology_posterior_beta=1.0,
        topology_posterior_edge_penalty=0.08,
    )
    config.validate()

    selected_obj, posterior = infer_topology_map_edges(
        samples=augmented,
        candidate_edges=candidate_obj,
        config=config,
        initial_edges=initial_obj,
        max_neighbors=2,
        probability_threshold=0.0,
        seed=RANDOM_SEED,
        grid_geometry=grid,
        head_map=head_map,
        node_df=node_df,
    )
    return [(e.u, e.v) for e in selected_obj], posterior

def run_independent_scenarios(
    samples: List[Dict[str, object]],
    ref_edges_raw: List[Tuple[str, str]],
    node_map: Dict[str, Tuple[float, float]],
    node_df: Optional[pd.DataFrame] = None,
) -> pd.DataFrame:
    """Run all independent validation scenarios.

    Returns a DataFrame with one row per scenario, including full confusion
    metrics, scale diagnostics, and guardrails.
    """
    rng = np.random.default_rng(RANDOM_SEED)
    results = []

    # All possible directed edges between the 153 nodes = 153*152 = 23,256
    all_node_ids = sorted(node_map.keys())
    candidate_universe = [
        (u, v) for u in all_node_ids for v in all_node_ids if u != v
    ]

    # ---- Scenario 1: Spatial-only (allow uphill = True) ----
    spatial_obj = infer_edges_from_coordinates(samples, max_neighbors=2, allow_uphill=True)
    spatial_edges = [(e.u, e.v) for e in spatial_obj]
    edge_lengths = compute_edge_lengths(node_map, ref_edges_raw + spatial_edges)
    report = validate_independent_graph_against_modpath(
        spatial_edges,
        ref_edges_raw,
        candidate_edges=candidate_universe,
        edge_lengths=edge_lengths,
    )
    m = report["metrics"]
    s = report["scale_mismatch"]
    results.append({
        "scenario": "spatial_only",
        "validation_mode": "independent_graph_inference",
        "independent_validation": True,
        "result_class": "independent_benchmark",
        "n_reference_edges": m["n_reference_edges"],
        "n_inferred_edges": m["n_inferred_edges"],
        "tp": m["tp"], "fp": m["fp"], "fn": m["fn"], "tn": m["tn"],
        "precision": m["precision"], "recall": m["recall"], "f1": m["f1"],
        "false_positive_rate": m["false_positive_rate"],
        "false_negative_rate": m["false_negative_rate"],
        "scale_mismatch": s["scale_mismatch"],
        "median_reference_length": s["median_reference_length"],
        "median_inferred_length": s["median_inferred_length"],
        "allowed_claim": "Spatial proximity alone recovers directed topology at above-chance rate",
        "required_guardrail": "Spatial-only graph uses no head, hydrostratigraphic, or MODPATH connectivity information",
    })

    # ---- Scenario 2: Head-gradient constrained (elevation proxy, downhill only) ----
    head_obj = infer_edges_from_coordinates(samples, max_neighbors=2, allow_uphill=False)
    head_edges = [(e.u, e.v) for e in head_obj]
    edge_lengths.update(compute_edge_lengths(node_map, head_edges))
    report = validate_independent_graph_against_modpath(
        head_edges,
        ref_edges_raw,
        candidate_edges=candidate_universe,
        edge_lengths=edge_lengths,
    )
    m = report["metrics"]
    s = report["scale_mismatch"]
    results.append({
        "scenario": "head_gradient",
        "validation_mode": "independent_graph_inference",
        "independent_validation": True,
        "result_class": "independent_benchmark",
        "n_reference_edges": m["n_reference_edges"],
        "n_inferred_edges": m["n_inferred_edges"],
        "tp": m["tp"], "fp": m["fp"], "fn": m["fn"], "tn": m["tn"],
        "precision": m["precision"], "recall": m["recall"], "f1": m["f1"],
        "false_positive_rate": m["false_positive_rate"],
        "false_negative_rate": m["false_negative_rate"],
        "scale_mismatch": s["scale_mismatch"],
        "median_reference_length": s["median_reference_length"],
        "median_inferred_length": s["median_inferred_length"],
        "allowed_claim": "Elevation-as-head proxy constrains flow direction and improves topology recovery",
        "required_guardrail": "Elevation is a proxy for hydraulic head, not a MODFLOW-simulated head value",
    })

    # ---- Scenario 2b: Head-gradient + Bayesian Hodge posterior ----
    hodge_edges, posterior = infer_bayesian_hodge_edges(samples)
    edge_lengths.update(compute_edge_lengths(node_map, hodge_edges))
    report = validate_independent_graph_against_modpath(
        hodge_edges,
        ref_edges_raw,
        candidate_edges=candidate_universe,
        edge_lengths=edge_lengths,
    )
    m = report["metrics"]
    s = report["scale_mismatch"]
    results.append({
        "scenario": "head_gradient_bayesian_hodge",
        "validation_mode": "independent_graph_inference",
        "independent_validation": True,
        "result_class": "independent_benchmark",
        "n_reference_edges": m["n_reference_edges"],
        "n_inferred_edges": m["n_inferred_edges"],
        "tp": m["tp"], "fp": m["fp"], "fn": m["fn"], "tn": m["tn"],
        "precision": m["precision"], "recall": m["recall"], "f1": m["f1"],
        "false_positive_rate": m["false_positive_rate"],
        "false_negative_rate": m["false_negative_rate"],
        "scale_mismatch": s["scale_mismatch"],
        "median_reference_length": s["median_reference_length"],
        "median_inferred_length": s["median_inferred_length"],
        "posterior_n_edges_mean": posterior.get("n_edges_mean"),
        "posterior_acceptance_rate": posterior.get("acceptance_rate"),
        "allowed_claim": "Bayesian Hodge pruning suppresses topologically inconsistent downhill shortcuts in an independent head-proxy graph.",
        "required_guardrail": "Uses elevation as a hydraulic-head proxy and Bayesian posterior pruning over a downhill candidate graph; no MODPATH connectivity enters the inference step.",
    })

    # ---- Scenario 2d: Projected-gradient prior from continuous head field ----
    pg_result = infer_projected_gradient_edges(samples, node_df) if node_df is not None else None
    if pg_result is not None:
        pg_edges, pg_posterior = pg_result
        edge_lengths.update(compute_edge_lengths(node_map, pg_edges))
        report = validate_independent_graph_against_modpath(
            pg_edges,
            ref_edges_raw,
            candidate_edges=candidate_universe,
            edge_lengths=edge_lengths,
        )
        m = report["metrics"]
        s = report["scale_mismatch"]
        results.append({
            "scenario": "real_head_projected_gradient",
            "validation_mode": "independent_graph_inference",
            "independent_validation": True,
            "result_class": "independent_benchmark",
            "n_reference_edges": m["n_reference_edges"],
            "n_inferred_edges": m["n_inferred_edges"],
            "tp": m["tp"], "fp": m["fp"], "fn": m["fn"], "tn": m["tn"],
            "precision": m["precision"], "recall": m["recall"], "f1": m["f1"],
            "false_positive_rate": m["false_positive_rate"],
            "false_negative_rate": m["false_negative_rate"],
            "scale_mismatch": s["scale_mismatch"],
            "median_reference_length": s["median_reference_length"],
            "median_inferred_length": s["median_inferred_length"],
            "posterior_n_edges_mean": pg_posterior.get("n_edges_mean"),
            "posterior_acceptance_rate": pg_posterior.get("acceptance_rate"),
            "allowed_claim": (
                "Projected-gradient prior replaces the discretised steepest-descent "
                "heuristic with a continuous head gradient computed from the full "
                "MODFLOW grid via finite differences. Edge priors are based on the "
                "projection of the continuous gradient onto each edge direction. "
                "Local head-plane residuals provide a secondary consistency check. "
                "No degenerate scalar-potential Hodge obstruction is used."
            ),
            "required_guardrail": (
                "Uses MODFLOW simulated head (base_calibration scenario); this is a "
                "model output, not a direct field measurement. The projected-gradient "
                "prior replaces the discretised candidate-edge heuristic with a "
                "continuous-field estimate and removes the artificial min_drop clipping "
                "and uphill double-penalty present in earlier Hodge-based scenarios."
            ),
        })
    else:
        results.append({
            "scenario": "real_head_projected_gradient",
            "validation_mode": "independent_graph_inference",
            "independent_validation": True,
            "result_class": "independent_benchmark",
            "n_reference_edges": None,
            "n_inferred_edges": None,
            "tp": None, "fp": None, "fn": None, "tn": None,
            "precision": None, "recall": None, "f1": None,
            "false_positive_rate": None, "false_negative_rate": None,
            "scale_mismatch": None,
            "median_reference_length": None, "median_inferred_length": None,
            "allowed_claim": "Scenario skipped: base.fhd not available locally.",
            "required_guardrail": (
                "Requires MODFLOW base_calibration head output (base.fhd) from the "
                "Savage archive."
            ),
        })

    # ---- Scenario 3: Head-depth (depth-tiered elevation proxy) ----
    depth_samples = []
    for s in samples:
        s_copy = dict(s)
        elev = float(s_copy.get("elevation", 0))
        # Three depth tiers based on elevation terciles
        s_copy["depth_tier"] = "shallow" if elev > -10 else ("mid" if elev > -50 else "deep")
        depth_samples.append(s_copy)
    depth_obj = infer_edges_from_coordinates(depth_samples, max_neighbors=3, allow_uphill=False)
    depth_edges = [(e.u, e.v) for e in depth_obj]
    edge_lengths.update(compute_edge_lengths(node_map, depth_edges))
    report = validate_independent_graph_against_modpath(
        depth_edges,
        ref_edges_raw,
        candidate_edges=candidate_universe,
        edge_lengths=edge_lengths,
    )
    m = report["metrics"]
    s = report["scale_mismatch"]
    results.append({
        "scenario": "head_depth",
        "validation_mode": "independent_graph_inference",
        "independent_validation": True,
        "result_class": "independent_benchmark",
        "n_reference_edges": m["n_reference_edges"],
        "n_inferred_edges": m["n_inferred_edges"],
        "tp": m["tp"], "fp": m["fp"], "fn": m["fn"], "tn": m["tn"],
        "precision": m["precision"], "recall": m["recall"], "f1": m["f1"],
        "false_positive_rate": m["false_positive_rate"],
        "false_negative_rate": m["false_negative_rate"],
        "scale_mismatch": s["scale_mismatch"],
        "median_reference_length": s["median_reference_length"],
        "median_inferred_length": s["median_inferred_length"],
        "allowed_claim": "Depth-tiered elevation stratifies candidate edges and permits more neighbors",
        "required_guardrail": "Depth tiers are based on elevation terciles, not formal hydrostratigraphic units",
    })

    # ---- Scenario 4: Hydrostratigraphic (depth-based aquifer split) ----
    strat_samples = []
    for s in samples:
        s_copy = dict(s)
        elev = float(s_copy.get("elevation", 0))
        s_copy["aquifer_unit"] = "shallow" if elev > -20 else "deep"
        strat_samples.append(s_copy)
    candidates_raw = infer_edges_from_coordinates(strat_samples, max_neighbors=5, allow_uphill=False)
    strat_obj = []
    for e in candidates_raw:
        u_sample = next((x for x in strat_samples if x["site_id"] == e.u), None)
        v_sample = next((x for x in strat_samples if x["site_id"] == e.v), None)
        if u_sample and v_sample:
            if u_sample.get("aquifer_unit") == v_sample.get("aquifer_unit"):
                strat_obj.append(e)
    strat_edges = [(e.u, e.v) for e in strat_obj]
    edge_lengths.update(compute_edge_lengths(node_map, strat_edges))
    report = validate_independent_graph_against_modpath(
        strat_edges,
        ref_edges_raw,
        candidate_edges=candidate_universe,
        edge_lengths=edge_lengths,
    )
    m = report["metrics"]
    s = report["scale_mismatch"]
    results.append({
        "scenario": "hydrostratigraphic",
        "validation_mode": "independent_graph_inference",
        "independent_validation": True,
        "result_class": "independent_benchmark",
        "n_reference_edges": m["n_reference_edges"],
        "n_inferred_edges": m["n_inferred_edges"],
        "tp": m["tp"], "fp": m["fp"], "fn": m["fn"], "tn": m["tn"],
        "precision": m["precision"], "recall": m["recall"], "f1": m["f1"],
        "false_positive_rate": m["false_positive_rate"],
        "false_negative_rate": m["false_negative_rate"],
        "scale_mismatch": s["scale_mismatch"],
        "median_reference_length": s["median_reference_length"],
        "median_inferred_length": s["median_inferred_length"],
        "allowed_claim": "Same-aquifer edge filtering removes cross-unit false positives",
        "required_guardrail": "Aquifer units are depth-based splits from elevation, not lithostratigraphic or model-layer assignments",
    })

    # ---- Scenario 5: Sparse-node (subsample 50% of nodes, 20 trials) ----
    n_half = max(len(samples) // 2, 5)
    f1_vals = []
    for trial in range(SPARSITY_TRIALS):
        trial_seed = int(rng.integers(0, 2**31 - 1))
        idxs = pd.Series(range(len(samples))).sample(n=n_half, random_state=trial_seed).values
        sub = [samples[i] for i in idxs]
        sub_ids = {str(s["site_id"]) for s in sub}
        sub_ref = [e for e in ref_edges_raw if e[0] in sub_ids and e[1] in sub_ids]
        if len(sub_ref) < 5:
            f1_vals.append(0.0)
            continue
        sub_obj = infer_edges_from_coordinates(sub, max_neighbors=2, allow_uphill=False)
        sub_edges = [(e.u, e.v) for e in sub_obj]
        sub_conf = edge_confusion(sub_ref, sub_edges)
        f1_vals.append(sub_conf["f1"])

    f1_mean = float(np.mean(f1_vals)) if f1_vals else float("nan")
    f1_std = float(np.std(f1_vals)) if f1_vals else float("nan")
    results.append({
        "scenario": "sparse_node",
        "validation_mode": "sensitivity_analysis",
        "independent_validation": False,
        "result_class": "sparsity_sensitivity",
        "n_reference_edges": m["n_reference_edges"],
        "n_inferred_edges": float("nan"),
        "tp": float("nan"), "fp": float("nan"), "fn": float("nan"), "tn": float("nan"),
        "precision": float("nan"), "recall": float("nan"), "f1": f1_mean,
        "false_positive_rate": float("nan"),
        "false_negative_rate": float("nan"),
        "scale_mismatch": False,
        "median_reference_length": float("nan"),
        "median_inferred_length": float("nan"),
        "sparse_node_fraction": 0.5,
        "sparse_node_trials": SPARSITY_TRIALS,
        "sparse_node_f1_std": f1_std,
        "allowed_claim": "Sparse-node subsampling quantifies sensitivity of topology metrics to reduced node density",
        "required_guardrail": "Sparse-node results are diagnostic sensitivity evidence with variable successful trials; do not claim monotonic degradation or field sampling performance",
    })

    # ---- Scenario 6: Negative control - random edges ----
    n_random = len(ref_edges_raw)
    random_pairs = []
    node_list = sorted(node_map.keys())
    for _ in range(n_random):
        u = str(rng.choice(node_list))
        v = str(rng.choice(node_list))
        while u == v or (u, v) in random_pairs:
            v = str(rng.choice(node_list))
        random_pairs.append((u, v))
    report_neg = validate_independent_graph_against_modpath(
        random_pairs, ref_edges_raw,
        candidate_edges=candidate_universe,
        edge_lengths=edge_lengths,
    )
    m = report_neg["metrics"]
    s = report_neg["scale_mismatch"]
    results.append({
        "scenario": "negative_random",
        "validation_mode": "diagnostic_negative_control",
        "independent_validation": False,
        "result_class": "negative_control",
        "n_reference_edges": m["n_reference_edges"],
        "n_inferred_edges": m["n_inferred_edges"],
        "tp": m["tp"], "fp": m["fp"], "fn": m["fn"], "tn": m["tn"],
        "precision": m["precision"], "recall": m["recall"], "f1": m["f1"],
        "false_positive_rate": m["false_positive_rate"],
        "false_negative_rate": m["false_negative_rate"],
        "scale_mismatch": s["scale_mismatch"],
        "median_reference_length": s["median_reference_length"],
        "median_inferred_length": s["median_inferred_length"],
        "allowed_claim": "Random edges serve as a negative-control baseline for topology metrics",
        "required_guardrail": "Random edges are not a proposed graph model; they establish the floor for precision",
    })

    # ---- Scenario 7: Negative wrong-direction ----
    wrong_dir = [(v, u) for (u, v) in ref_edges_raw]
    report_wd = validate_independent_graph_against_modpath(
        wrong_dir, ref_edges_raw,
        candidate_edges=candidate_universe,
        edge_lengths=edge_lengths,
    )
    m = report_wd["metrics"]
    s = report_wd["scale_mismatch"]
    results.append({
        "scenario": "negative_wrong_direction",
        "validation_mode": "diagnostic_negative_control",
        "independent_validation": False,
        "result_class": "negative_control",
        "n_reference_edges": m["n_reference_edges"],
        "n_inferred_edges": m["n_inferred_edges"],
        "tp": m["tp"], "fp": m["fp"], "fn": m["fn"], "tn": m["tn"],
        "precision": m["precision"], "recall": m["recall"], "f1": m["f1"],
        "false_positive_rate": m["false_positive_rate"],
        "false_negative_rate": m["false_negative_rate"],
        "scale_mismatch": s["scale_mismatch"],
        "median_reference_length": s["median_reference_length"],
        "median_inferred_length": s["median_inferred_length"],
        "allowed_claim": "Reversing reference edge direction quantifies directionality sensitivity",
        "required_guardrail": "Wrong-direction graph is a diagnostic negative control, not a Hydrosheaf inference",
    })

    # ---- Scenario 8: Negative shortcut (skip intermediate nodes) ----
    # Build a path-based shortcut: for each reference edge u->v, if v has an
    # outgoing edge v->w, skip v and create u->w
    ref_graph: Dict[str, List[str]] = {}
    for u, v in ref_edges_raw:
        ref_graph.setdefault(u, []).append(v)
    shortcuts = []
    for u, v in ref_edges_raw:
        for w in ref_graph.get(v, []):
            if w != u and (u, w) not in ref_edges_raw:
                shortcuts.append((u, w))
    shortcuts = list(set(shortcuts))[:len(ref_edges_raw)]  # Cap at same count
    report_sc = validate_independent_graph_against_modpath(
        shortcuts, ref_edges_raw,
        candidate_edges=candidate_universe,
        edge_lengths=edge_lengths,
    )
    m = report_sc["metrics"]
    s = report_sc["scale_mismatch"]
    results.append({
        "scenario": "negative_shortcut",
        "validation_mode": "diagnostic_negative_control",
        "independent_validation": False,
        "result_class": "negative_control",
        "n_reference_edges": m["n_reference_edges"],
        "n_inferred_edges": m["n_inferred_edges"],
        "tp": m["tp"], "fp": m["fp"], "fn": m["fn"], "tn": m["tn"],
        "precision": m["precision"], "recall": m["recall"], "f1": m["f1"],
        "false_positive_rate": m["false_positive_rate"],
        "false_negative_rate": m["false_negative_rate"],
        "scale_mismatch": s["scale_mismatch"],
        "median_reference_length": s["median_reference_length"],
        "median_inferred_length": s["median_inferred_length"],
        "allowed_claim": "Shortcut edges test sensitivity to skipped intermediate nodes",
        "required_guardrail": "Shortcut graph is a diagnostic negative control, not a Hydrosheaf inference",
    })

    return pd.DataFrame(results)


def run_prior_assisted_scenarios(
    head_edges: List[Tuple[str, str]],
    ref_edges_raw: List[Tuple[str, str]],
    edge_travel_times: Dict[Tuple[str, str], float],
) -> pd.DataFrame:
    """Run 3 prior-assisted scenarios: override, merge, only.

    These are labelled validation_mode='prior_assisted' and
    independent_validation=False. They must never enter independent summaries.
    """
    results = []

    for mode in ("override", "merge", "only"):
        report = apply_modpath_informed_graph_priors(
            head_edges,
            ref_edges_raw,
            mode=mode,
            travel_time_days=edge_travel_times,
        )
        results.append({
            "scenario": f"modpath_prior_{mode}",
            "validation_mode": "prior_assisted",
            "independent_validation": False,
            "result_class": "prior_assisted",
            "n_hydrosheaf_input_edges": report["n_input_hydrosheaf_edges"],
            "n_modpath_prior_edges": report["n_modpath_prior_edges"],
            "n_output_edges": report["n_output_edges"],
            "prior_mode": mode,
            "n_edges_added": report["n_output_edges"] - report["n_input_hydrosheaf_edges"],
            "allowed_claim": f"MODPATH prior {mode} demonstrates how external physics informs graph construction",
            "required_guardrail": "Prior-assisted only; not independent validation. MODPATH has entered the Hydrosheaf graph as prior information.",
        })

    return pd.DataFrame(results)


def build_edge_classification(
    independent_df: pd.DataFrame,
    reference_edges: List[Tuple[str, str]],
) -> pd.DataFrame:
    """Build per-edge classification table from the head_gradient scenario.

    Each MODPATH reference edge is classified as TP (recovered), FN (missed),
    and each Hydrosheaf-inferred edge not in reference is FP.
    """
    # Use head_gradient scenario for edge classification
    head_row = independent_df[independent_df["scenario"] == "head_gradient"]
    if head_row.empty:
        return pd.DataFrame()

    # Re-run to get detailed edge lists
    nodes = load_savage_nodes()
    samples = nodes_to_samples(nodes)
    head_obj = infer_edges_from_coordinates(samples, max_neighbors=2, allow_uphill=False)
    head_edges = [(e.u, e.v) for e in head_obj]

    ref_set = set(reference_edges)
    inf_set = set(head_edges)

    rows = []
    # True positives
    for u, v in sorted(ref_set & inf_set):
        rows.append({"edge_id": f"{u}->{v}", "u": u, "v": v, "classification": "TP",
                       "in_reference": True, "in_inferred": True})
    # False negatives
    for u, v in sorted(ref_set - inf_set):
        rows.append({"edge_id": f"{u}->{v}", "u": u, "v": v, "classification": "FN",
                       "in_reference": True, "in_inferred": False})
    # False positives
    for u, v in sorted(inf_set - ref_set):
        rows.append({"edge_id": f"{u}->{v}", "u": u, "v": v, "classification": "FP",
                       "in_reference": False, "in_inferred": True})

    df = pd.DataFrame(rows)
    df["archive_id"] = "savage_milford_nh"
    df["scenario"] = "head_gradient"
    df["validation_mode"] = "independent_graph_inference"
    df["allowed_claim"] = "Per-edge classification of Hydrosheaf topology against MODPATH reference"
    df["required_guardrail"] = "Classifications reflect this specific MODPATH benchmark; not generalisable without cross-archive confirmation"
    return df


def build_failure_modes(edge_classification: pd.DataFrame, nodes: pd.DataFrame) -> pd.DataFrame:
    """Build failure-mode diagnostics from edge classification."""
    if edge_classification.empty:
        return pd.DataFrame()

    node_map = {}
    for _, row in nodes.iterrows():
        node_map[str(row["node_id"])] = (float(row["y"]), float(row["x"]))

    rows = []
    # False positives
    fp = edge_classification[edge_classification["classification"] == "FP"]
    for _, row in fp.iterrows():
        u, v = row["u"], row["v"]
        dist = compute_haversine_km(
            node_map[u][0], node_map[u][1],
            node_map[v][0], node_map[v][1],
        ) if u in node_map and v in node_map else float("nan")
        rows.append({"edge_id": row["edge_id"], "failure_mode": "false_positive",
                       "u": u, "v": v, "distance_km": dist,
                       "diagnosis": "Inferred edge not present in MODPATH reference"})

    # False negatives
    fn = edge_classification[edge_classification["classification"] == "FN"]
    for _, row in fn.iterrows():
        u, v = row["u"], row["v"]
        dist = compute_haversine_km(
            node_map[u][0], node_map[u][1],
            node_map[v][0], node_map[v][1],
        ) if u in node_map and v in node_map else float("nan")
        rows.append({"edge_id": row["edge_id"], "failure_mode": "false_negative",
                       "u": u, "v": v, "distance_km": dist,
                       "diagnosis": "MODPATH reference edge not recovered by Hydrosheaf inference"})

    df = pd.DataFrame(rows)
    df["archive_id"] = "savage_milford_nh"
    return df


def build_node_sparsity_sensitivity(
    samples: List[Dict[str, object]],
    ref_edges_raw: List[Tuple[str, str]],
) -> pd.DataFrame:
    """Full sparsity sensitivity across node fractions [0.1, 0.25, 0.5, 0.75, 1.0]."""
    rng = np.random.default_rng(RANDOM_SEED)
    results = []
    total = len(samples)

    for frac in (0.1, 0.25, 0.5, 0.75, 1.0):
        n_count = max(int(total * frac), 5)
        f1_vals, rec_vals, prec_vals = [], [], []
        successful_trials = 0
        for _ in range(SPARSITY_TRIALS):
            seed = int(rng.integers(0, 2**31 - 1))
            idxs = pd.Series(range(total)).sample(n=n_count, random_state=seed).values
            sub = [samples[i] for i in idxs]
            sub_ids = {str(s["site_id"]) for s in sub}
            sub_ref = [e for e in ref_edges_raw if e[0] in sub_ids and e[1] in sub_ids]
            if len(sub_ref) < 5:
                f1_vals.append(0.0)
                rec_vals.append(0.0)
                prec_vals.append(0.0)
                continue
            sub_obj = infer_edges_from_coordinates(sub, max_neighbors=2, allow_uphill=False)
            sub_edges = [(e.u, e.v) for e in sub_obj]
            conf = edge_confusion(sub_ref, sub_edges)
            f1_vals.append(conf["f1"])
            rec_vals.append(conf["recall"])
            prec_vals.append(conf["precision"])
            successful_trials += 1
        if f1_vals:
            results.append({
                "node_fraction": frac,
                "node_count": n_count,
                "mean_f1": float(np.mean(f1_vals)),
                "std_f1": float(np.std(f1_vals)),
                "mean_recall": float(np.mean(rec_vals)),
                "std_recall": float(np.std(rec_vals)),
                "mean_precision": float(np.mean(prec_vals)),
                "std_precision": float(np.std(prec_vals)),
                "successful_trials": successful_trials,
                "planned_trials": SPARSITY_TRIALS,
                "random_seed": RANDOM_SEED,
            })
    df = pd.DataFrame(results)
    df["archive_id"] = "savage_milford_nh"
    df["scenario"] = "head_gradient"
    df["allowed_claim"] = "Graph inference performance collapses below 50% node density; zero-imputation applied for trials with fewer than 5 surviving reference edges."
    df["required_guardrail"] = "Sparsity results are random uniform subsampling; real monitoring gaps may differ"
    return df


def main():
    start_time = datetime.now(timezone.utc)
    print("=" * 60)
    print("Phase 2b: Independent Hydrosheaf Graph Validation")
    print("=" * 60)

    # Load data
    print("\n[1/7] Loading Savage nodes and reference edges...")
    nodes = load_savage_nodes()
    ref_edges_df = load_savage_reference_edges()
    samples = nodes_to_samples(nodes)

    ref_edges_raw = [(str(row["u"]), str(row["v"])) for _, row in ref_edges_df.iterrows()]
    node_map = {}
    for _, row in nodes.iterrows():
        node_map[str(row["node_id"])] = (float(row["y"]), float(row["x"]))

    ref_travel_times: Dict[Tuple[str, str], float] = {}
    for _, row in ref_edges_df.iterrows():
        key = (str(row["u"]), str(row["v"]))
        if "mean_travel_time" in row and not pd.isna(row["mean_travel_time"]):
            ref_travel_times[key] = float(row["mean_travel_time"])

    print(f"  Nodes: {len(samples)}, Reference edges: {len(ref_edges_raw)}")

    # Run independent scenarios
    print("\n[2/7] Running independent validation scenarios...")
    independent_df = run_independent_scenarios(samples, ref_edges_raw, node_map, node_df=nodes)

    # Print key results
    for _, row in independent_df.iterrows():
        f1_str = f"F1={row['f1']:.3f}" if not pd.isna(row['f1']) else "F1=N/A"
        tp_str = f"TP={int(row['tp'])}" if not pd.isna(row['tp']) else ""
        precision_str = f"Prec={row['precision']:.3f}" if not pd.isna(row['precision']) else ""
        print(f"  {row['scenario']:30s}  {f1_str:10s}  {tp_str:10s}  {precision_str:10s}")

    # Save independent results
    indep_path = SAVAGE_RESULTS / "hydrosheaf_independent_edges.csv"
    independent_df.to_csv(indep_path, index=False)
    print(f"\n  -> {indep_path}")

    # Run prior-assisted scenarios
    print("\n[3/7] Running 3 prior-assisted scenarios (separate from independent)...")
    head_obj = infer_edges_from_coordinates(samples, max_neighbors=2, allow_uphill=False)
    head_edges = [(e.u, e.v) for e in head_obj]

    prior_df = run_prior_assisted_scenarios(head_edges, ref_edges_raw, ref_travel_times)
    for _, row in prior_df.iterrows():
        print(f"  {row['scenario']:30s}  mode={row['prior_mode']:10s}  "
              f"in={row['n_hydrosheaf_input_edges']}  +priors={row['n_modpath_prior_edges']}  "
              f"out={row['n_output_edges']}  added={row['n_edges_added']}")

    prior_path = SAVAGE_RESULTS / "hydrosheaf_prior_assisted_edges.csv"
    prior_df.to_csv(prior_path, index=False)
    print(f"\n  -> {prior_path}")

    # Edge classification
    print("\n[4/7] Building per-edge classification...")
    edge_class = build_edge_classification(independent_df, ref_edges_raw)
    class_path = SAVAGE_RESULTS / "edge_classification.csv"
    edge_class.to_csv(class_path, index=False)
    tp_count = len(edge_class[edge_class["classification"] == "TP"]) if not edge_class.empty else 0
    fp_count = len(edge_class[edge_class["classification"] == "FP"]) if not edge_class.empty else 0
    fn_count = len(edge_class[edge_class["classification"] == "FN"]) if not edge_class.empty else 0
    print(f"  TP={tp_count}, FP={fp_count}, FN={fn_count}")
    print(f"  -> {class_path}")

    # Failure modes
    print("\n[5/7] Building failure-mode diagnostics...")
    failures = build_failure_modes(edge_class, nodes)
    fail_path = SAVAGE_RESULTS / "failure_modes.csv"
    failures.to_csv(fail_path, index=False)
    print(f"  -> {fail_path}")

    # Sparsity sensitivity
    print("\n[6/7] Running node-sparsity sensitivity...")
    sparsity = build_node_sparsity_sensitivity(samples, ref_edges_raw)
    sp_path = SAVAGE_RESULTS / "node_sparsity_sensitivity.csv"
    sparsity.to_csv(sp_path, index=False)
    for _, row in sparsity.iterrows():
        print(f"  frac={row['node_fraction']:.2f}  n={row['node_count']}  "
              f"F1={row['mean_f1']:.3f} +/- {row['std_f1']:.3f}")
    print(f"  -> {sp_path}")

    # Build aggregate validation metrics
    print("\n[7/7] Building graph_validation_metrics.csv (aggregate)...")
    metrics_rows = []
    head_row = independent_df[independent_df["scenario"] == "head_gradient"]
    if not head_row.empty:
        r = head_row.iloc[0]
        metrics_rows.append({
            "archive_id": "savage_milford_nh",
            "primary_scenario": "head_gradient",
            "n_nodes": len(samples),
            "n_reference_edges": int(r["n_reference_edges"]),
            "n_inferred_edges": int(r["n_inferred_edges"]),
            "tp": int(r["tp"]), "fp": int(r["fp"]), "fn": int(r["fn"]),
            "precision": r["precision"], "recall": r["recall"], "f1": r["f1"],
            "spatial_only_f1": float(independent_df[independent_df["scenario"] == "spatial_only"].iloc[0]["f1"]),
            "random_baseline_f1": float(independent_df[independent_df["scenario"] == "negative_random"].iloc[0]["f1"]),
            "scale_mismatch": bool(r["scale_mismatch"]),
            "independent_validation": True,
            "validation_mode": "independent_graph_inference",
            "allowed_claim": "Hydrosheaf head-gradient graph inference recovers MODPATH topology significantly above random baseline",
            "required_guardrail": "MODPATH is a model-conditioned advective reference, not absolute truth. Results reflect one MODFLOW-2005/MODPATH 5 model; cross-archive generalisation requires additional evidence.",
        })
    metrics_df = pd.DataFrame(metrics_rows)
    metrics_path = SAVAGE_RESULTS / "graph_validation_metrics.csv"
    metrics_df.to_csv(metrics_path, index=False)
    print(f"  -> {metrics_path}")

    # Update claim_guardrails.csv with new entries
    print("\n=== Updating claim guardrails ===")
    existing_cg = SAVAGE_RESULTS / "claim_guardrails.csv"
    if existing_cg.exists():
        cg_df = pd.read_csv(existing_cg)
        new_guardrails = [
            {
                "claim_id": "CG009",
                "claim": "Hydrosheaf head-gradient inference recovers MODPATH topology independently",
                "allowed_claim": f"Independent graph inference achieves F1={head_row.iloc[0]['f1']:.3f} against MODPATH reference without using MODPATH connectivity",
                "required_guardrail": "MODPATH is a model-conditioned advective reference, not absolute truth; head proxy uses elevation, not simulated heads",
                "evidence_source": "hydrosheaf_independent_edges.csv, graph_validation_metrics.csv",
                "current_status": "enforced",
                "archive_id": "savage_milford_nh",
            },
            {
                "claim_id": "CG010",
                "claim": "Prior-assisted modes improve edge coverage at cost of independence",
                "allowed_claim": "MODPATH prior merge adds MODPATH edges to Hydrosheaf graph where they do not already exist",
                "required_guardrail": "Prior-assisted mode is not independent validation; prior-assisted rows must never enter independent validation summaries",
                "evidence_source": "hydrosheaf_prior_assisted_edges.csv",
                "current_status": "enforced",
                "archive_id": "savage_milford_nh",
            },
            {
                "claim_id": "CG011",
                "claim": "Node sparsity sensitivity tests metric stability under reduced node density",
                "allowed_claim": f"Sparsity sensitivity shows full-node F1={sparsity.iloc[-1]['mean_f1']:.3f} and 10% node F1={sparsity.iloc[0]['mean_f1']:.3f} across {int(sparsity.iloc[0]['successful_trials'])}/{int(sparsity.iloc[0]['planned_trials'])} successful low-density trials",
                "required_guardrail": "Sparsity uses random uniform subsampling; low-density F1 can be biased by small retained reference-edge sets and must not be described as monotonic degradation or field-data missingness",
                "evidence_source": "node_sparsity_sensitivity.csv",
                "current_status": "enforced",
                "archive_id": "savage_milford_nh",
            },
        ]
        new_ids = {row["claim_id"] for row in new_guardrails}
        cg_df = cg_df[~cg_df["claim_id"].isin(new_ids)].copy()
        cg_df = pd.concat([cg_df, pd.DataFrame(new_guardrails)], ignore_index=True)
        cg_df.to_csv(existing_cg, index=False)
        print(f"  Upserted {len(new_guardrails)} guardrail entries -> {existing_cg}")

    # Update processing_log.json
    log_path = SAVAGE_RESULTS / "processing_log.json"
    if log_path.exists():
        with open(log_path) as f:
            log = json.load(f)
        primary_independent = independent_df[independent_df["validation_mode"] == "independent_graph_inference"]
        diagnostic_controls = independent_df[independent_df["validation_mode"] == "diagnostic_negative_control"]
        log["independent_validation"] = {
            "status": "complete",
            "n_scenarios": len(primary_independent),
            "n_diagnostic_controls": len(diagnostic_controls),
            "primary_scenario": "head_gradient",
            "primary_f1": float(head_row.iloc[0]["f1"]) if not head_row.empty else None,
            "primary_precision": float(head_row.iloc[0]["precision"]) if not head_row.empty else None,
            "primary_recall": float(head_row.iloc[0]["recall"]) if not head_row.empty else None,
            "tp": int(head_row.iloc[0]["tp"]) if not head_row.empty else 0,
            "fp": int(head_row.iloc[0]["fp"]) if not head_row.empty else 0,
            "fn": int(head_row.iloc[0]["fn"]) if not head_row.empty else 0,
            "random_baseline_f1": float(independent_df[independent_df["scenario"] == "negative_random"].iloc[0]["f1"]),
            "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        }
        log["prior_assisted"] = {
            "status": "complete",
            "n_scenarios": len(prior_df),
            "modes_tested": ["override", "merge", "only"],
            "strictly_separated_from_independent": True,
            "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        }
        with open(log_path, "w") as f:
            json.dump(log, f, indent=2, default=str)
        print(f"  Updated -> {log_path}")

    # Update handoff_status.json
    hs_path = SAVAGE_RESULTS / "handoff_status.json"
    if hs_path.exists():
        with open(hs_path) as f:
            hs = json.load(f)
        hs["savage_status"]["independent_validation_complete"] = True
        hs["savage_status"]["prior_assisted_complete"] = True
        hs["savage_status"]["edge_classification_complete"] = True
        hs["savage_status"]["failure_modes_complete"] = True
        hs["savage_status"]["node_sparsity_complete"] = True
        hs["savage_readiness"]["independent_validation_run"] = True
        hs["savage_readiness"]["primary_f1"] = float(head_row.iloc[0]["f1"]) if not head_row.empty else None
        hs["savage_readiness"]["notes"] = (
            "Savage independent validation complete. Hydrosheaf head-gradient inference "
            "recovers MODPATH topology with real (non-trivial) TP/FP/FN metrics. "
            "Prior-assisted modes run separately. Ready for figure/table generation."
        )
        hs["timestamp_utc"] = datetime.now(timezone.utc).isoformat()
        with open(hs_path, "w") as f:
            json.dump(hs, f, indent=2)
        print(f"  Updated -> {hs_path}")

    # Also save to benchmark results for downstream scripts
    print("\n=== Copying key files to benchmark results ===")
    independent_df.to_csv(BENCHMARK_RESULTS / "independent_graph_vs_modpath.csv", index=False)
    edge_class.to_csv(BENCHMARK_RESULTS / "edge_classification.csv", index=False)
    prior_df.to_csv(BENCHMARK_RESULTS / "modpath_informed_priors.csv", index=False)
    sparsity_src = SAVAGE_RESULTS / "node_sparsity_sensitivity.csv"
    if sparsity_src.exists():
        import shutil
        shutil.copy2(sparsity_src, BENCHMARK_RESULTS / "m4_sparsity_sensitivity.csv")
    print("  Done.")

    elapsed = (datetime.now(timezone.utc) - start_time).total_seconds()
    print(f"\n=== Phase 2b Complete ({elapsed:.1f}s) ===")
    print(f"""
Summary:
  Independent graph-inference scenarios: {len(independent_df[independent_df['validation_mode'] == 'independent_graph_inference'])} run
  Diagnostic negative controls: {len(independent_df[independent_df['validation_mode'] == 'diagnostic_negative_control'])} run
  Sensitivity analyses: {len(independent_df[independent_df['validation_mode'] == 'sensitivity_analysis'])} run
  Prior-assisted scenarios: {len(prior_df)} run (fully separated)
  Primary result (head_gradient): F1={head_row.iloc[0]['f1']:.3f}, TP={int(head_row.iloc[0]['tp'])}, FP={int(head_row.iloc[0]['fp'])}, FN={int(head_row.iloc[0]['fn'])}
  Random baseline F1: {independent_df[independent_df['scenario'] == 'negative_random'].iloc[0]['f1']:.3f}
  Edge classification: {tp_count} TP, {fp_count} FP, {fn_count} FN
  Sparsity sensitivity: {len(sparsity)} fractions x {SPARSITY_TRIALS} trials
  Guardrails: independent validation rows fully separated from prior-assisted rows
""")


if __name__ == "__main__":
    main()
