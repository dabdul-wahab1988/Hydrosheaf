"""Hydrosheaf package entry points."""

from .config import Config, DEFAULT_ION_ORDER
from .inference.edge_fit import EdgeResult, fit_edge
from .graph.build import infer_edges_probabilistic
from .graph3d.build_3d import build_network_3d, infer_edges_3d_probabilistic
from .inference.network_fit import (
    fit_edges,
    fit_network,
    edge_process_maps,
    infer_edges,
    predict_node_ec_tds,
    summarize_network,
)
from .models.ec_tds import calibrate_ec_tds, predict_ec_tds
from .models.reactions import build_reaction_dictionary
from .phreeqc.constraints import build_edge_bounds
from .phreeqc.runner import run_phreeqc
from .isotopes import compute_d_excess, evaporation_index, fit_lmwl
from .coda_sbp import ilr_from_sbp, robust_zscore
from .nitrate_source_v2 import infer_node_posteriors, NitrateSourceResult
from .api import (
    attach_temporal_results,
    auto_disable_missing_modules,
    build_vadose_priors,
    validate_required_inputs,
    fit_network_pipeline,
    fit_network_with_priors,
    fit_temporal_edges,
)
from .physics.priors import PhysicsPrior, apply_physics_priors
from .temporal import TemporalEdgeResult, TemporalNode, TimeSeriesSample
from .temporal.time_series import load_time_series_csv
from .vadose.contracts import (
    VadoseForcingSample,
    VadoseLinksRow,
    VadoseProfile,
    VadoseRunConfig,
)
from .tuning import TuningReport, tune_reaction_hyperparameters

__all__ = [
    "Config",
    "DEFAULT_ION_ORDER",
    "EdgeResult",
    "fit_edge",
    "fit_edges",
    "fit_network",
    "edge_process_maps",
    "infer_edges",
    "infer_edges_probabilistic",
    "infer_edges_3d_probabilistic",
    "build_network_3d",
    "predict_node_ec_tds",
    "summarize_network",
    "calibrate_ec_tds",
    "predict_ec_tds",
    "build_reaction_dictionary",
    "run_phreeqc",
    "build_edge_bounds",
    "compute_d_excess",
    "evaporation_index",
    "fit_lmwl",
    "ilr_from_sbp",
    "robust_zscore",
    "infer_node_posteriors",
    "NitrateSourceResult",
    "attach_temporal_results",
    "auto_disable_missing_modules",
    "build_vadose_priors",
    "validate_required_inputs",
    "fit_network_pipeline",
    "fit_network_with_priors",
    "fit_temporal_edges",
    "PhysicsPrior",
    "apply_physics_priors",
    "TemporalEdgeResult",
    "TemporalNode",
    "TimeSeriesSample",
    "load_time_series_csv",
    "VadoseForcingSample",
    "VadoseLinksRow",
    "VadoseProfile",
    "VadoseRunConfig",
    "TuningReport",
    "tune_reaction_hyperparameters",
]

__version__ = "0.1.0"
