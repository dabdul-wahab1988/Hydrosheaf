"""Hydrosheaf package entry points."""

from .config import Config, DEFAULT_ION_ORDER
from .inference.edge_fit import EdgeResult, fit_edge
from .graph.build import infer_edges_probabilistic
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
from .models.sheaf import edge_residual
from .phreeqc.constraints import build_edge_bounds
from .phreeqc.runner import run_phreeqc
from .isotopes import compute_d_excess, evaporation_index, fit_lmwl

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
    "predict_node_ec_tds",
    "summarize_network",
    "calibrate_ec_tds",
    "predict_ec_tds",
    "build_reaction_dictionary",
    "edge_residual",
    "run_phreeqc",
    "build_edge_bounds",
    "compute_d_excess",
    "evaporation_index",
    "fit_lmwl",
]

__version__ = "0.1.0"
