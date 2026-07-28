"""Hydrosheaf package entry points."""

from importlib import import_module
from typing import Any, Dict, List, Tuple

_EXPORTS: Dict[str, Tuple[str, str]] = {
    "Config": (".config", "Config"),
    "DEFAULT_ION_ORDER": (".config", "DEFAULT_ION_ORDER"),
    "EdgeResult": (".inference.edge_fit", "EdgeResult"),
    "fit_edge": (".inference.edge_fit", "fit_edge"),
    "fit_edges": (".inference.network_fit", "fit_edges"),
    "fit_network": (".inference.network_fit", "fit_network"),
    "edge_process_maps": (".inference.network_fit", "edge_process_maps"),
    "infer_edges": (".inference.network_fit", "infer_edges"),
    "predict_node_ec_tds": (".inference.network_fit", "predict_node_ec_tds"),
    "summarize_network": (".inference.network_fit", "summarize_network"),
    "infer_edges_probabilistic": (".graph.build", "infer_edges_probabilistic"),
    "infer_edges_3d_probabilistic": (
        ".graph3d.build_3d",
        "infer_edges_3d_probabilistic",
    ),
    "build_network_3d": (".graph3d.build_3d", "build_network_3d"),
    "calibrate_ec_tds": (".models.ec_tds", "calibrate_ec_tds"),
    "predict_ec_tds": (".models.ec_tds", "predict_ec_tds"),
    "build_reaction_dictionary": (".models.reactions", "build_reaction_dictionary"),
    "EvidenceLiftedResolution": (
        ".models.evidence_lifted",
        "EvidenceLiftedResolution",
    ),
    "evidence_lifted_resolution": (
        ".models.evidence_lifted",
        "evidence_lifted_resolution",
    ),
    "stoichiometric_equivalence_classes": (
        ".models.evidence_lifted",
        "stoichiometric_equivalence_classes",
    ),
    "run_phreeqc": (".phreeqc.runner", "run_phreeqc"),
    "build_edge_bounds": (".phreeqc.constraints", "build_edge_bounds"),
    "compute_d_excess": (".isotopes", "compute_d_excess"),
    "evaporation_index": (".isotopes", "evaporation_index"),
    "fit_lmwl": (".isotopes", "fit_lmwl"),
    "ilr_from_sbp": (".coda_sbp", "ilr_from_sbp"),
    "robust_zscore": (".coda_sbp", "robust_zscore"),
    "infer_node_posteriors": (".nitrate_source_v2", "infer_node_posteriors"),
    "NitrateSourceResult": (".nitrate_source_v2", "NitrateSourceResult"),
    "attach_temporal_results": (".api", "attach_temporal_results"),
    "auto_disable_missing_modules": (".api", "auto_disable_missing_modules"),
    "build_vadose_priors": (".api", "build_vadose_priors"),
    "validate_required_inputs": (".api", "validate_required_inputs"),
    "fit_network_pipeline": (".api", "fit_network_pipeline"),
    "fit_network_with_priors": (".api", "fit_network_with_priors"),
    "fit_temporal_edges": (".api", "fit_temporal_edges"),
    "PhysicsPrior": (".physics.priors", "PhysicsPrior"),
    "apply_physics_priors": (".physics.priors", "apply_physics_priors"),
    "TemporalEdgeResult": (".temporal", "TemporalEdgeResult"),
    "TemporalNode": (".temporal", "TemporalNode"),
    "TimeSeriesSample": (".temporal", "TimeSeriesSample"),
    "load_time_series_csv": (".temporal.time_series", "load_time_series_csv"),
    "VadoseForcingSample": (".vadose.contracts", "VadoseForcingSample"),
    "VadoseLinksRow": (".vadose.contracts", "VadoseLinksRow"),
    "VadoseProfile": (".vadose.contracts", "VadoseProfile"),
    "VadoseRunConfig": (".vadose.contracts", "VadoseRunConfig"),
    "TuningReport": (".tuning", "TuningReport"),
    "tune_reaction_hyperparameters": (".tuning", "tune_reaction_hyperparameters"),
    "AcquisitionConfig": (
        ".calibration.bayesian_active_learning",
        "AcquisitionConfig",
    ),
    "MeasurementOption": (
        ".calibration.bayesian_active_learning",
        "MeasurementOption",
    ),
    "PredictiveScenario": (
        ".calibration.bayesian_active_learning",
        "PredictiveScenario",
    ),
    "expected_information_gain": (
        ".calibration.bayesian_active_learning",
        "expected_information_gain",
    ),
    "expected_brier_risk_reduction": (
        ".calibration.bayesian_active_learning",
        "expected_brier_risk_reduction",
    ),
    "rank_measurement_options": (
        ".calibration.bayesian_active_learning",
        "rank_measurement_options",
    ),
    "select_measurement_batch": (
        ".calibration.bayesian_active_learning",
        "select_measurement_batch",
    ),
    "update_hypothesis_posterior": (
        ".calibration.bayesian_active_learning",
        "update_hypothesis_posterior",
    ),
    # NEW WORKFLOW
    "analyze_dataset": (".workflows.auto", "analyze_dataset"),
    "ClaimRecord": (".validation", "ClaimRecord"),
    "EvidenceLevel": (".validation", "EvidenceLevel"),
    "assess_claim_records": (".validation", "assess_claim_records"),
    "classification_metrics": (".validation", "classification_metrics"),
    "interval_coverage": (".validation", "interval_coverage"),
    "regression_metrics": (".validation", "regression_metrics"),
    "topology_metrics": (".validation", "topology_metrics"),
    "apply_modpath_informed_graph_priors": (
        ".validation",
        "apply_modpath_informed_graph_priors",
    ),
    "build_modpath_informed_graph_priors": (
        ".validation",
        "build_modpath_informed_graph_priors",
    ),
    "edge_confusion": (".validation", "edge_confusion"),
    "normalize_directed_edges": (".validation", "normalize_directed_edges"),
    "scale_mismatch_diagnostics": (".validation", "scale_mismatch_diagnostics"),
    "validate_independent_graph_against_modpath": (
        ".validation",
        "validate_independent_graph_against_modpath",
    ),
    "fit_sparse_reaction_once": (".validation", "fit_sparse_reaction_once"),
    "l1_penalty_sensitivity": (".validation", "l1_penalty_sensitivity"),
    "missing_ion_sensitivity": (".validation", "missing_ion_sensitivity"),
    "thermodynamic_bound_violations": (
        ".validation",
        "thermodynamic_bound_violations",
    ),
    "validate_sparse_inverse_reaction_model": (
        ".validation",
        "validate_sparse_inverse_reaction_model",
    ),
}

__all__ = list(_EXPORTS)
__version__ = "0.5.1"


def __getattr__(name: str) -> Any:
    if name not in _EXPORTS:
        raise AttributeError(f"module 'hydrosheaf' has no attribute '{name}'")
    module_name, attr_name = _EXPORTS[name]
    module = import_module(module_name, __name__)
    value = getattr(module, attr_name)
    globals()[name] = value
    return value


def __dir__() -> List[str]:
    return sorted(set(globals().keys()) | set(__all__))
