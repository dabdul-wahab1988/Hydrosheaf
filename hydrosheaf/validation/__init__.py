"""Validation helpers for manuscript-grade Hydrosheaf evidence.

The validation package keeps benchmark metrics and claim guardrails close to
the code that produces Hydrosheaf results. It is intentionally lightweight so
benchmark scripts and manuscript tables can share the same definitions.
"""

from .claims import (
    ClaimRecord,
    EvidenceLevel,
    assess_claim_records,
    evidence_level_allows,
)
from .evidence import (
    EdgeEvidenceAnnotation,
    EdgeEvidenceClass,
    classify_edge_evidence,
)
from .metrics import (
    classification_metrics,
    interval_coverage,
    regression_metrics,
    topology_metrics,
)
from .topology import (
    apply_modpath_informed_graph_priors,
    build_modpath_informed_graph_priors,
    edge_confusion,
    normalize_directed_edges,
    scale_mismatch_diagnostics,
    validate_independent_graph_against_modpath,
)
from .reaction import (
    fit_sparse_reaction_once,
    l1_penalty_sensitivity,
    missing_ion_sensitivity,
    thermodynamic_bound_violations,
    validate_sparse_inverse_reaction_model,
)
from .modpath_archive import (
    scan_modpath_archive,
    load_modpath_endpoints,
    load_modpath_pathlines,
    standardise_endpoint_columns,
    standardise_pathline_columns,
    build_node_mapping,
    build_reference_edges,
    summarise_archive_evidence,
)

__all__ = [
    "ClaimRecord",
    "EvidenceLevel",
    "EdgeEvidenceAnnotation",
    "EdgeEvidenceClass",
    "assess_claim_records",
    "classify_edge_evidence",
    "classification_metrics",
    "EdgeEvidenceAnnotation",
    "EdgeEvidenceClass",
    "classify_edge_evidence",
    "evidence_level_allows",
    "interval_coverage",
    "regression_metrics",
    "topology_metrics",
    "apply_modpath_informed_graph_priors",
    "build_modpath_informed_graph_priors",
    "edge_confusion",
    "normalize_directed_edges",
    "scale_mismatch_diagnostics",
    "validate_independent_graph_against_modpath",
    "fit_sparse_reaction_once",
    "l1_penalty_sensitivity",
    "missing_ion_sensitivity",
    "thermodynamic_bound_violations",
    "validate_sparse_inverse_reaction_model",
    "scan_modpath_archive",
    "load_modpath_endpoints",
    "load_modpath_pathlines",
    "standardise_endpoint_columns",
    "standardise_pathline_columns",
    "build_node_mapping",
    "build_reference_edges",
    "summarise_archive_evidence",
]
