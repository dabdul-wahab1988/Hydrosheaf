from __future__ import annotations

from hydrosheaf.config import Config
from hydrosheaf.inference.edge_fit import fit_edge


def _config(**updates: object) -> Config:
    values = {
        "phreeqc_enabled": False,
        "gibbs_enabled": False,
        "isotope_enabled": False,
        "sheaf_isotope_enabled": False,
        "sheaf_cl_enabled": False,
        "sheaf_age_enabled": False,
        "transport_models_enabled": ["evap"],
        "active_minerals": [],
        "exchange_enabled": False,
        "reaction_processes_enabled": [],
    }
    values.update(updates)
    return Config(**values)


def test_ratio_penalty_is_optional_and_does_not_change_matrix_schema() -> None:
    upstream = {"Ca": 1.0, "Mg": 1.0, "Na": 1.0, "HCO3": 2.0, "Cl": 1.0, "SO4": 1.0, "NO3": 1.0}
    downstream = {"Ca": 2.0, "Mg": 1.0, "Na": 1.0, "HCO3": 2.0, "Cl": 1.0, "SO4": 1.0, "NO3": 1.0}
    baseline = fit_edge(
        [upstream[ion] for ion in ("Ca", "Mg", "Na", "HCO3", "Cl", "SO4", "NO3")],
        [downstream[ion] for ion in ("Ca", "Mg", "Na", "HCO3", "Cl", "SO4", "NO3")],
        _config(),
        obs_u=upstream,
        obs_v=downstream,
    )
    ratio_enabled = fit_edge(
        [upstream[ion] for ion in ("Ca", "Mg", "Na", "HCO3", "Cl", "SO4", "NO3")],
        [downstream[ion] for ion in ("Ca", "Mg", "Na", "HCO3", "Cl", "SO4", "NO3")],
        _config(ratio_features_enabled=True, ratio_penalty_weight=1.0),
        obs_u=upstream,
        obs_v=downstream,
    )
    assert ratio_enabled.ratio_fit_pairs >= 2
    assert ratio_enabled.ratio_fit_penalty >= 0.0
    assert len(ratio_enabled.z_labels) == len(baseline.z_labels)
    assert len(ratio_enabled.z_extents) == len(ratio_enabled.z_labels)
