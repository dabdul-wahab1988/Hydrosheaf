from __future__ import annotations

from hydrosheaf.config import Config
from hydrosheaf.inference.network_fit import fit_network


def _record(node_id: str, site_id: str, *, cl: float) -> dict[str, object]:
    return {
        "node_id": node_id,
        "sample_id": node_id,
        "site_id": site_id,
        "Ca": 1.0,
        "Mg": 1.0,
        "Na": 1.0,
        "K": 1.0,
        "HCO3": 1.0,
        "Cl": cl,
        "SO4": 1.0,
        "NO3": 1.0,
        "F": 1.0,
        "Fe": None,
        "PO4": None,
        "EC": 10.0,
        "TDS": 10.0,
        "pH": 7.0,
    }


def test_repeated_site_records_are_addressable_by_node_id() -> None:
    samples = [
        _record("NG_W1_dry", "NG_W1", cl=1.0),
        _record("NG_W1_wet", "NG_W1", cl=1.2),
    ]
    config = Config(
        phreeqc_enabled=False,
        gibbs_enabled=False,
        isotope_enabled=False,
        sheaf_isotope_enabled=False,
        sheaf_cl_enabled=False,
        sheaf_age_enabled=False,
        transport_models_enabled=["evap"],
        active_minerals=[],
        exchange_enabled=False,
        reaction_processes_enabled=[],
        sparse_panel_enabled=True,
    )
    results = fit_network(samples, [("NG_W1_dry", "NG_W1_wet")], config)
    assert len(results) == 1
    assert results[0].u == "NG_W1_dry"
    assert results[0].v == "NG_W1_wet"
