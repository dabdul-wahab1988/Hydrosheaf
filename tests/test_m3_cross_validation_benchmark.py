from __future__ import annotations

import importlib.util
from pathlib import Path


SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "M3"
    / "m3_age_benchmark"
    / "scripts"
    / "run_m3_cross_validation_benchmark.py"
)


def _load_module():
    spec = importlib.util.spec_from_file_location("m3_cv_benchmark", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_cv_reportability_guard_requires_identifiable_positive_tau():
    module = _load_module()
    assert module._is_reportable_fit(
        {"fit_identifiable": True, "fit_mean_age_years": 12.0}
    )
    assert not module._is_reportable_fit(
        {"fit_identifiable": False, "fit_mean_age_years": 12.0}
    )
    assert not module._is_reportable_fit(
        {"fit_identifiable": True, "fit_mean_age_years": 0.0}
    )


def test_cv_graph_edges_cannot_include_full_fit_neighbours():
    module = _load_module()
    edges = [("withheld-a", "withheld-b"), ("withheld-a", "full-fit-c")]
    assert module._filter_edges_to_nodes(edges, {"withheld-a", "withheld-b"}) == [
        ("withheld-a", "withheld-b")
    ]


def test_mask_tracer_removes_only_requested_reported_input():
    module = _load_module()
    row = {"LPM_TracersMod": "3H, SF6, 14C"}
    masked = module.mask_tracer_in_row(row, "SF6")
    assert masked["LPM_TracersMod"] == "3H|14C"
    assert row["LPM_TracersMod"] == "3H, SF6, 14C"


def test_withheld_target_uses_the_same_correction_scale_as_the_fit():
    module = _load_module()
    c14, c14_source = module._withheld_observation(
        {"c14_pmc": 60.0, "corrected_c14_pmc": 92.0}, "14C"
    )
    assert (c14, c14_source) == (92.0, "selected_corrected_c14_pmc")

    sf6, sf6_source = module._withheld_observation(
        {"sf6_pptv": 4.0, "dgm_sf6_pptv": 5.0}, "SF6"
    )
    assert (sf6, sf6_source) == (5.0, "dgm_sf6_pptv")
