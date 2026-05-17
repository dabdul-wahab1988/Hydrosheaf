from __future__ import annotations

import importlib.util
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = REPO_ROOT / "M3" / "m3_age_benchmark" / "scripts" / "run_m3_real_usgs_graph_benchmark.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("run_m3_real_usgs_graph_benchmark", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_graph_regularize_enforces_monotonic_edge_order():
    module = _load_module()
    ages = {"A": 20.0, "B": 10.0}
    out = module.graph_regularize(ages, [("A", "B")], 0.85, iterations=8)
    assert out["B"] >= out["A"] - 1.0e-4


def test_build_graph_families_creates_candidate_and_control_edges():
    module = _load_module()
    df = module.pd.DataFrame(
        [
            {"site_id": "A", "StudyUnit": "U", "lat": 0.0, "lon": 0.0, "depth_m": 10.0, "ref_age": 10.0, "est_age_multi": 10.0},
            {"site_id": "B", "StudyUnit": "U", "lat": 0.0, "lon": 0.1, "depth_m": 20.0, "ref_age": 20.0, "est_age_multi": 20.0},
            {"site_id": "C", "StudyUnit": "U", "lat": 0.0, "lon": 0.2, "depth_m": 30.0, "ref_age": 30.0, "est_age_multi": 30.0},
            {"site_id": "D", "StudyUnit": "U", "lat": 0.0, "lon": 0.3, "depth_m": 40.0, "ref_age": 40.0, "est_age_multi": 40.0},
            {"site_id": "E", "StudyUnit": "U", "lat": 0.0, "lon": 0.4, "depth_m": 50.0, "ref_age": 50.0, "est_age_multi": 50.0},
        ]
    )
    families, edge_rows = module.build_graph_families(df, min_unit_size=5)
    assert "study_unit_coordinate" in families
    assert "wrong_direction_negative_control" in families
    assert "randomized_negative_control" in families
    assert edge_rows["family"].nunique() >= 3
