"""Scientific-contract tests for the M7.1 blind benchmark."""

from __future__ import annotations

import json
import sys
from copy import deepcopy
from pathlib import Path


SCRIPTS = (
    Path(__file__).resolve().parents[1]
    / "M7"
    / "m7_integration_benchmark"
    / "scripts"
)
sys.path.insert(0, str(SCRIPTS))

from m7_1_inference import run_blind_inference  # noqa: E402
from m7_1_truth_model import generate_blind_aquifer, poison_truth  # noqa: E402
from hydrosheaf.config import Config  # noqa: E402
from hydrosheaf.graph.types import Edge  # noqa: E402
from hydrosheaf.sheaf.topology_refine import refine_edges_with_sheaf  # noqa: E402


def test_truth_generator_is_reproducible_and_seed_sensitive():
    first = generate_blind_aquifer(10123)
    replay = generate_blind_aquifer(10123)
    other = generate_blind_aquifer(10124)
    assert first == replay
    assert first.observations != other.observations
    assert first.true_edges


def test_truth_poisoning_cannot_change_inference():
    case = generate_blind_aquifer(10100)
    poisoned = poison_truth(case)
    original = run_blind_inference(case.observations, threshold=0.2, seed=case.seed)
    attacked = run_blind_inference(
        poisoned.observations, threshold=0.2, seed=poisoned.seed
    )
    assert json.dumps(original.serializable(), sort_keys=True) == json.dumps(
        attacked.serializable(), sort_keys=True
    )


def test_fast_tier_calls_real_modules_and_preserves_candidate_recall():
    case = generate_blind_aquifer(20200)
    result = run_blind_inference(case.observations, threshold=0.2, seed=case.seed)
    assert result.module_calls["infer_edges"] == 1
    assert result.module_calls["refine_edges_with_sheaf"] == 2
    assert result.module_calls["fit_network"] == 1
    assert result.module_calls["run_phreeqc"] == 0
    truth_ids = {f"{u}->{v}" for u, v in case.true_edges}
    assert truth_ids <= set(result.candidate_edges)
    assert all(
        0.0 <= float(row["joint_probability"]) <= 1.0
        for row in result.edge_scores
    )


def test_held_out_ions_cannot_affect_inference_scores_or_extents():
    case = generate_blind_aquifer(20201)
    attacked_rows = deepcopy(list(case.observations))
    for index, row in enumerate(attacked_rows):
        row["Mg"] = 1000.0 + index
        row["SO4"] = 2000.0 + index
        row["Fe"] = 3000.0 + index
    original = run_blind_inference(case.observations, threshold=0.2, seed=case.seed)
    attacked = run_blind_inference(attacked_rows, threshold=0.2, seed=case.seed)
    assert original.candidate_edges == attacked.candidate_edges
    keys = (
        "hydraulic_probability",
        "age_only_probability",
        "sheaf_multievidence_probability",
        "chemistry_probability",
        "joint_probability",
        "reaction_extents",
    )
    original_by_id = {
        str(row["edge_id"]): row for row in original.edge_scores
    }
    attacked_by_id = {
        str(row["edge_id"]): row for row in attacked.edge_scores
    }
    for edge_id in original_by_id:
        assert {
            key: original_by_id[edge_id][key] for key in keys
        } == {
            key: attacked_by_id[edge_id][key] for key in keys
        }


def test_sheaf_weight_does_not_overwrite_hydraulic_prior():
    samples = [
        {"site_id": "A", "Cl": 0.2},
        {"site_id": "B", "Cl": 0.21},
    ]
    edge = Edge(
        edge_id="A->B",
        u="A",
        v="B",
        attrs={"edge_confidence": 0.83, "p_uv": 0.83},
    )
    selected = refine_edges_with_sheaf(
        samples,
        [edge],
        Config(
            phreeqc_enabled=False,
            sheaf_age_enabled=False,
            sheaf_isotope_enabled=False,
            sheaf_max_iter=0,
            edge_max_neighbors=1,
        ),
    )
    assert selected
    assert selected[0].attrs["edge_confidence"] == 0.83
    assert selected[0].attrs["prior_edge_probability"] == 0.83
    assert "sheaf_weight" in selected[0].attrs
