import pytest

from hydrosheaf.validation.baselines import (
    ABSTAIN,
    REJECT,
    SELECT,
    BaselineRegistry,
    default_baseline_registry,
    edge_local_topology_baseline_spec,
    hydraulic_only_baseline_spec,
    permutation_control_baseline,
)


def test_default_registry_exposes_auditable_metadata_contract() -> None:
    registry = default_baseline_registry()

    assert registry.names() == (
        "edge_local_topology_support",
        "hydraulic_only_directional_gradient",
    )
    for record in registry.audit_table():
        assert record["fingerprint"]
        assert record["input_channels"]
        assert set(("tuning", "uncertainty", "abstention", "cost")).issubset(
            record
        )
        assert record["abstention"]["select_threshold"] > record["abstention"][
            "reject_threshold"
        ]
        assert set(("false_positive", "false_negative", "abstain", "measurement")).issubset(
            record["cost"]
        )


def test_registry_rejects_duplicate_baseline_names() -> None:
    spec = hydraulic_only_baseline_spec()
    registry = BaselineRegistry([spec])

    with pytest.raises(ValueError, match="already registered"):
        registry.register(spec)


def test_hydraulic_baseline_scores_every_candidate_truth_blind() -> None:
    spec = hydraulic_only_baseline_spec()
    candidates = [("A", "B"), ("B", "A"), {"source": "A", "target": "C"}]
    observations = {
        "hydraulic": {
            "node_heads": {"A": 10.0, "B": 9.0, "C": 10.0},
            "edges": {("A", "C"): {"distance": 10.0}},
        },
        # Contradictory non-declared channel must not affect a hydraulic-only
        # baseline.
        "topology": {
            "edges": {
                ("A", "B"): {"probability": 0.01},
                ("B", "A"): {"probability": 0.99},
            }
        },
    }

    predictions = spec.score(candidates, observations)
    by_edge = {prediction.edge: prediction for prediction in predictions}

    assert set(by_edge) == {("A", "B"), ("B", "A"), ("A", "C")}
    assert by_edge[("A", "B")].decision == SELECT
    assert by_edge[("B", "A")].decision == REJECT
    assert by_edge[("A", "C")].decision == ABSTAIN
    assert all(pred.evidence_channel == "hydraulic" for pred in predictions)


def test_hydraulic_baseline_rejects_truth_fields_in_observations() -> None:
    spec = hydraulic_only_baseline_spec()

    with pytest.raises(ValueError, match="Truth/reference field is forbidden"):
        spec.score(
            [("A", "B")],
            {"hydraulic": {"edges": {("A", "B"): {"truth": True, "head_drop": 1.0}}}},
        )


def test_edge_local_topology_baseline_uses_only_declared_topology_channel() -> None:
    spec = edge_local_topology_baseline_spec()
    observations = {
        "hydraulic": {
            "edges": {
                ("A", "B"): {"head_drop": -10.0},
                ("B", "C"): {"head_drop": 10.0},
            }
        },
        "topology": {
            "edges": {
                "A->B": {"local_support": 0.85},
                "B->C": {"local_support": 0.15},
                "C->D": {"local_support": 0.50},
            }
        },
    }

    predictions = spec.score([("A", "B"), ("B", "C"), ("C", "D")], observations)
    by_edge = {prediction.edge: prediction for prediction in predictions}

    assert by_edge[("A", "B")].decision == SELECT
    assert by_edge[("B", "C")].decision == REJECT
    assert by_edge[("C", "D")].decision == ABSTAIN
    assert all(pred.evidence_channel == "topology" for pred in predictions)


def test_permutation_control_changes_only_declared_evidence_channel() -> None:
    base = hydraulic_only_baseline_spec()
    control = permutation_control_baseline(base, evidence_channel="permuted_hydraulic")

    base_audit = base.to_audit_record()
    control_audit = control.to_audit_record()
    for section in ("tuning", "uncertainty", "abstention", "cost"):
        assert control_audit[section] == base_audit[section]
    assert base_audit["input_channels"] == ["hydraulic"]
    assert control_audit["input_channels"] == ["permuted_hydraulic"]
    assert control_audit["control"] == {
        "type": "evidence_channel_permutation",
        "baseline": "hydraulic_only_directional_gradient",
        "original_channel": "hydraulic",
        "permuted_channel": "permuted_hydraulic",
    }

    candidates = [("A", "B"), ("B", "A")]
    permuted_channel_data = {"node_heads": {"A": 10.0, "B": 9.0}}
    observations = {
        "hydraulic": {"node_heads": {"A": 1.0, "B": 100.0}},
        "permuted_hydraulic": permuted_channel_data,
    }

    base_on_permuted_data = base.score(candidates, {"hydraulic": permuted_channel_data})
    control_predictions = control.score(candidates, observations)

    assert [pred.probability for pred in control_predictions] == [
        pred.probability for pred in base_on_permuted_data
    ]
    assert [pred.decision for pred in control_predictions] == [
        pred.decision for pred in base_on_permuted_data
    ]
    assert all(
        pred.evidence_channel == "permuted_hydraulic"
        for pred in control_predictions
    )
