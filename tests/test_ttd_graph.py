import numpy as np
import pytest

from hydrosheaf.nuclear.ttd_graph import (
    MassTransportMap,
    audit_ttd_graph_compatibility,
)
from hydrosheaf.nuclear.ttd_identified import (
    AgeFunctional,
    TracerConstraint,
    solve_ttd_identified_set,
)


def _fixed_report(mass):
    ages = [0.0, 10.0]
    return solve_ttd_identified_set(
        ages,
        [TracerConstraint("bin-0", [1.0, 0.0], mass[0], 1.0)],
        [AgeFunctional("young", [1.0, 0.0], 1.0e-10)],
        sigma_multiplier=0.0,
    )


def test_transport_map_must_be_nonnegative_and_mass_conserving():
    with pytest.raises(ValueError, match="columns must sum to one"):
        MassTransportMap("bad", [[0.5, 0.0], [0.0, 0.5]], "independent", "test")
    with pytest.raises(ValueError, match="non-negative"):
        MassTransportMap("bad", [[1.1, 0.0], [-0.1, 1.0]], "independent", "test")


def test_identity_map_reports_compatible_frozen_local_sets():
    source = _fixed_report([0.25, 0.75])
    target = _fixed_report([0.25, 0.75])
    transport = MassTransportMap(
        "source-to-target", np.eye(2), "independent", "hydraulic test"
    )
    audit = audit_ttd_graph_compatibility(source, target, transport)

    assert audit.status == "COMPATIBLE"
    assert audit.minimum_l1_gap == pytest.approx(0.0)
    assert audit.tightening_authorized is False
    assert audit.metadata["conditioning_performed"] is False


def test_incompatible_map_returns_obstruction_without_tightening_local_sets():
    source = _fixed_report([1.0, 0.0])
    target = _fixed_report([0.0, 1.0])
    transport = MassTransportMap(
        "source-to-target", np.eye(2), "independent", "hydraulic test"
    )
    source_bound_before = source.bounds["young"].lower
    target_bound_before = target.bounds["young"].lower
    audit = audit_ttd_graph_compatibility(source, target, transport)

    assert audit.status == "INCOMPATIBLE"
    assert audit.minimum_l1_gap == pytest.approx(2.0)
    assert source.bounds["young"].lower == source_bound_before
    assert target.bounds["young"].lower == target_bound_before
    assert audit.tightening_authorized is False
