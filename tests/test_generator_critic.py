from __future__ import annotations

import sys
from pathlib import Path

from hydrosheaf.validation import (
    CriticSeverity,
    OBSERVATION_STRESS_SCENARIOS,
    apply_observation_stress,
    audit_generator_case,
    audit_source_independence,
    consolidate_head_evidence,
)

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from independent_lattice_generator import generate_independent_lattice  # noqa: E402
from hydrosheaf.validation.generator_critic import _nonnegative_tracer  # noqa: E402


def test_critic_checks_source_independence_and_seed_behaviour() -> None:
    case = generate_independent_lattice(9301)
    repeat = generate_independent_lattice(9301)
    alternate = generate_independent_lattice(9302)
    report = audit_generator_case(
        "analytic_lattice_v1",
        case,
        source_path=SCRIPTS / "independent_lattice_generator.py",
        repeat_case=repeat,
        alternate_case=alternate,
    )

    assert report.blockers == ()
    assert report.summary["BLOCKER"] == 0
    assert report.summary["MAJOR"] >= 1
    finding_ids = {finding.check_id for finding in report.findings}
    assert "TOPOLOGY_COMPLEXITY" not in finding_ids
    assert "LAYER_HETEROGENEITY" not in finding_ids


def test_critic_rejects_inference_imports(tmp_path: Path) -> None:
    source = tmp_path / "bad_generator.py"
    source.write_text("from hydrosheaf.api import fit_network_pipeline\n", encoding="utf-8")

    findings = audit_source_independence(source)

    assert any(finding.severity is CriticSeverity.BLOCKER for finding in findings)
    assert any(finding.check_id == "SOURCE_INDEPENDENCE" for finding in findings)


def test_critic_accepts_declared_blind_observation_stress() -> None:
    case = generate_independent_lattice(9301)
    consolidated, head_audit = consolidate_head_evidence(
        case.observations,
        model=case.provenance["head_channel_covariance_model"],
    )
    scenarios = {
        scenario: apply_observation_stress(
            consolidated,
            scenario,
            seed=9301 + index * 7_919,
        ).observations
        for index, scenario in enumerate(OBSERVATION_STRESS_SCENARIOS)
    }
    provenance = dict(case.provenance)
    provenance["head_channel_covariance_model"] = dict(head_audit["model"])
    report = audit_generator_case(
        "analytic_lattice_v1",
        case,
        source_path=SCRIPTS / "independent_lattice_generator.py",
        repeat_case=generate_independent_lattice(9301),
        alternate_case=generate_independent_lattice(9302),
        observation_scenarios=scenarios,
        provenance_override=provenance,
        covariance_consumed=True,
    )

    assert "MISSINGNESS_STRESS" not in {
        finding.check_id for finding in report.findings
    }


def test_critic_rejects_negative_tracer_values() -> None:
    assert _nonnegative_tracer(-1.0) is False
    assert _nonnegative_tracer("-2.5") is False
    assert _nonnegative_tracer("<-0.1") is False
    assert _nonnegative_tracer("<0.1") is True
