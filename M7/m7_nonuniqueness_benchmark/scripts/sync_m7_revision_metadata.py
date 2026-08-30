"""Synchronise artifact metadata with the 2026-07-30 reviewer revision."""

from __future__ import annotations

import json
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]
MANUSCRIPT = PROJECT / "manuscript"


def load(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def save(path: Path, value) -> None:
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def by_id(items, key="id"):
    return {item[key]: item for item in items}


def main() -> int:
    registry_path = MANUSCRIPT / "artifact_registry.json"
    registry = load(registry_path)
    reg = by_id(registry["artifacts"])
    reg["FIG-5"].update(
        path="figures/publication/figure5_ghana_supportability_boundary_m7_only.png",
        caption=(
            "Northern Ghana evidence and claim boundary using evidence from this study only. "
            "Panels a and b map observed evidence to defensible claims; panel c "
            "reports 160 wells, 140 complete seasonal pairs and 320 seasonal "
            "observations; panel d reports the truth-free seasonal hold-forward "
            "comparison. No external field-transfer experiment is included."
        ),
        technical_description=(
            "Built only from results/m7_3_locked/ghana_data_scope_audit.json and "
            "results/supporting_validation/field_prequential outputs."
        ),
    )
    reg["TAB-1"].update(
        objective_ids=["OBJ-1", "OBJ-2", "OBJ-3", "OBJ-4", "OBJ-5"],
        caption="Design and claim map for all seven audits.",
        technical_description=(
            "Synthesis of the locked process-based integration, public-pipeline, "
            "representation-benchmark and estimator-diagnostic protocols plus the "
            "Ghana data-scope audit; introduces no fitted statistic."
        ),
        interpretation=(
            "Separates generator systems, samples, comparators, lock points, "
            "multiplicity families and permitted claims."
        ),
    )
    reg["FIG-6"].update(
        caption=(
            "Sheaf-versus-weighted-graph representation comparison with simultaneous 95% intervals "
            "controlling all 120 published contrasts; intended printed width "
            "7.08 in and minimum label size 8 pt."
        ),
        technical_description=(
            "Balanced 2 x 2 summary generated from immutable representation-benchmark case metrics and "
            "the post-review full-family max-z audit."
        ),
    )
    reg["FIG-7"].update(
        caption=(
            "Local-first/global-fallback estimator diagnostic with simultaneous 95% intervals controlling "
            "all 560 published contrasts. The selected nominal hybrid has local "
            "weight 1.0 and is local-first/global-fallback."
        ),
        interpretation=(
            "Native maps and conflict localisation remain supported, but overall "
            "edge-local ranking/calibration gains and scenario ranking gains do not "
            "survive full-family correction."
        ),
    )
    reg["TAB-8"]["caption"] = (
        "Representation-benchmark model means across 64 cases; full-family contrasts are in Table S10."
    )
    reg["TAB-9"]["caption"] = (
        "Estimator-diagnostic model means across 128 cases; selected local-first arms have local "
        "weight 1.0 and full-family contrasts are in Table S11."
    )
    new_registry = [
        {
            "id": "TAB-S10",
            "kind": "table",
            "path": "tables/publication/tableS10_m7_4_multiplicity_adjusted.csv",
            "objective_ids": ["OBJ-5"],
            "caption": "All 120 representation-benchmark contrasts with max-z simultaneous 95% intervals.",
            "technical_description": "10,000 shared case-block bootstrap resamples; one complete published family.",
            "interpretation": "Controls family-wise selection across every published representation-benchmark scenario, comparator and metric.",
            "validation_status": "PASS",
        },
        {
            "id": "TAB-S11",
            "kind": "table",
            "path": "tables/publication/tableS11_m7_5_multiplicity_adjusted.csv",
            "objective_ids": ["OBJ-5"],
            "caption": "All 560 estimator-diagnostic contrasts with max-z simultaneous 95% intervals.",
            "technical_description": "10,000 shared case-block bootstrap resamples; one complete published family.",
            "interpretation": "Removes unsupported scenario-ranking claims while retaining map and conflict signals that survive family correction.",
            "validation_status": "PASS",
        },
        {
            "id": "TAB-S12",
            "kind": "table",
            "path": "tables/publication/tableS12_precision_and_power.csv",
            "objective_ids": ["OBJ-1", "OBJ-2", "OBJ-3", "OBJ-5"],
            "caption": "Post-review precision and future-replication planning audit.",
            "technical_description": "20,000 empirical simulations with explicitly post-review, non-field-validated margins.",
            "interpretation": "Separates realised statistical precision from practical magnitude without retroactive power claims.",
            "validation_status": "PASS",
        },
        {
            "id": "TAB-S13",
            "kind": "table",
            "path": "tables/publication/tableS13_public_pipeline_selection.csv",
            "objective_ids": ["OBJ-5"],
            "caption": "Public-pipeline selection rule, thresholds and confusion-count audit.",
            "technical_description": "Deterministic audit of all recovered candidates across four pipeline arms.",
            "interpretation": "Explains identical conditional F1 and distinguishes conditional scoring from end-to-end candidate recall.",
            "validation_status": "PASS",
        },
    ]
    existing = by_id(registry["artifacts"])
    registry["artifacts"].extend(item for item in new_registry if item["id"] not in existing)
    save(registry_path, registry)

    plan_path = MANUSCRIPT / "analysis_plan.json"
    plan = load(plan_path)
    entries = by_id(plan["artifacts"])
    entries["FIG-5"].update(
        path="figures/publication/figure5_ghana_supportability_boundary_m7_only.png",
        inputs=[
            "results/m7_3_locked/ghana_data_scope_audit.json",
            "results/supporting_validation/field_prequential_summary.csv",
            "results/supporting_validation/field_prequential_audit.json",
        ],
        validation_checks=[
            "all four panels derive only from this study's Ghana scope and truth-free hold-forward audits",
            "no external field-transfer result is imported",
            "no panel implies field validation of topology, age or reaction truth",
        ],
    )
    entries["TAB-1"].update(
        objective_ids=["OBJ-1", "OBJ-2", "OBJ-3", "OBJ-4", "OBJ-5"],
        generator="scripts/assemble_m7_manuscript.py",
        inputs=[
            "docs/m7_3_protocol.md",
            "M7_SYSTEM_ACCEPTANCE_PROTOCOL.md",
            "M7_SHEAF_VS_GRAPH_PROTOCOL.md",
            "M7_ROBUST_HYBRID_SHEAF_PROTOCOL.md",
            "results/m7_3_locked/ghana_data_scope_audit.json",
        ],
        validation_checks=[
            "all seven audits are present once",
            "generator systems, samples, comparators, locks and multiplicity families match source protocols",
        ],
    )
    fig6_input = "tables/publication/tableS10_m7_4_multiplicity_adjusted.csv"
    if fig6_input not in entries["FIG-6"]["inputs"]:
        entries["FIG-6"]["inputs"].append(fig6_input)
    entries["FIG-6"]["validation_checks"] = [
        "identity-limit sheaf and graph predictions are exactly equal",
        "intervals control all 120 published representation-benchmark contrasts as one family",
        "printed width 7.08 in and minimum label size 8 pt are recorded",
    ]
    fig7_input = "tables/publication/tableS11_m7_5_multiplicity_adjusted.csv"
    if fig7_input not in entries["FIG-7"]["inputs"]:
        entries["FIG-7"]["inputs"].append(fig7_input)
    entries["FIG-7"]["validation_checks"] = [
        "all values trace to immutable locked-test tables",
        "intervals control all 560 published estimator-diagnostic contrasts as one family",
        "selected local weight 1.0 is disclosed",
        "printed width 7.08 in and minimum label size 8 pt are recorded",
    ]
    for item in new_registry:
        if item["id"] not in entries:
            plan["artifacts"].append(
                {
                    "id": item["id"],
                    "kind": "table",
                    "path": item["path"],
                    "objective_ids": item["objective_ids"],
                    "generator": "scripts/post_review_statistical_audit.py",
                    "inputs": ["results/post_review_audit_20260730/manifest.json"],
                    "validation_checks": [item["technical_description"]],
                }
            )
    save(plan_path, plan)

    validation_path = MANUSCRIPT / "artifact_validation.json"
    validation = load(validation_path)
    checks = by_id(validation, key="artifact_id")
    checks["FIG-5"]["checks"].update(
        scientific="All panels use this study's Ghana scope or truth-free hold-forward evidence; external field-transfer results are absent.",
        reproducibility="Deterministic from locked study-scope and field_prequential records.",
    )
    checks["FIG-6"]["checks"].update(
        statistical="Simultaneous intervals verified against the complete 120-contrast max-z audit.",
        presentation="Balanced layout visually checked at 7.08 in intended width with minimum labels of 8 pt.",
    )
    checks["FIG-7"]["checks"].update(
        statistical="Simultaneous intervals verified against the complete 560-contrast max-z audit.",
        scientific="Local weight 1.0, failed overall gate, supported native-map control and conflict localisation are all explicit.",
        presentation="Balanced layout visually checked at 7.08 in intended width with minimum labels of 8 pt.",
    )
    new_validation = {
        "TAB-S10": "Complete 120-contrast representation-benchmark family and max-z intervals verified against immutable inputs.",
        "TAB-S11": "Complete 560-contrast estimator-diagnostic family and max-z intervals verified against immutable inputs.",
        "TAB-S12": "Planning status and post-review margins are explicit on every row.",
        "TAB-S13": "Candidate, selection and confusion counts reproduce the public-pipeline case records.",
    }
    for artifact_id, note in new_validation.items():
        if artifact_id not in checks:
            validation.append(
                {
                    "artifact_id": artifact_id,
                    "status": "PASS",
                    "checks": {
                        "statistical": note,
                        "scientific": "Claim wording in the main text matches the audited table.",
                        "presentation": "Compact supplement view points to the complete CSV.",
                        "reproducibility": "Generated by scripts/post_review_statistical_audit.py from hash-verified immutable records.",
                    },
                    "reviewer_notes": "Added in the 2026-07-30 technical revision.",
                }
            )
    save(validation_path, validation)

    proposal_path = PROJECT / "proposal.normalized.json"
    proposal = load(proposal_path)
    proposal["title"] = (
        "Conditional evidence integration and the incremental contribution of sheaf "
        "structure in controlled-synthetic groundwater benchmarks"
    )
    proposal["scope"] = (
        "The study uses a MODFLOW/MODPATH flow-tracer-chemistry generator for the "
        "process-based integration and public-pipeline audits, a separate scalar "
        "affine generator for the representation benchmark and estimator diagnostic, and a "
        "descriptive Northern Ghana scope audit. Results are conditional capability "
        "tests, not whole-framework or field validation."
    )
    proposal["manuscript_requirements"]["tables"] = [
        "TAB-1", "TAB-7", "TAB-8", "TAB-9",
        *[f"TAB-S{i}" for i in range(1, 14)],
    ]
    save(proposal_path, proposal)

    methods_path = MANUSCRIPT / "methods" / "methods_manifest.json"
    methods = load(methods_path)
    required = methods["required_supplementary_sections"]
    for name in [
        "Generator construct validity and non-transfer boundary",
        "Post-review precision and practical-magnitude audit",
    ]:
        if name not in required:
            required.append(name)
    save(methods_path, methods)
    print("Synchronized M7 reviewer-revision metadata")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
