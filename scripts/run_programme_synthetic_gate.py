"""Run the first executable, synthetic-only programme validation gate.

This is an orchestration adapter, not a second inference implementation.  It
uses the independent M7 v3 multilayer generator and the public HydroSheaf pipeline, while
the shared truth-blind and scoring helpers live in ``hydrosheaf.validation``.
The output is a development execution record; it is not a locked superiority
benchmark and it does not perform field validation.
"""

from __future__ import annotations

import argparse
from copy import deepcopy
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import shutil
import subprocess
import sys
from typing import Any, Mapping, Sequence

import numpy as np

REPO = Path(__file__).resolve().parents[1]
M7_SCRIPTS = REPO / "M7" / "m7_nonuniqueness_benchmark" / "scripts"
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
if str(M7_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(M7_SCRIPTS))

from hydrosheaf.api import fit_network_pipeline  # noqa: E402
from hydrosheaf.inference.network_fit import infer_edges  # noqa: E402
from hydrosheaf.nuclear.input_history import (  # noqa: E402
    InputHistory,
    build_default_tritium_input,
)
from hydrosheaf.nuclear.tracer_inputs import extend_history_to_present  # noqa: E402
from hydrosheaf.validation import (  # noqa: E402
    ScenarioPrediction,
    DecisionKind,
    IdentifiabilityStatus,
    ProgrammeDecision,
    ProgrammeRun,
    ProgrammeStage,
    StageStatus,
    build_discrepancy_report,
    finite_numeric_mapping,
    prepare_truth_blind_rows,
    programme_stages_from_status,
    required_stages_completed,
    score_age_posteriors,
    score_reaction_families,
    score_topology,
)

from independent_modflow_generator_v3 import (  # noqa: E402
    IndependentAquiferV3,
    generate_independent_aquifer_v3,
)
from strong_inference import strong_config  # noqa: E402

RUN_ID = "RUN-PROGRAMME-SYNTHETIC-DEV-20260801-01"
QUICK_SEEDS = (9101,)
DEVELOPMENT_SEEDS = (9101, 9102, 9103)
CONDITIONS = ("hydraulic_only", "local_age", "full_sheaf", "age_permuted")
DISCREPANCY_SCENARIOS = ("age_history_tropical",)
RUN_CONDITIONS = CONDITIONS + DISCREPANCY_SCENARIOS
PROGRAMME_STAGES_NOT_IMPLEMENTED = (
    "next_measurement",
)
PROGRAMME_SOURCE_FILES = (
    "hydrosheaf/validation/__init__.py",
    "hydrosheaf/validation/programme_contract.py",
    "hydrosheaf/validation/programme_gate.py",
    "hydrosheaf/validation/discrepancy.py",
    "M7/m7_nonuniqueness_benchmark/scripts/independent_modflow_generator_v3.py",
    "scripts/run_programme_synthetic_gate.py",
)
DEFAULT_OUTPUT = REPO / ".codex_work" / "programme-validation" / RUN_ID
DEFAULT_SIMULATOR_WORKSPACE = REPO / ".codex_work" / "programme-validation-simulators"
DEFAULT_BIN_DIR = REPO / ".codex_work" / "modflow-bin"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _git_revision() -> str:
    try:
        return subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=REPO,
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
    except (FileNotFoundError, subprocess.CalledProcessError):
        return "UNAVAILABLE"


def _git_worktree_dirty() -> bool:
    try:
        result = subprocess.run(
            ["git", "status", "--porcelain", "--untracked-files=all"],
            cwd=REPO,
            check=True,
            capture_output=True,
            text=True,
        )
        return bool(result.stdout.strip())
    except (FileNotFoundError, subprocess.CalledProcessError):
        return True


def _json_default(value: object) -> object:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.generic):
        return value.item()
    raise TypeError(f"Object is not JSON serialisable: {type(value).__name__}")


def _safe_output(output: Path, *, overwrite: bool) -> Path:
    """Allow replacement only inside the generated-workspace directory."""

    generated_root = (REPO / ".codex_work").resolve()
    resolved = output.resolve()
    if not resolved.is_relative_to(generated_root):
        raise ValueError("Programme gate output must stay under .codex_work.")
    if resolved == generated_root:
        raise ValueError("Refusing to use .codex_work itself as an output directory.")
    if resolved.exists():
        if not overwrite:
            raise FileExistsError(
                f"Output exists: {resolved}; pass --overwrite to replace this generated run."
            )
        shutil.rmtree(resolved)
    resolved.mkdir(parents=True, exist_ok=False)
    return resolved


def _condition_config(condition: str, *, head_sigma_m: float | None = None):
    config = strong_config(phreeqc_enabled=False)
    config.residence_time_tracer = "3H"
    config.edge_max_neighbors = 3
    config.sheaf_soft_beta = 1.0
    config.hydraulic_hodge_head_key = "head_meas"
    if head_sigma_m is not None and float(head_sigma_m) > 0.0:
        config.edge_sigma_meas = float(head_sigma_m)
    if condition == "hydraulic_only":
        config.sheaf_age_enabled = False
    elif condition == "local_age":
        config.sheaf_weight_global = 0.0
        config.sheaf_max_iter = 1
    else:
        config.sheaf_weight_global = 1.0
        config.sheaf_max_iter = 3
    config.validate()
    return config


def _input_history_for_condition(condition: str) -> InputHistory:
    if condition == "age_history_tropical":
        return extend_history_to_present(
            build_default_tritium_input("tropical"),
            "tropical",
            present_year=2026.0,
        )
    return InputHistory(
        np.asarray([1850.0, 2035.0]),
        np.asarray([6.2, 6.2]),
    )


def _prepare_rows(
    case: IndependentAquiferV3,
    *,
    permute_age: bool,
    observations: Sequence[Mapping[str, object]] | None = None,
) -> list[dict[str, object]]:
    return prepare_truth_blind_rows(
        case.observations if observations is None else observations,
        tracer_source="tritium_TU",
        tracer_target="3H",
        permute_tracer_seed=case.seed + 880_301 if permute_age else None,
    )


def _run_condition(
    case: IndependentAquiferV3,
    condition: str,
    *,
    simulator_workspace: Path,
    mf6_executable: Path,
    mp7_executable: Path,
    age_samples: int,
    observations: Sequence[Mapping[str, object]] | None = None,
    observation_scenario: str = "complete",
    head_sigma_m: float | None = None,
    head_evidence_audit: Mapping[str, object] | None = None,
) -> dict[str, Any]:
    rows = _prepare_rows(
        case,
        permute_age=condition == "age_permuted",
        observations=observations,
    )
    config = _condition_config(condition, head_sigma_m=head_sigma_m)
    candidates = infer_edges(rows, method="probabilistic", config=config)
    sheaf_enabled = condition != "hydraulic_only"
    required = ["network_fit"]
    if sheaf_enabled:
        required = ["nuclear_age", "sheaf_refinement", "network_fit"]

    network_results, extras = fit_network_pipeline(
        rows,
        deepcopy(candidates),
        config,
        auto_disable_missing=False,
        sheaf_refinement_enabled=sheaf_enabled,
        nuclear_inference_options={
            "input_hist": _input_history_for_condition(condition),
            "sampler": "grid",
            "n_samples": int(age_samples),
            "n_chains": 4,
            "random_seed": int(case.seed) + 10_000 * RUN_CONDITIONS.index(condition),
            "max_age_years": 200.0,
            "root_age_prior_median_years": 20.0,
            "root_age_prior_log_sigma": 1.2,
        },
        strict_stage_completion=True,
        required_stages=required,
    )

    stage_status = extras["stage_status"]
    age_results = extras.get("nuclear_results")
    topology = score_topology(case.true_edges, candidates, extras["edges"])
    age = score_age_posteriors(case.true_ages_years, age_results)
    reaction = score_reaction_families(case.true_processes, network_results)
    metrics_finite = finite_numeric_mapping(topology)
    if age.get("status") == "scored":
        metrics_finite = metrics_finite and bool(age.get("all_finite", False))
    if reaction.get("status") == "scored":
        reaction_metrics = reaction.get("metrics", {})
        metrics_finite = metrics_finite and isinstance(reaction_metrics, Mapping)
        if isinstance(reaction_metrics, Mapping):
            metrics_finite = metrics_finite and finite_numeric_mapping(reaction_metrics)

    diagnostics: Mapping[str, object] = {}
    if isinstance(age_results, Mapping):
        raw_diagnostics = age_results.get("_diagnostics", {})
        if isinstance(raw_diagnostics, Mapping):
            diagnostics = raw_diagnostics

    age_predictions: dict[str, dict[str, float]] = {}
    if isinstance(age_results, Mapping):
        for node, result in age_results.items():
            if node == "_diagnostics" or not isinstance(result, Mapping):
                continue
            required_age_fields = ("mean_age_years", "age_95_low", "age_95_high")
            if all(field in result for field in required_age_fields):
                age_predictions[str(node)] = {
                    field: float(result[field]) for field in required_age_fields
                }

    return {
        "seed": int(case.seed),
        "condition": condition,
        "observation_scenario": str(observation_scenario),
        "n_observations": len(rows),
        "required_stages": required,
        "stage_status": stage_status,
        "stages_complete": required_stages_completed(stage_status, required),
        "nuclear_converged": (
            bool(diagnostics.get("converged")) if sheaf_enabled else None
        ),
        "truth_blind": True,
        "nonempty_outputs": bool(rows and candidates and extras["edges"]),
        "metrics_finite": bool(metrics_finite),
        "topology": topology,
        "hydrosheaf_candidate_edges": [
            [str(edge.u), str(edge.v)] for edge in candidates
        ],
        "hydrosheaf_selected_edges": [
            [str(edge.u), str(edge.v)] for edge in extras["edges"]
        ],
        "age": age,
        "age_predictions": age_predictions,
        "reaction": reaction,
        "reaction_stage": {
            "status": "completed",
            "detail": f"scored reaction outputs for {len(network_results)} fitted edges",
        },
        "generator_workspace": str(simulator_workspace),
        "simulator_executables": {
            "mf6": str(mf6_executable),
            "mp7": str(mp7_executable),
        },
        "head_evidence_audit": dict(head_evidence_audit or {}),
    }


def _build_discrepancy_reports(
    seed: int,
    condition_rows: Mapping[str, Mapping[str, object]],
) -> list[dict[str, object]]:
    """Compare nominal and alternative input-history scenario outputs."""

    nominal = condition_rows.get("full_sheaf", {}).get("age_predictions", {})
    alternative = condition_rows.get("age_history_tropical", {}).get(
        "age_predictions", {}
    )
    if not isinstance(nominal, Mapping) or not isinstance(alternative, Mapping):
        return []

    reports: list[dict[str, object]] = []
    for node in sorted(set(nominal) & set(alternative)):
        nominal_row = nominal[node]
        alternative_row = alternative[node]
        if not isinstance(nominal_row, Mapping) or not isinstance(alternative_row, Mapping):
            continue
        fields = ("mean_age_years", "age_95_low", "age_95_high")
        if not all(field in nominal_row and field in alternative_row for field in fields):
            continue
        report = build_discrepancy_report(
            f"age:{node}",
            [
                ScenarioPrediction(
                    "nominal_input_history",
                    estimate=float(nominal_row["mean_age_years"]),
                    lower=float(nominal_row["age_95_low"]),
                    upper=float(nominal_row["age_95_high"]),
                    scenario_family="age_forward_model",
                    discrepancy_tags=("steady_input_history",),
                    provenance={"seed": int(seed), "condition": "full_sheaf"},
                ),
                ScenarioPrediction(
                    "tropical_input_history",
                    estimate=float(alternative_row["mean_age_years"]),
                    lower=float(alternative_row["age_95_low"]),
                    upper=float(alternative_row["age_95_high"]),
                    scenario_family="age_forward_model",
                    discrepancy_tags=("tropical_input_history",),
                    provenance={
                        "seed": int(seed),
                        "condition": "age_history_tropical",
                    },
                ),
            ],
            disagreement_threshold=0.25,
            scale_floor=5.0,
            compatible_hypotheses=("nominal_input_history", "tropical_input_history"),
        )
        reports.append({"seed": int(seed), **report.to_dict()})
    return reports


def _write_json(path: Path, payload: object) -> None:
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True, default=_json_default),
        encoding="utf-8",
    )


def run_benchmark(
    *,
    output: Path,
    simulator_workspace: Path,
    mf6_executable: Path,
    mp7_executable: Path,
    seeds: Sequence[int],
    age_samples: int = 600,
    overwrite: bool = False,
) -> dict[str, object]:
    if int(age_samples) < 100:
        raise ValueError("age_samples must be at least 100 for a smoke gate.")
    output = _safe_output(Path(output), overwrite=overwrite)
    simulator_workspace = Path(simulator_workspace)
    mf6_executable = Path(mf6_executable).resolve()
    mp7_executable = Path(mp7_executable).resolve()
    if not mf6_executable.exists() or not mp7_executable.exists():
        raise FileNotFoundError("Both official mf6.exe and mp7.exe are required.")

    simulator_workspace.mkdir(parents=True, exist_ok=True)
    metric_rows: list[dict[str, Any]] = []
    programme_stages = []
    generator_records: dict[str, Mapping[str, object]] = {}
    errors: list[dict[str, object]] = []
    discrepancy_reports: list[dict[str, object]] = []
    discrepancy_stage_status: dict[str, bool] = {}

    for seed in map(int, seeds):
        try:
            case = generate_independent_aquifer_v3(
                seed,
                simulator_workspace / f"case_{seed}",
                mf6_executable,
                mp7_executable,
            )
            generator_records[str(seed)] = dict(case.provenance)
            if case.provenance.get("imports_hydrosheaf") is not False:
                raise RuntimeError("Independent generator import audit failed.")
            condition_rows: dict[str, dict[str, object]] = {}
            for condition in RUN_CONDITIONS:
                try:
                    result = _run_condition(
                        case,
                        condition,
                        simulator_workspace=simulator_workspace,
                        mf6_executable=mf6_executable,
                        mp7_executable=mp7_executable,
                        age_samples=age_samples,
                    )
                    metric_rows.append(result)
                    condition_rows[condition] = result
                    programme_stages.extend(
                        programme_stages_from_status(
                            result["stage_status"],
                            name_prefix=f"case_{seed}:{condition}:",
                        )
                    )
                    programme_stages.append(
                        ProgrammeStage(
                            name=f"case_{seed}:{condition}:chemistry_reaction",
                            status=StageStatus.COMPLETED,
                            detail=result["reaction_stage"]["detail"],
                        )
                    )
                except Exception as exc:  # record failed cases before continuing
                    detail = f"{type(exc).__name__}: {exc}"
                    errors.append(
                        {"seed": seed, "condition": condition, "error": detail}
                    )
                    programme_stages.append(
                        ProgrammeStage(
                            name=f"case_{seed}:{condition}:runner",
                            status=StageStatus.FAILED,
                            detail=detail,
                        )
                    )
            reports = _build_discrepancy_reports(seed, condition_rows)
            discrepancy_reports.extend(reports)
            discrepancy_ok = bool(reports)
            discrepancy_stage_status[str(seed)] = discrepancy_ok
            programme_stages.append(
                ProgrammeStage(
                    name=f"case_{seed}:model_discrepancy",
                    status=(StageStatus.COMPLETED if discrepancy_ok else StageStatus.FAILED),
                    detail=(
                        f"reported {len(reports)} age scenario comparisons"
                        if discrepancy_ok
                        else "No comparable age scenario outputs were produced."
                    ),
                )
            )
        except Exception as exc:
            detail = f"{type(exc).__name__}: {exc}"
            errors.append({"seed": seed, "condition": "generation", "error": detail})
            programme_stages.append(
                ProgrammeStage(
                    name=f"case_{seed}:generation",
                    status=StageStatus.FAILED,
                    detail=detail,
                )
            )
            discrepancy_stage_status[str(seed)] = False

    for stage_name in PROGRAMME_STAGES_NOT_IMPLEMENTED:
        programme_stages.append(
            ProgrammeStage(
                name=f"programme:{stage_name}",
                status=StageStatus.SKIPPED,
                detail="Not wired into the root-level synthetic runner yet.",
            )
        )

    independent_generator = bool(generator_records) and all(
        record.get("imports_hydrosheaf") is False
        for record in generator_records.values()
    )
    all_stages_complete = bool(metric_rows) and all(
        bool(row["stages_complete"]) for row in metric_rows
    )
    all_metrics_finite = bool(metric_rows) and all(
        bool(row["metrics_finite"]) for row in metric_rows
    )
    all_outputs_nonempty = bool(metric_rows) and all(
        bool(row["nonempty_outputs"]) for row in metric_rows
    )
    truth_blind = bool(metric_rows) and all(bool(row["truth_blind"]) for row in metric_rows)
    all_discrepancy_stages_complete = bool(discrepancy_stage_status) and all(
        discrepancy_stage_status.values()
    )
    discrepancy_reports_finite = bool(discrepancy_reports) and all(
        finite_numeric_mapping(
            {
                "estimate": report["estimate"],
                "lower": report["lower"],
                "upper": report["upper"],
                "scenario_range": report["scenario_range"],
                "relative_disagreement": report["relative_disagreement"],
            }
        )
        for report in discrepancy_reports
    )
    execution_gate = bool(
        not errors
        and independent_generator
        and truth_blind
        and all_stages_complete
        and all_metrics_finite
        and all_outputs_nonempty
        and all_discrepancy_stages_complete
        and discrepancy_reports_finite
    )

    # The module is fixed to the independent generator adapter.  The observed
    # provenance audit remains a separate execution-gate check, so a failed
    # generation can still be recorded as a failed programme run rather than
    # being lost while constructing the result contract.
    programme_run = ProgrammeRun(
        run_id=RUN_ID,
        generator="independent_multilayer_modflow_modpath_mixing_v3",
        generator_independent=True,
        stages=tuple(programme_stages),
        decisions=(
            ProgrammeDecision(
                decision_kind=DecisionKind.ABSTAIN,
                identifiability=IdentifiabilityStatus.UNKNOWN,
                reason=(
                    "The current synthetic gate executes inference and scoring, "
                    "but the joint next-measurement policy is not yet implemented."
                ),
                scenario_count=0,
                provenance={"field_validation": "deferred"},
            ),
        ),
        claim_boundary=(
            "development controlled-synthetic execution and conditional scoring; "
            "not field validation or universal superiority"
        ),
        truth_released_for_scoring=False,
    )

    _write_json(output / "case_metrics.json", metric_rows)
    _write_json(output / "discrepancy_reports.json", discrepancy_reports)
    _write_json(output / "generator_provenance.json", generator_records)
    _write_json(output / "stage_status.json", [stage.to_dict() for stage in programme_stages])
    _write_json(output / "programme_record.json", programme_run.to_dict())
    _write_json(
        output / "execution_gate.json",
        {
            "status": "PASS" if execution_gate else "FAIL",
            "gate_scope": (
                "public topology-age-reaction pipeline plus declared "
                "input-history discrepancy reporting"
            ),
            "checks": {
                "independent_generator": independent_generator,
                "truth_blind": truth_blind,
                "all_required_stages_completed": all_stages_complete,
                "all_metrics_finite": all_metrics_finite,
                "all_outputs_nonempty": all_outputs_nonempty,
                "all_discrepancy_stages_completed": all_discrepancy_stages_complete,
                "discrepancy_reports_finite": discrepancy_reports_finite,
                "no_recorded_errors": not errors,
            },
            "programme_workflow_status": "INCOMPLETE",
            "missing_programme_stages": list(PROGRAMME_STAGES_NOT_IMPLEMENTED),
            "performance_gate": "DEFERRED_UNTIL_LOCKED_SYNTHETIC_TEST",
            "field_validation": "DEFERRED_UNTIL_FIELD_DATA_READY",
        },
    )
    _write_json(output / "errors.json", errors)

    manifest: dict[str, object] = {
        "run_id": RUN_ID,
        "status": "PASS" if execution_gate else "FAIL",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "git_revision": _git_revision(),
        "git_worktree_dirty": _git_worktree_dirty(),
        "source_hashes": {
            relative: _sha256(REPO / relative) for relative in PROGRAMME_SOURCE_FILES
        },
        "python": sys.version,
        "platform": platform.platform(),
        "seeds": list(map(int, seeds)),
        "conditions": list(CONDITIONS),
        "discrepancy_scenarios": list(DISCREPANCY_SCENARIOS),
        "missing_programme_stages": list(PROGRAMME_STAGES_NOT_IMPLEMENTED),
        "age_samples": int(age_samples),
        "generator": "independent_multilayer_modflow_modpath_mixing_v3",
        "generator_independent": independent_generator,
        "discrepancy_report_count": len(discrepancy_reports),
        "model_disagreement_count": sum(
            report["identifiability"] == IdentifiabilityStatus.MODEL_DISAGREEMENT.value
            for report in discrepancy_reports
        ),
        "truth_blind_contract": (
            "observations were copied and audited before generator truth was used "
            "for post-inference scoring"
        ),
        "claim_boundary": (
            "development controlled-synthetic execution; no field-validation or "
            "universal-superiority claim"
        ),
        "errors": errors,
        "artifacts": {},
    }
    for path in sorted(output.iterdir()):
        if path.is_file() and path.name != "run_manifest.json":
            manifest["artifacts"][path.name] = _sha256(path)  # type: ignore[index]
    _write_json(output / "run_manifest.json", manifest)
    return {"manifest": manifest, "execution_gate": execution_gate}


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--simulator-workspace", type=Path, default=DEFAULT_SIMULATOR_WORKSPACE
    )
    parser.add_argument("--bin-dir", type=Path, default=DEFAULT_BIN_DIR)
    parser.add_argument("--seed", action="append", type=int)
    parser.add_argument("--age-samples", type=int, default=600)
    parser.add_argument("--quick", action="store_true", help="one non-claim-bearing smoke case")
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    seeds = tuple(args.seed) if args.seed else (QUICK_SEEDS if args.quick else DEVELOPMENT_SEEDS)
    output = args.output
    if args.quick and output == DEFAULT_OUTPUT:
        output = DEFAULT_OUTPUT.parent / f"{RUN_ID}-SMOKE"
    result = run_benchmark(
        output=output,
        simulator_workspace=args.simulator_workspace,
        mf6_executable=args.bin_dir / "mf6.exe",
        mp7_executable=args.bin_dir / "mp7.exe",
        seeds=seeds,
        age_samples=args.age_samples,
        overwrite=args.overwrite,
    )
    print(json.dumps(result, indent=2, sort_keys=True, default=_json_default))
    return 0 if bool(result["execution_gate"]) else 1


if __name__ == "__main__":
    raise SystemExit(main())
