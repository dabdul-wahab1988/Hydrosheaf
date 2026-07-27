"""Run the locked M7.3 conditional-integration benchmark."""

from __future__ import annotations

import argparse
from copy import deepcopy
import json
from pathlib import Path
import shutil
import subprocess
import sys
from typing import Dict, Mapping, Sequence

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from hydrosheaf.graph.types import Edge  # noqa: E402
from hydrosheaf.inference.network_fit import fit_network  # noqa: E402
from hydrosheaf.phreeqc.runner import run_phreeqc  # noqa: E402
from independent_modflow_generator import (  # noqa: E402
    ION_ORDER,
    IndependentAquifer,
    generate_independent_aquifer,
)
from run_supporting_validation import (  # noqa: E402
    _json_default,
    _label_rows,
    _truth_family,
    _write_case_artifacts,
)
from strong_inference import (  # noqa: E402
    StrongInferenceResult,
    _reaction_family,
    run_strong_inference,
    strong_config,
)

from m7_3_analysis import (  # noqa: E402
    CORE_IONS,
    ENHANCED_IONS,
    audit_ghana_workbook,
    bootstrap_evidence_contrasts,
    bootstrap_topology_age_contrasts,
    evaluate_evidence_conditions,
    fit_evidence_models,
    reaction_support_summary,
    topology_age_sensitivity,
)

DEV_SEEDS = tuple(range(5201, 5207))
TEST_SEEDS = tuple(range(5301, 5313))
SMOKE_DEV_SEEDS = (5901, 5902)
SMOKE_TEST_SEEDS = (5911, 5912, 5913)
DEFAULT_OUTPUT = (
    REPO_ROOT / "M7" / "m7_nonuniqueness_benchmark" / "results" / "m7_3_locked"
)
DEFAULT_SIMULATOR_WORKSPACE = REPO_ROOT / ".codex_work" / "m7_3_simulators"
DEFAULT_BIN_DIR = REPO_ROOT / ".codex_work" / "modflow-bin"
GHANA_WORKBOOK = REPO_ROOT / "data" / "FieldData" / "NorthenGhana" / "NorthernGhana.xlsx"


def _git_head() -> str:
    completed = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=REPO_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    return completed.stdout.strip()


def _truth_edges(case: IndependentAquifer) -> list[Edge]:
    return [
        Edge(
            edge_id=f"{u}->{v}",
            u=str(u),
            v=str(v),
            attrs={
                "edge_source_tier": "externally_generated_truth_for_reaction_scoring"
            },
        )
        for u, v in case.true_edges
    ]


def _perturb_chemistry(
    observations: Sequence[Mapping[str, object]],
    *,
    seed: int,
    bootstrap_index: int,
    sigma: float = 0.03,
) -> list[Dict[str, object]]:
    rng = np.random.default_rng(int(seed) * 100_003 + int(bootstrap_index) * 1009 + 37)
    rows = [dict(row) for row in observations]
    for row in rows:
        multipliers = rng.lognormal(0.0, float(sigma), len(ION_ORDER))
        for ion, multiplier in zip(ION_ORDER, multipliers):
            row[ion] = max(1.0e-8, float(row[ion]) * float(multiplier))
    return rows


def _reaction_bootstrap_case(
    case: IndependentAquifer,
    *,
    n_bootstrap: int,
) -> pd.DataFrame:
    """Bootstrap reaction support on externally true edges only."""

    original = [dict(row) for row in case.observations]
    true_edges = _truth_edges(case)
    phreeqc_config = strong_config(
        phreeqc_enabled=True,
        measured_ions=ENHANCED_IONS,
    )
    phreeqc_results = run_phreeqc(original, phreeqc_config)
    tiers = {
        "core": strong_config(
            phreeqc_enabled=True,
            measured_ions=CORE_IONS,
        ),
        "enhanced": phreeqc_config,
    }
    rows = []
    for bootstrap_index in range(int(n_bootstrap)):
        perturbed = _perturb_chemistry(
            original,
            seed=case.seed,
            bootstrap_index=bootstrap_index,
        )
        for tier, config in tiers.items():
            fitted = fit_network(
                perturbed,
                deepcopy(true_edges),
                config,
                phreeqc_results=phreeqc_results,
            )
            by_edge = {result.edge_id: result for result in fitted}
            for edge in true_edges:
                result = by_edge.get(edge.edge_id)
                if result is None or not result.z_labels:
                    dominant = None
                else:
                    extent_map = {
                        str(label): float(extent)
                        for label, extent in zip(
                            result.z_labels,
                            result.z_extents,
                        )
                    }
                    dominant = (
                        max(extent_map, key=lambda label: abs(extent_map[label]))
                        if extent_map
                        else None
                    )
                true_process = case.true_processes[edge.edge_id]
                rows.append(
                    {
                        "seed": int(case.seed),
                        "edge_id": edge.edge_id,
                        "tier": tier,
                        "bootstrap": int(bootstrap_index),
                        "predicted_reaction": dominant,
                        "predicted_family": _reaction_family(dominant),
                        "true_process": true_process,
                        "true_family": _truth_family(true_process),
                    }
                )
    return pd.DataFrame(rows)


def _run_case(
    *,
    seed: int,
    split: str,
    output: Path,
    simulator_workspace: Path,
    mf6_executable: Path,
    mp7_executable: Path,
    age_draws: int,
) -> tuple[IndependentAquifer, StrongInferenceResult, list[Dict[str, object]]]:
    case = generate_independent_aquifer(
        int(seed),
        simulator_workspace / f"{split}_{seed}",
        mf6_executable,
        mp7_executable,
    )
    result = run_strong_inference(
        case.observations,
        int(seed),
        run_posterior=False,
        age_draws=int(age_draws),
        age_chains=4,
        age_travel_cost_weight=0.0,
    )
    _write_case_artifacts(output, case, result, split)
    return case, result, _label_rows(case, result, split)


def _contrast_record(
    contrasts: pd.DataFrame,
    *,
    contrast: str,
    metric: str,
    tracer_regime: str | None = None,
) -> Dict[str, object]:
    selected = contrasts[
        (contrasts["contrast"] == contrast) & (contrasts["metric"] == metric)
    ]
    if tracer_regime is not None and "tracer_regime" in selected:
        selected = selected[selected["tracer_regime"] == tracer_regime]
    if len(selected) != 1:
        return {}
    return selected.iloc[0].to_dict()


def run_benchmark(
    *,
    output: Path,
    simulator_workspace: Path,
    mf6_executable: Path,
    mp7_executable: Path,
    dev_seeds: Sequence[int] = DEV_SEEDS,
    test_seeds: Sequence[int] = TEST_SEEDS,
    age_draws: int = 600,
    age_particles: int = 50_000,
    reaction_bootstrap: int = 64,
    paired_bootstrap: int = 10_000,
    overwrite: bool = False,
) -> Dict[str, object]:
    """Execute the Git-locked M7.3 analysis."""

    output = Path(output)
    if output.exists():
        if not overwrite:
            raise FileExistsError(
                f"Output exists: {output}. Pass overwrite=True explicitly."
            )
        shutil.rmtree(output)
    output.mkdir(parents=True)
    simulator_workspace.mkdir(parents=True, exist_ok=True)
    protocol_commit = _git_head()

    development_cases: Dict[int, IndependentAquifer] = {}
    development_results: Dict[int, StrongInferenceResult] = {}
    development_rows = []
    for seed in dev_seeds:
        case, result, rows = _run_case(
            seed=int(seed),
            split="development",
            output=output,
            simulator_workspace=simulator_workspace,
            mf6_executable=mf6_executable,
            mp7_executable=mp7_executable,
            age_draws=age_draws,
        )
        development_cases[int(seed)] = case
        development_results[int(seed)] = result
        development_rows.extend(rows)
    development = pd.DataFrame(development_rows)
    models = fit_evidence_models(development)

    test_cases: Dict[int, IndependentAquifer] = {}
    test_results: Dict[int, StrongInferenceResult] = {}
    test_rows = []
    age_sensitivity_rows = []
    reaction_bootstrap_rows = []
    for seed in test_seeds:
        case, result, rows = _run_case(
            seed=int(seed),
            split="locked_test",
            output=output,
            simulator_workspace=simulator_workspace,
            mf6_executable=mf6_executable,
            mp7_executable=mp7_executable,
            age_draws=age_draws,
        )
        test_cases[int(seed)] = case
        test_results[int(seed)] = result
        test_rows.extend(rows)
        for tracer_regime in ("informative", "tritium_only"):
            age_sensitivity_rows.append(
                topology_age_sensitivity(
                    observations=case.observations,
                    true_ages=case.true_ages_years,
                    true_edges=case.true_edges,
                    seed=int(seed),
                    regime=tracer_regime,
                    n_particles=int(age_particles),
                    order_scale_years=5.0,
                )
            )
        reaction_bootstrap_rows.append(
            _reaction_bootstrap_case(
                case,
                n_bootstrap=int(reaction_bootstrap),
            )
        )

    test = pd.DataFrame(test_rows)
    evidence_summary, evidence_case, conflict_summary = evaluate_evidence_conditions(
        test, models
    )
    evidence_contrasts = bootstrap_evidence_contrasts(
        evidence_case,
        n_bootstrap=int(paired_bootstrap),
    )

    topology_age = pd.concat(age_sensitivity_rows, ignore_index=True)
    topology_age_contrasts = bootstrap_topology_age_contrasts(
        topology_age,
        n_bootstrap=int(paired_bootstrap),
    )

    reaction_predictions = pd.concat(reaction_bootstrap_rows, ignore_index=True)
    (
        reaction_support,
        reaction_edge_summary,
        reaction_summary,
    ) = reaction_support_summary(
        reaction_predictions,
    )

    ghana_audit = audit_ghana_workbook(GHANA_WORKBOOK)

    development.to_csv(output / "development_edge_features.csv", index=False)
    test.to_csv(output / "locked_test_edge_features.csv", index=False)
    evidence_summary.to_csv(output / "evidence_panel_summary.csv", index=False)
    evidence_case.to_csv(output / "evidence_case_metrics.csv", index=False)
    conflict_summary.to_csv(output / "evidence_conflict_summary.csv", index=False)
    evidence_contrasts.to_csv(
        output / "evidence_case_bootstrap_contrasts.csv",
        index=False,
    )
    topology_age.to_csv(output / "topology_age_sensitivity.csv", index=False)
    topology_age_contrasts.to_csv(
        output / "topology_age_bootstrap_contrasts.csv",
        index=False,
    )
    reaction_predictions.to_csv(
        output / "reaction_bootstrap_predictions.csv",
        index=False,
    )
    reaction_support.to_csv(
        output / "reaction_family_support.csv",
        index=False,
    )
    reaction_edge_summary.to_csv(
        output / "reaction_edge_nonuniqueness.csv",
        index=False,
    )
    reaction_summary.to_csv(
        output / "reaction_nonuniqueness_summary.csv",
        index=False,
    )
    (output / "frozen_evidence_models.json").write_text(
        json.dumps(models, indent=2, default=_json_default),
        encoding="utf-8",
    )
    (output / "ghana_data_scope_audit.json").write_text(
        json.dumps(ghana_audit, indent=2, default=_json_default),
        encoding="utf-8",
    )

    candidate_recall = float(
        np.mean(
            [
                len(
                    {f"{u}->{v}" for u, v in test_cases[seed].true_edges}
                    & set(test_results[seed].candidate_edges)
                )
                / len(test_cases[seed].true_edges)
                for seed in test_seeds
            ]
        )
    )
    age_increment_pr = _contrast_record(
        evidence_contrasts,
        contrast="native_incremental_age",
        metric="pr_auc",
    )
    age_increment_log_loss = _contrast_record(
        evidence_contrasts,
        contrast="native_incremental_age",
        metric="log_loss",
    )
    correct_ambiguous_mae = _contrast_record(
        topology_age_contrasts,
        contrast="correct_minus_no_topology",
        metric="age_mae_years",
        tracer_regime="tritium_only",
    )
    reversed_ambiguous_mae = _contrast_record(
        topology_age_contrasts,
        contrast="reversed_minus_correct",
        metric="age_mae_years",
        tracer_regime="tritium_only",
    )
    reaction_index = reaction_summary.set_index(["tier", "true_process"])
    manifest: Dict[str, object] = {
        "benchmark": "M7.3 conditional integration and non-uniqueness",
        "protocol_commit_before_5301_series": protocol_commit,
        "development_seeds": list(map(int, dev_seeds)),
        "locked_test_seeds": list(map(int, test_seeds)),
        "n_development_cases": len(dev_seeds),
        "n_locked_test_cases": len(test_seeds),
        "age_draws_per_chain": int(age_draws),
        "age_importance_particles": int(age_particles),
        "reaction_bootstrap_per_case": int(reaction_bootstrap),
        "paired_case_bootstrap": int(paired_bootstrap),
        "candidate_recall": candidate_recall,
        "all_local_age_cases_converged": bool(
            all(
                bool(result.age_diagnostics.get("converged"))
                for result in test_results.values()
            )
        ),
        "all_topology_age_importance_runs_stable": bool(
            topology_age["importance_stable_ess_ge_400"].all()
        ),
        "native_incremental_age_pr_auc": age_increment_pr,
        "native_incremental_age_log_loss": age_increment_log_loss,
        "tritium_only_correct_topology_mae_contrast": correct_ambiguous_mae,
        "tritium_only_reversed_topology_mae_contrast": reversed_ambiguous_mae,
        "reaction_overall": {
            tier: reaction_index.loc[(tier, "ALL")].to_dict()
            for tier in ("core", "enhanced")
        },
        "carbonate_recovery": {
            tier: {
                process: reaction_index.loc[(tier, process)].to_dict()
                for process in (
                    "carbonate_weathering",
                    "carbonate_precipitation",
                )
            }
            for tier in ("core", "enhanced")
        },
        "ghana_scope": {
            "environmental_age_tracer_panel_available": ghana_audit[
                "environmental_age_tracer_panel_available"
            ],
            "screen_intervals_available": ghana_audit["screen_intervals_available"],
            "time_varying_head_series_available": ghana_audit[
                "time_varying_head_series_available"
            ],
            "single_occasion_head_proxy_possible": ghana_audit[
                "single_occasion_head_proxy_possible"
            ],
            "coordinates_masked": ghana_audit["coordinates_masked"],
            "independent_field_connectivity_truth_available": ghana_audit[
                "independent_field_connectivity_truth_available"
            ],
        },
        "claim_guardrail": (
            "Integration is evaluated as conditional uncertainty reduction. "
            "Synthetic truth is model-conditioned. Ghana evidence supports "
            "component diagnostics and non-identifiability mapping, not field "
            "age, exact topology, or reaction truth."
        ),
    }
    (output / "manifest.json").write_text(
        json.dumps(manifest, indent=2, default=_json_default),
        encoding="utf-8",
    )
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--simulator-workspace",
        type=Path,
        default=DEFAULT_SIMULATOR_WORKSPACE,
    )
    parser.add_argument("--bin-dir", type=Path, default=DEFAULT_BIN_DIR)
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    dev_seeds = SMOKE_DEV_SEEDS if args.quick else DEV_SEEDS
    test_seeds = SMOKE_TEST_SEEDS if args.quick else TEST_SEEDS
    output = (
        args.output
        if not args.quick or args.output != DEFAULT_OUTPUT
        else DEFAULT_OUTPUT.parent / "m7_3_quick"
    )
    manifest = run_benchmark(
        output=output,
        simulator_workspace=args.simulator_workspace,
        mf6_executable=args.bin_dir / "mf6.exe",
        mp7_executable=args.bin_dir / "mp7.exe",
        dev_seeds=dev_seeds,
        test_seeds=test_seeds,
        age_draws=500 if args.quick else 600,
        age_particles=8_000 if args.quick else 50_000,
        reaction_bootstrap=4 if args.quick else 64,
        paired_bootstrap=500 if args.quick else 10_000,
        overwrite=args.overwrite,
    )
    print(json.dumps(manifest, indent=2, default=_json_default))


if __name__ == "__main__":
    main()
