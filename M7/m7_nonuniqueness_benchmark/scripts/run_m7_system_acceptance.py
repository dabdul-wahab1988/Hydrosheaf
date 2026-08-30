"""Strict public-pipeline acceptance on independent M7 aquifer cases."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import sys
from copy import deepcopy
from datetime import datetime, timezone
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np
import pandas as pd
from sklearn.metrics import (
    average_precision_score,
    brier_score_loss,
    f1_score,
    log_loss,
)

PROJECT = Path(__file__).resolve().parents[1]
REPO = Path(__file__).resolve().parents[3]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from hydrosheaf.api import fit_network_pipeline  # noqa: E402
from hydrosheaf.inference.network_fit import infer_edges  # noqa: E402
from hydrosheaf.nuclear.input_history import InputHistory  # noqa: E402

from independent_modflow_generator import (  # noqa: E402
    IndependentAquifer,
    generate_independent_aquifer,
)
from strong_inference import strong_config  # noqa: E402

RUN_ID = "RUN-M7-SYSTEM-20260728-01"
LOCKED_SEEDS = (6301, 6302, 6303, 6304, 6305, 6306)
QUICK_SEEDS = (6391, 6392)
TRUTH_FIELDS = {
    "true_age",
    "true_edges",
    "true_process",
    "true_age_hidden_from_inference",
}
CONDITIONS = ("hydraulic_only", "local_age", "full_sheaf", "age_permuted")


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


def _verify_lock(expected_seeds: Sequence[int], bootstrap_samples: int) -> dict:
    lock_path = PROJECT / "m7_system_acceptance_protocol.lock.json"
    lock = json.loads(lock_path.read_text(encoding="utf-8"))
    if lock.get("run_id") != RUN_ID:
        raise RuntimeError("System-acceptance lock has the wrong run ID.")
    if list(map(int, expected_seeds)) != lock.get("locked_seeds"):
        raise RuntimeError("Requested seeds do not match the locked protocol.")
    if int(bootstrap_samples) != int(lock.get("bootstrap_samples", -1)):
        raise RuntimeError("Bootstrap count does not match the locked protocol.")
    for relative, expected in lock.get("sha256", {}).items():
        path = (REPO / relative).resolve()
        if not path.is_relative_to(REPO.resolve()):
            raise RuntimeError(f"Locked path escapes repository: {relative}")
        actual = _sha256(path)
        if actual != expected:
            raise RuntimeError(f"Locked file changed: {relative}: {actual} != {expected}")
    return lock


def _prepare_rows(case: IndependentAquifer, *, permute_age: bool) -> list[dict]:
    rows = [dict(row) for row in case.observations]
    leaked = TRUTH_FIELDS & set().union(*(set(row) for row in rows))
    if leaked:
        raise ValueError(f"Truth fields leaked into inference rows: {sorted(leaked)}")
    tritium = np.asarray([float(row["tritium_TU"]) for row in rows])
    if permute_age:
        tritium = np.random.default_rng(case.seed + 880_301).permutation(tritium)
    for row, value in zip(rows, tritium):
        row["3H"] = float(value)
        row["sample_date"] = "2025-07-01"
    return rows


def _condition_config(condition: str):
    config = strong_config(phreeqc_enabled=False)
    config.residence_time_tracer = "3H"
    config.edge_max_neighbors = 3
    config.sheaf_soft_beta = 1.0
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


def _probabilities(candidates, retained, condition: str) -> np.ndarray:
    retained_by_id = {edge.edge_id: edge for edge in retained}
    values: list[float] = []
    for candidate in candidates:
        if condition == "hydraulic_only":
            attrs = candidate.attrs or {}
            value = attrs.get(
                "prior_edge_probability",
                attrs.get("p_uv", attrs.get("edge_confidence", 0.5)),
            )
        else:
            edge = retained_by_id.get(candidate.edge_id)
            value = 1.0e-6 if edge is None else (edge.attrs or {}).get(
                "sheaf_weight", 1.0e-6
            )
        values.append(float(np.clip(value, 1.0e-6, 1.0 - 1.0e-6)))
    return np.asarray(values, dtype=float)


def _metric_row(
    case: IndependentAquifer,
    condition: str,
    candidates,
    retained,
) -> dict:
    truth_ids = {f"{u}->{v}" for u, v in case.true_edges}
    candidate_ids = [edge.edge_id for edge in candidates]
    truth = np.asarray([edge_id in truth_ids for edge_id in candidate_ids], dtype=int)
    probability = _probabilities(candidates, retained, condition)
    selected = np.asarray(
        [edge_id in {edge.edge_id for edge in retained} for edge_id in candidate_ids],
        dtype=int,
    )
    return {
        "seed": int(case.seed),
        "condition": condition,
        "n_candidates": len(candidates),
        "n_true_edges": len(truth_ids),
        "n_selected": int(selected.sum()),
        "candidate_recall": float(len(truth_ids & set(candidate_ids)) / len(truth_ids)),
        "pr_auc": float(average_precision_score(truth, probability)),
        "brier": float(brier_score_loss(truth, probability)),
        "log_loss": float(log_loss(truth, probability, labels=[0, 1])),
        "selected_f1": float(f1_score(truth, selected, zero_division=0)),
    }


def _bootstrap_contrast(
    metrics: pd.DataFrame,
    left: str,
    right: str,
    metric: str,
    *,
    samples: int,
    seed: int,
) -> dict:
    pivot = metrics.pivot(index="seed", columns="condition", values=metric)
    differences = (pivot[left] - pivot[right]).to_numpy(float)
    rng = np.random.default_rng(seed)
    boot = np.empty(int(samples), dtype=float)
    for index in range(int(samples)):
        boot[index] = float(
            np.mean(rng.choice(differences, size=len(differences), replace=True))
        )
    return {
        "left": left,
        "right": right,
        "metric": metric,
        "n_cases": len(differences),
        "mean_difference": float(np.mean(differences)),
        "ci95_low": float(np.quantile(boot, 0.025)),
        "ci95_high": float(np.quantile(boot, 0.975)),
    }


def _safe_output(output: Path, overwrite: bool) -> None:
    resolved = output.resolve()
    project_root = PROJECT.resolve()
    if not resolved.is_relative_to(project_root):
        raise ValueError("Output must stay within the M7 benchmark directory.")
    if resolved.exists():
        if not overwrite:
            raise FileExistsError(f"Output exists: {resolved}")
        if resolved == project_root:
            raise ValueError("Refusing to remove the benchmark root.")
        shutil.rmtree(resolved)
    resolved.mkdir(parents=True)


def _stable_frame(frame: pd.DataFrame) -> pd.DataFrame:
    """Round reported floats so cross-process artifact hashes are stable."""
    stable = frame.copy()
    float_columns = stable.select_dtypes(include=["float32", "float64"]).columns
    stable[float_columns] = stable[float_columns].round(12)
    return stable


def run_benchmark(
    *,
    output: Path,
    simulator_workspace: Path,
    mf6_executable: Path,
    mp7_executable: Path,
    seeds: Sequence[int] = LOCKED_SEEDS,
    bootstrap_samples: int = 5000,
    overwrite: bool = False,
) -> dict:
    if tuple(map(int, seeds)) == LOCKED_SEEDS:
        _verify_lock(seeds, bootstrap_samples)
    _safe_output(output, overwrite)
    simulator_workspace.mkdir(parents=True, exist_ok=True)
    metric_rows: list[dict] = []
    prediction_rows: list[dict] = []
    stage_records: dict[str, dict] = {}
    generator_records: dict[str, dict] = {}

    input_history = InputHistory(
        np.asarray([1850.0, 2035.0]),
        np.asarray([6.2, 6.2]),
    )
    for seed in seeds:
        case = generate_independent_aquifer(
            int(seed),
            simulator_workspace / f"case_{int(seed)}",
            mf6_executable,
            mp7_executable,
        )
        generator_records[str(seed)] = dict(case.provenance)
        for condition in CONDITIONS:
            rows = _prepare_rows(case, permute_age=condition == "age_permuted")
            config = _condition_config(condition)
            candidates = infer_edges(rows, method="probabilistic", config=config)
            sheaf_enabled = condition != "hydraulic_only"
            required = ["network_fit"]
            if sheaf_enabled:
                required = ["nuclear_age", "sheaf_refinement", "network_fit"]
            _, extras = fit_network_pipeline(
                rows,
                deepcopy(candidates),
                config,
                auto_disable_missing=False,
                sheaf_refinement_enabled=sheaf_enabled,
                nuclear_inference_options={
                    "input_hist": input_history,
                    "sampler": "grid",
                    "n_samples": 600,
                    "n_chains": 4,
                    "random_seed": int(seed) + 10_000 * CONDITIONS.index(condition),
                    "max_age_years": 200.0,
                    "root_age_prior_median_years": 20.0,
                    "root_age_prior_log_sigma": 1.2,
                },
                strict_stage_completion=True,
                required_stages=required,
            )
            retained = extras["edges"]
            metric_rows.append(
                _metric_row(case, condition, candidates, retained)
            )
            stage_records[f"{seed}:{condition}"] = extras["stage_status"]
            truth_ids = {f"{u}->{v}" for u, v in case.true_edges}
            probabilities = _probabilities(candidates, retained, condition)
            selected_ids = {edge.edge_id for edge in retained}
            for edge, probability in zip(candidates, probabilities):
                prediction_rows.append(
                    {
                        "seed": int(seed),
                        "condition": condition,
                        "edge_id": edge.edge_id,
                        "u": edge.u,
                        "v": edge.v,
                        "is_true_edge": int(edge.edge_id in truth_ids),
                        "probability": float(probability),
                        "selected": int(edge.edge_id in selected_ids),
                    }
                )

    metrics = _stable_frame(pd.DataFrame(metric_rows))
    predictions = _stable_frame(pd.DataFrame(prediction_rows))
    summary = metrics.groupby("condition", sort=True).agg(
        n_cases=("seed", "count"),
        candidate_recall_mean=("candidate_recall", "mean"),
        pr_auc_mean=("pr_auc", "mean"),
        brier_mean=("brier", "mean"),
        log_loss_mean=("log_loss", "mean"),
        selected_f1_mean=("selected_f1", "mean"),
    ).reset_index()
    summary = _stable_frame(summary)

    contrast_specs = (
        ("full_sheaf", "hydraulic_only"),
        ("full_sheaf", "local_age"),
        ("full_sheaf", "age_permuted"),
    )
    contrasts = _stable_frame(pd.DataFrame(
        [
            _bootstrap_contrast(
                metrics,
                left,
                right,
                metric,
                samples=bootstrap_samples,
                seed=2026072800 + 100 * index + metric_index,
            )
            for index, (left, right) in enumerate(contrast_specs)
            for metric_index, metric in enumerate(("pr_auc", "brier", "selected_f1"))
        ]
    ))

    all_completed = all(
        record[stage]["status"] == "completed"
        for record in stage_records.values()
        for stage in ("network_fit",)
    ) and all(
        record[stage]["status"] == "completed"
        for key, record in stage_records.items()
        if not key.endswith(":hydraulic_only")
        for stage in ("nuclear_age", "sheaf_refinement")
    )
    finite_metrics = bool(
        np.isfinite(metrics[["pr_auc", "brier", "log_loss", "selected_f1"]]).all().all()
    )
    independent_generator = all(
        record.get("imports_hydrosheaf") is False
        for record in generator_records.values()
    )
    candidate_recall = float(
        metrics[metrics["condition"] == "hydraulic_only"]["candidate_recall"].mean()
    )
    system_gate = bool(
        all_completed and finite_metrics and independent_generator and candidate_recall >= 0.80
    )

    def contrast(left: str, right: str, metric: str) -> Mapping[str, object]:
        row = contrasts[
            (contrasts["left"] == left)
            & (contrasts["right"] == right)
            & (contrasts["metric"] == metric)
        ]
        return row.iloc[0].to_dict()

    full_h_pr = contrast("full_sheaf", "hydraulic_only", "pr_auc")
    full_h_brier = contrast("full_sheaf", "hydraulic_only", "brier")
    full_perm_pr = contrast("full_sheaf", "age_permuted", "pr_auc")
    full_h_f1 = contrast("full_sheaf", "hydraulic_only", "selected_f1")
    incremental_claim = bool(
        system_gate
        and float(full_h_pr["ci95_low"]) > 0.0
        and float(full_h_brier["ci95_high"]) < 0.0
        and float(full_perm_pr["ci95_low"]) > 0.0
        and float(full_h_f1["mean_difference"]) >= 0.0
    )
    claim_decision = {
        "system_acceptance_passed": system_gate,
        "incremental_full_sheaf_claim_allowed": incremental_claim,
        "global_section_contribution": contrast(
            "full_sheaf", "local_age", "pr_auc"
        ),
        "allowed_claim": (
            "On fresh code-independent controlled-synthetic cases, the strict public "
            "pipeline completed and the full sheaf improved topology discrimination "
            "under the locked native and adverse-control criteria."
            if incremental_claim
            else "The strict public pipeline is evaluated as a controlled-synthetic "
            "system execution and ablation result; no incremental topology-superiority "
            "claim is allowed."
        ),
        "guardrail": (
            "Controlled-synthetic system acceptance is not field validation and does "
            "not establish general HydroSheaf superiority."
        ),
    }

    metrics.to_csv(output / "case_metrics.csv", index=False, lineterminator="\n")
    predictions.to_csv(output / "edge_predictions.csv", index=False, lineterminator="\n")
    summary.to_csv(output / "condition_summary.csv", index=False, lineterminator="\n")
    contrasts.to_csv(output / "paired_bootstrap_contrasts.csv", index=False, lineterminator="\n")
    (output / "stage_status.json").write_text(
        json.dumps(stage_records, indent=2), encoding="utf-8"
    )
    (output / "generator_provenance.json").write_text(
        json.dumps(generator_records, indent=2), encoding="utf-8"
    )
    (output / "claim_decision.json").write_text(
        json.dumps(claim_decision, indent=2), encoding="utf-8"
    )
    manifest = {
        "run_id": RUN_ID,
        "status": "PASS" if system_gate else "FAIL",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "git_revision": _git_revision(),
        "seeds": list(map(int, seeds)),
        "bootstrap_samples": int(bootstrap_samples),
        "mean_candidate_recall": candidate_recall,
        "all_requested_stages_completed": all_completed,
        "all_metrics_finite": finite_metrics,
        "generator_independent": independent_generator,
        "excluded_capabilities": ["temporal_series", "graph_3d", "vadose"],
        "artifacts": {},
    }
    for path in sorted(output.iterdir()):
        if path.is_file() and path.name != "run_manifest.json":
            manifest["artifacts"][path.name] = _sha256(path)
    (output / "run_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )
    return {"manifest": manifest, "claim_decision": claim_decision}


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=PROJECT / "results" / RUN_ID,
    )
    parser.add_argument(
        "--simulator-workspace",
        type=Path,
        default=REPO / ".codex_work" / "m7-system-simulators",
    )
    parser.add_argument(
        "--bin-dir", type=Path, default=REPO / ".codex_work" / "modflow-bin"
    )
    parser.add_argument("--bootstrap-samples", type=int, default=5000)
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    if args.bootstrap_samples < 100:
        raise ValueError("bootstrap-samples must be at least 100")
    result = run_benchmark(
        output=args.output,
        simulator_workspace=args.simulator_workspace,
        mf6_executable=args.bin_dir / "mf6.exe",
        mp7_executable=args.bin_dir / "mp7.exe",
        seeds=QUICK_SEEDS if args.quick else LOCKED_SEEDS,
        bootstrap_samples=args.bootstrap_samples,
        overwrite=args.overwrite,
    )
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
