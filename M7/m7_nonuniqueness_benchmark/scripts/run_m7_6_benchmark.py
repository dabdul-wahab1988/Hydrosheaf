"""Run the locked M7.6 auxiliary diagnostic for the unresolved M3 mechanism.

The runner has an explicit truth-blind boundary.  Development cases are fit
and scored first.  For locked test cases, observations and hidden truth are
written separately; all predictions are completed before the truth files are
read for scoring and mechanism stratification.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
import subprocess
import sys
from typing import Mapping, Sequence

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from independent_modflow_generator_v2 import (  # noqa: E402
    IndependentAquiferV2,
    generate_independent_aquifer_v2,
)
from m7_6_analysis import (  # noqa: E402
    REDOX_PROCESSES,
    TARGETS,
    build_edge_features,
    build_node_features,
    bootstrap_case_contrasts,
    evaluate_case_panels,
    fit_models,
    model_summary,
    predict_case_panels,
    score_case_predictions,
    summarise_m3_mechanism,
    ttd_feasibility_diagnostics,
)

DEV_SEEDS = tuple(range(5401, 5407))
TEST_SEEDS = tuple(range(5501, 5513))
SMOKE_DEV_SEEDS = (5951, 5952)
SMOKE_TEST_SEEDS = (5961, 5962, 5963)
NUISANCE_SCALES = {"none": 0.0, "mild": 0.5, "severe": 1.0}
NUISANCE_LEVELS = tuple(NUISANCE_SCALES)
DEFAULT_OUTPUT = REPO_ROOT / "M7" / "m7_nonuniqueness_benchmark" / "results" / "RUN-M7-6-M3-MECHANISM-20260731-01"
DEFAULT_SIMULATOR_WORKSPACE = REPO_ROOT / ".codex_work" / "m7_6_simulators"
DEFAULT_BIN_DIR = REPO_ROOT / ".codex_work" / "modflow-bin"


@dataclass(frozen=True)
class BlindCase:
    seed: int
    nuisance_level: str
    observations: tuple[dict[str, object], ...]
    observation_path: Path
    truth_path: Path
    provenance: Mapping[str, object]


def _git_head() -> str:
    completed = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=REPO_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    return completed.stdout.strip()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _json_default(value: object) -> object:
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, (np.bool_,)):
        return bool(value)
    if isinstance(value, Path):
        return str(value)
    raise TypeError(f"Cannot JSON serialise {type(value)!r}")


def _write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, default=_json_default), encoding="utf-8")


def _case_truth_payload(case: IndependentAquiferV2) -> dict[str, object]:
    return {
        "seed": int(case.seed),
        "true_edges": [list(edge) for edge in case.true_edges],
        "true_ages_years": dict(case.true_ages_years),
        "true_processes": dict(case.true_processes),
        "pathline_rows": list(case.pathline_rows),
        "true_nuisance": dict(case.true_nuisance),
        "true_recharge_isotopes": dict(case.true_recharge_isotopes),
        "true_reaction_isotopes": dict(case.true_reaction_isotopes),
    }


def _write_case_files(case: IndependentAquiferV2, case_dir: Path) -> tuple[Path, Path]:
    case_dir.mkdir(parents=True, exist_ok=True)
    observation_path = case_dir / "observations.json"
    truth_path = case_dir / "truth.json"
    _write_json(observation_path, {"seed": int(case.seed), "observations": list(case.observations), "provenance": case.provenance})
    _write_json(truth_path, _case_truth_payload(case))
    return observation_path, truth_path


def _load_truth(path: Path) -> dict[str, object]:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def _redox_labels(case: IndependentAquiferV2) -> dict[str, str]:
    process_by_node: dict[str, list[str]] = {}
    for edge_id, process in case.true_processes.items():
        _, downstream = str(edge_id).split("->", 1)
        process_by_node.setdefault(downstream, []).append(str(process))
    labels: dict[str, str] = {}
    for row in case.observations:
        site = str(row["site_id"])
        particle = site.rsplit("_M", 1)[0]
        milestone = int(site.rsplit("_M", 1)[1])
        processes = []
        for index in range(1, milestone + 1):
            processes.extend(process_by_node.get(f"{particle}_M{index}", ()))
        if any(process in REDOX_PROCESSES for process in processes):
            labels[site] = "reducing"
        elif "denitrification" in processes:
            labels[site] = "suboxic"
        else:
            labels[site] = "nonreducing"
    return labels


def _mixing_fraction(case: IndependentAquiferV2) -> dict[str, float]:
    """Define the synthetic T3 source fraction from recharge x-position.

    The generator's environmental-isotope field is a deterministic recharge
    gradient.  The x-position fraction is therefore a declared source-mixture
    target, not a claim about a real aquifer endmember fraction.
    """

    recharge_x = {
        particle: float(values["true_recharge_x_m"])
        for particle, values in case.true_recharge_isotopes.items()
    }
    maximum = max(recharge_x.values(), default=1.0)
    output: dict[str, float] = {}
    for row in case.observations:
        site = str(row["site_id"])
        output[site] = float(np.clip(recharge_x.get(site.rsplit("_M", 1)[0], 0.0) / maximum, 0.0, 1.0))
    return output


def _development_case_frames(case: IndependentAquiferV2) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, str]]:
    nodes = build_node_features(
        case.observations,
        true_mixing_fraction=_mixing_fraction(case),
        true_ages=case.true_ages_years,
    )
    edges = build_edge_features(case.observations, true_edges=case.true_edges)
    return nodes, edges, _redox_labels(case)


def _truth_frames(truth: Mapping[str, object], observations: Sequence[Mapping[str, object]]) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, str]]:
    true_edges = [tuple(edge) for edge in truth["true_edges"]]
    true_ages = {str(key): float(value) for key, value in truth["true_ages_years"].items()}
    recharge = truth.get("true_recharge_isotopes", {})
    maximum = max((float(value["true_recharge_x_m"]) for value in recharge.values()), default=1.0)
    fractions = {
        str(row["site_id"]): float(np.clip(float(recharge[str(row["site_id"]).rsplit("_M", 1)[0]]["true_recharge_x_m"]) / maximum, 0.0, 1.0))
        for row in observations
    }
    nodes = build_node_features(observations, true_mixing_fraction=fractions, true_ages=true_ages)
    edges = build_edge_features(observations, true_edges=true_edges)
    # Reconstruct the redox labels from truth only after test predictions.
    process_by_node: dict[str, list[str]] = {}
    for edge_id, process in truth["true_processes"].items():
        _, downstream = str(edge_id).split("->", 1)
        process_by_node.setdefault(downstream, []).append(str(process))
    labels: dict[str, str] = {}
    for row in observations:
        site = str(row["site_id"])
        particle = site.rsplit("_M", 1)[0]
        milestone = int(site.rsplit("_M", 1)[1])
        processes = []
        for index in range(1, milestone + 1):
            processes.extend(process_by_node.get(f"{particle}_M{index}", ()))
        labels[site] = "reducing" if any(process in REDOX_PROCESSES for process in processes) else "suboxic" if "denitrification" in processes else "nonreducing"
    return nodes, edges, labels


def _generate_development(
    *,
    seeds: Sequence[int],
    output: Path,
    simulator_workspace: Path,
    mf6: Path,
    mp7: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    node_frames = []
    edge_frames = []
    diagnostic_nodes = []
    diagnostic_pairs = []
    for level, scale in NUISANCE_SCALES.items():
        for seed in seeds:
            case = generate_independent_aquifer_v2(
                int(seed),
                simulator_workspace / f"development_{level}_{seed}",
                mf6,
                mp7,
                shared_nuisance=level != "none",
                nuisance_scale=float(scale),
            )
            case_dir = output / "development" / level / str(seed)
            _write_case_files(case, case_dir)
            nodes, edges, redox = _development_case_frames(case)
            nodes["nuisance_level"] = level
            edges["nuisance_level"] = level
            node_frames.append(nodes)
            edge_frames.append(edges)
            node_diag, pair_diag = ttd_feasibility_diagnostics(case.observations, redox_by_node=redox)
            node_diag["seed"] = int(seed)
            node_diag["nuisance_level"] = level
            pair_diag["seed"] = int(seed)
            pair_diag["nuisance_level"] = level
            diagnostic_nodes.append(node_diag)
            diagnostic_pairs.append(pair_diag)
    return (
        pd.concat(node_frames, ignore_index=True),
        pd.concat(edge_frames, ignore_index=True),
        pd.concat(diagnostic_nodes, ignore_index=True),
        pd.concat(diagnostic_pairs, ignore_index=True),
    )


def _incremental_comparisons() -> list[tuple[str, str, str]]:
    comparisons = []
    for panel in sorted(SUBSETS_FOR_COMPARISON):
        for stream in panel:
            baseline = "".join(char for char in panel if char != stream)
            if baseline:
                comparisons.append((panel, baseline, f"{stream}_added_to_{baseline}"))
    return comparisons


SUBSETS_FOR_COMPARISON = tuple(
    "".join(stream for index, stream in enumerate(("H", "C", "N", "E")) if mask & (1 << index))
    for mask in range(1, 1 << 4)
)


def _blind_test_predictions(
    *,
    seeds: Sequence[int],
    output: Path,
    simulator_workspace: Path,
    mf6: Path,
    mp7: Path,
    models: Mapping[str, Mapping[str, Mapping[str, object]]],
) -> tuple[pd.DataFrame, pd.DataFrame, dict[tuple[str, int], BlindCase]]:
    edge_predictions = []
    node_predictions = []
    cases: dict[tuple[str, int], BlindCase] = {}
    # Generation and all inference occur before any truth file is opened.
    for level, scale in NUISANCE_SCALES.items():
        for seed in seeds:
            case = generate_independent_aquifer_v2(
                int(seed),
                simulator_workspace / f"locked_test_{level}_{seed}",
                mf6,
                mp7,
                shared_nuisance=level != "none",
                nuisance_scale=float(scale),
            )
            observation_path, truth_path = _write_case_files(case, output / "locked_test" / level / str(seed))
            blind = BlindCase(
                seed=int(seed),
                nuisance_level=level,
                observations=tuple(dict(row) for row in case.observations),
                observation_path=observation_path,
                truth_path=truth_path,
                provenance=case.provenance,
            )
            cases[(level, int(seed))] = blind
            nodes = build_node_features(blind.observations)
            edges = build_edge_features(blind.observations)
            for condition in ("native", "permuted", "identity_relation"):
                edge_pred, node_pred = predict_case_panels(nodes, edges, models, condition=condition)
                edge_pred["nuisance_level"] = level
                node_pred["nuisance_level"] = level
                edge_predictions.append(edge_pred)
                node_predictions.append(node_pred)
    return pd.concat(edge_predictions, ignore_index=True), pd.concat(node_predictions, ignore_index=True), cases


def _score_locked_test(
    edge_predictions: pd.DataFrame,
    node_predictions: pd.DataFrame,
    cases: Mapping[tuple[str, int], BlindCase],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    edge_metrics = []
    node_metrics = []
    diagnostic_nodes = []
    diagnostic_pairs = []
    # This is the first point where locked-test truth is opened.
    for (level, seed), case in cases.items():
        truth = _load_truth(case.truth_path)
        nodes, edges, redox = _truth_frames(truth, case.observations)
        selected_edge_predictions = edge_predictions[(edge_predictions["nuisance_level"] == level) & (edge_predictions["seed"] == seed)]
        selected_node_predictions = node_predictions[(node_predictions["nuisance_level"] == level) & (node_predictions["seed"] == seed)]
        edge_truth, node_truth = score_case_predictions(nodes, edges, selected_edge_predictions, selected_node_predictions)
        edge_truth["nuisance_level"] = level
        node_truth["nuisance_level"] = level
        edge_metrics.append(edge_truth)
        node_metrics.append(node_truth)
        node_diag, pair_diag = ttd_feasibility_diagnostics(case.observations, redox_by_node=redox)
        node_diag["seed"] = int(seed)
        node_diag["nuisance_level"] = level
        pair_diag["seed"] = int(seed)
        pair_diag["nuisance_level"] = level
        diagnostic_nodes.append(node_diag)
        diagnostic_pairs.append(pair_diag)
    return (
        pd.concat(edge_metrics, ignore_index=True),
        pd.concat(node_metrics, ignore_index=True),
        pd.concat(diagnostic_nodes, ignore_index=True),
        pd.concat(diagnostic_pairs, ignore_index=True),
    )


def _bootstrap_group_difference(
    frame: pd.DataFrame,
    *,
    value: str,
    left_filter: Mapping[str, object],
    right_filter: Mapping[str, object],
    label: str,
    n_bootstrap: int,
    random_seed: int,
) -> dict[str, object]:
    left = frame.copy()
    right = frame.copy()
    for key, expected in left_filter.items():
        left = left[left[key] == expected]
    for key, expected in right_filter.items():
        right = right[right[key] == expected]
    left_case = left.groupby("seed")[value].mean()
    right_case = right.groupby("seed")[value].mean()
    common = sorted(set(left_case.index) & set(right_case.index))
    if not common:
        return {"contrast": label, "n_cases": 0, "mean_difference": float("nan"), "ci95_low": float("nan"), "ci95_high": float("nan")}
    paired = left_case.loc[common].to_numpy(float) - right_case.loc[common].to_numpy(float)
    rng = np.random.default_rng(int(random_seed))
    sampled = rng.choice(paired, size=(int(n_bootstrap), len(paired)), replace=True).mean(axis=1)
    return {
        "contrast": label,
        "mean_difference": float(np.mean(paired)),
        "ci95_low": float(np.quantile(sampled, 0.025)),
        "ci95_high": float(np.quantile(sampled, 0.975)),
        "n_cases": int(len(common)),
        "n_bootstrap": int(n_bootstrap),
        "resampling_unit": "independent_MODFLOW_case",
    }


def _m3_claim_decision(
    *,
    ttd_summary: pd.DataFrame,
    primary_contrasts: pd.DataFrame,
    node_metrics: pd.DataFrame,
    node_contrasts: pd.DataFrame,
) -> dict[str, object]:
    def row(label: str) -> dict[str, object]:
        selected = primary_contrasts[primary_contrasts["contrast"] == label]
        return selected.iloc[0].to_dict() if len(selected) == 1 else {}

    redox_cfc = row("reducing_minus_nonreducing_cfc11_pair_rate_severe")
    redox_cfc12 = row("reducing_minus_nonreducing_cfc12_pair_rate_severe")
    severe_full = row("severe_minus_none_full_infeasibility")
    cfc_supported = bool(
        redox_cfc
        and float(redox_cfc.get("ci95_low", 0.0)) > 0.0
        and (not redox_cfc12 or float(redox_cfc12.get("ci95_high", 0.0)) <= float(redox_cfc.get("ci95_low", 0.0)))
    )
    isotope = node_contrasts[
        (node_contrasts["target"] == "T2")
        & (node_contrasts["condition"] == "native")
        & (node_contrasts["nuisance_level"] == "none")
        & (node_contrasts["contrast"] == "E_added_to_N")
        & (node_contrasts["metric"] == "mae")
    ]
    e_gate = bool(
        len(isotope) == 1
        and abs(float(isotope.iloc[0]["mean_difference"])) <= 1.0e-12
        and abs(float(isotope.iloc[0]["ci95_low"])) <= 1.0e-12
        and abs(float(isotope.iloc[0]["ci95_high"])) <= 1.0e-12
    )
    return {
        "decision_type": "auxiliary_controlled_synthetic_M3_mechanism_diagnostic",
        "execution_status": "PASS" if e_gate else "FAIL",
        "isotope_age_control_passed": e_gate,
        "isotope_age_control": isotope.iloc[0].to_dict() if len(isotope) == 1 else {},
        "redox_stratified_cfc11_mechanism_supported": cfc_supported,
        "shared_nuisance_full_infeasibility_contrast": severe_full,
        "interpretation": (
            "The result tests whether the declared synthetic nuisance and CFC-degradation mechanisms can produce an M3-like signature. "
            "It is not field validation and does not resolve the cause of the USGS infeasibility."
        ),
        "claim_if_supported": (
            "Under the tested controlled-synthetic generator, reducing-path CFC degradation can produce a family-specific infeasibility signature."
        ),
        "claim_prohibited": [
            "the mechanism caused the USGS or Ghana observations",
            "field age recovery was validated",
            "the M7 framework is generally superior",
        ],
        "ttd_summary_rows": int(len(ttd_summary)),
    }


def run_benchmark(
    *,
    output: Path = DEFAULT_OUTPUT,
    simulator_workspace: Path = DEFAULT_SIMULATOR_WORKSPACE,
    mf6_executable: Path = DEFAULT_BIN_DIR / "mf6.exe",
    mp7_executable: Path = DEFAULT_BIN_DIR / "mp7.exe",
    dev_seeds: Sequence[int] = DEV_SEEDS,
    test_seeds: Sequence[int] = TEST_SEEDS,
    paired_bootstrap: int = 10_000,
    overwrite: bool = False,
) -> dict[str, object]:
    output = Path(output)
    if output.exists():
        if not overwrite:
            raise FileExistsError(f"Output exists: {output}. Pass --overwrite explicitly.")
        # Only the run-specific output directory is removed.
        import shutil

        shutil.rmtree(output)
    output.mkdir(parents=True, exist_ok=True)
    mf6_executable = Path(mf6_executable).resolve()
    mp7_executable = Path(mp7_executable).resolve()
    if not mf6_executable.exists() or not mp7_executable.exists():
        raise FileNotFoundError("Both official mf6.exe and mp7.exe are required.")

    protocol_commit = _git_head()
    dev_nodes, dev_edges, dev_ttd_nodes, dev_ttd_pairs = _generate_development(
        seeds=dev_seeds,
        output=output,
        simulator_workspace=Path(simulator_workspace),
        mf6=mf6_executable,
        mp7=mp7_executable,
    )
    models = fit_models(dev_nodes, dev_edges)
    dev_edge_metrics = []
    dev_node_metrics = []
    for condition in ("native", "permuted", "identity_relation"):
        for seed in sorted(dev_edges["seed"].unique()):
            for level in NUISANCE_LEVELS:
                nodes = dev_nodes[(dev_nodes["seed"] == seed) & (dev_nodes["nuisance_level"] == level)]
                edges = dev_edges[(dev_edges["seed"] == seed) & (dev_edges["nuisance_level"] == level)]
                edge_metrics, node_metrics = evaluate_case_panels(nodes, edges, models, condition=condition)
                edge_metrics["nuisance_level"] = level
                node_metrics["nuisance_level"] = level
                dev_edge_metrics.append(edge_metrics)
                dev_node_metrics.append(node_metrics)
    dev_edge_metrics_frame = pd.concat(dev_edge_metrics, ignore_index=True)
    dev_node_metrics_frame = pd.concat(dev_node_metrics, ignore_index=True)

    # The locked-test truth files are not opened by this call or by any model
    # fitting code below it.
    edge_predictions, node_predictions, blind_cases = _blind_test_predictions(
        seeds=test_seeds,
        output=output,
        simulator_workspace=Path(simulator_workspace),
        mf6=mf6_executable,
        mp7=mp7_executable,
        models=models,
    )
    edge_metrics, node_metrics, test_ttd_nodes, test_ttd_pairs = _score_locked_test(edge_predictions, node_predictions, blind_cases)

    comparisons = _incremental_comparisons()
    edge_contrasts = bootstrap_case_contrasts(
        edge_metrics,
        metric_columns=("pr_auc", "brier", "log_loss", "mean_edge_entropy", "overconfident_error_fraction"),
        comparisons=comparisons,
        n_bootstrap=int(paired_bootstrap),
        random_seed=7608,
    )
    node_contrasts = bootstrap_case_contrasts(
        node_metrics,
        metric_columns=("mae", "mean_interval_width", "coverage"),
        comparisons=comparisons,
        n_bootstrap=int(paired_bootstrap),
        random_seed=7609,
    )
    ttd_summary = summarise_m3_mechanism(test_ttd_nodes, test_ttd_pairs, nuisance_level="all")

    # Primary M3 mechanism contrasts are deliberately a small, explicit family.
    ttd_case = (
        test_ttd_nodes.groupby(["seed", "nuisance_level", "redox_class"], as_index=False)
        .agg(
            full_infeasibility_rate=("full_panel_infeasible", "mean"),
            cfc_pair_rate=("cfc_cfc_infeasible_pair_count", lambda values: float(np.mean(np.asarray(values) > 0))),
            tritium_pair_rate=("tritium_pair_infeasible", "mean"),
        )
    )
    ttd_pairs_case = (
        test_ttd_pairs[test_ttd_pairs["pair_family"].isin(["cfc_cfc", "tritium_pair"])]
        .groupby(["seed", "nuisance_level", "redox_class", "pair_family", "cfc11_in_pair", "cfc12_in_pair"], as_index=False)
        .agg(pair_rate=("pair_infeasible", "mean"))
    )
    primary_rows = [
        _bootstrap_group_difference(ttd_case, value="full_infeasibility_rate", left_filter={"nuisance_level": "severe"}, right_filter={"nuisance_level": "none"}, label="severe_minus_none_full_infeasibility", n_bootstrap=int(paired_bootstrap), random_seed=7610),
        _bootstrap_group_difference(ttd_case, value="cfc_pair_rate", left_filter={"nuisance_level": "severe", "redox_class": "reducing"}, right_filter={"nuisance_level": "severe", "redox_class": "nonreducing"}, label="reducing_minus_nonreducing_cfc_pair_rate_severe", n_bootstrap=int(paired_bootstrap), random_seed=7611),
        _bootstrap_group_difference(ttd_case, value="tritium_pair_rate", left_filter={"nuisance_level": "severe", "redox_class": "reducing"}, right_filter={"nuisance_level": "severe", "redox_class": "nonreducing"}, label="reducing_minus_nonreducing_tritium_pair_rate_severe", n_bootstrap=int(paired_bootstrap), random_seed=7612),
    ]
    cfc11_case = ttd_pairs_case[(ttd_pairs_case["pair_family"] == "cfc_cfc") & (ttd_pairs_case["cfc11_in_pair"] == 1)].copy()
    cfc12_case = ttd_pairs_case[(ttd_pairs_case["pair_family"] == "cfc_cfc") & (ttd_pairs_case["cfc12_in_pair"] == 1)].copy()
    primary_rows.append(
        _bootstrap_group_difference(cfc11_case, value="pair_rate", left_filter={"nuisance_level": "severe", "redox_class": "reducing"}, right_filter={"nuisance_level": "severe", "redox_class": "nonreducing"}, label="reducing_minus_nonreducing_cfc11_pair_rate_severe", n_bootstrap=int(paired_bootstrap), random_seed=7613)
    )
    primary_rows.append(
        _bootstrap_group_difference(cfc12_case, value="pair_rate", left_filter={"nuisance_level": "severe", "redox_class": "reducing"}, right_filter={"nuisance_level": "severe", "redox_class": "nonreducing"}, label="reducing_minus_nonreducing_cfc12_pair_rate_severe", n_bootstrap=int(paired_bootstrap), random_seed=7614)
    )
    primary_contrasts = pd.DataFrame(primary_rows)
    claim_decision = _m3_claim_decision(ttd_summary=ttd_summary, primary_contrasts=primary_contrasts, node_metrics=node_metrics, node_contrasts=node_contrasts)

    dev_nodes.to_csv(output / "development_node_features.csv", index=False)
    dev_edges.to_csv(output / "development_edge_features.csv", index=False)
    dev_edge_metrics_frame.to_csv(output / "development_edge_metrics.csv", index=False)
    dev_node_metrics_frame.to_csv(output / "development_node_metrics.csv", index=False)
    dev_ttd_nodes.to_csv(output / "development_ttd_node_diagnostics.csv", index=False)
    dev_ttd_pairs.to_csv(output / "development_ttd_pair_diagnostics.csv", index=False)
    edge_predictions.to_csv(output / "locked_test_edge_predictions_blind.csv", index=False)
    node_predictions.to_csv(output / "locked_test_node_predictions_blind.csv", index=False)
    edge_metrics.to_csv(output / "locked_test_edge_metrics.csv", index=False)
    node_metrics.to_csv(output / "locked_test_node_metrics.csv", index=False)
    edge_contrasts.to_csv(output / "locked_test_edge_bootstrap_contrasts.csv", index=False)
    node_contrasts.to_csv(output / "locked_test_node_bootstrap_contrasts.csv", index=False)
    test_ttd_nodes.to_csv(output / "locked_test_ttd_node_diagnostics.csv", index=False)
    test_ttd_pairs.to_csv(output / "locked_test_ttd_pair_diagnostics.csv", index=False)
    ttd_summary.to_csv(output / "m3_mechanism_ttd_summary.csv", index=False)
    primary_contrasts.to_csv(output / "m3_mechanism_primary_contrasts.csv", index=False)
    _write_json(output / "frozen_fusion_models.json", model_summary(models))
    _write_json(output / "m7_6_claim_decision.json", claim_decision)

    manifest = {
        "benchmark": "M7.6 auxiliary controlled-synthetic M3 mechanism diagnostic",
        "protocol_commit_before_locked_test": protocol_commit,
        "runner_commit": _git_head(),
        "development_seeds": list(map(int, dev_seeds)),
        "locked_test_seeds": list(map(int, test_seeds)),
        "nuisance_scales": dict(NUISANCE_SCALES),
        "paired_case_bootstrap": int(paired_bootstrap),
        "mf6": {"path": str(mf6_executable), "sha256": _sha256(mf6_executable)},
        "mp7": {"path": str(mp7_executable), "sha256": _sha256(mp7_executable)},
        "truth_blind_contract": "all locked-test predictions were written before truth.json files were opened",
        "mixing_target_definition": "recharge x-position divided by the maximum recharge x-position in the synthetic case",
        "claim_decision": claim_decision,
        "claim_guardrail": "controlled-synthetic mechanism diagnostic only; no field validation and no resolution of the USGS cause",
    }
    _write_json(output / "manifest.json", manifest)
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--simulator-workspace", type=Path, default=DEFAULT_SIMULATOR_WORKSPACE)
    parser.add_argument("--bin-dir", type=Path, default=DEFAULT_BIN_DIR)
    parser.add_argument("--quick", action="store_true", help="non-claim-bearing smoke run")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--paired-bootstrap", type=int, default=10_000)
    args = parser.parse_args()
    dev_seeds = SMOKE_DEV_SEEDS if args.quick else DEV_SEEDS
    test_seeds = SMOKE_TEST_SEEDS if args.quick else TEST_SEEDS
    output = args.output
    if args.quick and output == DEFAULT_OUTPUT:
        output = DEFAULT_OUTPUT.parent / "RUN-M7-6-M3-MECHANISM-SMOKE-20260731-01"
    manifest = run_benchmark(
        output=output,
        simulator_workspace=args.simulator_workspace,
        mf6_executable=args.bin_dir / "mf6.exe",
        mp7_executable=args.bin_dir / "mp7.exe",
        dev_seeds=dev_seeds,
        test_seeds=test_seeds,
        paired_bootstrap=500 if args.quick else args.paired_bootstrap,
        overwrite=args.overwrite,
    )
    print(json.dumps(manifest, indent=2, default=_json_default))


if __name__ == "__main__":
    main()
