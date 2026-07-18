"""Run the optional M7 operational synthetic-twin extension.

The locked M7 integration CSVs are treated as immutable.  This command writes
only beneath ``results/digital_twin`` and ``figures/digital_twin`` plus the
generated results note ``docs/m7_digital_twin_results.md``.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[1]
RESULTS = HERE / "results" / "digital_twin"
FIGURES = HERE / "figures" / "digital_twin"
DOCS = HERE / "docs"

from m7_digital_twin import (  # noqa: E402
    TwinConfig,
    paired_skill,
    run_replicate,
    summarize_metrics,
)


LOCKED_BASELINE_RESULTS = (
    "m7_age_coherence_demo.csv",
    "m7_age_gain.csv",
    "m7_age_recovery.csv",
    "m7_edge_classification.csv",
    "m7_integration_gain.csv",
    "m7_trap_rejection.csv",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def locked_hashes() -> dict[str, str]:
    result_root = HERE / "results"
    return {
        name: sha256(result_root / name)
        for name in LOCKED_BASELINE_RESULTS
        if (result_root / name).exists()
    }


def source_hashes() -> dict[str, str]:
    return {
        path.name: sha256(path)
        for path in (
            Path(__file__).resolve(),
            Path(__file__).with_name("m7_digital_twin.py"),
        )
    }


def make_figures(
    timeseries: pd.DataFrame,
    summary: pd.DataFrame,
    skill: pd.DataFrame,
    figure_dir: Path,
) -> list[str]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure_dir.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.2), constrained_layout=True)

    ts = timeseries.query("node == 'N23' and variable == 'head_m'")
    axes[0, 0].plot(ts.time, ts.truth, color="#111827", lw=2.0, label="hidden truth")
    axes[0, 0].plot(
        ts.time, ts.updated_mean, color="#0072B2", lw=1.7, label="updated twin"
    )
    axes[0, 0].fill_between(
        ts.time,
        ts.updated_q05,
        ts.updated_q95,
        color="#56B4E9",
        alpha=0.24,
        label="90% interval",
    )
    axes[0, 0].plot(
        ts.time,
        ts.open_loop_mean,
        color="#D55E00",
        lw=1.2,
        ls="--",
        label="open loop",
    )
    observed = ts.dropna(subset=["observation"])
    axes[0, 0].scatter(
        observed.time,
        observed.observation,
        color="#009E73",
        s=11,
        alpha=0.75,
        label="observations",
    )
    evaluation_start = int(ts.loc[ts.forecast_evaluation_period == 1, "time"].min())
    axes[0, 0].axvline(evaluation_start, color="0.35", ls=":", lw=1)
    axes[0, 0].set(
        title="a  Sequential updating at downstream node N23",
        xlabel="Month",
        ylabel="Head (m)",
    )
    axes[0, 0].legend(frameon=False, fontsize=8, ncol=2)

    methods = ["updated_twin", "open_loop", "oracle_persistence"]
    colors = {"updated_twin": "#0072B2", "open_loop": "#D55E00", "oracle_persistence": "#555555"}
    for method in methods:
        sub = summary.query(
            "domain == 'all' and variable == 'head_m' and method == @method"
        ).sort_values("horizon")
        axes[0, 1].plot(
            sub.horizon,
            sub.rmse_mean,
            marker="o",
            lw=1.8,
            color=colors[method],
            label=method.replace("_", " "),
        )
    axes[0, 1].set(
        title="b  Prequential head forecast error",
        xlabel="Forecast horizon (months)",
        ylabel="RMSE (m)",
        xticks=[1, 3, 6],
    )
    axes[0, 1].legend(frameon=False, fontsize=8)

    coverage = summary.query(
        "domain == 'all' and method == 'updated_twin'"
    )
    for variable, group in coverage.groupby("variable"):
        axes[1, 0].plot(
            group.sort_values("horizon").horizon,
            group.sort_values("horizon").coverage90_mean,
            marker="o",
            lw=1.7,
            label=variable.replace("_", " "),
        )
    axes[1, 0].axhline(0.90, color="0.25", ls=":", lw=1)
    axes[1, 0].set(
        title="c  Predictive uncertainty calibration",
        xlabel="Forecast horizon (months)",
        ylabel="Empirical 90% coverage",
        xticks=[1, 3, 6],
        ylim=(0.0, 1.02),
    )
    axes[1, 0].legend(frameon=False, fontsize=8)

    ablation = skill.query(
        "horizon == 3 and variable == 'head_m' and "
        "comparator in ['open_loop', 'wrong_topology_updated', "
        "'shuffled_observation_updated', 'sensor_dropout_updated']"
    ).copy()
    ablation["label"] = ablation.comparator.str.replace("_", " ", regex=False)
    axes[1, 1].barh(
        ablation.label,
        100.0 * ablation.mean_rmse_skill,
        color=["#D55E00", "#CC79A7", "#E69F00", "#009E73"][: len(ablation)],
    )
    axes[1, 1].axvline(0.0, color="0.2", lw=0.8)
    axes[1, 1].set(
        title="d  Value of updating and information quality",
        xlabel="Updated-twin RMSE skill (%)",
    )

    for axis in axes.flat:
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)

    png = figure_dir / "figure_m7_digital_twin_extension.png"
    pdf = figure_dir / "figure_m7_digital_twin_extension.pdf"
    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    return [str(png.relative_to(HERE)), str(pdf.relative_to(HERE))]


def results_note(config: TwinConfig, summary: pd.DataFrame, skill: pd.DataFrame) -> str:
    headline = skill.query(
        "comparator == 'open_loop' and horizon == 3 and variable == 'head_m'"
    ).iloc[0]
    coverage = summary.query(
        "method == 'updated_twin' and horizon == 3 and "
        "variable == 'head_m' and domain == 'all'"
    ).iloc[0]
    oracle = skill.query(
        "comparator == 'oracle_persistence' and horizon == 3 and variable == 'head_m'"
    ).iloc[0]
    wrong = skill.query(
        "comparator == 'wrong_topology_updated' and horizon == 3 and variable == 'head_m'"
    ).iloc[0]
    return f"""# M7 Operational Synthetic-Twin Extension — Results

Generated deterministically from {config.n_replicates} independent process/observation
realisations with {config.ensemble_size} ensemble members. Forecasts are prequential:
the forecast issued at month *t* cannot use observations after month *t*.

## Locked headline diagnostics

- Three-month head RMSE skill versus the open-loop twin:
  **{100 * headline.mean_rmse_skill:.1f}%** (replicate percentile interval
  {100 * headline.skill_q025:.1f}% to {100 * headline.skill_q975:.1f}%);
  the updated twin was better in **{100 * headline.fraction_replicates_better:.1f}%**
  of independent realisations.
- Three-month head RMSE skill versus oracle persistence:
  **{100 * oracle.mean_rmse_skill:.1f}%**. Oracle persistence is intentionally
  conservative because it knows the exact current synthetic state.
- Empirical coverage of the nominal 90% three-month head interval:
  **{100 * coverage.coverage90_mean:.1f}%**. Raw ensemble coverage was
  **{100 * coverage.raw_coverage90_mean:.1f}%**; the mean rolling
  monitoring-residual spread multiplier available at issue time was
  **{coverage.calibration_factor_mean:.2f}**. The correction is reported because
  raw EnKF uncertainty was substantially under-dispersed.
- Three-month head RMSE skill versus the updated wrong-topology control:
  **{100 * wrong.mean_rmse_skill:.1f}%**.

## Interpretation

This extension demonstrates an operational mechanism: sparse observations update a
graph-constrained state ensemble, and that updated ensemble can be evaluated on
future states that were unavailable at issue time. It does **not** establish a field
digital twin, aquifer-specific fidelity, management safety, or transfer to Ghana.

The hidden truth and operational model intentionally differ in topology, coefficients,
nonlinearity, stochastic forcing, and an unmodelled chemistry pulse. This reduces the
inverse crime but does not reproduce all structural errors present in a real aquifer.
See `m7_digital_twin_protocol.md` for the claim boundary and required field extension.
"""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replicates", type=int, default=24)
    parser.add_argument("--ensemble-size", type=int, default=80)
    parser.add_argument("--seed", type=int, default=20260718)
    parser.add_argument("--quick", action="store_true", help="4 replicates, 32 members")
    parser.add_argument("--no-figures", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    config = TwinConfig(
        n_replicates=4 if args.quick else args.replicates,
        ensemble_size=32 if args.quick else args.ensemble_size,
        seed=args.seed,
    )
    config.validate()
    RESULTS.mkdir(parents=True, exist_ok=True)
    DOCS.mkdir(parents=True, exist_ok=True)
    before = locked_hashes()

    metric_frames: list[pd.DataFrame] = []
    diagnostic_frames: list[pd.DataFrame] = []
    timeseries: pd.DataFrame | None = None
    for replicate in range(config.n_replicates):
        metrics, diagnostics, representative = run_replicate(config, replicate)
        metric_frames.append(metrics)
        diagnostic_frames.append(diagnostics)
        if representative is not None:
            timeseries = representative
        print(f"replicate {replicate + 1}/{config.n_replicates} complete")

    all_metrics = pd.concat(metric_frames, ignore_index=True)
    all_diagnostics = pd.concat(diagnostic_frames, ignore_index=True)
    if timeseries is None:
        raise RuntimeError("Representative time series was not generated.")
    summary = summarize_metrics(all_metrics)
    skill = paired_skill(all_metrics)

    all_metrics.to_csv(RESULTS / "m7_dt_prequential_metrics.csv", index=False)
    all_diagnostics.to_csv(RESULTS / "m7_dt_assimilation_diagnostics.csv", index=False)
    timeseries.to_csv(RESULTS / "m7_dt_representative_timeseries.csv", index=False)
    summary.to_csv(RESULTS / "m7_dt_metric_summary.csv", index=False)
    skill.to_csv(RESULTS / "m7_dt_paired_skill.csv", index=False)

    figure_paths: list[str] = []
    if not args.no_figures:
        figure_paths = make_figures(timeseries, summary, skill, FIGURES)

    after = locked_hashes()
    if before != after:
        raise RuntimeError("Locked M7 baseline results changed during the extension run.")

    metadata = {
        "benchmark_name": "M7 operational synthetic-twin extension",
        "claim_class": "controlled synthetic operational capability",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "configuration": {
            **config.__dict__,
            "horizons": list(config.horizons),
        },
        "python": sys.version,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "source_sha256": source_hashes(),
        "locked_m7_baseline_sha256_before": before,
        "locked_m7_baseline_sha256_after": after,
        "locked_baseline_unchanged": before == after,
        "figures": figure_paths,
        "limitations": [
            "synthetic truth rather than a named field aquifer",
            "graph state-space surrogate rather than MODFLOW/PHREEQC co-simulation",
            "future seasonal climatology and pumping schedule treated as known inputs",
            "EnKF Gaussian update can underrepresent multimodality",
            "management decisions and safety constraints are not evaluated",
        ],
    }
    (RESULTS / "m7_dt_run_metadata.json").write_text(
        json.dumps(metadata, indent=2, default=str),
        encoding="utf-8",
    )
    (DOCS / "m7_digital_twin_results.md").write_text(
        results_note(config, summary, skill),
        encoding="utf-8",
    )
    print(f"M7 operational synthetic-twin extension complete -> {RESULTS}")


if __name__ == "__main__":
    main()
