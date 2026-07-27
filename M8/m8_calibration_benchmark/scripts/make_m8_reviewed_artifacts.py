"""Create scientifically reviewed presentation artifacts from the locked M8 run."""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PROJECT = Path(__file__).resolve().parents[1]
RUN_ID = "RUN-M8-CONFIRM-20260728-01"
ARTIFACTS = PROJECT / "manuscript" / "artifacts"
RUN_ROOT = PROJECT / "provenance" / "runs" / RUN_ID
PARAMETERS = ("dispersivity", "decay")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    manifest = json.loads((RUN_ROOT / "run_manifest.json").read_text(encoding="utf-8"))
    if manifest.get("status") != "PASS":
        raise RuntimeError(f"Source run is not PASS: {RUN_ID}")

    design = pd.read_csv(ARTIFACTS / "m8_transport_parameter_summary.csv")
    oed = pd.read_csv(ARTIFACTS / "m8_transport_oed_summary.csv")
    kinetics = pd.read_csv(ARTIFACTS / "m8_kinetics_structural_summary.csv")
    candidates = pd.read_csv(ARTIFACTS / "m8_transport_candidate_scores.csv")

    # Random uses a different candidate in each paired replicate; no single time
    # belongs in a strategy-level summary.
    oed.loc[oed["strategy"] == "random", "candidate_time_days"] = np.nan
    reviewed_table = ARTIFACTS / "m8_transport_oed_summary_reviewed.csv"
    oed.to_csv(reviewed_table, index=False, lineterminator="\n")

    colours = {"dispersivity": "#1f77b4", "decay": "#d62728"}
    plt.rcParams.update({"font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9})
    fig, axes = plt.subplots(2, 2, figsize=(12, 9), constrained_layout=True)

    ax = axes[0, 0]
    for parameter in PARAMETERS:
        frame = design.loc[design["parameter"] == parameter]
        ax.scatter(
            np.log10(frame["condition"]),
            frame["median_abs_log10_error"],
            color=colours[parameter],
            label=parameter,
            alpha=0.85,
        )
    ax.set_xlabel("log10 condition number (whitened log-parameter FIM)")
    ax.set_ylabel("Median absolute log10 error")
    ax.set_title("A  Fixed-design recovery varies by parameter")
    ax.legend(frameon=False)
    ax.grid(alpha=0.2)

    ax = axes[0, 1]
    for parameter, marker in (("dispersivity", "o"), ("decay", "s")):
        ax.plot(
            candidates["candidate_time_days"],
            candidates[f"sd_log10_{parameter}"],
            marker + "-",
            color=colours[parameter],
            label=parameter,
        )
    ax.axvline(50.0, color="black", linestyle="--", linewidth=0.8)
    ax.annotate("common optimum: 50 d", (53, 1.18), fontsize=8)
    ax.set_xlabel("Added sampling time (d)")
    ax.set_ylabel("Predicted marginal SD (log10 units)")
    ax.set_title("B  All local criteria select the front observation")
    ax.legend(frameon=False)
    ax.grid(alpha=0.2)

    ax = axes[1, 0]
    local = oed.loc[oed["start_regime"] == "local"].copy()
    strategies = ["balanced", "random", "worst_joint", "no_new_measurement"]
    labels = ["common optimum\n(50 d)", "random\n(varied)", "redundant\n(240 d)", "no new"]
    x = np.arange(len(strategies), dtype=float)
    width = 0.36
    for offset, parameter in ((-width / 2, "dispersivity"), (width / 2, "decay")):
        frame = (
            local.loc[local["parameter"] == parameter]
            .set_index("strategy")
            .reindex(strategies)
        )
        ax.bar(
            x + offset,
            frame["median_abs_log10_error"],
            width,
            color=colours[parameter],
            label=parameter,
        )
    ax.set_xticks(x, labels)
    ax.set_ylabel("Median absolute log10 error")
    ax.set_title("C  One front sample improves both parameters")
    ax.legend(frameon=False)
    ax.grid(axis="y", alpha=0.2)

    ax = axes[1, 1]
    kinetic_labels = kinetics["design"].str.replace("_", " ").tolist()
    ratios = np.maximum(kinetics["eigen_ratio"].to_numpy(float), 1e-18)
    bars = ax.bar(
        np.arange(len(kinetic_labels)), ratios, color=["#999999", "#666666", "#2ca02c"]
    )
    ax.set_yscale("log")
    ax.set_xticks(np.arange(len(kinetic_labels)), kinetic_labels, rotation=15)
    ax.set_ylabel("Smallest/largest FIM eigenvalue")
    ax.set_title("D  Kinetic k*A confounding needs external A evidence")
    ax.axhline(1e-10, color="black", linestyle="--", linewidth=0.8, label="rank tolerance")
    ax.legend(frameon=False)
    for bar, rank in zip(bars, kinetics["rank"]):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height(),
            f"rank {int(rank)}",
            ha="center",
            va="bottom",
        )
    ax.grid(axis="y", alpha=0.2)

    fig.suptitle(
        "M8 confirmatory controlled-synthetic benchmark",
        fontsize=13,
        fontweight="bold",
    )
    figure = ARTIFACTS / "m8_confirmatory_figure_reviewed.png"
    fig.savefig(figure, dpi=300)
    fig.savefig(figure.with_suffix(".pdf"))
    plt.close(fig)

    outputs = [reviewed_table, figure, figure.with_suffix(".pdf")]
    reviewed_manifest = {
        "schema_version": "1.0",
        "source_run_id": RUN_ID,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "amendment_class": "presentation_only",
        "numerical_results_changed": False,
        "producer": "scripts/make_m8_reviewed_artifacts.py",
        "producer_sha256": sha256_file(Path(__file__)),
        "outputs": [
            {
                "path": path.relative_to(PROJECT).as_posix(),
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
            for path in outputs
        ],
    }
    (PROJECT / "provenance" / "reviewed_artifact_manifest.json").write_text(
        json.dumps(reviewed_manifest, indent=2), encoding="utf-8"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
