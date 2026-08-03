"""TAB-1: common stress/signal taxonomy applied to every retained M6, M7 and
M8 experiment cited in the manuscript.

Classification only (stress axis; internal signal reported; external signal
reported; whether a negative/structural control exists), applied by direct
inspection of each component's own locked DECISIONS.md, README.md and result
files. No experiment is re-labelled in a way that contradicts its own
component's claim ledger (DECISIONS.md D2, C1).

Run:  .venv/Scripts/python.exe O4/analysis/python/derive_taxonomy.py
"""

from __future__ import annotations

import pandas as pd

from _common import write

ROWS = [
    dict(
        component="M6", experiment="Tier ablation T4->T0 (E3)",
        stress_axis="data limitation (evidence tier)",
        internal_signal="mean mechanism-resolution score (MRS)",
        external_signal="true/gated identifiability class (% non-identifiable)",
        has_negative_control=True, control_type="evidence-gate-off (ungated) ablation",
        has_synthetic_ground_truth=True,
        source="results/m6_tier_ablation.csv; results/m6_field_gate_structural.csv",
    ),
    dict(
        component="M6", experiment="Edge-set perturbation (E4)",
        stress_axis="model-form (assumed connectivity)",
        internal_signal="per-edge MRS / identifiability class",
        external_signal="network-level process-composition shift (total-variation distance)",
        has_negative_control=False, control_type="",
        has_synthetic_ground_truth=False,
        source="results/m6_edge_sensitivity.csv",
    ),
    dict(
        component="M6", experiment="External sparse transfer, Talensi/Lower Anayari (E5)",
        stress_axis="data limitation + domain shift",
        internal_signal="mean MRS",
        external_signal="% non-identifiable",
        has_negative_control=False, control_type="",
        has_synthetic_ground_truth=False,
        source="results/m6_external_summary.csv",
    ),
    dict(
        component="M6", experiment="Synthetic validation with known truth (ED2)",
        stress_axis="model-form (independent extended model)",
        internal_signal="exact-mineral F1 (fit-quality proxy)",
        external_signal="true dominant-process recovery by process class",
        has_negative_control=False, control_type="",
        has_synthetic_ground_truth=True,
        source="results/m6_synthetic_recovery_by_tier.csv",
    ),
    dict(
        component="M7", experiment="Native evidence-panel integration (age/chem/hydraulics increments)",
        stress_axis="evidence combination",
        internal_signal="posterior mean edge-entropy reduction",
        external_signal="PR-AUC / Brier / log-loss vs independent MODFLOW6/MODPATH7 truth",
        has_negative_control=True, control_type="permuted-stream adverse controls (same experiment)",
        has_synthetic_ground_truth=True,
        source="results/m7_3_locked/evidence_case_bootstrap_contrasts.csv",
    ),
    dict(
        component="M7", experiment="Adverse controls: permuted age / permuted hydraulics / joint",
        stress_axis="evidence misspecification",
        internal_signal="posterior mean edge-entropy reduction",
        external_signal="PR-AUC / Brier / log-loss vs independent MODFLOW6/MODPATH7 truth",
        has_negative_control=True, control_type="self (each adverse condition IS the control)",
        has_synthetic_ground_truth=True,
        source="results/m7_3_locked/evidence_case_bootstrap_contrasts.csv",
    ),
    dict(
        component="M7", experiment="M7.4 sheaf vs graph, heterogeneous-affine stratum",
        stress_axis="model-form (restriction-map heterogeneity)",
        internal_signal="log-loss / expected calibration error",
        external_signal="PR-AUC / selected-edge F1 vs independent MODFLOW6/MODPATH7 truth",
        has_negative_control=True, control_type="permuted-map adverse control (same run)",
        has_synthetic_ground_truth=True,
        source="results/RUN-M7-SHEAF-VS-GRAPH-20260729-01/paired_bootstrap_contrasts.csv",
    ),
    dict(
        component="M7", experiment="M7.5 robust/local-first-global-fallback hybrid",
        stress_axis="evidence misspecification (false-edge downweighting)",
        internal_signal="Brier score / log-loss",
        external_signal="PR-AUC vs independent MODFLOW6/MODPATH7 truth",
        has_negative_control=True, control_type="permuted-map adverse control (same run)",
        has_synthetic_ground_truth=True,
        source="results/RUN-M7-ROBUST-HYBRID-SHEAF-20260729-01/locked_test/paired_bootstrap_contrasts.csv",
    ),
    dict(
        component="M8", experiment="Fixed 16-design sweep, matched analytical model",
        stress_axis="data limitation (sampling-design placement)",
        internal_signal="optimiser success rate; nominal 95% coverage",
        external_signal="realised median absolute log10 parameter-recovery error vs known truth",
        has_negative_control=False, control_type="",
        has_synthetic_ground_truth=True,
        source="manuscript/artifacts/m8_transport_parameter_summary.csv",
    ),
    dict(
        component="M8", experiment="Optimal-design experiment, independent numerical truth",
        stress_axis="model-form (calibration model vs independent numerical solver)",
        internal_signal="linearised 95% interval coverage",
        external_signal="realised median absolute log10 parameter-recovery error vs known truth",
        has_negative_control=True, control_type="no-new-measurement / random-candidate baselines",
        has_synthetic_ground_truth=True,
        source="results/RUN-M8-INDEPENDENT-20260728-01/strategy_summary.csv",
    ),
    dict(
        component="M8", experiment="Kinetic rate-law structural confound",
        stress_axis="data limitation (measurement type: residence time only vs +surface area)",
        internal_signal="optimiser convergence; objective value",
        external_signal="Fisher-information numerical rank / condition number",
        has_negative_control=True, control_type="doubled-product off-ridge sensitivity check",
        has_synthetic_ground_truth=True,
        source="manuscript/artifacts/m8_kinetics_structural_summary.csv",
    ),
    dict(
        component="M8", experiment="Frontier topology active-learning qualification",
        stress_axis="evidence-value / sequential decision design",
        internal_signal="joint-hypothesis entropy reduction per unit cost",
        external_signal="edge Brier score / PR-AUC vs independent MODFLOW6/MODPATH7 truth",
        has_negative_control=True, control_type="random-feasible-acquisition baseline; realised oracle ceiling",
        has_synthetic_ground_truth=True,
        source="provenance/runs/RUN-M8-FRONTIER-AL-20260728-01/strategy_summary.csv",
    ),
]


def main() -> None:
    df = pd.DataFrame(ROWS)
    write(df, "taxonomy.csv")


if __name__ == "__main__":
    main()
