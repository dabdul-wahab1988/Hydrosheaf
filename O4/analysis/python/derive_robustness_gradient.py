"""TAB-2/FIG-3 (M6 rows): the tier-ablation internal-vs-external divergence.

Reads M6's per-well, per-tier locked result table (800 rows: 160 wells x 5
tiers) and recomputes the tier summary directly from those rows, rather than
transcribing `docs/m6_results_summary.md`'s own prose table. This reproduces
that prose table exactly (cross-checked as part of DECISIONS.md D3) and
additionally emits the negative control (evidence-gate-off) and the external
sparse-transfer (Talensi/Lower Anayari) rows on the same schema so FIG-3 can
plot all of them together.

Run:  .venv/Scripts/python.exe O4/analysis/python/derive_robustness_gradient.py
"""

from __future__ import annotations

import pandas as pd

from _common import M6, wilson_ci, write

TIER_ORDER = [
    "tier4_full_metadata",
    "tier3_sr_sio2",
    "tier2_fluoride",
    "tier1_isotopes",
    "tier0_majors",
]
TIER_LABEL = {
    "tier4_full_metadata": "T4 (full)",
    "tier3_sr_sio2": "T3 (+Sr/SiO2)",
    "tier2_fluoride": "T2 (+F)",
    "tier1_isotopes": "T1 (+isotopes)",
    "tier0_majors": "T0 (majors)",
}


def tier_ablation() -> pd.DataFrame:
    df = pd.read_csv(M6 / "results" / "m6_tier_ablation.csv")
    rows = []
    for tier in TIER_ORDER:
        g = df.loc[df["tier"] == tier]
        n = len(g)
        k_non_id = int((g["resolution_class"] == "non_identifiable").sum())
        pct_non_id = float(k_non_id / n * 100)
        pct_partial = float((g["resolution_class"] == "partially_identifiable").mean() * 100)
        ci_lo, ci_hi = wilson_ci(k_non_id, n)
        rows.append(
            dict(
                component="M6 robustness",
                axis_x=TIER_LABEL[tier],
                tier_order=TIER_ORDER.index(tier),
                n=n,
                internal_signal_name="mean_mrs",
                internal_signal_value=float(g["mrs"].mean()),
                external_signal_name="pct_non_identifiable",
                external_signal_value=pct_non_id,
                external_ci_low=round(ci_lo * 100, 1),
                external_ci_high=round(ci_hi * 100, 1),
                external_signal_pct_partial=pct_partial,
                condition="field_dataset_tier_ablation",
                source="results/m6_tier_ablation.csv, recomputed directly",
            )
        )
    return pd.DataFrame(rows)


def negative_control() -> pd.DataFrame:
    """Evidence-gate-off ablation and structural-rank diagnostic: proves the
    tier collapse is the conservative prior, not a classifier artefact
    (DECISIONS.md C4), and gives the independent structural-rank/silicate-
    coherence direction check referenced alongside it."""
    g = pd.read_csv(M6 / "results" / "m6_field_gate_structural.csv").set_index("tier")
    g = g.reindex(TIER_ORDER)
    rows = []
    for tier in TIER_ORDER:
        row = g.loc[tier]
        rows.append(
            dict(
                component="M6 robustness",
                axis_x=TIER_LABEL[tier],
                tier_order=TIER_ORDER.index(tier),
                n=160,
                internal_signal_name="mean_mrs",
                internal_signal_value=float("nan"),
                external_signal_name="pct_non_identifiable_gate_off",
                external_signal_value=float(row["ungated_frac_non_identifiable"]) * 100,
                external_signal_pct_partial=float("nan"),
                structural_rank=int(row["rank"]),
                structural_nullity=int(row["nullity"]),
                silicate_coherence=float(row["silicate_coherence"]),
                condition="evidence_gate_disabled_negative_control",
                source="results/m6_field_gate_structural.csv",
            )
        )
    return pd.DataFrame(rows)


def external_transfer() -> pd.DataFrame:
    """Talensi / Lower Anayari external sparse transfer (E5): the internal
    signal (mean MRS) barely separates a 37%-non-identifiable dataset from a
    95%-non-identifiable one."""
    df = pd.read_csv(M6 / "results" / "m6_external_summary.csv")
    df = df.rename(
        columns={
            "mean_mrs": "internal_signal_value",
            "frac_non_identifiable": "external_signal_value",
            "frac_partial": "external_signal_pct_partial",
        }
    )
    df["external_signal_value"] = df["external_signal_value"] * 100
    df["external_signal_pct_partial"] = df["external_signal_pct_partial"] * 100
    cis = df.apply(
        lambda r: wilson_ci(round(r["external_signal_value"] / 100 * r["n_edges"]), int(r["n_edges"])),
        axis=1,
    )
    df["external_ci_low"] = [round(c[0] * 100, 1) for c in cis]
    df["external_ci_high"] = [round(c[1] * 100, 1) for c in cis]
    df["internal_signal_name"] = "mean_mrs"
    df["external_signal_name"] = "pct_non_identifiable"
    df["component"] = "M6 robustness"
    df["condition"] = "external_field_transfer"
    df["source"] = "results/m6_external_summary.csv"
    return df


def main() -> None:
    parts = [tier_ablation(), negative_control()]
    df = pd.concat([p for p in parts if len(p)], ignore_index=True)
    write(df, "robustness_gradient.csv")

    ext = external_transfer()
    write(ext, "robustness_external_transfer.csv")


if __name__ == "__main__":
    main()
