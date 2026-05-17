"""Audit raw versus USGS dissolved-gas-corrected tracer effects for M3."""

from __future__ import annotations

import argparse
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[3]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

SCRIPT_DIR = Path(__file__).resolve().parent
BENCHMARK_DIR = SCRIPT_DIR.parents[0]
RESULT_DIR = BENCHMARK_DIR / "results"
DOCS_DIR = BENCHMARK_DIR / "docs"
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import run_m3_usgs_benchmark as usgs


DEFAULT_DESIGN_RESULTS = RESULT_DIR / "m3_design_matrix_results.csv"
DEFAULT_AUDIT_OUTPUT = RESULT_DIR / "m3_gas_correction_audit.csv"
DEFAULT_SUMMARY_OUTPUT = RESULT_DIR / "m3_gas_correction_audit_summary.csv"

GAS_FIELDS = (
    ("tritium", "tritium_TU"),
    ("he3_trit", "he3_trit_TU"),
    ("sf6", "sf6_pptv"),
    ("cfc11", "cfc11_pptv"),
    ("cfc12", "cfc12_pptv"),
    ("cfc113", "cfc113_pptv"),
)


def _numeric(value) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return float("nan")
    return number if math.isfinite(number) else float("nan")


def _changed(raw: float, corrected: float, *, rel_tol: float = 1.0e-6, abs_tol: float = 1.0e-12) -> bool:
    if not math.isfinite(raw) or not math.isfinite(corrected):
        return False
    return abs(raw - corrected) > max(abs_tol, rel_tol * max(abs(raw), abs(corrected), 1.0))


def _safe_ratio(corrected: float, raw: float) -> float:
    if not math.isfinite(corrected) or not math.isfinite(raw) or raw == 0:
        return float("nan")
    return corrected / raw


def _apparent_3h3he_age(tritium: float, he3_trit: float) -> float:
    tritium = _numeric(tritium)
    he3_trit = _numeric(he3_trit)
    if not math.isfinite(tritium) or not math.isfinite(he3_trit) or tritium <= 0 or he3_trit < 0:
        return float("nan")
    lambda_tritium = math.log(2.0) / 12.32
    return math.log1p(he3_trit / tritium) / lambda_tritium


def _paired_gas_effects(design_results: pd.DataFrame) -> pd.DataFrame:
    baseline = design_results[design_results["scenario_id"] == "parity_reported_corrected"][
        ["site_id", "log10_error", "est_age_multi", "within_factor_2"]
    ].rename(
        columns={
            "log10_error": "parity_log10_error",
            "est_age_multi": "parity_est_age_multi",
            "within_factor_2": "parity_within_factor_2",
        }
    )
    raw = design_results[design_results["scenario_id"] == "ablation_raw_gases"][
        ["site_id", "log10_error", "est_age_multi", "within_factor_2"]
    ].rename(
        columns={
            "log10_error": "raw_gas_log10_error",
            "est_age_multi": "raw_gas_est_age_multi",
            "within_factor_2": "raw_gas_within_factor_2",
        }
    )
    paired = baseline.merge(raw, on="site_id", how="inner")
    paired["delta_log10_error_raw_minus_corrected"] = paired["raw_gas_log10_error"] - paired["parity_log10_error"]
    paired["raw_gases_improved"] = paired["delta_log10_error_raw_minus_corrected"] < 0
    paired["raw_gases_degraded"] = paired["delta_log10_error_raw_minus_corrected"] > 0
    paired["raw_gases_gained_factor_2"] = ~paired["parity_within_factor_2"].astype(bool) & paired["raw_gas_within_factor_2"].astype(bool)
    paired["raw_gases_lost_factor_2"] = paired["parity_within_factor_2"].astype(bool) & ~paired["raw_gas_within_factor_2"].astype(bool)
    return paired


def build_gas_audit(design_results: pd.DataFrame, source_df: pd.DataFrame) -> pd.DataFrame:
    paired = _paired_gas_effects(design_results)
    source_cols = [
        "site_id", "reference_age_years", "LPM_TracersMod", "dissolved_gas_correction",
        "dgm_name", "dgm_temp_c", "dgm_excess_air_cckg", "dgm_fractionation", "dgm_excess_n2_mgl",
    ]
    for _, field in GAS_FIELDS:
        source_cols.extend([field, f"raw_{field}"])
    source_cols = [col for col in source_cols if col in source_df.columns]
    source = source_df[source_cols].copy()
    audit = paired.merge(source, on="site_id", how="left")
    audit["ref_age"] = audit["reference_age_years"].map(usgs._parse_age)
    audit["age_class"] = audit["ref_age"].map(usgs._age_class)

    changed_columns = []
    ratio_columns = []
    for label, field in GAS_FIELDS:
        raw_col = f"raw_{field}"
        corrected_col = field
        if raw_col not in audit.columns or corrected_col not in audit.columns:
            continue
        raw_values = audit[raw_col].map(_numeric)
        corrected_values = audit[corrected_col].map(_numeric)
        delta_col = f"{label}_corrected_minus_raw"
        ratio_col = f"{label}_corrected_raw_ratio"
        changed_col = f"{label}_changed_by_dgm"
        audit[delta_col] = corrected_values - raw_values
        audit[ratio_col] = [
            _safe_ratio(corrected, raw)
            for corrected, raw in zip(corrected_values, raw_values)
        ]
        audit[changed_col] = [
            _changed(raw, corrected)
            for raw, corrected in zip(raw_values, corrected_values)
        ]
        changed_columns.append(changed_col)
        ratio_columns.append(ratio_col)

    audit["n_changed_gas_fields"] = audit[changed_columns].sum(axis=1) if changed_columns else 0
    audit["changed_gas_fields"] = audit.apply(
        lambda row: ";".join(col.replace("_changed_by_dgm", "") for col in changed_columns if bool(row.get(col))),
        axis=1,
    )
    audit["any_gas_value_changed"] = audit["n_changed_gas_fields"] > 0

    if {"tritium_TU", "he3_trit_TU", "raw_tritium_TU", "raw_he3_trit_TU"}.issubset(audit.columns):
        audit["corrected_3h3he_apparent_age"] = [
            _apparent_3h3he_age(t, h)
            for t, h in zip(audit["tritium_TU"], audit["he3_trit_TU"])
        ]
        audit["raw_3h3he_apparent_age"] = [
            _apparent_3h3he_age(t, h)
            for t, h in zip(audit["raw_tritium_TU"], audit["raw_he3_trit_TU"])
        ]
        audit["delta_3h3he_apparent_age_raw_minus_corrected"] = (
            audit["raw_3h3he_apparent_age"] - audit["corrected_3h3he_apparent_age"]
        )
    return audit


def summarize_audit(audit: pd.DataFrame) -> pd.DataFrame:
    if audit.empty:
        return pd.DataFrame()
    return (
        audit.groupby("age_class")
        .agg(
            n=("site_id", "count"),
            median_delta_log10_error=("delta_log10_error_raw_minus_corrected", "median"),
            mean_delta_log10_error=("delta_log10_error_raw_minus_corrected", "mean"),
            raw_improved_fraction=("raw_gases_improved", "mean"),
            raw_degraded_fraction=("raw_gases_degraded", "mean"),
            raw_gained_factor_2=("raw_gases_gained_factor_2", "sum"),
            raw_lost_factor_2=("raw_gases_lost_factor_2", "sum"),
            any_gas_value_changed_fraction=("any_gas_value_changed", "mean"),
            median_changed_gas_fields=("n_changed_gas_fields", "median"),
        )
        .reset_index()
    )


def write_report(audit: pd.DataFrame, summary: pd.DataFrame, output: Path) -> None:
    DOCS_DIR.mkdir(parents=True, exist_ok=True)
    improved = audit.sort_values("delta_log10_error_raw_minus_corrected").head(10)
    degraded = audit.sort_values("delta_log10_error_raw_minus_corrected", ascending=False).head(10)
    lines = [
        "# M3 Gas-Correction Audit",
        "",
        f"- Generated: {datetime.now(timezone.utc).isoformat()}",
        f"- Audit output: `{output}`",
        f"- Paired rows: {len(audit)}",
        f"- Rows where raw gases improved log10 error: {int(audit['raw_gases_improved'].sum())}",
        f"- Rows where raw gases degraded log10 error: {int(audit['raw_gases_degraded'].sum())}",
        f"- Rows with any changed gas value: {int(audit['any_gas_value_changed'].sum())}",
        "",
        "## Summary By Age Class",
        "",
        "```text",
        summary.to_string(index=False) if not summary.empty else "No rows.",
        "```",
        "",
        "## Largest Raw-Gas Improvements",
        "",
        "```text",
        improved[[
            "site_id", "age_class", "delta_log10_error_raw_minus_corrected",
            "parity_est_age_multi", "raw_gas_est_age_multi", "ref_age",
            "changed_gas_fields", "dissolved_gas_correction",
        ]].to_string(index=False),
        "```",
        "",
        "## Largest Raw-Gas Degradations",
        "",
        "```text",
        degraded[[
            "site_id", "age_class", "delta_log10_error_raw_minus_corrected",
            "parity_est_age_multi", "raw_gas_est_age_multi", "ref_age",
            "changed_gas_fields", "dissolved_gas_correction",
        ]].to_string(index=False),
        "```",
        "",
        "## Interpretation Guardrail",
        "",
        "A negative delta means the raw-gas ablation was closer to the USGS reference age than the corrected-gas parity run for the same site. This does not prove raw gases are generally superior; it identifies rows where the dissolved-gas correction pathway, tracer-mode masking, or Hydrosheaf atmospheric-equivalent assumptions need targeted review.",
        "",
        "Rows with changed `tritium;he3_trit` should be inspected through `raw_3h3he_apparent_age`, `corrected_3h3he_apparent_age`, and `delta_3h3he_apparent_age_raw_minus_corrected`, because small concentration changes can produce large apparent-age changes when tritium is low.",
        "",
    ]
    (DOCS_DIR / "m3_gas_correction_audit.md").write_text("\n".join(lines), encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Audit M3 raw versus corrected gas-tracer effects.")
    parser.add_argument("--design-results", type=Path, default=DEFAULT_DESIGN_RESULTS)
    parser.add_argument("--output", type=Path, default=DEFAULT_AUDIT_OUTPUT)
    args = parser.parse_args(argv)

    design_results = pd.read_csv(args.design_results)
    source_df = usgs.load_usgs_national_dataset()
    audit = build_gas_audit(design_results, source_df)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    audit.to_csv(args.output, index=False)
    summary = summarize_audit(audit)
    summary.to_csv(DEFAULT_SUMMARY_OUTPUT, index=False)
    write_report(audit, summary, args.output)
    manifest = {
        "run_utc": datetime.now(timezone.utc).isoformat(),
        "design_results": str(args.design_results),
        "output": str(args.output),
        "summary_output": str(DEFAULT_SUMMARY_OUTPUT),
        "paired_rows": int(len(audit)),
    }
    (RESULT_DIR / "m3_gas_correction_audit_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"Wrote {len(audit)} gas-correction audit rows to {args.output}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
