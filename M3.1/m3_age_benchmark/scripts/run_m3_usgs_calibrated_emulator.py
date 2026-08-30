"""USGS-calibrated parity emulator for Hydrosheaf M3.

This script calibrates strict-parity Hydrosheaf output against USGS-reported
ages using leave-study-unit-out cross-validation. It is explicitly labelled
as calibrated emulation, not independent validation.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any, Dict, List, Mapping, Sequence

import numpy as np
import pandas as pd
from sklearn.linear_model import Ridge
from sklearn.isotonic import IsotonicRegression

REPO_ROOT = Path(__file__).resolve().parents[3]
RESULTS_DIR = REPO_ROOT / "M3" / "m3_age_benchmark" / "results"


def _finite_float(val: Any) -> float | None:
    try:
        f = float(val)
        return f if math.isfinite(f) else None
    except (TypeError, ValueError):
        return None


def _log10(val: Any) -> float | None:
    f = _finite_float(val)
    if f is not None and f > 0:
        return math.log10(f)
    return None


def _build_features(df: pd.DataFrame) -> pd.DataFrame:
    """Return feature matrix from strict-parity output."""
    feats = pd.DataFrame()
    age_source = df["est_age_total_years"] if "est_age_total_years" in df.columns else df["est_age_multi"]
    feats["log10_est_age"] = age_source.apply(lambda x: _log10(x))
    feats["n_tracers"] = pd.to_numeric(df.get("n_tracers", 0), errors="coerce")
    feats["reported_model_dm"] = (df.get("reported_model", "") == "DM").astype(float)
    feats["reported_model_em"] = (df.get("reported_model", "") == "EM").astype(float)
    feats["reported_model_pfm"] = (df.get("reported_model", "") == "PFM").astype(float)
    feats["reported_model_ga"] = (df.get("reported_model", "") == "GA").astype(float)
    feats["age_class_intermediate"] = (df.get("age_class", "") == "intermediate_50_1000").astype(float)
    feats["age_class_old"] = (df.get("age_class", "") == "old_1000_30000").astype(float)
    feats["age_class_very_old"] = (df.get("age_class", "") == "very_old_gt_30000").astype(float)
    feats["log10_rpt_chisquare"] = pd.to_numeric(df.get("Rpt_ChiSquare", np.nan), errors="coerce").apply(lambda x: math.log10(max(x, 1e-6)) if _finite_float(x) else np.nan)
    feats["frac_anthropocene"] = pd.to_numeric(df.get("FracAnthropocene", np.nan), errors="coerce")
    feats["frac_holocene"] = pd.to_numeric(df.get("FracHolocene", np.nan), errors="coerce")
    feats["frac_pleistocene"] = pd.to_numeric(df.get("FracPleistocene", np.nan), errors="coerce")
    feats["depth_m"] = pd.to_numeric(df.get("Depth_m", np.nan), errors="coerce")
    feats = feats.fillna(0.0)
    return feats


def _target_log10_age(df: pd.DataFrame) -> np.ndarray:
    """Return log10 of USGS reported total age."""
    if "Rpt_TotAge_yrs" in df.columns:
        source = df["Rpt_TotAge_yrs"]
    elif "ref_age" in df.columns:
        source = df["ref_age"]
    else:
        raise ValueError("Input must contain `Rpt_TotAge_yrs` or `ref_age`.")
    return np.array([_log10(v) for v in source], dtype=object)


def leave_study_unit_out_calibration(
    df: pd.DataFrame,
    *,
    study_unit_col: str = "StudyUnit",
    method: str = "ridge",
) -> pd.DataFrame:
    """Run leave-study-unit-out calibration and return predictions."""
    df = df.copy()
    df["_target"] = _target_log10_age(df)
    df["_fold"] = ""
    df["_calibrated"] = np.nan
    df["_raw_pred"] = np.nan
    df["_residual"] = np.nan
    df["_held_out_study_unit"] = ""

    if study_unit_col not in df.columns:
        raise ValueError(f"Input must contain `{study_unit_col}` for leave-study-unit-out calibration.")

    study_units = df[study_unit_col].dropna().unique()
    features = _build_features(df)
    feature_cols = list(features.columns)
    X = features.values
    y = pd.to_numeric(df["_target"], errors="coerce").values

    for su in study_units:
        test_mask = df[study_unit_col] == su
        train_mask = (~test_mask) & np.isfinite(y) & np.isfinite(X).all(axis=1)
        test_mask = test_mask & np.isfinite(X).all(axis=1)
        if train_mask.sum() < 10 or test_mask.sum() < 1:
            continue

        X_train, y_train = X[train_mask], y[train_mask]
        X_test = X[test_mask]

        if method == "ridge":
            model = Ridge(alpha=1.0)
        elif method == "isotonic":
            # Isotonic regression works on 1D; use primary feature
            model = IsotonicRegression(out_of_bounds="clip")
            X_train = X_train[:, 0]
            X_test = X_test[:, 0]
        else:
            raise ValueError(f"Unknown calibration method: {method}")

        model.fit(X_train, y_train)
        raw_pred = X_test if method == "isotonic" else X_test[:, 0]
        calibrated = model.predict(X_test)

        df.loc[test_mask, "_fold"] = su
        df.loc[test_mask, "_calibrated"] = calibrated
        df.loc[test_mask, "_raw_pred"] = raw_pred
        df.loc[test_mask, "_residual"] = calibrated - y[test_mask]
        df.loc[test_mask, "_held_out_study_unit"] = su

    return df


def summarize_calibration(result: pd.DataFrame) -> pd.DataFrame:
    """Return headline metrics for calibrated-emulation predictions."""
    target = pd.to_numeric(result["_target"], errors="coerce")
    raw = pd.to_numeric(result["_raw_pred"], errors="coerce")
    calibrated = pd.to_numeric(result["_calibrated"], errors="coerce")
    valid_raw = target.notna() & raw.notna()
    valid_calibrated = target.notna() & calibrated.notna()

    def metrics(pred: pd.Series, valid: pd.Series, label: str) -> dict[str, float | str]:
        residual = pred[valid] - target[valid]
        y = target[valid]
        ss_tot = float(((y - y.mean()) ** 2).sum()) if len(y) else float("nan")
        ss_res = float((residual**2).sum()) if len(residual) else float("nan")
        return {
            "mode": label,
            "N": int(valid.sum()),
            "median_abs_log10_error": float(residual.abs().median()) if len(residual) else float("nan"),
            "log10_rmse": float(np.sqrt(np.mean(residual**2))) if len(residual) else float("nan"),
            "log10_r2": 1.0 - ss_res / ss_tot if ss_tot and math.isfinite(ss_tot) else float("nan"),
            "within_factor_2": float((residual.abs() <= math.log10(2.0)).mean()) if len(residual) else float("nan"),
            "within_factor_10": float((residual.abs() <= 1.0).mean()) if len(residual) else float("nan"),
        }

    return pd.DataFrame(
        [
            metrics(raw, valid_raw, "uncalibrated_strict_parity_on_heldout_folds"),
            metrics(calibrated, valid_calibrated, "usgs_calibrated_benchmark_emulation"),
        ]
    )


def main() -> int:
    parser = argparse.ArgumentParser(description="USGS-calibrated parity emulator")
    parser.add_argument("--input", type=Path, required=True, help="Strict-parity CSV output")
    parser.add_argument("--output", type=Path, default=RESULTS_DIR / "m3_usgs_calibrated_parity.csv")
    parser.add_argument("--method", choices=["ridge", "isotonic"], default="ridge")
    args = parser.parse_args()

    df = pd.read_csv(args.input)
    result = leave_study_unit_out_calibration(df, method=args.method)

    out = pd.DataFrame({
        "site_id": result["site_id"],
        "calibration_method": args.method,
        "fold_id": result["_fold"],
        "held_out_study_unit": result["_held_out_study_unit"],
        "log10_est_age": result["_raw_pred"],
        "log10_calibrated_age": result["_calibrated"],
        "log10_residual": result["_residual"],
        "log10_reported_age": result["_target"],
    })

    args.output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, index=False)
    summary = summarize_calibration(result)
    summary_path = args.output.with_name(f"{args.output.stem}_summary.csv")
    summary.to_csv(summary_path, index=False)
    manifest_path = args.output.with_name(f"{args.output.stem}_manifest.json")
    manifest_path.write_text(
        json.dumps(
            {
                "input": str(args.input),
                "output": str(args.output),
                "summary": str(summary_path),
                "method": args.method,
                "mode": "USGS-calibrated benchmark emulation",
                "n_rows": int(len(out)),
                "n_calibrated_rows": int(pd.to_numeric(out["log10_calibrated_age"], errors="coerce").notna().sum()),
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    print(f"Wrote {len(out)} rows to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
