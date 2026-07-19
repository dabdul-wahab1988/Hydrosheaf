"""Leakage-audited one-step field prequential evaluation.

The Northern Ghana workbook has one dry and one wet observation per well.  It
cannot support a multi-year operational digital twin.  It can, however, support
a genuine within-campaign prequential experiment: before each August sampling
batch is revealed, predict that batch using all March observations and only
earlier August batches, then update.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Tuple

import numpy as np
import pandas as pd


ION_COLUMNS = (
    "Ca_mg_L",
    "Mg_mg_L",
    "Na_mg_L",
    "K_mg_L",
    "HCO3_mg_L",
    "Cl_mg_L",
    "SO4_mg_L",
    "NO3_mg_L",
    "F_mg_L",
    "Sr_mg_L",
    "SiO2_mg_L",
)


@dataclass(frozen=True)
class FieldPrequentialResult:
    predictions: pd.DataFrame
    summary: pd.DataFrame
    audit: Dict[str, object]


def _normalise_date(values: Iterable[object]) -> pd.Series:
    return pd.to_datetime(pd.Series(values), errors="raise").dt.normalize()


def _build_upstream_features(
    wells: pd.DataFrame,
    dry: pd.DataFrame,
) -> pd.DataFrame:
    """Build a transparent head-directed k-nearest graph from March data."""

    merged = wells.merge(
        dry[["Well_ID", *ION_COLUMNS]],
        on="Well_ID",
        how="inner",
        validate="one_to_one",
    ).copy()
    merged["head_proxy_m"] = (
        merged["Elevation_m"] - merged["Static_Water_Level_m"]
    )
    coordinates = merged[["Latitude", "Longitude"]].to_numpy(float)
    heads = merged["head_proxy_m"].to_numpy(float)
    chemistry = np.log1p(merged[list(ION_COLUMNS)].to_numpy(float))
    graph_features = np.zeros_like(chemistry)
    parent_counts = np.zeros(len(merged), dtype=int)
    for target in range(len(merged)):
        delta = coordinates - coordinates[target]
        distances = np.sqrt(np.sum(delta**2, axis=1))
        eligible = np.flatnonzero(
            (heads > heads[target])
            & (np.arange(len(merged)) != target)
        )
        same_aquifer = eligible[
            merged.iloc[eligible]["Aquifer_Type"].to_numpy()
            == merged.iloc[target]["Aquifer_Type"]
        ]
        if same_aquifer.size:
            eligible = same_aquifer
        if eligible.size == 0:
            graph_features[target] = chemistry[target]
            continue
        parents = eligible[np.argsort(distances[eligible])[:3]]
        weights = 1.0 / np.maximum(distances[parents], 1e-6)
        weights /= weights.sum()
        graph_features[target] = weights @ chemistry[parents]
        parent_counts[target] = len(parents)
    frame = pd.DataFrame(
        graph_features,
        columns=[f"upstream_log_{ion}" for ion in ION_COLUMNS],
    )
    frame.insert(0, "Well_ID", merged["Well_ID"].to_numpy())
    frame["graph_parent_count"] = parent_counts
    return frame


def _design_matrix(
    wells: pd.DataFrame,
    dry: pd.DataFrame,
) -> Tuple[pd.DataFrame, np.ndarray]:
    numeric_static = [
        "Elevation_m",
        "Borehole_Depth_m",
        "Static_Water_Level_m",
        "Distance_River_km",
        "Distance_Farm_km",
        "Distance_Settlement_km",
        "Recharge_Zone_Score_0_1",
        "Aridity_Stress_0_1",
        "Fluoride_Geogenic_Risk_0_1",
        "Anthropogenic_Risk_0_1",
    ]
    categories = ["Region", "Aquifer_Type", "Geology_Group", "Land_Use"]
    dry_features = [
        "pH",
        "EC_uS_cm",
        "TDS_mg_L",
        "Temperature_C",
        *ION_COLUMNS,
        "d18O_permil",
        "d2H_permil",
    ]
    dry_for_design = dry.drop(
        columns=[column for column in categories if column in dry.columns]
    )
    merged = dry_for_design.merge(
        wells[["Well_ID", *numeric_static, *categories]],
        on="Well_ID",
        how="inner",
        validate="one_to_one",
    )
    upstream = _build_upstream_features(wells, dry)
    merged = merged.merge(
        upstream,
        on="Well_ID",
        how="inner",
        validate="one_to_one",
    )
    numeric = merged[
        [
            *dry_features,
            *numeric_static,
            *[f"upstream_log_{ion}" for ion in ION_COLUMNS],
            "graph_parent_count",
        ]
    ].astype(float)
    numeric = numeric.fillna(numeric.median())
    for column in dry_features:
        if column in ION_COLUMNS or column in {"EC_uS_cm", "TDS_mg_L"}:
            numeric[column] = np.log1p(np.maximum(numeric[column], 0.0))
    categorical = pd.get_dummies(
        merged[categories].fillna("missing").astype(str),
        prefix=categories,
        dtype=float,
    )
    features = pd.concat(
        [merged[["Well_ID"]].reset_index(drop=True), numeric.reset_index(drop=True), categorical],
        axis=1,
    )
    matrix = features.drop(columns=["Well_ID"]).to_numpy(float)
    means = matrix.mean(axis=0)
    scales = matrix.std(axis=0)
    scales[scales < 1e-10] = 1.0
    matrix = (matrix - means) / scales
    matrix = np.column_stack([np.ones(len(matrix)), matrix])
    return features[["Well_ID"]], matrix


def _ridge_delta(
    x_train: np.ndarray,
    y_delta_train: np.ndarray,
    x_test: np.ndarray,
    alpha: float,
) -> np.ndarray:
    penalty = np.eye(x_train.shape[1]) * float(alpha)
    penalty[0, 0] = 0.0
    coefficients = np.linalg.solve(
        x_train.T @ x_train + penalty,
        x_train.T @ y_delta_train,
    )
    return x_test @ coefficients


def run_prequential_frames(
    wells: pd.DataFrame,
    hydrochemistry: pd.DataFrame,
    *,
    ridge_alpha: float = 10.0,
    minimum_ridge_samples: int = 16,
) -> FieldPrequentialResult:
    """Issue predictions before each wet-season batch is assimilated."""

    required_wells = {
        "Well_ID",
        "Latitude",
        "Longitude",
        "Elevation_m",
        "Static_Water_Level_m",
    }
    required_hydro = {
        "Well_ID",
        "Season",
        "Sampling_Date",
        "Data_Class",
        *ION_COLUMNS,
    }
    missing_wells = required_wells - set(wells.columns)
    missing_hydro = required_hydro - set(hydrochemistry.columns)
    if missing_wells or missing_hydro:
        raise ValueError(
            f"Missing field columns: wells={sorted(missing_wells)}, "
            f"hydrochemistry={sorted(missing_hydro)}"
        )

    hydro = hydrochemistry.copy()
    hydro["Sampling_Date"] = _normalise_date(hydro["Sampling_Date"]).to_numpy()
    dry = hydro[hydro["Season"].astype(str).str.lower() == "dry"].copy()
    wet = hydro[hydro["Season"].astype(str).str.lower() == "wet"].copy()
    quantitative = (
        dry["Data_Class"].astype(str).str.contains("Quantitative", case=False)
        & dry[list(ION_COLUMNS)].notna().all(axis=1)
    )
    dry = dry.loc[quantitative].copy()
    eligible_wells = set(dry["Well_ID"])
    wet = wet[
        wet["Well_ID"].isin(eligible_wells)
        & wet["Data_Class"].astype(str).str.contains("Quantitative", case=False)
        & wet[list(ION_COLUMNS)].notna().all(axis=1)
    ].copy()
    eligible_wells &= set(wet["Well_ID"])
    dry = dry[dry["Well_ID"].isin(eligible_wells)].copy()
    wet = wet[wet["Well_ID"].isin(eligible_wells)].copy()
    if dry["Well_ID"].duplicated().any() or wet["Well_ID"].duplicated().any():
        raise ValueError("Expected exactly one dry and one wet row per well.")
    if not eligible_wells:
        raise ValueError("No complete quantitative dry/wet field pairs.")

    ids, x = _design_matrix(
        wells[wells["Well_ID"].isin(eligible_wells)].copy(),
        dry,
    )
    row_index = {well_id: index for index, well_id in enumerate(ids["Well_ID"])}
    dry_by_id = dry.set_index("Well_ID")
    wet_by_id = wet.set_index("Well_ID")
    dry_log = np.log1p(
        dry_by_id.loc[ids["Well_ID"], list(ION_COLUMNS)].to_numpy(float)
    )
    wet_log = np.log1p(
        wet_by_id.loc[ids["Well_ID"], list(ION_COLUMNS)].to_numpy(float)
    )
    wet_dates = wet_by_id.loc[ids["Well_ID"], "Sampling_Date"].to_numpy()

    prediction_rows = []
    past_prequential_residuals: list[np.ndarray] = []
    unique_dates = sorted(pd.unique(wet_dates))
    for issue_index, issue_date in enumerate(unique_dates):
        test_indices = np.flatnonzero(wet_dates == issue_date)
        train_indices = np.flatnonzero(wet_dates < issue_date)
        if train_indices.size:
            mean_delta = (wet_log[train_indices] - dry_log[train_indices]).mean(
                axis=0
            )
        else:
            mean_delta = np.zeros(len(ION_COLUMNS), dtype=float)

        predictions = {
            "persistence": dry_log[test_indices],
            "expanding_mean_delta": dry_log[test_indices] + mean_delta,
        }
        if train_indices.size >= int(minimum_ridge_samples):
            ridge_delta = _ridge_delta(
                x[train_indices],
                wet_log[train_indices] - dry_log[train_indices],
                x[test_indices],
                ridge_alpha,
            )
            predictions["hydrosheaf_graph_ridge"] = (
                dry_log[test_indices] + ridge_delta
            )
        else:
            predictions["hydrosheaf_graph_ridge"] = (
                dry_log[test_indices] + mean_delta
            )

        if past_prequential_residuals:
            residual_history = np.vstack(past_prequential_residuals)
            half_width = np.maximum(
                0.20,
                np.quantile(np.abs(residual_history), 0.90, axis=0),
            )
        else:
            half_width = np.full(len(ION_COLUMNS), 0.35)

        for local_index, global_index in enumerate(test_indices):
            well_id = str(ids.iloc[global_index]["Well_ID"])
            for method, prediction in predictions.items():
                for ion_index, ion in enumerate(ION_COLUMNS):
                    observed = wet_log[global_index, ion_index]
                    predicted = prediction[local_index, ion_index]
                    prediction_rows.append(
                        {
                            "issue_batch_index": issue_index,
                            "issue_date": pd.Timestamp(issue_date),
                            "well_id": well_id,
                            "ion": ion,
                            "method": method,
                            "n_prior_wet_observations": int(train_indices.size),
                            "prediction_log1p": float(predicted),
                            "observed_log1p": float(observed),
                            "error_log1p": float(predicted - observed),
                            "interval90_low_log1p": float(
                                predicted - half_width[ion_index]
                            ),
                            "interval90_high_log1p": float(
                                predicted + half_width[ion_index]
                            ),
                            "strict_prequential": True,
                        }
                    )
        past_prequential_residuals.extend(
            wet_log[test_indices] - predictions["hydrosheaf_graph_ridge"]
        )

    predictions_frame = pd.DataFrame(prediction_rows)
    predictions_frame["covered90"] = (
        (
            predictions_frame["observed_log1p"]
            >= predictions_frame["interval90_low_log1p"]
        )
        & (
            predictions_frame["observed_log1p"]
            <= predictions_frame["interval90_high_log1p"]
        )
    )
    summary = (
        predictions_frame.groupby(["method", "ion"], as_index=False)
        .agg(
            n=("error_log1p", "size"),
            mae_log1p=("error_log1p", lambda value: float(np.mean(np.abs(value)))),
            rmse_log1p=(
                "error_log1p",
                lambda value: float(np.sqrt(np.mean(np.asarray(value) ** 2))),
            ),
            bias_log1p=("error_log1p", "mean"),
            coverage90=("covered90", "mean"),
        )
    )
    overall = (
        predictions_frame.groupby("method", as_index=False)
        .agg(
            n=("error_log1p", "size"),
            mae_log1p=("error_log1p", lambda value: float(np.mean(np.abs(value)))),
            rmse_log1p=(
                "error_log1p",
                lambda value: float(np.sqrt(np.mean(np.asarray(value) ** 2))),
            ),
            bias_log1p=("error_log1p", "mean"),
            coverage90=("covered90", "mean"),
        )
        .assign(ion="ALL")
    )
    summary = pd.concat([overall, summary], ignore_index=True)

    per_well = (
        predictions_frame.assign(
            absolute_error=lambda frame: frame["error_log1p"].abs()
        )
        .groupby(["well_id", "method"], as_index=False)["absolute_error"]
        .mean()
        .pivot(index="well_id", columns="method", values="absolute_error")
        .dropna()
    )
    rng = np.random.default_rng(2047)
    paired_bootstrap: Dict[str, Dict[str, float]] = {}
    for baseline in ("persistence", "expanding_mean_delta"):
        difference = (
            per_well["hydrosheaf_graph_ridge"] - per_well[baseline]
        ).to_numpy(float)
        draws = np.empty(5000, dtype=float)
        for draw in range(draws.size):
            indices = rng.integers(0, difference.size, difference.size)
            draws[draw] = float(np.mean(difference[indices]))
        paired_bootstrap[f"graph_ridge_minus_{baseline}"] = {
            "mean_paired_mae_difference_log1p": float(np.mean(difference)),
            "ci95_low": float(np.quantile(draws, 0.025)),
            "ci95_high": float(np.quantile(draws, 0.975)),
            "negative_favours_graph_ridge": True,
            "resampling_unit": "well",
        }
    audit: Dict[str, object] = {
        "n_complete_quantitative_pairs": int(len(eligible_wells)),
        "n_issue_batches": int(len(unique_dates)),
        "first_issue_date": str(pd.Timestamp(unique_dates[0]).date()),
        "last_issue_date": str(pd.Timestamp(unique_dates[-1]).date()),
        "future_wet_labels_used_at_issue": False,
        "march_features_known_before_august": True,
        "coordinates_spatially_masked": True,
        "field_claim": (
            "within-campaign one-step hydrochemistry hold-forward; not "
            "independent topology, age, or reaction-truth validation"
        ),
        "ridge_alpha_preregistered": float(ridge_alpha),
        "minimum_ridge_samples": int(minimum_ridge_samples),
        "paired_block_bootstrap": paired_bootstrap,
    }
    return FieldPrequentialResult(predictions_frame, summary, audit)


def run_field_prequential(
    workbook: Path,
) -> FieldPrequentialResult:
    wells = pd.read_excel(workbook, sheet_name="Wells_Nodes")
    hydro = pd.read_excel(workbook, sheet_name="Hydrochemistry_Seasonal")
    return run_prequential_frames(wells, hydro)
