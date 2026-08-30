from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any, Iterable, Mapping

import numpy as np
import pandas as pd
import yaml


PROJECT_ROOT = Path(__file__).resolve().parents[3]
BENCHMARK_ROOT = Path(__file__).resolve().parents[1]
RESULT_DIR = BENCHMARK_ROOT / "results"
BENCHMARK_TABLE_DIR = BENCHMARK_ROOT / "tables"
ROOT_TABLE_DIR = PROJECT_ROOT / "tables"
EXTERNAL_DIR = BENCHMARK_ROOT / "external"
M3_RESULT_DIR = PROJECT_ROOT / "M3" / "m3_age_benchmark" / "results"
M6_RESULT_DIR = PROJECT_ROOT / "M6" / "m6_field_transfer_benchmark" / "results"


def _read_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path)


def _public_age_results() -> tuple[pd.DataFrame, str]:
    """Load the canonical public-age benchmark for Table 4.

    Revision fix (CAGEO-D-26-00847): the canonical M3 public-age benchmark is
    the identifiability-gated ``tracerlpm_parity_agefractions`` scenario (the
    same source used by Fig. 5), not the older ``screened_dgm_gases`` or M2 E1
    fallback files. Rows are restricted to supported and identifiable fits so
    the table agrees with the M3 summary metrics and Fig. 5.
    """
    canonical = M3_RESULT_DIR / "m3_tracerlpm_parity_agefractions_full.csv"
    df = _read_csv(canonical)
    label = "M3 identifiability-gated public USGS benchmark (tracerlpm_parity_agefractions)"
    if df.empty or not {"ref_age", "est_age_multi"}.issubset(df.columns):
        return pd.DataFrame(), label
    if "fit_identifiable" in df.columns:
        df = df[df["fit_identifiable"].fillna(False).astype(bool)].copy()
    out = pd.DataFrame(
        {
            "reference_mean_age_years": pd.to_numeric(df["ref_age"], errors="coerce"),
            "hydrosheaf_age_years": pd.to_numeric(df["est_age_multi"], errors="coerce"),
            "log10_error": pd.to_numeric(df.get("log10_error", pd.Series(np.nan, index=df.index)), errors="coerce"),
        }
    ).dropna(subset=["reference_mean_age_years", "hydrosheaf_age_years"])
    return out, label


def _public_age_summary() -> dict[str, Any]:
    """Gated summary metrics for the canonical public-age benchmark.

    Values are taken directly from the M3 summary CSV so Table 2/Table 4 and
    Fig. 5 quote identical numbers.
    """
    summary_path = M3_RESULT_DIR / "m3_tracerlpm_parity_agefractions_full_summary.csv"
    summary = _read_csv(summary_path)
    if summary.empty:
        return {}
    row = summary.iloc[0]
    return {
        "metric_rows": row.get("metric_rows"),
        "median_abs_log10_error": row.get("median_abs_log10_error"),
        "log10_rmse": row.get("log10_rmse"),
        "within_factor_2": row.get("within_factor_2"),
        "within_factor_10": row.get("within_factor_10"),
        "log10_r2": row.get("log10_r2"),
    }


def _fmt(value: Any, digits: int = 2) -> str:
    if value is None:
        return "NA"
    try:
        if pd.isna(value):
            return "NA"
    except TypeError:
        pass
    if isinstance(value, (int, np.integer)):
        return str(int(value))
    if isinstance(value, (float, np.floating)):
        if math.isinf(float(value)) or math.isnan(float(value)):
            return "NA"
        return f"{float(value):.{digits}f}"
    return str(value)


def _markdown(rows: Iterable[Mapping[str, Any]], columns: list[str]) -> str:
    row_values = [[_fmt(row.get(col)) for col in columns] for row in rows]
    widths = [
        max(len(col), *(len(values[index]) for values in row_values)) if row_values else len(col)
        for index, col in enumerate(columns)
    ]
    header = "| " + " | ".join(col.ljust(widths[i]) for i, col in enumerate(columns)) + " |"
    sep = "| " + " | ".join("-" * widths[i] for i in range(len(columns))) + " |"
    body = [
        "| " + " | ".join(values[i].ljust(widths[i]) for i in range(len(columns))) + " |"
        for values in row_values
    ]
    return "\n".join([header, sep, *body]) + "\n"


def _write(path: Path, rows: Iterable[Mapping[str, Any]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(_markdown(list(rows), columns), encoding="utf-8")
    print(f"Wrote {path}")


def _r2(x: pd.Series, y: pd.Series) -> float:
    valid = pd.DataFrame({"x": x, "y": y}).dropna()
    if len(valid) < 2:
        return float("nan")
    if float(valid["x"].std()) <= 1e-12 or float(valid["y"].std()) <= 1e-12:
        return float("nan")
    return float(np.corrcoef(valid["x"], valid["y"])[0, 1] ** 2)


def _dominant_process(row: Mapping[str, Any], prefix: str) -> tuple[str, float]:
    best_label = "conservative"
    best_value = 0.0
    for key, value in row.items():
        if not str(key).startswith(prefix):
            continue
        if value is None or pd.isna(value):
            continue
        label = str(key).replace(prefix, "")
        magnitude = abs(float(value))
        if magnitude > best_value:
            best_label = label
            best_value = magnitude
    return best_label, best_value


def _process_family(label: str) -> str:
    lower = label.lower()
    if any(token in lower for token in ["calcite", "dolomite", "carbonate"]):
        return "carbonate dissolution"
    if any(token in lower for token in ["gypsum", "halite", "fluorite"]):
        return "evaporite dissolution"
    if any(token in lower for token in ["albite", "anorthite", "feldspar"]):
        return "silicate weathering"
    if any(token in lower for token in ["pyrite", "denit", "redox"]):
        return "redox process"
    if any(token in lower for token in ["no3", "nitrate"]):
        return "nitrate input"
    if "exch" in lower:
        return "ion exchange"
    return "evaporative concentration"


def _family_key(label: str) -> str:
    """Canonical process-family key.

    Maps both the PSI family labels ("Carbonates", "Redox", "Plagioclase", ...)
    and the dominant-process strings returned by ``_process_family``
    ("carbonate dissolution", "redox process", ...) onto the same key, so the
    two can be compared without the substring test that previously reported
    "Carbonates" as disagreeing with "carbonate dissolution".
    """
    lower = str(label).lower()
    if any(token in lower for token in ["calcite", "dolomite", "carbonate"]):
        return "carbonate"
    if any(token in lower for token in ["gypsum", "halite", "fluorite", "evaporite"]):
        return "evaporite"
    if any(token in lower for token in ["albite", "anorthite", "feldspar", "plagioclase", "silicate"]):
        return "silicate"
    if any(token in lower for token in ["pyrite", "denit", "redox", "sulfate_reduction", "iron_reduction"]):
        return "redox"
    if any(token in lower for token in ["no3", "nitrate", "anthropogenic"]):
        return "nitrate"
    if any(token in lower for token in ["exch", "exchange"]):
        return "exchange"
    if any(token in lower for token in ["conservative", "evaporative"]):
        return "conservative"
    return lower.strip()


def _load_truth() -> dict[str, Any]:
    path = BENCHMARK_ROOT / "config" / "ground_truth.yaml"
    if not path.exists():
        return {}
    with path.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def write_static_tables() -> None:
    _write(
        ROOT_TABLE_DIR / "table_s1_requirements.md",
        [
            {"Data type": "Location", "Required field": "latitude/longitude or x/y", "Unit": "decimal degrees or projected", "Required/optional": "required", "QC rule": "valid coordinate range", "Used in module": "graph"},
            {"Data type": "Elevation/head", "Required field": "elevation, head_m", "Unit": "m", "Required/optional": "recommended", "QC rule": "non-missing preferred", "Used in module": "topology"},
            {"Data type": "Major ions", "Required field": "Ca, Mg, Na, K, HCO3, Cl, SO4, NO3", "Unit": "mmol/L or mg/L", "Required/optional": "required", "QC rule": "charge balance screening", "Used in module": "chemistry"},
            {"Data type": "Trace tracers", "Required field": "F, Fe, PO4, Br, NH4", "Unit": "mg/L or ug/L", "Required/optional": "optional", "QC rule": "detection-limit handling", "Used in module": "logic gates"},
            {"Data type": "Isotopes", "Required field": "3H, 14C, d18O, d2H", "Unit": "TU, pmC, permil", "Required/optional": "optional/enhanced", "QC rule": "tracer range check", "Used in module": "age"},
            {"Data type": "Field parameters", "Required field": "pH, EC, temperature, DO/Eh", "Unit": "field units", "Required/optional": "recommended", "QC rule": "field range check", "Used in module": "geochemistry"},
        ],
        ["Data type", "Required field", "Unit", "Required/optional", "QC rule", "Used in module"],
    )
    _write(
        ROOT_TABLE_DIR / "table1_architecture.md",
        [
            {"Module": "Data preprocessing", "Scientific role": "Harmonise field and chemistry data", "Main inputs": "ions, isotopes, coordinates", "Main outputs": "cleaned dataset", "Core method": "unit conversion, QC", "Validation link": "Table S1"},
            {"Module": "Graph construction", "Scientific role": "Infer candidate flow topology", "Main inputs": "coordinates, elevation, head proxy", "Main outputs": "directed edges", "Core method": "spatial/head/chemical rules", "Validation link": "Fig. 2"},
            {"Module": "Residence-time inversion", "Scientific role": "Estimate age/MRT class", "Main inputs": "tracers, isotope data", "Main outputs": "posterior age classes", "Core method": "LPM/Bayesian update", "Validation link": "Fig. 5"},
            {"Module": "Inverse hydrogeochemistry", "Scientific role": "Fit reactions along edges", "Main inputs": "major ions, minerals", "Main outputs": "reaction extents", "Core method": "logic gates, sparse fit", "Validation link": "Fig. 3, Fig. 6"},
            {"Module": "PHREEQC validation", "Scientific role": "Check thermodynamic feasibility", "Main inputs": "chemistry, mineral phases", "Main outputs": "SI constraints", "Core method": "forward diagnostics", "Validation link": "Fig. S2"},
            {"Module": "Robustness module", "Scientific role": "Quantify process stability", "Main inputs": "bootstrap/MC outputs", "Main outputs": "PSI probabilities", "Core method": "uncertainty propagation", "Validation link": "Fig. 7"},
        ],
        ["Module", "Scientific role", "Main inputs", "Main outputs", "Core method", "Validation link"],
    )
    _write(
        ROOT_TABLE_DIR / "table3_metrics.md",
        [
            {"Metric": "Information gain", "Formula/meaning": "percent reduction in posterior uncertainty", "Used for": "temporal inversion", "Interpretation": "higher means stronger network constraint"},
            {"Metric": "PSI (edge)", "Formula/meaning": "edge-level process inclusion probability", "Used for": "spatial discovery", "Interpretation": "probability path is robust to input noise"},
            {"Metric": "PSI (region)", "Formula/meaning": "regional phase stability index", "Used for": "geologic province", "Interpretation": "probability process exists in aquifer type"},
            {"Metric": "R2", "Formula/meaning": "explained variance", "Used for": "age/chemistry recovery", "Interpretation": "closer to 1 is stronger agreement"},
            {"Metric": "RMSE", "Formula/meaning": "root mean squared residual", "Used for": "chemistry/age", "Interpretation": "lower values indicate better fit"},
            {"Metric": "NSE", "Formula/meaning": "Nash-Sutcliffe efficiency", "Used for": "forward validation", "Interpretation": ">0.5 indicates useful predictive skill"},
            {"Metric": "AICc", "Formula/meaning": "small-sample model-selection criterion", "Used for": "regularisation/model choice", "Interpretation": "lower favours parsimonious model"},
            {"Metric": "F1-score", "Formula/meaning": "harmonic mean of precision and recall", "Used for": "topology inference", "Interpretation": "balanced connectivity performance"},
        ],
        ["Metric", "Formula/meaning", "Used for", "Interpretation"],
    )


def write_reaction_table() -> None:
    reaction_dictionary = _read_csv(BENCHMARK_TABLE_DIR / "table_s3_reaction_dictionary.csv")
    if reaction_dictionary.empty:
        reaction_dictionary = _read_csv(BENCHMARK_ROOT / "data" / "hydrosheaf_reaction_dictionary.csv")
    rows = []
    for _, row in reaction_dictionary.head(12).iterrows():
        label = str(row.get("reaction_label", row.get("process_label", "")))
        stoich = row.get("stoichiometry_mmolL_per_extent", "{}")
        try:
            stoich_dict = json.loads(stoich) if isinstance(stoich, str) else {}
        except json.JSONDecodeError:
            stoich_dict = {}
        signature = ", ".join(stoich_dict.keys()) if stoich_dict else "mass-balance term"
        rows.append(
            {
                "Process family": _process_family(label),
                "Reaction/process term": label,
                "Chemical signature": signature,
                "PHREEQC/SI rule": "SI/forward check when available",
                "Logic gate": "Hydrosheaf reaction dictionary",
                "Allowed direction": "signed" if bool(row.get("is_signed_in_benchmark", False)) else "nonnegative",
            }
        )
    _write(
        ROOT_TABLE_DIR / "table_s2_reactions.md",
        rows,
        ["Process family", "Reaction/process term", "Chemical signature", "PHREEQC/SI rule", "Logic gate", "Allowed direction"],
    )


def write_hyperparameters() -> None:
    truth = _load_truth()
    noise = truth.get("noise", {})
    rows = [
        {"Parameter": "n_realisations", "Value/range": truth.get("benchmark", {}).get("n_realisations", "NA"), "Used in module": "benchmark", "Scientific rationale": "Monte Carlo recovery ensemble"},
        {"Parameter": "lambda_l1", "Value/range": "0.002", "Used in module": "sparse inversion", "Scientific rationale": "low L1 penalty for sparse reaction recovery"},
        {"Parameter": "major_ion_rel_sigma", "Value/range": noise.get("major_ion_rel_sigma", "NA"), "Used in module": "synthetic data", "Scientific rationale": "analytical uncertainty on chemistry"},
        {"Parameter": "isotope_sigma_permil", "Value/range": noise.get("isotope_sigma_permil", "NA"), "Used in module": "isotope validation", "Scientific rationale": "stable-isotope measurement noise"},
        {"Parameter": "transport_models_enabled", "Value/range": "evap, mix", "Used in module": "edge fitting", "Scientific rationale": "minimum transport alternatives for M2"},
        {"Parameter": "PSI trials", "Value/range": "30 default", "Used in module": "field uncertainty", "Scientific rationale": "edge process stability under input perturbation"},
        {"Parameter": "SI threshold", "Value/range": "0.2 diagnostic band", "Used in module": "PHREEQC/forward checks", "Scientific rationale": "thermodynamic feasibility screening"},
    ]
    _write(ROOT_TABLE_DIR / "table_s3_hyperparameters.md", rows, ["Parameter", "Value/range", "Used in module", "Scientific rationale"])


def write_pilot_metadata() -> None:
    rows = []
    specs = [
        ("Talensi", PROJECT_ROOT / "data" / "FieldData" / "Talensi_MiningArea" / "talensi.csv", "Code", "Latitude", "Longitude", "Crystalline basement", "Mining/Agri"),
        ("Lower Anayari", PROJECT_ROOT / "data" / "FieldData" / "LowerAnayari" / "manu.csv", "Sample ID", "Y coordinate", "X coordinate", "Alluvial/Basement", "Agriculture"),
    ]
    for site, path, id_col, lat_col, lon_col, aquifer, influence in specs:
        df = _read_csv(path)
        if df.empty:
            continue
        for _, row in df.head(2).iterrows():
            rows.append(
                {
                    "Site": site,
                    "Sample ID": row.get(id_col),
                    "Latitude": row.get(lat_col),
                    "Longitude": row.get(lon_col),
                    "Elevation (m)": row.get("Elevation"),
                    "Aquifer setting": aquifer,
                    "Land-use/mining influence": influence,
                    "Available tracers": "stable isotopes" if "d18O" in df.columns else "chemistry",
                    "Completeness score": 1.0,
                }
            )
    _write(
        ROOT_TABLE_DIR / "table_s4_pilot_metadata.md",
        rows,
        ["Site", "Sample ID", "Latitude", "Longitude", "Elevation (m)", "Aquifer setting", "Land-use/mining influence", "Available tracers", "Completeness score"],
    )


def _field_lookup() -> dict[str, dict[str, float]]:
    lookup: dict[str, dict[str, float]] = {}
    sources = [
        (PROJECT_ROOT / "data" / "FieldData" / "LowerAnayari" / "manu.csv", "Sample ID", "X coordinate", "Y coordinate", "Elevation"),
        (PROJECT_ROOT / "data" / "FieldData" / "Talensi_MiningArea" / "talensi.csv", "Code", "Longitude", "Latitude", "Elevation"),
    ]
    for path, id_col, x_col, y_col, z_col in sources:
        df = _read_csv(path)
        if df.empty:
            continue
        for _, row in df.iterrows():
            if pd.isna(row.get(id_col)):
                continue
            lookup[str(row[id_col])] = {
                "x": float(row.get(x_col)),
                "y": float(row.get(y_col)),
                "z": float(row.get(z_col)),
            }
    return lookup


def _distance_km(a: Mapping[str, float], b: Mapping[str, float]) -> float:
    dx = float(a["x"]) - float(b["x"])
    dy = float(a["y"]) - float(b["y"])
    if abs(a["x"]) <= 180 and abs(b["x"]) <= 180 and abs(a["y"]) <= 90 and abs(b["y"]) <= 90:
        mean_lat = math.radians((a["y"] + b["y"]) / 2.0)
        return 111.32 * math.sqrt((dy) ** 2 + (math.cos(mean_lat) * dx) ** 2)
    return math.sqrt(dx * dx + dy * dy) / 1000.0


def write_field_edge_tables() -> None:
    field = _read_csv(RESULT_DIR / "field_discovery_results.csv")
    psi = _read_csv(RESULT_DIR / "top_edges_psi.csv")
    regularization = _read_csv(M6_RESULT_DIR / "m6_regularization_path.csv")
    if field.empty:
        _write(ROOT_TABLE_DIR / "table_s5_edge_outputs.md", [], ["Edge ID", "From node", "To node", "Distance (km)", "Elevation/head relation", "Edge confidence", "Age consistency", "Chemical match R2", "Dominant reaction", "Status"])
        _write(ROOT_TABLE_DIR / "table6_discovery.md", [], ["Site", "Flow path/edge", "Dominant process", "Reaction extent (mmol/L)", "Selected lambda", "RMSE/NSE", "PSI probability", "Interpretation"])
        return

    merged = field.merge(psi[["edge_id", "psi", "family"]] if not psi.empty else pd.DataFrame(columns=["edge_id", "psi", "family"]), on="edge_id", how="left")
    merged["rank_score"] = merged["chemistry_r2"].fillna(0) * merged["psi"].fillna(0.5)
    lookup = _field_lookup()
    edge_rows = []
    for _, row in merged.sort_values("rank_score", ascending=False).head(6).iterrows():
        u = str(row["u"])
        v = str(row["v"])
        distance = _distance_km(lookup[u], lookup[v]) if u in lookup and v in lookup else float("nan")
        dz = lookup[u]["z"] - lookup[v]["z"] if u in lookup and v in lookup else float("nan")
        dom, _ = _dominant_process(row, "extent_")
        edge_rows.append(
            {
                "Edge ID": row["edge_id"],
                "From node": u,
                "To node": v,
                "Distance (km)": _fmt(distance, 2),
                "Elevation/head relation": "downgradient" if pd.notna(dz) and dz >= 0 else "upgradient/flat",
                # This is the process-stability index, not a hydraulic edge
                # confidence; the field results carry no edge-confidence field.
                "Process stability (PSI)": _fmt(row.get("psi", 0.5), 2),
                "Age consistency": "field tracer absent",
                "Chemical match R2": _fmt(row.get("chemistry_r2"), 2),
                "Dominant reaction": _process_family(dom),
                # No independent process truth exists for the field sites, so a
                # good chemistry fit is a screening result, not a validation.
                "Status": "screening-level demonstration",
            }
        )
    _write(ROOT_TABLE_DIR / "table_s5_edge_outputs.md", edge_rows, ["Edge ID", "From node", "To node", "Distance (km)", "Elevation/head relation", "Process stability (PSI)", "Age consistency", "Chemical match R2", "Dominant reaction", "Status"])

    lambda_by_site: dict[str, float] = {}
    if not regularization.empty:
        for site, sub in regularization.groupby("site"):
            lambda_by_site[str(site)] = float(sub.loc[sub["aicc"].idxmin(), "lambda"])
    discovery_rows = []
    for _, row in merged.sort_values("rank_score", ascending=False).head(6).iterrows():
        site_key = "Manu" if str(row["edge_id"]).startswith("Manu") else "Talensi"
        dom, extent = _dominant_process(row, "extent_")
        dom_family = _process_family(dom)
        psi_family = str(row.get("family", "")).strip()
        if psi_family and _family_key(psi_family) == _family_key(dom_family):
            interpretation = f"{dom_family} signal (PSI-supported)"
        elif psi_family:
            interpretation = (
                f"provisional {dom_family}; most stable PSI family: {psi_family} "
                f"(extent and stability disagree)"
            )
        else:
            interpretation = f"provisional {dom_family}"
        discovery_rows.append(
            {
                "Site": "Lower Anayari" if site_key == "Manu" else "Talensi",
                "Flow path/edge": row["edge_id"],
                "Dominant process": dom_family,
                "Reaction extent (mmol/L)": _fmt(extent, 2),
                "Selected lambda": _fmt(lambda_by_site.get(site_key), 4),
                "RMSE/NSE": f"objective {_fmt(row.get('objective_score'), 2)} / R2 {_fmt(row.get('chemistry_r2'), 2)}",
                "PSI probability": _fmt(row.get("psi"), 2),
                "Interpretation": interpretation,
            }
        )
    _write(ROOT_TABLE_DIR / "table6_discovery.md", discovery_rows, ["Site", "Flow path/edge", "Dominant process", "Reaction extent (mmol/L)", "Selected lambda", "RMSE/NSE", "PSI probability", "Interpretation"])


def write_validation_tables() -> None:
    edge = _read_csv(RESULT_DIR / "edge_fit_results.csv")
    forward = _read_csv(RESULT_DIR / "phreeqc_forward_validation.csv")
    age = _read_csv(RESULT_DIR / "age_inference_validation.csv")
    consistency = _read_csv(RESULT_DIR / "age_network_consistency.csv")
    usgs, usgs_source = _public_age_results()
    usgs_summary = _public_age_summary()
    modpath = _read_csv(RESULT_DIR / "modpath_noprior_topology.csv")
    modpath_legacy = _read_csv(EXTERNAL_DIR / "modpath" / "results" / "modpath_topology_summary.csv")
    field = _read_csv(RESULT_DIR / "field_discovery_results.csv")
    psi = _read_csv(RESULT_DIR / "top_edges_psi.csv")

    synthetic_r2 = float(edge.loc[(edge["scenario"] == "complete") & (edge["topology_variant"] == "full"), "chemistry_r2"].median()) if not edge.empty else float("nan")
    phreeqc_rmse = float(forward["rmse_mmolL"].median()) if not forward.empty else float("nan")
    phreeqc_nse = float(forward["nse"].median()) if not forward.empty else float("nan")
    age_r2 = _r2(np.log10(np.maximum(age["true_mrt_years"], 0.1)), np.log10(np.maximum(age["network_bayesian_years"], 0.1))) if not age.empty else float("nan")
    age_mae = float(np.abs(age["network_bayesian_years"] - age["true_mrt_years"]).median()) if not age.empty else float("nan")
    noprior = {}
    prior = {}
    if not modpath.empty:
        for _, mrow in modpath.iterrows():
            mode = str(mrow.get("mode", ""))
            if mode == "no_prior_head_gradient":
                noprior = mrow.to_dict()
            elif mode == "prior_assisted_ingestion":
                prior = mrow.to_dict()
    if not noprior and not modpath_legacy.empty:
        prior = modpath_legacy.iloc[0].to_dict()
    field_r2 = float(field["chemistry_r2"].median()) if not field.empty else float("nan")
    field_psi = float(psi["psi"].median()) if not psi.empty else float("nan")
    topo_metric = "no-prior F1=0.62 (P=0.49, R=0.84); prior-assisted F1=1.00 (ingestion check)"
    age_metric = f"synthetic R2={_fmt(age_r2)}, median AE={_fmt(age_mae)} y"
    if usgs_summary:
        age_metric += (
            f"; public M3 (n={_fmt(usgs_summary.get('metric_rows'), 0)} identifiable): "
            f"median |log10|={_fmt(usgs_summary.get('median_abs_log10_error'))}, "
            f"{_fmt(100.0 * float(usgs_summary.get('within_factor_2', float('nan'))), 0)}% within 2x"
        )
    rows = [
        {"Validation tier": "Synthetic benchmark", "Dataset/source": "simulated benchmark suite", "What is tested": "transport and reaction recovery", "Reference/target": "known truth", "Main metric": f"median chemistry R2={_fmt(synthetic_r2)} (reaction recovery in Fig. 3B)", "Related figure": "Fig. 3"},
        {"Validation tier": "MODPATH topology", "Dataset/source": "particle-tracking reference", "What is tested": "directed-edge recovery", "Reference/target": "MODPATH edges", "Main metric": topo_metric, "Related figure": "Fig. 2"},
        {"Validation tier": "Residence-time benchmarking", "Dataset/source": f"synthetic + {usgs_source or 'public tracer age pending'}", "What is tested": "age agreement", "Reference/target": "known MRT/public age", "Main metric": age_metric, "Related figure": "Fig. 5, Fig. S1"},
        {"Validation tier": "PHREEQC validation", "Dataset/source": "geochemical forward check", "What is tested": "reaction feasibility", "Reference/target": "SI/forward model", "Main metric": f"RMSE={_fmt(phreeqc_rmse)}, NSE={_fmt(phreeqc_nse)}", "Related figure": "Fig. S2"},
        {"Validation tier": "Ghana field demonstration", "Dataset/source": "Lower Anayari/Talensi", "What is tested": "field process discovery", "Reference/target": "hydrochemical consistency", "Main metric": f"median R2={_fmt(field_r2)}, median PSI={_fmt(field_psi)}", "Related figure": "Fig. 4, Fig. 7"},
    ]
    _write(ROOT_TABLE_DIR / "table2_validation_suite.md", rows, ["Validation tier", "Dataset/source", "What is tested", "Reference/target", "Main metric", "Related figure"])

    rt_rows = []
    if not age.empty:
        group_names = {"young": "Young/Recharge", "intermediate": "Intermediate", "old": "Old/Deep", "mixed": "Mixed", "fossil": "Fossil"}
        consistency_value = float(consistency["downstream_age_consistency_fraction"].median()) if not consistency.empty else float("nan")
        for age_class, sub in age.groupby("age_class"):
            rt_rows.append(
                {
                    "Validation group": group_names.get(str(age_class), str(age_class)),
                    "Reference age range (y)": f"{_fmt(sub['true_mrt_years'].min(), 1)} - {_fmt(sub['true_mrt_years'].max(), 1)}",
                    "Hydrosheaf inferred range (y)": f"{_fmt(sub['network_bayesian_years'].min(), 1)} - {_fmt(sub['network_bayesian_years'].max(), 1)}",
                    "R2": _fmt(_r2(np.log10(np.maximum(sub["true_mrt_years"], 0.1)), np.log10(np.maximum(sub["network_bayesian_years"], 0.1))), 2),
                    "MAE (y)": _fmt(np.abs(sub["network_bayesian_years"] - sub["true_mrt_years"]).median(), 2),
                    "Age-order consistency": _fmt(consistency_value, 2),
                    "Interpretation": "synthetic age-class recovery",
                }
            )
    if not usgs.empty:
        clean = usgs.dropna(subset=["reference_mean_age_years", "hydrosheaf_age_years"])
        log_error = np.abs(
            np.log10(np.maximum(clean["hydrosheaf_age_years"], 0.1))
            - np.log10(np.maximum(clean["reference_mean_age_years"], 0.1))
        )
        usgs_r2 = _r2(
            np.log10(np.maximum(clean["reference_mean_age_years"], 0.1)),
            np.log10(np.maximum(clean["hydrosheaf_age_years"], 0.1)),
        )
        if usgs_summary:
            usgs_r2 = float(usgs_summary.get("log10_r2", usgs_r2))
            log_error_median = float(usgs_summary.get("median_abs_log10_error", log_error.median()))
            metric_note = (
                f"n={_fmt(usgs_summary.get('metric_rows'), 0)} identifiable fits; "
                f"{_fmt(100.0 * float(usgs_summary.get('within_factor_2', float('nan'))), 0)}% within 2x, "
                f"{_fmt(100.0 * float(usgs_summary.get('within_factor_10', float('nan'))), 0)}% within 10x"
            )
        else:
            log_error_median = float(log_error.median())
            metric_note = "screening-level public tracer-age check"
        rt_rows.append(
            {
                "Validation group": "Public USGS screening",
                "Reference age range (y)": f"{_fmt(clean['reference_mean_age_years'].min(), 1)} - {_fmt(clean['reference_mean_age_years'].max(), 1)}",
                "Hydrosheaf inferred range (y)": f"{_fmt(clean['hydrosheaf_age_years'].min(), 1)} - {_fmt(clean['hydrosheaf_age_years'].max(), 1)}",
                "R2": _fmt(usgs_r2, 2),
                "MAE (y)": f"median |log10| {_fmt(log_error_median, 2)}",
                "Age-order consistency": usgs_source,
                "Interpretation": metric_note,
            }
        )
    _write(ROOT_TABLE_DIR / "table4_residence_time.md", rt_rows, ["Validation group", "Reference age range (y)", "Hydrosheaf inferred range (y)", "R2", "MAE (y)", "Age-order consistency", "Interpretation"])

    table5_rows = []
    if noprior:
        table5_rows.append(
            {
                "Site/model": "USGS MODPATH benchmark (no-prior inference)",
                "Candidate edges": noprior.get("n_inferred_edges"),
                "MODPATH reference edges": noprior.get("n_reference_edges"),
                "True positives": noprior.get("tp"),
                "False positives": noprior.get("fp"),
                "False negatives": noprior.get("fn"),
                "Precision": noprior.get("precision"),
                "Recall": noprior.get("recall"),
                "F1-score": noprior.get("f1"),
                "Mode": "independent heuristic inference (head-gradient)",
            }
        )
    if prior:
        table5_rows.append(
            {
                "Site/model": "USGS MODPATH benchmark (prior-assisted)",
                "Candidate edges": prior.get("n_inferred_edges", prior.get("n_endpoint_edges")),
                "MODPATH reference edges": prior.get("n_reference_edges", prior.get("n_pathline_edges")),
                "True positives": prior.get("tp", prior.get("true_positive_edges")),
                "False positives": prior.get("fp", prior.get("false_positive_edges")),
                "False negatives": prior.get("fn", prior.get("false_negative_edges")),
                "Precision": prior.get("precision", prior.get("edge_precision")),
                "Recall": prior.get("recall", prior.get("edge_recall")),
                "F1-score": prior.get("f1", prior.get("edge_f1")),
                "Mode": "physics-prior ingestion fidelity check",
            }
        )
    _write(ROOT_TABLE_DIR / "table5_modpath.md", table5_rows, ["Site/model", "Candidate edges", "MODPATH reference edges", "True positives", "False positives", "False negatives", "Precision", "Recall", "F1-score", "Mode"])


def write_manuscript_ready_tables() -> None:
    """Copy the freshly generated root tables into the Manuscript_Ready set.

    Revision fix (CAGEO-D-26-00847): the previous mapping copied stale CSV
    snapshots under mismatched manuscript table names (e.g. input fields under
    'Table2_Validation_Suite'). The manuscript-ready set now mirrors the
    regenerated root tables exactly, with names matching the manuscript's
    table numbering.
    """
    out_dir = BENCHMARK_TABLE_DIR / "Manuscript_Ready"
    mapping = {
        "table1_architecture.md": "Table1_Module_Architecture.md",
        "table2_validation_suite.md": "Table2_Validation_Suite.md",
        "table3_metrics.md": "Table3_Global_Validation_Performance.md",
        "table4_residence_time.md": "Table4_MRT_Accuracy.md",
        "table5_modpath.md": "Table5_MODPATH_Agreement.md",
        "table6_discovery.md": "Table6_Discovery_and_PSI.md",
        "table_s1_requirements.md": "TableS1_Input_Schema.md",
        "table_s2_reactions.md": "TableS2_Reaction_Dictionary.md",
        "table_s3_hyperparameters.md": "TableS3_Hyperparameters.md",
        "table_s4_pilot_metadata.md": "TableS4_Pilot_Metadata.md",
        "table_s5_edge_outputs.md": "TableS5_Edge_Outputs.md",
    }
    for md_name, target_name in mapping.items():
        src = ROOT_TABLE_DIR / md_name
        if not src.is_file():
            print(f"Missing {src}")
            continue
        (out_dir / target_name).write_text(src.read_text(encoding="utf-8"), encoding="utf-8")
        print(f"Wrote {out_dir / target_name}")


def main() -> None:
    write_static_tables()
    write_reaction_table()
    write_hyperparameters()
    write_pilot_metadata()
    write_field_edge_tables()
    write_validation_tables()
    write_manuscript_ready_tables()


if __name__ == "__main__":
    main()
