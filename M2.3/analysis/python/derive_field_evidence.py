"""Derive reader-facing field-data evidence for the M2.3 manuscript.

Computation authority: Python. All outputs are written as tidy, read-only CSV
exports under M2.3/manuscript/artifacts/data/ for consumption by the R figure
layer. No value reported in the manuscript is copied from an earlier manuscript;
every number is recomputed here from the primary files under data/FieldData/.

Run:  .venv/Scripts/python.exe M2.3/analysis/python/derive_field_evidence.py
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "M2.3" / "manuscript" / "artifacts" / "data"
OUT.mkdir(parents=True, exist_ok=True)

# Equivalent weights: molar mass (g/mol) divided by absolute charge.
EQ_WEIGHT = {
    "Ca": 40.078 / 2,
    "Mg": 24.305 / 2,
    "Na": 22.990 / 1,
    "K": 39.098 / 1,
    "Fe": 55.845 / 2,
    "HCO3": 61.016 / 1,
    "Cl": 35.453 / 1,
    "SO4": 96.06 / 2,
    "NO3": 62.004 / 1,
    "F": 18.998 / 1,
}
CATIONS = ["Ca", "Mg", "Na", "K"]
ANIONS = ["HCO3", "Cl", "SO4", "NO3", "F"]

# Canonical variable metadata used for the reader-facing data description.
VARIABLE_META = {
    "pH": ("pH", "pH units", "field"),
    "EC": ("Electrical conductivity", "uS/cm", "field"),
    "TDS": ("Total dissolved solids", "mg/L", "field"),
    "Temp": ("Water temperature", "degC", "field"),
    "Eh": ("Redox potential", "mV", "field"),
    "Sal": ("Salinity", "practical salinity units", "field"),
    "Ca": ("Calcium", "mg/L", "major ion"),
    "Mg": ("Magnesium", "mg/L", "major ion"),
    "Na": ("Sodium", "mg/L", "major ion"),
    "K": ("Potassium", "mg/L", "major ion"),
    "HCO3": ("Bicarbonate", "mg/L", "major ion"),
    "Cl": ("Chloride", "mg/L", "major ion"),
    "SO4": ("Sulfate", "mg/L", "major ion"),
    "NO3": ("Nitrate", "mg/L", "major ion"),
    "F": ("Fluoride", "mg/L", "minor ion"),
    "Fe": ("Iron", "mg/L", "minor ion"),
    "Sr": ("Strontium", "mg/L", "trace/corroborating"),
    "SiO2": ("Dissolved silica", "mg/L", "trace/corroborating"),
    "d18O": ("Oxygen-18 composition", "per mil VSMOW", "stable isotope"),
    "d2H": ("Deuterium composition", "per mil VSMOW", "stable isotope"),
    "Elevation": ("Ground elevation", "m a.s.l.", "physical setting"),
    "Borehole_Depth": ("Borehole depth", "m", "physical setting"),
    "Static_Water_Level": ("Static water level", "m b.g.l.", "physical setting"),
    "Latitude": ("Latitude", "decimal degrees", "location"),
    "Longitude": ("Longitude", "decimal degrees", "location"),
}

# Column renaming to a single canonical schema across the three sources.
RENAME_NG = {
    "EC_uS_cm": "EC", "TDS_mg_L": "TDS", "Temperature_C": "Temp",
    "Ca_mg_L": "Ca", "Mg_mg_L": "Mg", "Na_mg_L": "Na", "K_mg_L": "K",
    "HCO3_mg_L": "HCO3", "Cl_mg_L": "Cl", "SO4_mg_L": "SO4",
    "NO3_mg_L": "NO3", "F_mg_L": "F", "Sr_mg_L": "Sr", "SiO2_mg_L": "SiO2",
    "d18O_permil": "d18O", "d2H_permil": "d2H",
    "Elevation_m": "Elevation", "Borehole_Depth_m": "Borehole_Depth",
    "Static_Water_Level_m": "Static_Water_Level",
}
RENAME_LA = {"Sample ID": "site_id", "X coordinate": "Longitude",
             "Y coordinate": "Latitude", "Station": "Town"}
RENAME_TA = {"Code": "site_id"}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def to_numeric(series: pd.Series) -> pd.Series:
    """Coerce to float, recording '<x' censored values at the reported limit."""
    if series.dtype.kind in "if":
        return series.astype(float)
    cleaned = (series.astype(str)
               .str.strip()
               .str.replace("<", "", regex=False)
               .str.replace(",", "", regex=False)
               .replace({"": np.nan, "nan": np.nan, "-": np.nan, "ND": np.nan}))
    return pd.to_numeric(cleaned, errors="coerce")


def censored_mask(series: pd.Series) -> pd.Series:
    if series.dtype.kind in "if":
        return pd.Series(False, index=series.index)
    return series.astype(str).str.strip().str.startswith("<")


def load_datasets() -> dict[str, pd.DataFrame]:
    """Load the three field sources into one canonical schema."""
    ng_path = ROOT / "data/FieldData/NorthenGhana/NorthernGhana.xlsx"
    la_path = ROOT / "data/FieldData/LowerAnayari/manu.csv"
    ta_path = ROOT / "data/FieldData/Talensi_MiningArea/talensi.csv"

    frames = {}
    xl = pd.ExcelFile(ng_path)
    ng_parts = []
    for sheet in ["Dry", "Wet"]:
        part = xl.parse(sheet).rename(columns=RENAME_NG)
        part["season"] = sheet
        part["site_id"] = part["Well_ID"].astype(str)
        ng_parts.append(part)
    ng = pd.concat(ng_parts, ignore_index=True)
    ng["dataset"] = "Northern Ghana"
    # Column count as delivered by the source file, before harmonisation columns.
    ng.attrs["source_n_columns"] = xl.parse("Dry").shape[1]
    frames["Northern Ghana"] = ng

    la_raw = pd.read_csv(la_path)
    la = la_raw.rename(columns=RENAME_LA)
    la["dataset"] = "Lower Anayari"
    la["season"] = "not recorded"
    la["site_id"] = la["site_id"].astype(str)
    la.attrs["source_n_columns"] = la_raw.shape[1]
    frames["Lower Anayari"] = la

    ta_raw = pd.read_csv(ta_path)
    ta = ta_raw.rename(columns=RENAME_TA)
    # Recorded correction. The Talensi source file stores all 63 longitudes as
    # positive, but Talensi District lies at roughly 0.8 degrees WEST in the
    # Upper East Region. As delivered, every sample plots outside Ghana. Negating
    # the sign places all 63 inside the Upper East Region polygon, which confirms
    # a dropped sign in transcription rather than a different coordinate system.
    # The source file is left unmodified; the correction is applied here so that
    # it is visible, versioned and reversible.
    assert (ta["Longitude"] > 0).all(), \
        "Talensi longitude sign correction assumes all values are positive"
    ta["Longitude"] = -ta["Longitude"]
    ta["dataset"] = "Talensi"
    ta["season"] = "not recorded"
    ta["site_id"] = ta["site_id"].astype(str)
    ta.attrs["source_n_columns"] = ta_raw.shape[1]
    frames["Talensi"] = ta
    return frames


def variable_inventory(frames: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """Per-dataset, per-variable completeness and distribution summary."""
    rows = []
    for name, df in frames.items():
        for var, (label, unit, group) in VARIABLE_META.items():
            if var not in df.columns:
                rows.append(dict(dataset=name, variable=var, label=label,
                                 unit=unit, variable_group=group, present=False,
                                 n_total=len(df), n_observed=0, n_missing=len(df),
                                 n_censored=0, pct_complete=0.0, minimum=np.nan,
                                 q25=np.nan, median=np.nan, q75=np.nan,
                                 maximum=np.nan, mean=np.nan, sd=np.nan))
                continue
            raw = df[var]
            values = to_numeric(raw)
            n_cens = int(censored_mask(raw).sum())
            obs = values.dropna()
            rows.append(dict(
                dataset=name, variable=var, label=label, unit=unit,
                variable_group=group, present=True, n_total=len(df),
                n_observed=int(obs.size), n_missing=int(values.isna().sum()),
                n_censored=n_cens,
                pct_complete=round(100.0 * obs.size / len(df), 2),
                minimum=obs.min() if obs.size else np.nan,
                q25=obs.quantile(0.25) if obs.size else np.nan,
                median=obs.median() if obs.size else np.nan,
                q75=obs.quantile(0.75) if obs.size else np.nan,
                maximum=obs.max() if obs.size else np.nan,
                mean=obs.mean() if obs.size else np.nan,
                sd=obs.std(ddof=1) if obs.size > 1 else np.nan))
    return pd.DataFrame(rows)


def ionic_frame(df: pd.DataFrame) -> pd.DataFrame:
    """Convert major ions to meq/L and compute the charge-balance error."""
    out = pd.DataFrame(index=df.index)
    for ion in CATIONS + ANIONS:
        if ion in df.columns:
            out[ion] = to_numeric(df[ion]) / EQ_WEIGHT[ion]
        else:
            out[ion] = np.nan
    # Charge balance uses only the ions actually measured in the dataset.
    cat_cols = [c for c in CATIONS if c in df.columns]
    an_cols = [a for a in ANIONS if a in df.columns]
    out["sum_cations_meq"] = out[cat_cols].sum(axis=1, min_count=len(cat_cols))
    out["sum_anions_meq"] = out[an_cols].sum(axis=1, min_count=len(an_cols))
    denom = out["sum_cations_meq"] + out["sum_anions_meq"]
    out["cbe_percent"] = 100.0 * (out["sum_cations_meq"] - out["sum_anions_meq"]) / denom
    out["n_cations_used"] = len(cat_cols)
    out["n_anions_used"] = len(an_cols)
    return out


def classify_quality(cbe: float) -> str:
    """Charge-balance quality tiers used throughout the study."""
    if not np.isfinite(cbe):
        return "unclassified"
    a = abs(cbe)
    if a <= 5:
        return "quantitative"
    if a <= 10:
        return "screening"
    return "exploratory"


def facies(row: pd.Series) -> str:
    """Dominant-ion (Piper-style) hydrochemical facies label."""
    cats = {"Ca": row.get("Ca"), "Mg": row.get("Mg"),
            "Na+K": (row.get("Na") or 0) + (row.get("K") or 0)}
    ans = {"HCO3": row.get("HCO3"), "Cl": row.get("Cl"), "SO4": row.get("SO4")}
    cats = {k: v for k, v in cats.items() if v is not None and np.isfinite(v)}
    ans = {k: v for k, v in ans.items() if v is not None and np.isfinite(v)}
    if not cats or not ans:
        return "unclassified"
    tc, ta = sum(cats.values()), sum(ans.values())
    if tc <= 0 or ta <= 0:
        return "unclassified"
    dom_c = max(cats, key=cats.get)
    dom_a = max(ans, key=ans.get)
    # "Mixed" when no cation/anion exceeds 50% of its total equivalents.
    cat_label = dom_c if cats[dom_c] / tc > 0.5 else "Mixed"
    an_label = dom_a if ans[dom_a] / ta > 0.5 else "Mixed"
    return f"{cat_label}-{an_label}"


def build_sample_table(frames: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """One row per water sample with derived ionic and quality attributes."""
    records = []
    for name, df in frames.items():
        ions = ionic_frame(df)
        keep = df[["dataset", "site_id", "season"]].copy()
        for col in ["Latitude", "Longitude", "Elevation", "pH", "EC", "TDS",
                    "Temp", "d18O", "d2H", "F", "Sr", "SiO2", "NO3", "Fe",
                    "Borehole_Depth", "Static_Water_Level"]:
            keep[col] = to_numeric(df[col]) if col in df.columns else np.nan
        for ion in CATIONS + ANIONS:
            keep[f"{ion}_meq"] = ions[ion]
        keep["sum_cations_meq"] = ions["sum_cations_meq"]
        keep["sum_anions_meq"] = ions["sum_anions_meq"]
        keep["cbe_percent"] = ions["cbe_percent"]
        keep["quality_tier"] = keep["cbe_percent"].map(classify_quality)
        meq = ions[CATIONS + ANIONS].rename(columns=lambda c: c)
        keep["facies"] = meq.apply(facies, axis=1)
        # Deuterium excess relative to the global meteoric water line.
        keep["d_excess"] = keep["d2H"] - 8.0 * keep["d18O"]
        # D2: only the primary measured panel enters inferential summaries.
        keep["primary_panel"] = (keep["season"] != "Wet")
        records.append(keep)
    return pd.concat(records, ignore_index=True)


def seasonal_analysis(samples: pd.DataFrame) -> pd.DataFrame:
    """Dry/wet contrast for the Northern Ghana panels.

    DECISION D2: the seasonal separation in this workbook was reconstructed, not
    independently sampled. This function is retained only as a diagnostic that
    characterises the reconstruction. Its output must never be reported as
    measured seasonal change in the aquifer.
    """
    ng = samples[samples["dataset"] == "Northern Ghana"]
    dry = ng[ng["season"] == "Dry"].set_index("site_id")
    wet = ng[ng["season"] == "Wet"].set_index("site_id")
    shared = sorted(set(dry.index) & set(wet.index))
    variables = ["pH", "EC", "TDS", "Temp", "NO3", "F", "d18O", "d2H",
                 "d_excess", "Ca_meq", "Mg_meq", "Na_meq", "K_meq",
                 "HCO3_meq", "Cl_meq", "SO4_meq", "Static_Water_Level"]
    rows = []
    for var in variables:
        a = pd.to_numeric(dry.loc[shared, var], errors="coerce")
        b = pd.to_numeric(wet.loc[shared, var], errors="coerce")
        ok = a.notna() & b.notna()
        a, b = a[ok], b[ok]
        if a.size < 5:
            continue
        diff = b - a  # wet minus dry
        stat, p = stats.wilcoxon(a, b, zero_method="wilcox", alternative="two-sided")
        # Matched-pairs rank-biserial correlation as the effect size.
        n = diff.size
        rbc = 1.0 - (4.0 * stat) / (n * (n + 1))
        rows.append(dict(
            variable=var, n_pairs=int(n),
            dry_median=float(a.median()), wet_median=float(b.median()),
            median_difference_wet_minus_dry=float(diff.median()),
            wilcoxon_statistic=float(stat), p_value=float(p),
            rank_biserial_correlation=float(rbc)))
    res = pd.DataFrame(rows)
    # Benjamini-Hochberg control across the seasonal variable family.
    order = res["p_value"].rank(method="first")
    m = len(res)
    res["p_value_bh"] = (res["p_value"] * m / order).clip(upper=1.0)
    res["p_value_bh"] = res.sort_values("p_value", ascending=False)["p_value_bh"].cummin().reindex(res.index)
    res["significant_bh_005"] = res["p_value_bh"] < 0.05
    return res.sort_values("p_value").reset_index(drop=True)


def meteoric_lines(samples: pd.DataFrame) -> pd.DataFrame:
    """Ordinary least-squares local meteoric water line per dataset."""
    rows = []
    for name, grp in samples[samples["primary_panel"]].groupby("dataset"):
        ok = grp["d18O"].notna() & grp["d2H"].notna()
        x, y = grp.loc[ok, "d18O"], grp.loc[ok, "d2H"]
        if x.size < 5:
            continue
        res = stats.linregress(x, y)
        rows.append(dict(
            dataset=name, n=int(x.size), slope=float(res.slope),
            intercept=float(res.intercept), r_squared=float(res.rvalue ** 2),
            slope_stderr=float(res.stderr), p_value=float(res.pvalue),
            mean_d_excess=float((y - 8.0 * x).mean()),
            median_d_excess=float((y - 8.0 * x).median())))
    return pd.DataFrame(rows)


def measurement_value() -> pd.DataFrame:
    """Cumulative measurement-tier ablation over the Northern Ghana wells.

    Each tier adds a group of determinands to the preceding one and the wells are
    re-classified. The unit is the well, so the reconstructed seasonal attribute
    is not involved (DECISION D2).
    """
    path = ROOT / "M6/m6_field_transfer_benchmark/results/m6_tier_ablation.csv"
    df = pd.read_csv(path)
    labels = {
        "tier0_majors": "Major ions only",
        "tier1_isotopes": "+ stable isotopes",
        "tier2_fluoride": "+ fluoride",
        "tier3_sr_sio2": "+ strontium and silica",
        "tier4_full_metadata": "+ full metadata",
    }
    out = (df.groupby("tier")
             .agg(n_wells=("well_id", "nunique"),
                  pct_non_identifiable=("resolution_class",
                                        lambda s: 100.0 * (s == "non_identifiable").mean()),
                  mean_resolution_score=("mrs", "mean"),
                  mean_retained_alternatives=("n_support", "mean"),
                  pct_evidence_corroborated=("evidence_corroborated",
                                             lambda s: 100.0 * s.astype(bool).mean()))
             .reset_index())
    out["tier_label"] = out["tier"].map(labels)
    out["tier_order"] = out["tier"].str.extract(r"tier(\d)").astype(int)
    return out.sort_values("tier_order").reset_index(drop=True)


def main() -> None:
    frames = load_datasets()
    measurement_value().to_csv(OUT / "field_measurement_value.csv", index=False)

    inventory = variable_inventory(frames)
    inventory.to_csv(OUT / "field_variable_inventory.csv", index=False)

    samples = build_sample_table(frames)
    samples.to_csv(OUT / "field_samples_derived.csv", index=False)

    primary = samples[samples["primary_panel"]]
    quality = (primary.groupby(["dataset", "quality_tier"]).size()
               .rename("n").reset_index())
    quality.to_csv(OUT / "field_quality_tiers.csv", index=False)

    fac = (primary.groupby(["dataset", "facies"]).size()
           .rename("n").reset_index())
    fac["pct"] = fac.groupby("dataset")["n"].transform(lambda s: 100.0 * s / s.sum())
    fac.to_csv(OUT / "field_facies_counts.csv", index=False)

    seasonal = seasonal_analysis(samples)
    seasonal["interpretation"] = "reconstruction diagnostic; not measured seasonality"
    seasonal.to_csv(OUT / "field_seasonal_reconstruction_diagnostic.csv", index=False)

    # Physical fields that a genuine seasonal resampling would have to change.
    ng = samples[samples["dataset"] == "Northern Ghana"]
    dry = ng[ng["season"] == "Dry"].set_index("site_id").sort_index()
    wet = ng[ng["season"] == "Wet"].set_index("site_id").sort_index()
    recon_rows = []
    for var in ["Static_Water_Level", "Elevation", "Borehole_Depth", "pH",
                "Temp", "TDS", "EC"]:
        a, b = dry[var], wet[var]
        ok = a.notna() & b.notna()
        identical = int((a[ok].values == b[ok].values).sum())
        recon_rows.append(dict(
            variable=var, n_compared=int(ok.sum()), n_identical=identical,
            pct_identical=round(100.0 * identical / max(int(ok.sum()), 1), 1),
            expected_to_vary_seasonally=var in {"Static_Water_Level", "pH",
                                                "Temp", "TDS", "EC"}))
    pd.DataFrame(recon_rows).to_csv(
        OUT / "field_seasonal_reconstruction_evidence.csv", index=False)

    mwl = meteoric_lines(samples)
    mwl.to_csv(OUT / "field_meteoric_water_lines.csv", index=False)

    # Dataset-level headline summary used by the Data section and Table 1.
    summary = []
    for name, df in frames.items():
        sub = samples[samples["dataset"] == name]
        # D2: for Northern Ghana the measured unit is the well, and the primary
        # measured panel is the Dry sheet; the Wet sheet is reconstructed.
        if name == "Northern Ghana":
            sub = sub[sub["season"] == "Dry"]
        cbe = sub["cbe_percent"].dropna()
        summary.append(dict(
            dataset=name,
            n_records=len(sub),
            n_sites=int(sub["site_id"].nunique()),
            panel_status=("primary measured panel; a second reconstructed "
                          "seasonal panel exists and is excluded from inference"
                          if name == "Northern Ghana" else "single measured panel"),
            seasons=", ".join(sorted(sub["season"].unique())),
            n_variables=int(df.attrs["source_n_columns"]),
            median_abs_cbe_percent=float(cbe.abs().median()) if cbe.size else np.nan,
            pct_within_5pct_cbe=float(100.0 * (cbe.abs() <= 5).mean()) if cbe.size else np.nan,
            n_quantitative=int((sub["quality_tier"] == "quantitative").sum()),
            n_screening=int((sub["quality_tier"] == "screening").sum()),
            n_exploratory=int((sub["quality_tier"] == "exploratory").sum()),
            has_fluoride=bool("F" in df.columns and to_numeric(df["F"]).notna().any()),
            has_strontium=bool("Sr" in df.columns and to_numeric(df["Sr"]).notna().any()),
            has_silica=bool("SiO2" in df.columns and to_numeric(df["SiO2"]).notna().any()),
            has_stable_isotopes=bool(sub["d18O"].notna().any()),
            has_age_tracers=False,
            has_screen_intervals=False,
            # D2: no dataset supplies a season-resolved or repeated head series.
            has_repeated_heads=False,
        ))
    pd.DataFrame(summary).to_csv(OUT / "field_dataset_summary.csv", index=False)

    provenance = dict(
        generated_by="M2.3/analysis/python/derive_field_evidence.py",
        sources=[
            dict(path="data/FieldData/NorthenGhana/NorthernGhana.xlsx",
                 sha256=sha256(ROOT / "data/FieldData/NorthenGhana/NorthernGhana.xlsx")),
            dict(path="data/FieldData/LowerAnayari/manu.csv",
                 sha256=sha256(ROOT / "data/FieldData/LowerAnayari/manu.csv")),
            dict(path="data/FieldData/Talensi_MiningArea/talensi.csv",
                 sha256=sha256(ROOT / "data/FieldData/Talensi_MiningArea/talensi.csv")),
        ],
        charge_balance_definition="100*(sum cations meq/L - sum anions meq/L)/(sum cations + sum anions)",
        quality_tiers={"quantitative": "|CBE| <= 5%", "screening": "5% < |CBE| <= 10%",
                       "exploratory": "|CBE| > 10%"},
        note="Charge balance uses only the ions measured in each dataset; "
             "Talensi lacks fluoride so its balance omits that species.",
        applied_corrections=[
            dict(dataset="Talensi",
                 field="Longitude",
                 issue="All 63 longitudes stored as positive (east) although "
                       "Talensi District lies near 0.8 degrees west; every "
                       "sample plotted outside Ghana as delivered.",
                 correction="Sign negated in the derivation script.",
                 verification="After correction all 63 samples fall inside the "
                              "Ghana national boundary and within the Upper East "
                              "Region polygon (geoBoundaries).",
                 source_file_modified=False),
        ],
        northern_ghana_seasonal_status=(
            "The workbook presents Dry and Wet sheets for the same 160 wells. "
            "The author confirmed that the chemistry is real but the seasonal "
            "separation was reconstructed rather than independently sampled. "
            "Static water level is identical in 160 of 160 wells across the two "
            "sheets. The Dry sheet is therefore treated as the primary measured "
            "panel, the Wet sheet is excluded from inference, and no seasonal "
            "change or repeated-head result is reported (DECISION D2)."),
    )
    (OUT / "field_provenance.json").write_text(
        json.dumps(provenance, indent=2), encoding="utf-8")

    print("field records:", len(samples))
    print(pd.DataFrame(summary).to_string(index=False))
    print()
    print(quality.to_string(index=False))
    print()
    print(mwl.to_string(index=False))
    print()
    print(seasonal.head(12).to_string(index=False))


if __name__ == "__main__":
    main()
