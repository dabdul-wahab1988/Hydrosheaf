"""Harmonised loaders for the three public groundwater-age releases used by M3.

M3 Section 2.2 declares three data sources.  Until now only the first was
actually loaded; this module adds the other two and combines them into the
common analysis table.

Releases
--------
1. **National public-supply** (Jurgens et al., 2022a) - the original source.
   Full coverage: coordinates, depth, atmospheric-equivalent tracers, reported
   LPM family/tracer set, reported total age, age fractions.
2. **Western Principal Aquifers** (Faulkner & Jurgens, 2019).  Tracer values are
   taken from the ``LPM_Meas_Tracer{i}`` / ``LPM_Tracer{i}_Name`` pairs, which
   are already in the units the reported LPM was fitted in (TU, pptv, pmC,
   ccpg), so no solubility conversion is applied.
3. **MRVA alluvial** (Gratzer et al., 2025a).  Well attributes come from the
   ``Table1_Wells`` shapefile DBF, which ScienceBase attaches as a shapefile
   *facet* rather than as a listed file.

Capability flags - not every release supports every benchmark scenario
----------------------------------------------------------------------
The three releases are not interchangeable, and pretending otherwise would
repeat the unequal-support problem Section 2.8 already warns about.  Two boolean
columns record what each row can actually support:

``has_coordinates``
    False for the entire Western release - it publishes no latitude or longitude
    (verified against all five tables and the column-definitions table).  Rows
    without coordinates can never enter graph construction, so the graph
    benchmark support is unchanged by adding Western.

``has_reported_lpm``
    False for MRVA - the release's ``Table8_LPMs.csv`` returns HTTP 404 from
    ScienceBase, so the reported LPM family, tracer set and unsaturated-zone
    travel time are unavailable and strict-TracerLPM-parity scenarios cannot be
    reproduced for those rows.

Deduplication
-------------
332 Western ``SampleID`` values also appear in the national release.  The
national row is kept, because it carries coordinates and the reported-LPM
fields; the duplicate Western row is dropped.
"""
from __future__ import annotations

import re
import struct
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[3]
INPUT = REPO / "M2" / "m2_benchmark" / "external" / "usgs_age" / "input"
WESTERN_DIR = INPUT / "WesternPrincipalAquifers_2004_2018" / "Tables as txt"
MRVA_DIR = INPUT / "MRVA_GroundwaterAge_2018_20"

NODATA = -9999.0


def _num(s):
    return pd.to_numeric(s, errors="coerce")


def _decimal_year(value):
    ts = pd.to_datetime(value, errors="coerce")
    if pd.isna(ts):
        return np.nan
    start = pd.Timestamp(year=ts.year, month=1, day=1)
    end = pd.Timestamp(year=ts.year + 1, month=1, day=1)
    return ts.year + (ts - start).total_seconds() / (end - start).total_seconds()


# ---------------------------------------------------------------------------
# Western Principal Aquifers
# ---------------------------------------------------------------------------

_WESTERN_TRACER_COLUMN = {
    "3H": "tritium_TU",
    "3He(trit)": "he3_trit_TU",
    "SF6": "sf6_pptv",
    "CFC-11": "cfc11_pptv",
    "CFC-12": "cfc12_pptv",
    "CFC-113": "cfc113_pptv",
    "14C": "c14_pmc",
    "4He": "he4_ccpg",
}


def load_western_dataset() -> pd.DataFrame:
    path = WESTERN_DIR / "Table_1_WPAS_AgeInterpretations.txt"
    if not path.exists():
        return pd.DataFrame()
    d = pd.read_csv(path, sep="\t", low_memory=False, encoding="latin-1")

    out = pd.DataFrame({
        "site_id": d["SampleID"].astype(str),
        "StudyUnit": d["StudyUnit"],
        "reference_age_years": _num(d.get("LPM_MeanAgeFinal_yrs")),
        "reported_model_name": d.get("LPM_Name"),
        "reported_uztt_years": _num(d.get("LPM_UZtt_yrs")),
    })
    out["sample_year"] = d["SampleDate"].map(_decimal_year)

    # Tracer values in the units the reported LPM used.
    for col in set(_WESTERN_TRACER_COLUMN.values()):
        out[col] = np.nan
    for i in range(1, 11):
        name_col, val_col, err_col = (
            f"LPM_Tracer{i}_Name", f"LPM_Meas_Tracer{i}", f"LPM_Meas_Tracer{i}_Err")
        if name_col not in d.columns or val_col not in d.columns:
            continue
        values, errs = _num(d[val_col]), _num(d.get(err_col))
        for tracer, target in _WESTERN_TRACER_COLUMN.items():
            mask = (d[name_col] == tracer) & values.notna()
            if not mask.any():
                continue
            out.loc[mask, target] = values[mask]
            sigma_col = (target.replace("_TU", "_sigma_TU")
                               .replace("_pptv", "_sigma_pptv")
                               .replace("_pmc", "_sigma_pmc")
                               .replace("_ccpg", "_sigma_ccpg"))
            if sigma_col not in out.columns:
                out[sigma_col] = np.nan
            out.loc[mask, sigma_col] = errs[mask] if errs is not None else np.nan

    # Perforation depth is reported in feet below land surface.
    bottom_ft = _num(d.get("BottomPerforationsWell_ftbls"))
    out["depth_m"] = bottom_ft * 0.3048

    # This release publishes no well locations.
    out["lat"] = np.nan
    out["lon"] = np.nan
    out["has_coordinates"] = False
    out["has_reported_lpm"] = out["reported_model_name"].notna()
    out["source_release"] = "western_principal_aquifers"
    return out


# ---------------------------------------------------------------------------
# MRVA alluvial aquifer
# ---------------------------------------------------------------------------

def read_dbf(path: Path) -> pd.DataFrame:
    """Minimal dBase-III reader (the MRVA wells layer ships only as a shapefile)."""
    b = path.read_bytes()
    nrec, hlen, rlen = struct.unpack("<IHH", b[4:12])
    fields, off = [], 32
    while b[off] != 0x0D:
        fields.append((b[off:off + 11].split(b"\x00")[0].decode("latin-1"), b[off + 16]))
        off += 32
    rows = []
    for i in range(nrec):
        rec = b[hlen + i * rlen: hlen + (i + 1) * rlen]
        if not rec or rec[:1] == b"*":
            continue
        pos, item = 1, {}
        for name, flen in fields:
            item[name] = rec[pos:pos + flen].decode("latin-1").strip()
            pos += flen
        rows.append(item)
    return pd.DataFrame(rows)


def load_mrva_dataset() -> pd.DataFrame:
    attrs = MRVA_DIR / "Table1_Wells_attributes.csv"
    dbf = MRVA_DIR / "Table1_Wells.dbf"
    if attrs.exists():
        wells = pd.read_csv(attrs)
    elif dbf.exists():
        wells = read_dbf(dbf)
        wells.to_csv(attrs, index=False)
    else:
        return pd.DataFrame()

    ages_path = MRVA_DIR / "Table2_MeanAgeSummary.csv"
    trit_path = MRVA_DIR / "Table3_tritium.csv"
    if not ages_path.exists() or not trit_path.exists():
        return pd.DataFrame()

    for col in ("ALT_VA", "WELL_DEPTH", "DEC_LAT_VA", "DEC_LONG_V"):
        wells[col] = _num(wells[col])
    wells = wells[(wells.WELL_DEPTH > 0) & (wells.ALT_VA > NODATA + 1)]

    ages = pd.read_csv(ages_path, encoding="utf-8-sig")
    ages["avg_mean_age"] = _num(ages["avg_mean_age"])

    trit = pd.read_csv(trit_path, encoding="utf-8-sig")
    trit["Trit_TU"] = _num(trit["Trit_TU"])
    trit["year"] = pd.to_datetime(trit["date"], errors="coerce").dt.year
    trit = (trit.dropna(subset=["Trit_TU"])
                .groupby("siteag", as_index=False)
                .agg(tritium_TU=("Trit_TU", "median"), sample_year=("year", "median")))

    df = (wells.merge(ages[["siteag", "avg_mean_age"]], on="siteag", how="inner")
               .merge(trit, on="siteag", how="inner"))

    out = pd.DataFrame({
        "site_id": df["siteag"].astype(str),
        "lat": df["DEC_LAT_VA"],
        "lon": df["DEC_LONG_V"],
        "depth_m": df["WELL_DEPTH"],
        "tritium_TU": df["tritium_TU"],
        "reference_age_years": df["avg_mean_age"],
        "sample_year": df["sample_year"].fillna(2019.0),
    })
    out["StudyUnit"] = "Mississippi River Valley alluvial aquifer (MRVA)"
    out["has_coordinates"] = True
    # Table8_LPMs.csv is unavailable (HTTP 404 from ScienceBase), so the reported
    # LPM family, tracer set and UZ travel time cannot be reconstructed.
    out["has_reported_lpm"] = False
    out["reported_model_name"] = np.nan
    out["source_release"] = "mrva_alluvial"
    return out


# ---------------------------------------------------------------------------
# Combined
# ---------------------------------------------------------------------------

def load_combined_dataset(national: pd.DataFrame) -> pd.DataFrame:
    """Union the three releases, de-duplicating on site_id.

    `national` is passed in rather than imported to avoid a circular import with
    run_m3_usgs_benchmark.
    """
    nat = national.copy()
    nat["source_release"] = "national_public_supply"
    nat["has_coordinates"] = nat["lat"].notna() & nat["lon"].notna()
    nat["has_reported_lpm"] = True

    frames = [nat]
    known = set(nat["site_id"].astype(str))
    for loader in (load_western_dataset, load_mrva_dataset):
        extra = loader()
        if extra.empty:
            continue
        extra = extra[~extra["site_id"].astype(str).isin(known)]
        extra = extra[_num(extra["reference_age_years"]) > 0]
        known |= set(extra["site_id"].astype(str))
        frames.append(extra)

    combined = pd.concat(frames, ignore_index=True, sort=False)

    # Fields the scenario runner expects to exist for every row.
    if "he4_accumulation_rate_ccpg_per_year" in combined.columns:
        combined["he4_accumulation_rate_ccpg_per_year"] = combined[
            "he4_accumulation_rate_ccpg_per_year"].fillna(1e-11)
    else:
        combined["he4_accumulation_rate_ccpg_per_year"] = 1e-11
    combined["dissolved_gas_correction"] = combined.get(
        "dissolved_gas_correction", pd.Series(index=combined.index, dtype=object)).fillna("dgm_sf6")
    combined["he4_source"] = combined.get(
        "he4_source", pd.Series(index=combined.index, dtype=object)).fillna("calibrated")

    for col in ("tritium_TU", "he3_trit_TU", "sf6_pptv", "cfc11_pptv", "cfc12_pptv", "cfc113_pptv"):
        if col not in combined.columns:
            combined[col] = np.nan
        raw = f"raw_{col}"
        if raw not in combined.columns:
            combined[raw] = np.nan

    raw_gas_columns = [f"raw_{gas}_pptv" for gas in ("sf6", "cfc11", "cfc12", "cfc113")]
    combined["raw_gas_atmospheric_equivalent_available"] = combined[raw_gas_columns].notna().any(axis=1)

    return combined.dropna(subset=["sample_year"]).reset_index(drop=True)


def coverage_report(df: pd.DataFrame) -> pd.DataFrame:
    """Per-release coverage of the capabilities each scenario family requires."""
    rows = []
    for release, g in df.groupby("source_release"):
        rows.append({
            "source_release": release,
            "n_rows": len(g),
            "with_reference_age": int(_num(g["reference_age_years"]).notna().sum()),
            "with_coordinates": int(g["has_coordinates"].fillna(False).sum()),
            "with_reported_lpm": int(g["has_reported_lpm"].fillna(False).sum()),
            "with_tritium": int(_num(g.get("tritium_TU")).notna().sum()),
            "n_study_units": int(g["StudyUnit"].nunique()),
        })
    return pd.DataFrame(rows).sort_values("n_rows", ascending=False)
