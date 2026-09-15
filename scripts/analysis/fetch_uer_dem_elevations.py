"""Fetch and validate public DEM elevations for the Upper East Region (UER) dataset.

Extracts SRTM 30m and ASTER 30m elevations from OpenTopoData API (api.opentopodata.org)
for all 237 sampling points in Northern Ghana New.
Validates against the surveyed elevations recorded in the native table and outputs
data/FieldData/derived/uer_elevations_dem.csv.
"""

from __future__ import annotations

import json
from pathlib import Path
import time
import urllib.parse
import urllib.request

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
JOIN_CSV = ROOT / "outputs" / "2026-09-05_northern_ghana_geology_join" / "NorthernGhanaNew_geology_join.csv"
OUT_DIR = ROOT / "data" / "FieldData" / "derived"
OUT_CSV = OUT_DIR / "uer_elevations_dem.csv"

API_TEMPLATE = "https://api.opentopodata.org/v1/{ds}"
DATASETS = ("srtm30m", "aster30m")
BATCH_SIZE = 100


def fetch_dem_batch(ds: str, pts: list[tuple[float, float]]) -> list[float | None]:
    """Fetch elevation values for a list of (lat, lon) tuples in batches."""
    results: list[float | None] = []
    for i in range(0, len(pts), BATCH_SIZE):
        chunk = pts[i : i + BATCH_SIZE]
        locs = "|".join(f"{lat:.6f},{lon:.6f}" for lat, lon in chunk)
        url = API_TEMPLATE.format(ds=ds) + "?" + urllib.parse.urlencode({"locations": locs})
        req = urllib.request.Request(
            url,
            headers={"User-Agent": "Hydrosheaf-UER-Elevation/1.0 (academic research)"},
        )
        print(f"  Querying {ds} for points {i + 1}..{i + len(chunk)} of {len(pts)}...")
        with urllib.request.urlopen(req, timeout=60) as resp:
            data = json.loads(resp.read().decode())
        if data.get("status") != "OK":
            raise RuntimeError(f"OpenTopoData error for {ds}: {data.get('status')} {data.get('error')}")
        results.extend(item.get("elevation") for item in data["results"])
        time.sleep(1.1)  # Respect API rate limit (<= 1 req/sec)
    return results


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    if not JOIN_CSV.exists():
        raise FileNotFoundError(f"Geology join CSV not found: {JOIN_CSV}")

    df = pd.read_csv(JOIN_CSV)
    print(f"Loaded {len(df)} records from {JOIN_CSV.name}")

    pts = list(
        zip(
            df["latitude_candidate_dd"].astype(float),
            df["longitude_candidate_dd"].astype(float),
        )
    )

    rec = pd.DataFrame(
        {
            "sample_no": df["sample_no"].astype(int),
            "site_id": [f"NGN_{int(no):03d}" for no in df["sample_no"]],
            "community": df["community"],
            "sample_type": df["sample_type"],
            "latitude_dd": [lat for lat, _ in pts],
            "longitude_dd": [lon for _, lon in pts],
            "elev_surveyed": df["elevation_m"].astype(float),
        }
    )

    for ds in DATASETS:
        print(f"\nFetching {ds} ...")
        rec[f"elev_{ds}"] = fetch_dem_batch(ds, pts)

    # Primary completed elevation:
    # Use surveyed where available, otherwise SRTM 30m
    rec["elevation_completed_m"] = np.where(
        rec["elev_surveyed"].notna(),
        rec["elev_surveyed"],
        rec["elev_srtm30m"],
    )
    rec["elevation_source"] = np.where(
        rec["elev_surveyed"].notna(),
        "surveyed",
        "srtm30m",
    )
    rec["dem_difference_srtm_minus_aster"] = rec["elev_srtm30m"] - rec["elev_aster30m"]

    # Reconcile against surveyed control points
    surveyed_mask = rec["elev_surveyed"].notna()
    n_surveyed = int(surveyed_mask.sum())
    print(f"\n=== Reconciling DEM with {n_surveyed} surveyed points ===")
    for ds in DATASETS:
        diff = rec.loc[surveyed_mask, f"elev_{ds}"] - rec.loc[surveyed_mask, "elev_surveyed"]
        for _, row in rec[surveyed_mask].iterrows():
            print(
                f"  Sample {row['sample_no']:03d} ({row['community']}): "
                f"surveyed = {row['elev_surveyed']:.1f} m, {ds} = {row[f'elev_{ds}']:.1f} m "
                f"(diff = {row[f'elev_{ds}'] - row['elev_surveyed']:+.1f} m)"
            )
        print(f"  {ds} mean bias: {diff.mean():+.2f} m, RMSE: {(diff ** 2).mean() ** 0.5:.2f} m")

    print("\n=== DEM Elevation Summary across all 237 points ===")
    print(
        rec[["elev_srtm30m", "elev_aster30m", "elevation_completed_m", "dem_difference_srtm_minus_aster"]].describe()
    )

    rec.to_csv(OUT_CSV, index=False)
    print(f"\nSuccessfully wrote {OUT_CSV} ({len(rec)} rows)")


if __name__ == "__main__":
    main()
