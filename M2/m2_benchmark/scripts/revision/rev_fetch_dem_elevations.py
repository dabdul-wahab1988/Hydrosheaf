"""Sample public DEMs at the surveyed well coordinates for both field sites.

Lower Anayari elevations are recorded per village (41 wells, 6 unique values,
exactly constant within station), so no within-village pair carries a direction.
This script queries per-well elevations by coordinate.

Talensi is fetched as a CONTROL: it has genuine per-well surveyed elevations, so
comparing DEM against survey there measures the DEM's real vertical error in
this exact terrain, which is what sets the error budget for using DEM
elevations at Lower Anayari.

Source: OpenTopoData public API (api.opentopodata.org), datasets srtm30m and
aster30m. Coordinates are sent to that service; no other data leaves the
machine. Results are cached to disk so the API is queried only once.
"""

from __future__ import annotations

import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[4]
OUT = ROOT / "data" / "FieldData" / "derived"
API = "https://api.opentopodata.org/v1/{ds}"
DATASETS = ("srtm30m", "aster30m")
BATCH = 100


def fetch(ds: str, pts: list[tuple[float, float]]) -> list[float | None]:
    out: list[float | None] = []
    for i in range(0, len(pts), BATCH):
        chunk = pts[i:i + BATCH]
        locs = "|".join(f"{la:.6f},{lo:.6f}" for la, lo in chunk)
        url = API.format(ds=ds) + "?" + urllib.parse.urlencode({"locations": locs})
        with urllib.request.urlopen(url, timeout=60) as r:
            data = json.loads(r.read().decode())
        if data.get("status") != "OK":
            raise RuntimeError(f"{ds}: {data.get('status')} {data.get('error')}")
        out.extend(x["elevation"] for x in data["results"])
        time.sleep(1.1)          # public API asks for <= 1 call/second
    return out


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    man = pd.read_csv(ROOT / "data" / "FieldData" / "LowerAnayari" / "manu.csv")
    tal = pd.read_csv(ROOT / "data" / "FieldData" / "Talensi_MiningArea" / "talensi.csv")

    frames = []
    for site, df, idc, latc, lonc in [
            ("LowerAnayari", man, "Sample ID", "Y coordinate", "X coordinate"),
            ("Talensi", tal, "Code", "Latitude", "Longitude")]:
        sub = df.dropna(subset=[latc, lonc]).copy()
        pts = list(zip(sub[latc].astype(float), sub[lonc].astype(float)))
        rec = pd.DataFrame({
            "site": site,
            "site_id": sub[idc].astype(str).to_numpy(),
            "lat": [a for a, _ in pts],
            "lon": [b for _, b in pts],
            "elev_recorded": sub["Elevation"].astype(float).to_numpy(),
        })
        for ds in DATASETS:
            print(f"  {site}: querying {ds} for {len(pts)} wells ...")
            rec[f"elev_{ds}"] = fetch(ds, pts)
        frames.append(rec)

    out = pd.concat(frames, ignore_index=True)
    path = OUT / "well_elevations_dem.csv"
    out.to_csv(path, index=False)
    print(f"\nwrote {path}  ({len(out)} wells)")

    print("\n=== CONTROL: DEM vs surveyed per-well elevation at Talensi ===")
    t = out[out.site == "Talensi"]
    for ds in DATASETS:
        d = t[f"elev_{ds}"] - t["elev_recorded"]
        print(f"  {ds:<10} bias {d.mean():+6.2f} m   RMSE {(d ** 2).mean() ** 0.5:6.2f} m   "
              f"MAE {d.abs().mean():5.2f} m   |d|<2 m: {(d.abs() < 2).mean():.0%}")

    print("\n=== Lower Anayari: recorded (village) vs DEM (per well) ===")
    m = out[out.site == "LowerAnayari"]
    for ds in DATASETS:
        print(f"  {ds:<10} unique values {m[f'elev_{ds}'].nunique():>3} of {len(m)} "
              f"(recorded: {m['elev_recorded'].nunique()})   "
              f"median offset {(m[f'elev_{ds}'] - m['elev_recorded']).median():+.1f} m")


if __name__ == "__main__":
    main()
