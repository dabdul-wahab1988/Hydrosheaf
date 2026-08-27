"""Decide whether DEM elevations can direct the Lower Anayari lateral edges.

Uses the measured DEM error, not an assumed one. Talensi has genuine per-well
surveyed elevations, so DEM-minus-survey there gives the vertical error budget
for this terrain; that budget is then applied to the Lower Anayari within-village
pairs, whose recorded elevations are village constants and therefore carry no
direction.

The test is whether the DEM elevation DIFFERENCE across a pair is large enough
to be distinguishable from the DEM's own noise -- not merely large enough to
clear the p_ij threshold at the optimistic sigma currently configured.
"""

from __future__ import annotations

from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm

ROOT = Path(__file__).resolve().parents[4]
EDGE_P_MIN = 0.75
GRADIENT_MIN = 1e-4
SIGMA_CONFIGURED = 1.0     # Table S3 elevation-tier sigma


def km(la1, lo1, la2, lo2):
    p = np.pi / 180
    a = (np.sin((la2 - la1) * p / 2) ** 2
         + np.cos(la1 * p) * np.cos(la2 * p) * np.sin((lo2 - lo1) * p / 2) ** 2)
    return 2 * 6371 * np.arcsin(np.sqrt(a))


def main() -> None:
    dem = pd.read_csv(ROOT / "data" / "FieldData" / "derived" / "well_elevations_dem.csv")
    man_raw = pd.read_csv(ROOT / "data" / "FieldData" / "LowerAnayari" / "manu.csv")
    station = {str(r["Sample ID"]): str(r["Station"]) for _, r in man_raw.iterrows()}

    # ---- error budget, measured at Talensi ---------------------------------
    tal = dem[dem.site == "Talensi"]
    err = tal["elev_srtm30m"] - tal["elev_recorded"]
    rmse = float((err ** 2).mean() ** 0.5)
    sigma_pair = rmse * np.sqrt(2.0)
    print("=== measured DEM error budget (Talensi control, n=%d) ===" % len(tal))
    print(f"  SRTM bias {err.mean():+.2f} m, RMSE {rmse:.2f} m")
    print(f"  => per-PAIR elevation-difference noise sigma = {sigma_pair:.2f} m")
    dz_naive = norm.ppf(EDGE_P_MIN) * SIGMA_CONFIGURED * np.sqrt(2)
    dz_honest = norm.ppf(EDGE_P_MIN) * sigma_pair
    print(f"  dz needed for p_ij>=0.75 at the CONFIGURED sigma ({SIGMA_CONFIGURED} m): "
          f"{dz_naive:.2f} m")
    print(f"  dz needed for p_ij>=0.75 at the MEASURED sigma  ({rmse:.2f} m): "
          f"{dz_honest:.2f} m")

    # ---- do the recorded village elevations agree with the DEM? ------------
    m = dem[dem.site == "LowerAnayari"].copy()
    m["station"] = m["site_id"].map(station)
    print("\n=== are the recorded village elevations wrong? ===")
    g = m.groupby("station").agg(n=("site_id", "size"),
                                 recorded=("elev_recorded", "first"),
                                 dem_mean=("elev_srtm30m", "mean"),
                                 dem_sd=("elev_srtm30m", "std"))
    g["offset"] = g["recorded"] - g["dem_mean"]
    print(g.round(1).to_string())
    print(f"  village-mean agreement: median |offset| = {g['offset'].abs().median():.1f} m, "
          f"correlation r = {np.corrcoef(g['recorded'], g['dem_mean'])[0, 1]:.3f}")

    # ---- the 100 within-village pairs --------------------------------------
    fd = pd.read_csv(ROOT / "M2" / "m2_benchmark" / "results" / "field_discovery_results.csv")
    mn = fd[fd["edge_id"].str.startswith("Manu")]
    z = dict(zip(m["site_id"], m["elev_srtm30m"]))
    lat = dict(zip(m["site_id"], m["lat"]))
    lon = dict(zip(m["site_id"], m["lon"]))

    rows = []
    for u, v in zip(mn["u"], mn["v"]):
        u, v = str(u), str(v)
        if station.get(u) != station.get(v):
            continue
        if u not in z or v not in z:
            continue
        d = km(lat[u], lon[u], lat[v], lon[v])
        rows.append({"u": u, "v": v, "d_km": d, "dz": z[u] - z[v]})
    p = pd.DataFrame(rows)
    adz = p["dz"].abs()
    print(f"\n=== the {len(p)} within-village (lateral) pairs, using DEM elevations ===")
    print(f"  |dz| from DEM: median {adz.median():.1f} m, "
          f"IQR {adz.quantile(.25):.1f}-{adz.quantile(.75):.1f} m, max {adz.max():.1f} m")
    print(f"  identical DEM elevation (dz = 0): {(adz == 0).sum()} pairs")

    grad_ok = adz / (p["d_km"] * 1000) >= GRADIENT_MIN
    naive = (adz >= dz_naive) & grad_ok
    honest = (adz >= dz_honest) & grad_ok
    print(f"\n  recovered at the CONFIGURED sigma (1.0 m)  : {int(naive.sum()):>3} of {len(p)}")
    print(f"  recovered at the MEASURED  sigma ({rmse:.1f} m) : {int(honest.sum()):>3} of {len(p)}")

    # signal-to-noise: is the DEM difference bigger than the DEM's own noise?
    snr = adz / sigma_pair
    print(f"\n  signal-to-noise |dz| / {sigma_pair:.1f} m:")
    for lo, hi, lab in [(0, 1, "below noise      (direction meaningless)"),
                        (1, 2, "1-2x noise       (direction unreliable)"),
                        (2, np.inf, ">2x noise        (direction usable)")]:
        k = int(((snr >= lo) & (snr < hi)).sum())
        print(f"    {lab}: {k:>3} pairs ({k / len(p):.0%})")

    print("\n=== verdict ===")
    usable = int((snr >= 2).sum())
    print(f"  Of {len(p)} lateral pairs, {usable} ({usable / len(p):.0%}) have a DEM elevation")
    print(f"  difference exceeding twice the DEM's own pairwise noise.")
    if usable / len(p) < 0.5:
        print("  => DEM elevations cannot direct the majority of these edges.")
        print("     Adopting them would assign directions that are predominantly noise.")
    else:
        print("  => DEM elevations can direct most of these edges.")


if __name__ == "__main__":
    main()
