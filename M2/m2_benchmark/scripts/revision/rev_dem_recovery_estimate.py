"""Estimate how many Lower Anayari lateral edges a per-well DEM would recover.

Lower Anayari elevations are recorded per village, not per well: 41 wells carry
6 unique values, exactly constant within each of the 6 stations. Every one of
the 100 within-village edges therefore has dz = 0 and is retained as a lateral
mixing candidate rather than as a directed edge. Sampling a DEM at the surveyed
well coordinates would give per-well elevations. This script estimates, before
committing to that, how many of the 100 laterals would actually be recovered.

An edge is recovered when it clears BOTH graph filters:
  p_ij  = Phi(dz / sigma_dh) >= edge_p_min (0.75), with sigma_dh = sqrt(2) *
          sigma_elev = 1.414 m, hence dz >= 0.954 m
  slope = |dz| / d >= edge_gradient_min (1e-4)
At the observed within-village separations the gradient floor is not binding, so
recovery is governed by the ~0.95 m confidence threshold.

The required dz distribution is estimated three ways, all from data in hand:
  A. Lower Anayari's own between-village relief, giving a local slope that is
     extrapolated down to within-village separations
  B. Talensi's per-well elevations (same region, same survey campaign) as a
     direct empirical analogue of |dz| versus separation
  C. B rescaled by the relief ratio between the two catchments, as a
     conservative bound (Talensi is the hillier site)
"""

from __future__ import annotations

import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm

ROOT = Path(__file__).resolve().parents[4]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

SIGMA_ELEV = 1.0          # Table S3, elevation-tier head sigma (m)
EDGE_P_MIN = 0.75         # Table S3
GRADIENT_MIN = 1e-4       # Table S3
RNG = np.random.default_rng(20260820)
N_DRAWS = 20000


def km(a: np.ndarray, b: np.ndarray) -> float:
    """Great-circle distance in km for (lon, lat) degree pairs."""
    lon1, lat1 = np.radians(a)
    lon2, lat2 = np.radians(b)
    dlat, dlon = lat2 - lat1, lon2 - lon1
    h = np.sin(dlat / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2) ** 2
    return float(2 * 6371.0 * np.arcsin(np.sqrt(h)))


def main() -> None:
    sigma_dh = np.sqrt(2.0) * SIGMA_ELEV
    dz_needed = norm.ppf(EDGE_P_MIN) * sigma_dh
    print(f"sigma_dh = sqrt(2) x {SIGMA_ELEV} = {sigma_dh:.3f} m")
    print(f"dz required for p_ij >= {EDGE_P_MIN}: {dz_needed:.3f} m\n")

    man = pd.read_csv(ROOT / "data" / "FieldData" / "LowerAnayari" / "manu.csv")
    tal = pd.read_csv(ROOT / "data" / "FieldData" / "Talensi_MiningArea" / "talensi.csv")

    coord = {str(r["Sample ID"]): np.array([float(r["X coordinate"]),
                                            float(r["Y coordinate"])])
             for _, r in man.iterrows()}
    station = {str(r["Sample ID"]): str(r["Station"]) for _, r in man.iterrows()}
    velev = {str(r["Station"]): float(r["Elevation"]) for _, r in man.iterrows()}

    # the 100 within-village edges currently retained as laterals
    fd = pd.read_csv(ROOT / "M2" / "m2_benchmark" / "results" / "field_discovery_results.csv")
    mn = fd[fd["edge_id"].str.startswith("Manu")]
    lat_pairs = [(str(u), str(v)) for u, v in zip(mn["u"], mn["v"])
                 if station.get(str(u)) == station.get(str(v))]
    dists = np.array([km(coord[u], coord[v]) for u, v in lat_pairs])
    print(f"lateral (within-village) edges: {len(lat_pairs)}")
    print(f"separation km: median {np.median(dists):.2f}, "
          f"IQR {np.percentile(dists, 25):.2f}-{np.percentile(dists, 75):.2f}, "
          f"max {dists.max():.2f}")
    print(f"gradient floor binds only below dz = {GRADIENT_MIN * np.median(dists) * 1000:.3f} m "
          f"at the median separation -- not binding\n")

    # ---- A. Lower Anayari's own between-village slope -----------------------
    vil = sorted({station[s] for s in station})
    vcoord = {v: np.mean([coord[s] for s in station if station[s] == v], axis=0)
              for v in vil}
    slopes = []
    for a, b in combinations(vil, 2):
        d = km(vcoord[a], vcoord[b])
        if d > 0:
            slopes.append(abs(velev[a] - velev[b]) / d)
    slope_A = float(np.median(slopes))
    print(f"A. Lower Anayari between-village relief: median slope "
          f"{slope_A:.2f} m/km over {len(slopes)} village pairs")

    # ---- B/C. Talensi as an empirical analogue ------------------------------
    tcoord = [(np.array([float(r["Longitude"]), float(r["Latitude"])]),
               float(r["Elevation"]))
              for _, r in tal.iterrows()
              if pd.notna(r["Elevation"]) and pd.notna(r["Longitude"])]
    t_d, t_dz = [], []
    for (ca, za), (cb, zb) in combinations(tcoord, 2):
        d = km(ca, cb)
        if 0 < d <= dists.max() * 1.5:
            t_d.append(d)
            t_dz.append(abs(za - zb))
    t_d, t_dz = np.array(t_d), np.array(t_dz)
    relief_ratio = ((man["Elevation"].max() - man["Elevation"].min()) /
                    (tal["Elevation"].max() - tal["Elevation"].min()))
    print(f"B. Talensi analogue: {len(t_d)} well pairs within the same separation "
          f"range; median |dz| = {np.median(t_dz):.1f} m")
    print(f"C. relief ratio Lower Anayari / Talensi = {relief_ratio:.3f} "
          f"(conservative rescaling)\n")

    def recovered(dz_draws: np.ndarray, d: np.ndarray) -> np.ndarray:
        return (dz_draws >= dz_needed) & (dz_draws / (d * 1000.0) >= GRADIENT_MIN)

    print(f"{'estimate':<52}{'recovered of 100':>18}")
    print("-" * 70)

    # A: local slope x separation, with lognormal scatter about it
    exp_dz = slope_A * dists
    counts = []
    for _ in range(N_DRAWS // 200):
        draw = exp_dz * RNG.lognormal(mean=0.0, sigma=0.75, size=len(dists))
        counts.append(recovered(draw, dists).sum())
    a_lo, a_med, a_hi = np.percentile(counts, [5, 50, 95])
    print(f"{'A. Lower Anayari own between-village slope':<52}"
          f"{f'{a_med:.0f}  (90% CI {a_lo:.0f}-{a_hi:.0f})':>18}")

    # B: resample Talensi |dz| from pairs at comparable separation
    for label, scale in [("B. Talensi empirical analogue (raw)", 1.0),
                         ("C. Talensi analogue, relief-scaled (conservative)",
                          float(relief_ratio))]:
        counts = []
        for _ in range(N_DRAWS // 200):
            draw = np.empty(len(dists))
            for i, d in enumerate(dists):
                near = np.abs(t_d - d) <= max(0.15, 0.25 * d)
                pool = t_dz[near] if near.sum() >= 5 else t_dz
                draw[i] = RNG.choice(pool) * scale
            counts.append(recovered(draw, dists).sum())
        lo, med, hi = np.percentile(counts, [5, 50, 95])
        print(f"{label:<52}{f'{med:.0f}  (90% CI {lo:.0f}-{hi:.0f})':>18}")

    print("\nAll three estimates assume a DEM resolves per-well elevation to better")
    print(f"than the {dz_needed:.2f} m threshold. A 30 m SRTM/Copernicus DEM has a")
    print("vertical RMSE of roughly 2-4 m in low-relief terrain, which is LARGER")
    print("than the threshold, so the recovered count below is an upper bound:")
    print("DEM noise would both create and destroy edges near the cut-off.")

    # DEM vertical error sensitivity, run on the CONSERVATIVE (relief-scaled)
    # distribution -- the optimistic pool would understate the flip rate
    print(f"\nDEM error sensitivity, conservative relief-scaled dz distribution")
    print(f"{'DEM vertical RMSE (m)':<26}{'recovered':>12}{'direction flipped':>20}"
          f"{'|dz| < 2x RMSE':>18}")
    print("-" * 76)
    for rmse in (0.0, 1.0, 2.0, 4.0):
        rec, flip, marginal = [], [], []
        for _ in range(N_DRAWS // 400):
            true_dz = np.empty(len(dists))
            for i, d in enumerate(dists):
                near = np.abs(t_d - d) <= max(0.15, 0.25 * d)
                pool = t_dz[near] if near.sum() >= 5 else t_dz
                true_dz[i] = RNG.choice(pool) * float(relief_ratio)
            noise = RNG.normal(0.0, rmse * np.sqrt(2.0), size=len(dists)) if rmse else 0.0
            obs = true_dz + noise
            keep = recovered(np.abs(obs), dists)
            rec.append(int(keep.sum()))
            # DEM noise reversed the true downhill direction on a kept edge
            flip.append(int(((obs < 0) & keep).sum()))
            # kept, but |dz| is within the DEM's own noise band -> direction
            # is not trustworthy even where the filters pass
            marginal.append(int((keep & (np.abs(obs) < 2 * rmse)).sum()) if rmse else 0)
        print(f"{rmse:<26.1f}{np.median(rec):>12.0f}{np.median(flip):>20.0f}"
              f"{np.median(marginal):>18.0f}")


if __name__ == "__main__":
    main()
