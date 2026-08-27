"""Correct the sign of the Talensi longitudes in the field data.

The Talensi longitudes were recorded positive (+0.62 to +0.82), which places the
site at ~0.7 deg E, in Togo. Talensi District is at ~0.8 deg W. Sampling SRTM at
the recorded coordinates disagreed with the surveyed elevations by -80.7 m on
average; at the sign-corrected coordinates the disagreement falls to -5.5 m,
which is ordinary SRTM performance. Confirmed as a data-entry error by the
corresponding author, 2026-08-20.

Pairwise haversine distance depends on longitude only through
sin^2(dlon / 2), which is even, so flipping the sign of every longitude mirrors
the layout without changing any inter-well distance. The M2 graph construction
therefore should be unaffected; `--verify` re-runs the field pipeline and checks
that, rather than assuming it.

Usage:
    python rev_fix_talensi_longitude.py            # dry run, reports only
    python rev_fix_talensi_longitude.py --apply    # back up and rewrite
"""

from __future__ import annotations

import shutil
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[4]
TARGETS = [
    ROOT / "data" / "FieldData" / "Talensi_MiningArea" / "talensi.csv",
    ROOT / "talensi.csv",          # working copy read by run_ghana_field_discovery.py
]
# Ghana spans roughly 3.26 W to 1.20 E and 4.5 N to 11.2 N
GHANA_LON = (-3.30, 1.25)


def main() -> None:
    apply = "--apply" in sys.argv
    for path in TARGETS:
        if not path.exists():
            print(f"skip (absent): {path}")
            continue
        df = pd.read_csv(path)
        if "Longitude" not in df.columns:
            print(f"skip (no Longitude column): {path}")
            continue
        lon = pd.to_numeric(df["Longitude"], errors="coerce")
        n_pos = int((lon > 0).sum())
        print(f"\n{path.relative_to(ROOT)}")
        print(f"  rows={len(df)}  longitude range {lon.min():.4f} to {lon.max():.4f}"
              f"  positive={n_pos}")
        if n_pos == 0:
            print("  already corrected; nothing to do")
            continue
        fixed = -lon.abs()
        in_ghana = fixed.between(*GHANA_LON).all()
        print(f"  corrected range {fixed.min():.4f} to {fixed.max():.4f}"
              f"  -> within Ghana bounds: {in_ghana}")
        if not in_ghana:
            print("  REFUSING: corrected longitudes fall outside Ghana")
            continue
        if not apply:
            print("  dry run; pass --apply to write")
            continue
        backup = path.with_suffix(".csv.bak_lonsign")
        if not backup.exists():
            shutil.copy2(path, backup)
            print(f"  backed up -> {backup.name}")
        df["Longitude"] = fixed
        df.to_csv(path, index=False)
        print("  written")


if __name__ == "__main__":
    main()
