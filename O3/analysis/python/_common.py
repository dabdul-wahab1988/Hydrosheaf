"""Shared paths and helpers for the O3 cross-component synthesis layer.

Computation authority: Python. O3 performs no PHREEQC, MODFLOW/MODPATH or
TracerLPM re-runs (DECISION D2). Every function below reads an already-locked
result file under M3/M4/M5 and either passes a value through unchanged or
applies a disclosed, transparent arithmetic operation (a ratio, a sum, a
group-by) to already-published counts. Where an operation could disagree with
a value the source component already reports as its own headline number, both
values are kept and the discrepancy is written to
`manuscript/artifacts/data/evidence_discrepancies.csv`, mirroring M2.3's D5.
"""

from __future__ import annotations

import math
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "O3" / "manuscript" / "artifacts" / "data"
OUT.mkdir(parents=True, exist_ok=True)

M3 = ROOT / "M3" / "m3_age_benchmark"
M4 = ROOT / "M4" / "m4_topology_benchmark"
M5 = ROOT / "M5" / "m5_inverse_reaction_benchmark"


def wilson_ci(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    """Wilson score 95% CI for a binomial proportion k/n.

    Used only where the underlying quantity is a genuine count out of a known
    denominator (edges, fits, fitted rows); added after review flagged the
    absence of any uncertainty on point-estimate proportions (Methodology
    Issue 3, O3/manuscript/review/O3_manuscript_review.md).
    """
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    denom = 1 + z**2 / n
    centre = p + z**2 / (2 * n)
    adj = z * math.sqrt(p * (1 - p) / n + z**2 / (4 * n**2))
    return (round((centre - adj) / denom, 3), round((centre + adj) / denom, 3))


def write(df: pd.DataFrame, name: str) -> Path:
    path = OUT / name
    df.to_csv(path, index=False)
    print(f"wrote {path.relative_to(ROOT)} ({len(df)} rows)")
    return path
