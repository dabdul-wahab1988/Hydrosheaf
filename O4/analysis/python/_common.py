"""Shared paths and helpers for the O4 cross-component synthesis layer.

Computation authority: Python. O4 performs no PHREEQC, MODFLOW/MODPATH,
PEST-GLM or Bayesian active-learning re-run (DECISION D2). Every function
below reads an already-locked result file under M6/M7/M8 and either passes a
value through unchanged or applies a disclosed, transparent arithmetic
operation (a ratio, a re-grouping, a pass-through of an already-computed
paired bootstrap contrast) to already-published values. Where an operation
could disagree with a value the source component already reports as its own
headline number, both values are kept and the discrepancy is written to
`manuscript/artifacts/data/evidence_discrepancies.csv`, mirroring O3's D2/D5
precedent (itself mirroring M2.3's D5).

M7's public-pipeline system-acceptance run (RUN-M7-SYSTEM-20260728-01) is
never read here: it imports a hydrosheaf module (api.py) changed after that
run was locked (DECISIONS.md D3) and is excluded from every O4 claim.
"""

from __future__ import annotations

import math
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "O4" / "manuscript" / "artifacts" / "data"
OUT.mkdir(parents=True, exist_ok=True)

M6 = ROOT / "M6" / "m6_field_transfer_benchmark"
M7 = ROOT / "M7" / "m7_nonuniqueness_benchmark"
M8 = ROOT / "M8" / "m8_calibration_benchmark"


def wilson_ci(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    """Wilson score 95% CI for a binomial proportion k/n.

    Added after internal adversarial review flagged bare-point-estimate
    proportions with a known denominator (manuscript/review/
    O4_manuscript_review.md, Methodology finding), mirroring O3's own D7
    precedent for the identical gap.
    """
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    denom = 1 + z**2 / n
    centre = p + z**2 / (2 * n)
    adj = z * math.sqrt(p * (1 - p) / n + z**2 / (4 * n**2))
    return (round((centre - adj) / denom, 4), round((centre + adj) / denom, 4))


def write(df: pd.DataFrame, name: str) -> Path:
    path = OUT / name
    df.to_csv(path, index=False)
    print(f"wrote {path.relative_to(ROOT)} ({len(df)} rows)")
    return path
