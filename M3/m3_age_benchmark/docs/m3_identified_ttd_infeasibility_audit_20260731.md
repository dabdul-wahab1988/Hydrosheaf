# M3 identified-TTD infeasibility audit: why 27.85% of eligible folds admit no feasible TTD

**Run date:** 2026-07-31
**Baseline run:** `m3_identified_ttd_development_full_20260731.csv` (protocol `m3-identified-ttd-development-20260730`)
**Freeze commit:** `2e73d51`
**Status:** development diagnostic; no new scientific claim authorized

## Question

The frozen local baseline reported that 975 of 3,501 eligible site-tracer folds
(27.85%) admitted no non-negative unit-mass TTD satisfying the declared tracer
intervals. Before any graph edge can be audited against that baseline, the cause
of that infeasibility must be identified. An obstruction reported by
`audit_ttd_graph_compatibility` is not attributable to an edge if the local model
already fails at this rate for its own reasons.

## Correction to the baseline count

The QA document reports 975 inconsistent folds. The correct decomposition is:

| Reason | Count |
|---|---:|
| HiGHS status 8 / 15 / 3, model infeasible | 974 |
| HiGHS status 4, solve error | 1 |
| **Total reported as inconsistent** | **975** |

One row is a solver failure, not a demonstrated infeasibility, and should not be
counted as model-data tension.

## Finding 1: infeasibility is not tracer-specific

Marginal rates suggested CFCs were responsible (+42.9pp for CFC11, +42.4pp for
CFC113, +42.2pp for CFC12). That association is entirely confounded with the
number of conditioning tracers. Controlling for conditioning-set size, the CFC
effect disappears:

| Conditioning set size | With a CFC | Without a CFC | Lift |
|---|---:|---:|---:|
| 2 tracers | 14.3% (n=7) | 20.1% (n=1,619) | −5.9pp |
| 3 tracers | 45.6% (n=68) | 45.1% (n=700) | +0.4pp |
| 4–5 tracers | — | no CFC-free sets exist | no contrast |

14C, which an earlier reading of `held_out_tracer` counts appeared to implicate,
is in fact protective (−8.5pp). `held_out_tracer` records the tracer *withheld*,
which takes no part in constructing the feasible set; it must not be used to
attribute infeasibility.

The actual driver is the number of simultaneous constraints:

| Conditioning tracers | Infeasible | Rate |
|---:|---|---:|
| 1 | 9 / 594 | 1.5% |
| 2 | 327 / 1,626 | 20.1% |
| 3 | 347 / 768 | 45.2% |
| 4 | 188 / 303 | 62.0% |
| 5 | 104 / 138 | 75.4% |

## Finding 2: multiplicity and interval width do not explain it

`sigma_multiplier: 1.96` imposes an independent 95% interval per tracer with no
multiplicity adjustment, so some infeasibility is expected purely from geometry.
The protocol was rerun over the full national dataset at k = 2.5, 3.0, 4.0, 5.0,
and 6.0 and compared against the multiplicity null `1 − (2Φ(k) − 1)^n`.

Observed infeasibility % / multiplicity-null %:

| k | Overall | n=2 | n=3 | n=4 | n=5 |
|---|---:|---|---|---|---|
| 1.96 | 27.8% | 20.1 / 9.75 | 45.2 / 14.26 | 62.0 / 18.55 | 75.4 / 22.62 |
| 2.5 | 25.3% | 16.1 / 2.47 | 42.4 / 3.68 | 62.0 / 4.88 | 74.6 / 6.06 |
| 3.0 | 22.1% | 13.4 / 0.54 | 37.1 / 0.81 | 57.4 / 1.08 | 68.8 / 1.34 |
| 4.0 | 18.3% | 9.1 / 0.01 | 30.5 / 0.02 | 55.1 / 0.03 | 65.9 / 0.03 |
| 5.0 | 15.5% | 6.6 / 0.00 | 24.6 / 0.00 | 52.5 / 0.00 | 62.3 / 0.00 |
| 6.0 | 14.1% | 5.7 / 0.00 | 21.4 / 0.00 | 51.2 / 0.00 | 56.5 / 0.00 |

Two hypotheses are excluded:

- **Multiplicity.** At k ≥ 4 the null predicts effectively zero infeasibility,
  yet 18.3% of folds remain infeasible.
- **Uniform interval under-dispersion.** Widening every interval by a factor of
  three (1.96σ → 6σ) reduces infeasibility only from 27.8% to 14.1%, and the
  rise with conditioning-set size persists undiminished (5.7% → 56.5% at k=6).
  Had the sigmas simply been too tight by a constant factor, some k would have
  collapsed the observed curve onto the null at all n. None does.

## Finding 3: the observations are individually reproducible

A third hypothesis is that individual observations lie outside the forward
model's achievable range. For tracer *i*, `A[i,:] @ w` over the simplex is
bounded by `[min_j A[i,j], max_j A[i,j]]`; an observation whose interval misses
that range is unreproducible by any TTD, at any k.

Envelope violations across 3,429 reconstructible folds and 8,052 constraints:

| k | Folds with ≥1 violation | Violating constraints |
|---|---:|---:|
| 1.96 | 158 (4.6%) | 158 / 8,052 (2.0%) |
| 6.0 | 41 (1.2%) | 41 / 8,052 (0.5%) |

By tracer at k=6.0: 3H/3He 31/998 (3.1%), 14C 8/2,209 (0.4%), 3H 2/1,998 (0.1%).
No CFC or SF6 observation violates its envelope at any tested k.

Envelope violation therefore accounts for 1.2% of folds against 14.1% infeasible
at k=6. It is not the cause.

## Finding 4: the inconsistency is pairwise and within tracer families

Findings 1–3 establish only that infeasibility is joint rather than individual.
They do not locate it. An earlier draft of this document concluded from that
elimination that the single-TTD assumption was failing at these sites. **That
conclusion was wrong and is retracted.** The test below refutes it.

For every site with at least three constraints, all singletons and all pairs
were tested for feasibility at k=6.0, and for infeasible full panels the minimal
infeasible subset (MIS) was found by exhaustive search over increasing size.

Minimal infeasible subset size, 226 infeasible panels:

| MIS size | Count | Share |
|---:|---:|---:|
| 1 | 19 | 8.4% |
| 2 | 188 | **83.2%** |
| 3 | 17 | 7.5% |
| 4 | 2 | 0.9% |

Infeasibility is overwhelmingly **pairwise**. Had no common TTD existed because
a sample integrates multiple flow paths, pairs would mostly have been satisfiable
and infeasibility would have emerged only at higher order. It does not. The
size-1 cases are exactly the envelope violations of Finding 3 — the 3H/3He rate
here (15/488 = 3.1%) reproduces the independent envelope diagnostic exactly,
which cross-validates both.

Pairwise infeasibility rates among sites measuring both tracers:

| Pair | Rate | | Pair | Rate |
|---|---:|---|---|---:|
| CFC11 + CFC12 | 32.7% | | 14C + CFC12 | 10.3% |
| CFC11 + CFC113 | 29.5% | | 3H/3He + SF6 | 8.4% |
| 3H/3He + CFC11 | 21.1% | | 3H + SF6 | 6.0% |
| CFC113 + CFC12 | 19.4% | | 14C + 3H/3He | 5.5% |
| 3H + 3H/3He | 17.9% | | 14C + SF6 | 3.4% |
| 3H/3He + CFC113 | 15.0% | | **14C + 3H** | **0.8%** |

The structure is not diffuse. It concentrates on pairs that ought to agree most
strongly:

- **The three CFCs disagree with each other** at 19–33%, the highest rates in
  the panel, while each is individually feasible at every site measured (0/106,
  0/107, 0/98 infeasible alone).
- **3H disagrees with 3H/3He** at 17.9%. These share an element: 3H/3He is
  derived from tritium plus tritiogenic helium.
- **14C and 3H, the most independent pair, almost never conflict** (0.8%).

Tracers constraining *different* parts of the age grid stay consistent. Tracers
constraining the *same* part, and sharing a correction, conflict. That is the
signature of error in the per-site correction parameters applied across a
family — not of a failure of the single-TTD assumption.

## Conclusion

Infeasibility is dominated by **within-family tracer disagreement driven by
shared per-site correction parameters**:

1. **The CFC family.** Recharge temperature and excess-air corrections are
   estimated per site and applied to CFC11, CFC12, and CFC113 together. An error
   propagates to all three, leaving each individually in-range but mutually
   inconsistent. CFC-11 degradation under anoxic conditions affects one member
   only and would produce the same signature.
2. **The tritium family.** Separating tritiogenic from terrigenic helium is the
   prime suspect for the 3H versus 3H/3He conflict.

This is more specific and more tractable than the retracted conclusion. It is a
defect in correction provenance, not in the transit-time model.

Note that Finding 1 (no tracer-specific effect on infeasibility at fixed
conditioning-set size) and Finding 4 (CFC-CFC pairs are the most inconsistent)
are not in tension. Finding 1 concerns the *marginal* effect of a tracer's
presence, which is confounded with set size; Finding 4 conditions on which
specific pair is present.

## Consequences for the graph stage

The graph compatibility audit must not proceed against this baseline.
`audit_ttd_graph_compatibility` returns ABSTAIN when either local set is
infeasible, so a quarter of edge audits would abstain outright. More seriously,
for edges returning INCOMPATIBLE, a graph obstruction cannot be distinguished
from the local inconsistency documented here, which is present at ±6σ and rises
with panel size. Any obstruction reported now would be unattributable.

The remedy is now specific: audit the per-site recharge-temperature, excess-air,
and tritiogenic-helium separation parameters against the CFC-CFC and 3H/3H-3He
conflicts identified above. Because these are correction-provenance defects
rather than transit-time-model defects, they are correctable without abandoning
the single-TTD representation.

## Reproduction

```powershell
# k sweep (protocol variants differ only in sigma_multiplier and protocol_id)
foreach ($k in "2p5","3p0","4p0","5p0","6p0") {
  .\.venv\Scripts\python.exe M3\m3_age_benchmark\scripts\run_m3_identified_ttd_benchmark.py `
    --sources national --withhold-tracer all `
    --protocol M3\m3_age_benchmark\configs\sigma_sweep_20260731\protocol_k$k.yaml `
    --output M3\m3_age_benchmark\results\sigma_sweep_20260731\sweep_k$k.csv
}
```

Artifacts:

- `configs/sigma_sweep_20260731/protocol_k{2p5,3p0,4p0,5p0,6p0}.yaml`
- `results/sigma_sweep_20260731/sweep_k{2p5,3p0,4p0,5p0,6p0}.csv` and manifests

The envelope diagnostic (Finding 3) was run as a scratch script reusing
`build_tracer_constraints`; it is not yet a committed, tested artifact. It must
be promoted to `M3/m3_age_benchmark/scripts/` with tests before any of Finding 3
is cited outside this document.

## Claim boundary

This audit authorizes no new scientific claim. It establishes that the declared
observation-error model and the single-TTD forward model are jointly
inconsistent with the measured multi-tracer panels at a rate that invalidates
downstream graph attribution. It does not identify which of the two candidate
explanations is correct, and it does not revise any previously reported M3
metric.
