# Step 1: can the CFC inconsistency be explained? A negative result

**Run date:** 2026-07-31
**Predecessor:** `m3_identified_ttd_infeasibility_audit_20260731.md` (Finding 4)
**Status:** development diagnostic; **negative result — no explanation identified**

## Question

Finding 4 established that infeasibility in the local TTD baseline is pairwise
(83.2% of infeasible panels have a minimal infeasible subset of size 2) and
concentrated on tracers that share a correction: the three CFCs conflict with
each other at 19–33% while each is individually feasible at every site measured.

That finding attributed the conflict to error in shared per-site correction
parameters — recharge temperature and excess air. This document tests that
attribution. **It does not survive.**

## Test 1: refitting recharge temperature and excess air — premise invalid

The intended test was to scan the (recharge temperature, excess air) grid for
values reconciling the CFCs, and compare against the pipeline's fitted value.

The premise is wrong. Under the frozen protocol's `gas_correction_mode:
usgs_dgm`, `_apply_design_factors` reads `dgm_cfc11_pptv`, `dgm_cfc12_pptv`, and
`dgm_cfc113_pptv` **directly from the USGS release**. HydroSheaf's own
`fit_dissolved_gas_model` is not invoked in this configuration. Of 1,272 sites,
107 have at least two measured CFCs and **none** carries a finite locally-fitted
recharge temperature.

The shared correction is therefore applied upstream by USGS DGMETA, not by
HydroSheaf. It is not HydroSheaf's parameter to refit, and the attribution in
Finding 4 — that these are "shared per-site correction parameters" HydroSheaf
controls — is incorrect as stated. What HydroSheaf does do is consume the
corrected values as exact, attaching only the analytical measurement sigma.

## Test 2: common-mode shared scale — refuted

If the upstream correction carries unmodelled uncertainty, that error is common
to all three CFCs at a site. A shared multiplicative shift is the natural
first-order form. Unlike the k sweep — which widens each tracer's interval
*independently* — a shared scale moves all CFCs together, which is what a
common-mode correction error does.

Scanning a single shared scale factor s ∈ [0.50, 2.00] over 107 sites with ≥2
CFCs, at k = 1.96:

| Outcome | Count | Share |
|---|---:|---:|
| CFCs already feasible under pipeline | 39 / 107 | 36.4% |
| Infeasible, reconciled by some shared s | 19 / 68 | 27.9% |
| Infeasible, **not** reconciled at any s | 49 / 68 | **72.1%** |

Median reconciling scale 0.500 — the grid boundary — with only 2 of 19 within
±20% of unity. The apparent successes are degenerate: shrinking all CFCs toward
zero makes the constraints trivially satisfiable by an old-water TTD. They are
not evidence of a plausible correction error.

**A common-mode multiplicative error does not explain the conflict.**

## Test 3: single culprit tracer — refuted

Among 97 sites measuring all three CFCs, 64 are jointly infeasible. Dropping one
tracer and retesting the surviving pair:

| Dropped | Reconciles remaining pair |
|---|---:|
| CFC11 | 40 / 64 (62.5%) |
| CFC113 | 31 / 64 (48.4%) |
| CFC12 | 29 / 64 (45.3%) |

CFC-11 is modestly implicated, consistent with its known degradation under
anoxic conditions — but the lift over the other two is only about 15 points, and
**37.5% remain infeasible even after removing the single most helpful tracer**.
No single tracer accounts for the conflict.

## Status of candidate explanations

| Explanation | Verdict |
|---|---|
| Multiplicity of independent intervals | Excluded (audit, Finding 2) |
| Uniform interval under-dispersion | Excluded (audit, Finding 2; k=6 insufficient) |
| Observations outside achievable envelope | Excluded (audit, Finding 3; 1.2%) |
| No common TTD / multi-path sampling | Refuted (audit, Finding 4; 83% pairwise) |
| Recharge-temperature / excess-air refit | **Premise invalid** (correction is upstream) |
| Common-mode shared scale | **Refuted** (72% unreconcilable) |
| Single culprit CFC | **Refuted** (37.5% survive best drop) |

Seven candidate explanations have now been tested and none accounts for the
observed infeasibility.

## What this means

Step 1 did not confirm a diagnosis. It narrowed the space and eliminated the
leading candidate, which is a real result, but it did not identify a fix.

Consequently the downstream plan cannot proceed as scheduled:

- **Step 2 (correct the error model and re-run)** has no identified correction
  to apply. Applying one anyway — for instance inflating CFC sigmas until
  infeasibility falls to a comfortable level — would be tuning the uncertainty
  model to remove an inconvenient signal, with no independent justification. It
  must not be done.
- **Step 3 (re-freeze the baseline)** depends on step 2.
- **Step 5 (graph benefit experiment)** depends on step 3, and independently
  requires candidate edges with provenance independent of every evaluated
  tracer. No such edge set exists.

## Recommended next work

The remaining hypotheses require information not currently in the constraint
model:

1. **Redox-stratified analysis.** CFC-11 degradation is anoxic. Test whether
   CFC-11-implicated infeasibility concentrates in reducing conditions, using
   dissolved oxygen, iron, manganese, or sulphide where the release carries
   them. This is the strongest surviving lead and is testable with existing data.
2. **Propagate the upstream correction uncertainty.** The DGMETA correction has
   its own error, currently discarded. Obtain or bound it from the USGS release
   and carry it as a shared nuisance parameter across the CFC family rather than
   as independent per-tracer sigma. Note this is *not* the same as test 2: a
   nuisance parameter admits non-multiplicative, tracer-specific sensitivity to
   the same underlying (T, A) error, which a shared scale cannot represent.
3. **Compare raw against corrected.** Test 1's raw-versus-corrected comparison
   could not run because the raw pathway does not populate the pptv fields
   required. Establishing whether raw CFC values are more mutually consistent
   than DGM-corrected ones would directly indict or exonerate the upstream
   correction. This needs a small amount of plumbing.
4. **Contamination screening.** Urban and industrial recharge elevates specific
   CFCs non-uniformly. Any land-use covariate in the release would test this.

## Claim boundary

This document authorizes no new scientific claim and revises no reported metric.
It records that the correction-provenance explanation offered in the preceding
audit is not supported, and that the cause of the 27.85% local infeasibility
remains **unidentified**.

## Reproduction

Diagnostics are scratch scripts, not committed artifacts, and their numbers
should be treated as provisional until promoted with tests in the manner of
`scripts/run_m3_infeasibility_diagnostics.py`:

- `cfc_reconcile.py` — test 1, drives the real correction pipeline
  (`correct_environmental_tracers` → `build_lpm_tracer_observations` →
  `build_tracer_constraints`)
- `cfc_commonmode.py` — test 2
- `cfc_loo.py` — test 3
