# Redox-based CFC-11 exclusion rule: predeclared protocol and data-availability blocker

**Date:** 2026-07-31
**Predecessors:** `m3_identified_ttd_infeasibility_audit_20260731.md`, `m3_cfc_reconciliation_step1_20260731.md`
**Status:** protocol declared; **test cannot be executed on the current release**

## Purpose

Seven candidate explanations for the 27.85% local infeasibility have been tested
and rejected. The strongest surviving lead is CFC-11 degradation: dropping
CFC-11 reconciles 62.5% of infeasible CFC triples, against 45–48% for CFC-12 and
CFC-113. CFC-11 is known to degrade under anoxic, particularly methanogenic and
sulphate-reducing, conditions, while CFC-12 is comparatively persistent.

If that lead holds, it yields a **mechanism-based exclusion rule** — drop CFC-11
where conditions are reducing — which removes constraint conflicts for a stated
physical reason. This is categorically different from dropping folds because
they were infeasible, which is selection on the outcome and would invalidate the
graph experiment by deleting the very signal the graph audit measures.

This document declares the rule **before** any outcome is observed, so that it
cannot later be tuned to a result.

## Predeclared protocol

### Redox classification

Primary, following the USGS redox framework (McMahon and Chapelle, 2008):

- `DG_Meas_O2_mgL >= 0.5` → **oxic**
- `DG_Meas_O2_mgL < 0.5` → **anoxic**
- `DG_Meas_CH4_mgL` above detection → **methanogenic**, the most reducing class
  and where CFC-11 degradation is most expected

### Hypothesis

CFC-11's involvement in pairwise infeasibility is higher in anoxic than in oxic
samples.

### Primary statistics

1. Rate at which dropping CFC-11 reconciles an infeasible CFC set, stratified by
   redox class.
2. Pairwise infeasibility rate of CFC-11-containing pairs, oxic versus anoxic.

### Decision rule

The exclusion rule is authorized only if **both** hold:

- **(a) Effect.** CFC-11-implicated infeasibility is at least 15 percentage
  points higher in anoxic than oxic samples, with the 95% confidence interval of
  the difference excluding zero.
- **(b) Specificity control.** CFC-12, the persistent member of the family, does
  **not** show the same redox association. If all three CFCs track redox
  equally, the pattern is not selective degradation and the mechanistic
  justification fails.

Control (b) is essential and is the analogue of the reversed and randomized
graph controls: without it, any covariate correlated with sample depth or
aquifer type would appear to justify the exclusion.

### Claim boundary

Passing both gates authorizes only: *"CFC-11 observations under declared
reducing conditions are excluded as unreliable, on the basis of a predeclared
redox criterion and a specificity control."* It does not authorize any claim
about age recovery, and the exclusion must be applied identically in every
downstream comparison including all graph controls.

## Blocker: the required measurements do not co-occur

The test cannot be executed. The redox indicators and the CFC measurements are
in **disjoint subsets** of the release.

Coverage across all 1,272 national sites:

| Indicator | Sites with it | Also ≥2 CFCs | Also all 3 CFCs |
|---|---:|---:|---:|
| `DG_Meas_O2_mgL` (dissolved oxygen) | 286 | **0** | **0** |
| `DG_Meas_CH4_mgL` (methane) | 285 | **0** | **0** |
| `DG_Meas_H2_mgL` (hydrogen) | 0 | 0 | 0 |
| `DG_Meas_N2O_mgL` (nitrous oxide) | 24 | 23 | 21 |

Localising it further: CFC-11 appears in exactly one data release.

| DataRelease 6 (498 sites) | Count |
|---|---:|
| with dissolved oxygen | 46 |
| with methane | 46 |
| with CFC-11 | 54 |
| **with oxygen AND CFC-11** | **0** |
| **with methane AND CFC-11** | **0** |

Within the only release containing CFCs, the dissolved-gas redox suite and the
CFC suite were measured on entirely different wells.

### Why nitrous oxide is not an adequate substitute

N2O is the only overlapping indicator (23 sites), and it is **non-monotonic in
redox**. It accumulates during denitrification under suboxic and
nitrate-reducing conditions, then is itself reduced to N2 under strongly
reducing conditions. Low N2O is therefore ambiguous between fully oxic and
strongly reducing — which is precisely the discrimination the hypothesis
requires, since CFC-11 degradation occurs at the strongly reducing end. With
n=23 and an ambiguous proxy, the test would not be credible at any result, and
running it would risk manufacturing a justification for an exclusion rule.

## Recommended data acquisition

The blocker is data availability, not method. Dissolved oxygen is a routine
field parameter and USGS NWIS records it for most sampled wells; this particular
release simply did not carry it alongside the CFC panel.

1. Retrieve field dissolved oxygen (NWIS parameter code 00300) and, where
   available, dissolved methane, by `site_id` for the 54 DataRelease 6 sites
   carrying CFC-11.
2. Join on site and sample date, with a declared maximum date offset, since
   redox conditions vary in time.
3. Re-run this protocol unchanged. The rule above is frozen as of this document
   and must not be revised after the data arrives.

If NWIS coverage proves insufficient, the honest outcome is that this release
cannot support the exclusion rule, and the CFC conflict remains unexplained.

## Standing note

This is the third consecutive blocker in the M3 graph track that resolves to a
data-availability limit rather than a modelling limit:

1. Candidate graph edges with tracer-independent provenance — do not exist.
2. The upstream CFC gas correction — applied by USGS DGMETA, not refittable
   within HydroSheaf.
3. Redox indicators co-located with CFC measurements — not present.

The implementation is ahead of the evidence the public release can supply. That
should be stated plainly in any write-up rather than presented as work pending.
