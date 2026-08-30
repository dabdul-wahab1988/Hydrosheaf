# M8 decisions

## 2026-07-27 — Separate exploration from confirmation

The earlier 16-design transport sweep, 30-replicate OED output, pooled
two-parameter errors, and ad hoc parameter-specific correlations were treated
as exploratory. They motivated the confirmatory design but are not inputs to its
reported estimates.

## 2026-07-27 — Use log-scaled, whitened information

Raw `J.T @ J` was rejected as the confirmatory information matrix because the
two parameters have different units and scales. The confirmatory matrix uses
sensitivities to log10 parameters divided by the declared observation standard
deviations.

## 2026-07-27 — Remove start-centred regularisation

The production adapters construct prior means from their base values. When the
base value is also a random starting value, this changes the prior between
replicates. The primary confirmatory experiment therefore retains the production
forward adapter and PESTGLM optimiser but disables these priors.

## 2026-07-27 — Treat kinetics as a structural negative control

The PHREEQC rate laws contain `k*A` multiplicatively. Residence-time sampling
alone cannot separate the two constants. Kinetics is therefore used to test the
limit of sampling design and the effect of an independent surface-area
constraint, not to claim a second transport-like Pareto frontier.

## 2026-07-27 — Presentation-only amendment after B1 review

The immutable confirmatory run passed, but human review rejected two display
choices. The first figure title implied that target-specific criteria selected
different times even though every criterion selected 50 d. The OED table also
showed the first random draw as if it were the strategy's fixed time. The locked
raw outputs and all numerical estimates were preserved. A checksummed
postprocessor produced reviewed artifacts with an explicit common optimum and a
blank candidate time for the varying random strategy.

## 2026-07-31 — Manuscript draft from locked confirmatory results (r2m Track B)

The M8 manuscript was drafted from the locked confirmatory artifacts only.
The five exploratory pilot claims (C1, C2, C3, C3b, C3c) are excluded from the
manuscript, per M8_CLAIM_DECISION.md; the manuscript presents the common-50-d
optimal-design result, the kinetic k*A rank-one structural limit, the
independent-model parameter-specific robustness failure, the
NOT_ACTIONABLE_FOR_TRANSPORT_TIME_SELECTION portability result, and the
bounded topology active-learning qualification.

## 2026-07-31 — Presentation-only title amendment

The proposal title was changed from the internal project name ("M8 confirmatory
optimal-design and Bayesian active-learning benchmark") to the journal-facing
working title from M8_OUTLINE.md. This is a presentation-only amendment: no
objective, acceptance criterion, input, risk or evidence was changed. The
pre-amendment file is preserved byte-for-byte at
provenance/runs/RUN-M8-CONFIRM-20260728-01/proposal.normalized.locked.json,
mirroring the existing analysis_plan.locked.json precedent. The lock hashes in
m8_confirmatory_protocol.lock.json refer to the pre-amendment file and are
unchanged.

## 2026-07-31 — Independent verification pass on the manuscript draft

A read-only cross-check of every quantitative statement in the draft against
the locked artifacts found and fixed four issues before final assembly:
1. Results 3.1 misattributed the largest dispersivity error to `very_late`
   (0.2404); the true maximum is `c120_s5` (0.3818), with `very_late`
   second-worst (0.2404) and overall worst-conditioned (cond 7074.8).
2. The Results opening undercounted matched-model calibrations: 4,000 fixed-
   design plus 4,500 optimal-design = 8,500 (independent-model 1,000 unchanged).
3. The abstract's "in well-conditioned designs" qualifier on the 0.75-0.91
   coverage range was removed because the upper endpoint comes from the most
   ill-conditioned design (`very_late`).
4. Kinetic parameter correlation sign corrected to -0.957 in Results 3.3.
No other mismatches were found; all medians, bootstrap intervals, ranks,
condition numbers, counts, run IDs, seeds, dates and versions match the locked
artifacts.
