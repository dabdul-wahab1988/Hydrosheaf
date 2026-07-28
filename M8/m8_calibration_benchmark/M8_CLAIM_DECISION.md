# M8 confirmatory claim decision

Run: `RUN-M8-CONFIRM-20260728-01`
Status: PASS
Scope: controlled synthetic evaluation through production HydroSheaf adapters;
not field validation.

## Decision

The proposed headline that scalar identifiability metrics conceal opposing
parameter-specific optimal sampling times is **not confirmed**. After using a
noise-whitened Fisher matrix in log-parameter space, removing random
start-centred priors, and applying the locked decision rule, every
parameter-specific and joint local criterion selected the same 50 d front
observation from the three-point late-time base.

This is not a null result. It gives two supported findings:

1. **A front observation repaired the late-only transport design for both
   parameters.** Under local starts, the common 50 d choice reduced median
   dispersivity absolute log10 error from 0.2154 to 0.0173. The paired median
   difference was -0.1952 (95% bootstrap interval -0.2874 to -0.1656). Decay
   error decreased from 0.0182 to 0.0146, with a paired difference of -0.0040
   (-0.0071 to -0.00003). All calibrations succeeded, and the distant-start
   sensitivity produced essentially the same median errors.
2. **The kinetic adapter demonstrated a structural limit of sampling design.**
   One and four residence-time experiments both produced rank-one information
   for calcite `k` and `A`. Predictions along the constant-`k*A` ridge differed
   by at most 1.33e-6 declared observation standard deviations, whereas doubling
   the product changed predictions by 6.42 standard deviations. An independent
   surface-area observation increased the information rank to two.

The fixed four-point transport sweep still showed parameter-specific placement:
`c90_s45` had the lowest dispersivity error (0.0172), whereas `very_late` had the
lowest decay error (0.0140). That supports reporting per-parameter resolution,
but it does not rescue the stronger claim of opposed conditioning responses or
different optimal next observations. Replicate-wise Spearman correlations with
the corrected condition number were positive for both parameters (median 0.587
for dispersivity and 0.490 for decay); the 95% replicate interval excluded zero
only for dispersivity.

## Manuscript consequence

Withdraw these claims:

- scalar condition, D-, A- and E-optimality were actively misleading in the
  transport confirmation;
- dispersivity- and decay-targeted criteria selected different next samples;
- the earlier negative decay-versus-conditioning correlation was a robust
  physical result.

Retain this bounded claim:

> In a controlled 1-D transport benchmark, a noise-whitened log-parameter
> Fisher analysis selected one front observation that improved recovery of both
> dispersivity and decay from a late-only base. In the PHREEQC kinetic adapter,
> residence-time sampling could not separate a rate constant from reactive
> surface area because the forward law identified only their product; an
> independent surface-area constraint restored local identifiability.

The analytic toy benchmark may remain a positive control, but its 60-fold result
must not be presented as the real-adapter effect.

## Independent-model and active-learning amendment

Run: `RUN-M8-INDEPENDENT-20260728-01`
Status: PASS execution and numerical-convergence gates; the dual-parameter
robustness claim was not supported.

Truth was generated with a separately implemented implicit finite-volume
advection-dispersion-decay solver and calibrated with the production analytical
adapter. The 240-cell truth differed from the 480-cell reference by at most
0.203 declared observation standard deviations, passing the locked 0.25 gate.

The analytical E-optimal criterion and the independent development oracle both
selected 50 d. On 250 locked test replicates, that observation reduced median
dispersivity absolute log10 error from 0.8262 to 0.1674; the paired difference
was -0.6690 (95% bootstrap interval -0.7374 to -0.5757). However, decay error
increased from 0.1367 to 0.1541; the paired difference was +0.0210 (+0.0092 to
+0.0276). Thus the matched-model conclusion that 50 d improved both parameters
did not survive model-form discrepancy. The defensible result is
parameter-specific: the 50 d front observation was robustly useful for
dispersivity and was the equal-weight development oracle, but it degraded decay
recovery under the independent numerical truth.

The existing `rank_next_measurements` routine could not emit a non-tied
transport concentration sampling-time action. It consumes categorical topology
evidence, not candidate-time Jacobians. Its M8 transport-design claim is
therefore **not supported**; report it as
`NOT_ACTIONABLE_FOR_TRANSPORT_TIME_SELECTION` unless a separate Jacobian-aware
interface is implemented and prospectively validated.

## Frontier Bayesian active-learning qualification

Run: `RUN-M8-FRONTIER-AL-20260728-01`

Status: PASS for a bounded controlled-synthetic scientific workflow.

A separate explicit Bayesian value-of-information interface was implemented for
topology-measurement design; it does not retrospectively make the legacy
categorical routine actionable for transport-time OED. Across 24 untouched
independent MODFLOW 6/MODPATH 7 cases, the new policy improved Brier score over
random acquisition by a median -0.02483 (95% paired interval -0.02897 to
-0.01919) and improved joint-entropy reduction per cost by 0.03581 (0.02206 to
0.04218). It was noninferior, not superior, to the strong legacy
uncertainty-chemistry policy under the frozen margins. All seven registered
outputs regenerated byte-identically.

This qualifies the method as a controlled-synthetic scientific workflow with
frontier-level methodological safeguards. It does not establish field
effectiveness or general superiority. All 120 proposed-policy actions were
chemistry panels, so the usefulness of age and connectivity modalities remains
undemonstrated. See `M8_FRONTIER_ACTIVE_LEARNING_CLAIM_DECISION.md` for the full
claim boundary.

## Pilot-only claim disposition

The canonical-design claims previously labelled C1, C2, C3, C3b and C3c were
generated in exploratory pilot/conditioning runs and were not covered by either
locked confirmatory protocol. They are removed from the current contribution
set. Their numerical outputs may remain as hypothesis-generating diagnostics,
but no headline dissociation, coverage, threshold, detector, or manuscript
figure claim should be based on them without a fresh prospective confirmation.

## Updated bounded claim

> In controlled synthetic experiments, a 50 d front observation selected by a
> local Fisher criterion improved both transport parameters in a matched
> analytical benchmark, but under independently generated numerical truth it
> improved dispersivity while degrading decay recovery. In the production
> PHREEQC kinetic adapter, residence-time sampling could not separate a rate
> constant from reactive surface area because the forward law identified only
> their product; independent surface-area information restored local rank.

These results remain model-class-bounded and are not field validation.
