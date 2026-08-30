# M7.5 prospectively locked robust/hybrid sheaf diagnostic

**Run family:** `RUN-M7-ROBUST-HYBRID-SHEAF-20260729-01`  
**Protocol frozen:** 2026-07-29, before implementation, development fitting, or locked-test generation  
**Status at freeze:** analysis specified; outcomes unknown

## Question

M7.4 showed that the affine sheaf encoded non-identity relations and improved
planted conflict localisation, but did not outperform the stronger edge-local
weighted graph overall. Its principal adverse result was heterogeneous-affine
log loss. M7.5 asks whether that loss is caused by (i) self-influence of a
candidate edge on the global section used to score that edge, (ii) loss of
direct endpoint evidence when a global residual replaces a local residual, or
(iii) calibration fitted to the wrong representation. If none of these
prospectively specified corrections succeeds, the result will support the
bounded conclusion that the generator contains no additional globally useful
information for that stratum.

## Frozen data design

- The independent M7.4 generator is unchanged and remains the sole case
  generator. It does not import HydroSheaf.
- Scenarios remain `identity_limit`, `heterogeneous_affine`,
  `incompatible_cycles`, and `noisy_missing`.
- Development seeds are 8401--8464 (64 cases; 16 per scenario).
- Locked-test seeds are 8501--8628 (128 cases; 32 per scenario).
- Development and test seeds are disjoint from each other and from M7.4.
- Temporal, three-dimensional, and vadose-zone capabilities remain explicitly
  excluded.
- Locked-test truth is unavailable to representation selection, calibration,
  threshold selection, and model fitting; it is joined only for final scoring.

## Frozen estimator arms

All arms receive the same prior logit, one scalar residual feature, and the
local-missing indicator. This equal-dimensional design prevents a hybrid arm
from winning solely because it has more fitted coefficients.

1. `edge_local`: the observed endpoint affine residual; missing values are
   median-imputed within the development fit.
2. `identity_graph`: the original iteratively reweighted global section with
   identity restrictions.
3. `original_affine_global`: the M7.4 iteratively reweighted affine-section
   residual.
4. `original_hybrid`: a development-selected convex blend of log-transformed
   local and original-global residuals; when an endpoint is missing, the
   global residual is used alone.
5. `robust_affine_global`: a leave-one-edge-out (LOO) affine compatibility
   residual. At each of three fixed iterations, each candidate is evaluated
   against a section solved without that candidate, so it cannot make itself
   appear compatible. Weights are updated as
   `prior * exp(-0.5 * (LOO residual / median LOO residual)^2)` and clipped to
   `[1e-4, 1]`.
6. `robust_hybrid`: the same convex blend using the robust LOO residual.
7. `robust_hybrid_permuted`: a negative control using deterministically
   permuted affine maps, with the native robust-hybrid blend and calibrator
   applied unchanged.

The hybrid local weight is selected from `{0, 0.25, 0.50, 0.75, 1.0}` using
development cases only. `1.0` is the local boundary and `0.0` the global
boundary. Ties are resolved toward 0.50, then toward the smaller value.

## Frozen calibration analysis

The primary, deployable regime is arm-specific group-cross-fitted logistic
calibration. Candidate regularisation values are `C in {0.1, 1, 10}`. Eight
deterministic group folds are formed by seed, and the selection objective is
mean per-case log loss. For hybrid arms, `C` and the blend are selected jointly.
The decision threshold is selected from out-of-fold development probabilities
to maximise F1, with ties resolved nearest 0.50. The final calibrator is then
refitted to all development cases.

A secondary shared-calibrator regime stacks the development representation
blocks after arm-specific development standardisation and fits one common
logistic calibrator. It diagnoses whether calibration alone explains a
representation difference; it cannot establish the primary claim.

## Outcomes and uncertainty

Primary outcomes are per-case PR-AUC, Brier score, and log loss. Secondary
outcomes are ROC-AUC, selected F1, expected calibration error, false-confidence
rate, abstention coverage/accuracy, and conflict-localisation PR-AUC. Overall
paired differences are estimated by resampling whole cases with 10,000 fixed
bootstrap replicates and 95% percentile intervals. Scenario-specific and
mechanism contrasts are diagnostic and are not promoted to an unadjusted
family of confirmatory claims.

## Prespecified claim gates

The workflow/provenance gate requires exact file hashes, fresh disjoint seeds,
finite primary metrics, generator independence, immutable frozen settings, and
no truth-bearing inference features.

An **incremental robust-hybrid sheaf claim** is allowed only if, under the
primary separate-calibration regime:

- robust hybrid minus edge local has PR-AUC CI lower bound greater than zero;
- robust hybrid minus edge local has Brier and log-loss CI upper bounds below
  zero (lower scores are better);
- in the identity-limit stratum, PR-AUC degradation is no worse than 0.02 and
  Brier/log-loss degradation is no worse than 0.01/0.02, respectively; and
- native robust hybrid beats the frozen permuted-map control in PR-AUC and
  Brier score with intervals excluding zero in the favourable direction.

If the full gate fails, no overall superiority claim is allowed. A conditional
mechanistic claim may be made only where its paired interval excludes zero and
must name the scenario and comparator. Failure is retained and reported; no
post-test tuning, seed substitution, outcome switching, or silent rerun is
permitted.

## Frozen interpretation rules

- `original_hybrid > original_affine_global` implicates discarded local
  endpoint evidence.
- `robust_affine_global > original_affine_global` implicates candidate
  self-influence/false-edge contamination.
- `robust_hybrid > original_hybrid` implicates the robust global estimator.
- improvement confined to separate rather than shared calibration implicates
  calibration mismatch.
- no material improvement supports absence of additional usable global
  information in this controlled generator, not absence in aquifers generally.

## Reproducibility sequence

1. Freeze this protocol and hashes of the unchanged generator and production
   section solver.
2. Implement and unit-test the runner without locked-test execution.
3. Execute development only; write the selected blends, calibration constants,
   thresholds, and development artifact hashes.
4. Freeze the runner and development settings in a confirmatory lock.
5. Execute the 128 locked-test cases once, score, issue the claim decision, and
   write immutable result hashes.
6. Revise the manuscript from the locked outputs without deleting the M7.4
   adverse result.

## Claim boundary

This is a controlled-synthetic, scalar, static, two-dimensional graph study.
It can discriminate estimator and calibration mechanisms and can support a
bounded capability claim. It is not field validation and cannot establish
general HydroSheaf superiority, temporal performance, three-dimensional
performance, vadose-zone performance, or universal active-learning quality.
