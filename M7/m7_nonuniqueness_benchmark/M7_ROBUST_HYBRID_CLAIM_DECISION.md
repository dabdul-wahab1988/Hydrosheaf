# M7.5 robust/hybrid sheaf claim decision

**Run:** `RUN-M7-ROBUST-HYBRID-SHEAF-20260729-01`  
**Decision date:** 2026-07-29  
**Development cases:** 64 (seeds 8401--8464)  
**Locked-test cases:** 128 (seeds 8501--8628)  
**Locked-test executions:** one

## Decision

The execution, provenance, freshness, truth-separation and identity checks
passed. The study therefore has **scientific-workflow status** as a valid,
prospectively locked estimator comparison. The prespecified **incremental
robust-hybrid superiority claim failed** because all three primary outcomes
did not improve simultaneously.

The robust hybrid improved overall PR-AUC over the edge-local weighted graph
by 0.0200 (95% case-block bootstrap CI 0.0073 to 0.0324). Its Brier difference
was -0.00151 (-0.00419 to 0.00105), and its log-loss difference was +0.00333
(-0.00341 to 0.01009); the calibration intervals crossed zero and the mean
log-loss direction was adverse. It therefore cannot be described as generally
superior to the edge-local graph.

## What was learned

Development-only selection set both hybrid local weights to 1.0. Consequently,
the selected hybrids retained the local affine residual wherever both
endpoints were observed and used the global section only as a fallback for
missing endpoints. This rejects a broad local/global blend in the present
generator.

Restoring local evidence was beneficial. The original hybrid relative to the
original global residual changed PR-AUC by +0.00287 (0.00097 to 0.00470),
Brier score by -0.000641 (-0.001071 to -0.000220), and log loss by -0.00151
(-0.00250 to -0.00050). Thus, part of the M7.4 loss was caused by replacing
strong observed-endpoint evidence with a global residual.

The stronger LOO false-edge downweighting did not solve the remaining problem.
Robust global minus original global changed PR-AUC by -0.00480 (-0.00973 to
0.00048), worsened Brier score by +0.00290 (0.00195 to 0.00385), and worsened
log loss by +0.00734 (0.00492 to 0.00973). Robust hybrid did not improve
PR-AUC over original hybrid (-0.00078, -0.00560 to 0.00439) and significantly
worsened Brier score (+0.00167, 0.00080 to 0.00254) and log loss (+0.00426,
0.00210 to 0.00645). Candidate self-influence was therefore not the main
source of the M7.4 calibration loss; the LOO estimator was too destructive in
this dense false-candidate setting.

## Conditional contribution beyond the weighted graph

The selected hybrid improved PR-AUC in the planted incompatible-cycle stratum
by +0.0437 (0.0125 to 0.0751) and in the noisy/missing stratum by +0.0335
(0.0110 to 0.0559). No PR-AUC gain was established in the heterogeneous-affine
stratum (-0.0100, -0.0299 to 0.0121), where log loss remained significantly
worse by +0.0251 (0.0121 to 0.0373). In the noisy/missing stratum, Brier score
also improved by -0.00715 (-0.01269 to -0.00170).

Native robust-hybrid maps decisively beat the frozen permuted-map control:
PR-AUC +0.0441 (0.0321 to 0.0568), Brier -0.01230 (-0.01491 to -0.00976), and
log loss -0.03114 (-0.03764 to -0.02497). The map-edge correspondence is
therefore informative rather than decorative.

The secondary shared-calibrator regime increased the robust-hybrid PR-AUC
contrast against the edge-local graph to +0.0564 (0.0441 to 0.0690) and
improved Brier score by -0.00353 (-0.00643 to -0.00064), while the log-loss
interval still crossed zero (-0.00281, -0.00922 to 0.00381). Calibration
choice affects the magnitude of the apparent gain but does not rescue the
strict three-outcome gate.

## Allowed manuscript claim

In this controlled-synthetic scalar benchmark, affine restriction maps carry
edge-specific information beyond identity connectivity, and a hybrid that
retains direct local residuals while using global compatibility when endpoint
evidence is missing improves ranking in incompatible-cycle and noisy/missing
settings. The additional LOO robustification did not improve the hybrid and
worsened calibration. No claim of general predictive superiority over a
strong edge-local weighted graph is allowed.

## Prohibited claims

- general or universal sheaf superiority;
- field validation of topology, age or reactions;
- superiority in the heterogeneous-affine scenario;
- excellence of the LOO robust estimator;
- temporal, three-dimensional, vadose-zone or vector-stalk performance;
- active-learning superiority (active learning was not tested in M7.5).

The authoritative numerical record is `results/RUN-M7-ROBUST-HYBRID-SHEAF-
20260729-01/locked_test/claim_decision.json`; complete paired contrasts are in
`paired_bootstrap_contrasts.csv`.
