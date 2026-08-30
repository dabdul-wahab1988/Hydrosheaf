# M7.4 Protocol: What the sheaf contributes beyond a weighted graph

## Status and scope

This protocol defines a prospectively locked, controlled-synthetic comparator
experiment. It supplements, and does not alter, the locked M7.3 conditional-
integration benchmark or the strict M7 public-pipeline system-acceptance run.
The target question is:

> Does an affine sheaf section supply incremental topology-discrimination,
> calibration, or inconsistency-localisation value beyond a competence-matched
> ordinary weighted graph?

The experiment is static and scalar on two-dimensional abstract graph layouts.
Temporal-series, spatial three-dimensional, and vadose-zone capabilities are
excluded by instruction and cannot be inferred from its results.

## Scientific contrast

Every comparator receives the same directed candidate edges, weak hydraulic
prior probability, noisy node observations, missingness mask, development/test
split, calibration model, regularisation strength, optimisation routine,
iteration count, threshold-selection budget, and test cases. The two global
models also share one calibration link fitted to stacked identity and affine
development residuals, so only the global relation represented on an edge
changes in their direct comparison.

1. **Edge-local weighted graph (`weighted_graph`).** A calibrated edge model
   combines the common prior with the locally observed affine residual. It has
   no global section and cannot impute a missing endpoint from neighbouring
   constraints.
2. **Ordinary graph Laplacian (`graph_laplacian`).** The same global weighted
   least-squares routine is used with identity restrictions
   $x_v=x_u$. This is the scalar weighted-graph/Laplacian comparator.
3. **Affine sheaf (`affine_sheaf`).** The production HydroSheaf directed-section
   solver is used with edge-specific restrictions
   $x_v=\alpha_e x_u+b_e$.
4. **Permuted-map adverse control (`affine_sheaf_permuted`).** Affine maps are
   permuted among edges within each locked case and the native sheaf calibrator
   is applied without refitting. Marginal map distributions are preserved while
   edge-specific meaning is destroyed.

The graph-Laplacian and sheaf solutions are iteratively reweighted three times
using the same common prior and residual rule. The identity comparator is
implemented through the same numerical solver as the sheaf so a difference
cannot be attributed to a different optimiser. Each calibrated comparator has
the same three inputs: common prior logit, one log-transformed residual, and the
common local-missingness indicator. All use unweighted L2 logistic regression
with `C=1`, `lbfgs`, and 2,000 maximum iterations. The local model has one
calibrator; the identity-Laplacian and affine-sheaf development rows are stacked
to fit one shared global calibrator and threshold. Threshold selection uses the
same 81-point grid from 0.10 to 0.90. The permuted-map control reuses the native
global calibrator without refitting.

## Independent generator and locked split

The generator is implemented in `scripts/independent_sheaf_graph_generator.py`
and imports no HydroSheaf code. It creates rooted directed trees with planted
candidate alternatives, edge maps, scalar latent sections, noisy observations,
and controlled missingness. Thirty-two development cases (seeds 7401-7432)
fit all calibrators and thresholds. Sixty-four locked test cases (seeds
7501-7564) are scored without refitting. Consecutive seeds are balanced across
four predeclared strata:

- `identity_limit`: every restriction is the identity; graph and sheaf raw
  section residuals must coincide numerically.
- `heterogeneous_affine`: valid true edges carry non-identity gains and offsets.
- `incompatible_cycles`: reverse/cross edges plant globally incompatible cycles,
  with at least one endpoint hidden where possible so the test is not reducible
  to a local residual.
- `noisy_missing`: affine restrictions are evaluated with greater analytical
  noise and missingness.

The generator truth labels are excluded from every feature computation. Truth
is used for development calibration and is joined to locked-test predictions
only for scoring.

## Outcomes

Primary outcomes are per-case PR-AUC and Brier score. Secondary outcomes are
ROC-AUC, log loss, ten-bin expected calibration error, development-threshold
selected F1, false-confidence rate, confident-prediction coverage and accuracy,
planted-conflict localisation PR-AUC, and runtime. All contrasts are paired by
independent case. Percentile 95% confidence intervals use 10,000 case-block
bootstrap replicates.

## Predeclared gates and claim rule

The execution/equivalence gate requires finite primary metrics, a generator
that imports no HydroSheaf code, disjoint development and test seeds, and a
maximum raw residual difference of at most $10^{-10}$ between the identity
graph and affine sheaf in the identity-limit stratum.

A controlled-synthetic incremental sheaf claim is allowed only if all of the
following hold:

1. the execution/equivalence gate passes;
2. sheaf-minus-edge-local-graph PR-AUC has CI lower bound above zero;
3. sheaf-minus-edge-local-graph Brier score has CI upper bound below zero;
4. sheaf-minus-graph-Laplacian PR-AUC has CI lower bound above zero;
5. sheaf-minus-graph-Laplacian Brier score has CI upper bound below zero;
6. native-sheaf-minus-permuted-map PR-AUC has CI lower bound above zero;
7. in the identity-limit stratum, sheaf-minus-graph PR-AUC has CI lower bound
   at least -0.02 and Brier-score CI upper bound at most 0.01.

If the full rule fails but sheaf-minus-graph conflict-localisation PR-AUC has a
CI lower bound above zero in the incompatible-cycle stratum, only a conditional
conflict-localisation claim is allowed. Otherwise no incremental sheaf claim is
allowed.

No outcome can establish field validity, universal superiority, temporal or
three-dimensional performance, vadose performance, or superiority over every
possible graph model. A graph endowed with the same non-identity restriction
maps is mathematically approaching the sheaf model; the claimed increment is
therefore specifically over ordinary scalar weighted graphs with identity
coupling and edge-local weighted scoring.
