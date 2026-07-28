# M8 frontier Bayesian active-learning protocol

Planned locked run: `RUN-M8-FRONTIER-AL-20260728-01`

## Scientific question

Can an explicit, cost-aware and scenario-robust Bayesian value-of-information
policy select groundwater measurements that reduce uncertainty about a joint
flow-topology ensemble more efficiently than random sampling while remaining
noninferior to the existing uncertainty-priority heuristic?

## Method contract

The production method shall:

1. accept named joint hypotheses and normalised posterior weights;
2. represent an action by measurement type, target, feasibility, relative cost,
   and conditional predictive distributions under every hypothesis;
3. calculate scalar expected information gain by deterministic Gauss--Hermite
   quadrature;
4. combine model-discrepancy scenarios using a declared mean/worst-case robust
   utility and report both components;
5. combine normalised expected information gain (weight 0.05) with expected
   posterior Brier-risk reduction (weight 0.95), then divide the robust
   combined utility by relative cost;
6. update the posterior with scenario-averaged likelihoods;
7. select batches by marginal joint information gain using a reproducible Sobol
   estimator so redundant measurements are not double counted; and
8. abstain when the maximum robust information gain is below the locked floor.

The method follows the Bayesian experimental-design principle of maximising
mutual information between a future observation and uncertain hypotheses. Its
scenario-robust objective is deliberately described as a finite-scenario
robustness analysis, not as an implementation of a particular ambiguity-set
relaxation. Greedy joint information, rather than independent top-k scores, is
used for batch selection because independent ranking can select redundant
actions.

Primary methodological anchors:

- Go and Isaac (2022), *Robust expected information gain for optimal Bayesian
  experimental design using ambiguity sets*, Proceedings of Machine Learning
  Research 180:728--737,
  https://proceedings.mlr.press/v180/go22a.html.
- Kirsch, van Amersfoort and Gal (2019), *BatchBALD: Efficient and diverse
  batch acquisition for deep Bayesian active learning*, NeurIPS 2019,
  https://papers.neurips.cc/paper/8925-batchbald-efficient-and-diverse-batch-acquisition-for-deep-bayesian-active-learning.

## Independent closed-loop benchmark

- Truth source: the existing code-independent MODFLOW 6/MODPATH 7 heterogeneous
  aquifer generator. It imports no HydroSheaf package code.
- Initial information: HydroSheaf hydraulic candidate graph only. Generator
  truth and hidden measurement outcomes are unavailable to acquisition.
- Joint hypothesis ensemble: 256 topology particles. Raw hydraulic candidate
  probabilities are calibrated by a class-balanced logistic model fitted only
  on development cases; calibrated odds define a source-conditioned
  categorical distribution with an explicit unit-mass no-edge state. The same
  frozen prior calibration is supplied to every strategy.
- Hidden action types and relative resource costs:
  - major-ion chemistry panel: 2 units;
  - groundwater-age tracer: 5 units;
  - directed connectivity tracer test: 9 units.
- Budget: 10 relative units; maximum five sequential actions.
- Predictive likelihoods: class-conditional ridge regressions fitted only on
  development cases using distance, absolute hydraulic-head difference and
  hydraulic candidate probability. To control development overfit, effective
  coefficients equal 0.25 times the fitted deviation from a fully pooled
  class-mean model; residual variances are recalculated after shrinkage. Floors
  prevent zero residual variance.
- Robust scenarios:
  - nominal fitted likelihood, weight 0.50;
  - separation stress, with class separation multiplied by 0.65 and residual
    standard deviation multiplied by 1.25, weight 0.25;
  - noise stress, with residual standard deviation multiplied by 1.60, weight
    0.25.
- Robust acquisition weight: 0.75 on the worst scenario and 0.25 on the
  scenario-weighted mean.
- Decision utility weight: 0.95 on expected Brier-risk reduction and 0.05 on
  expected joint-hypothesis information, after normalising each by its current
  reducible risk.
- Expected-information abstention floor: 0.002 nats.

Development calibration seeds are 7401--7408. Development tuning seeds are
7451--7458. They may be inspected while finalising the method. The frozen
calibration model and every source listed in the lock must be hashed before any
locked-test seed is executed.

Locked-test seeds are 7601--7624. They must not be executed before the protocol
lock exists and passes source-hash verification.

## Comparators

All strategies receive the same initial topology particles, hidden outcomes,
budget, update rule and observation likelihoods.

- `robust_information_decision_per_cost`: proposed method.
- `mean_information_decision_per_cost`: the same information/decision utility
  without worst-case scenario weighting.
- `legacy_uncertainty_chemistry`: selects the edge with greatest marginal
  Bernoulli entropy and requests a fixed chemistry panel. This makes the legacy
  topology-priority concept executable without claiming it can choose a
  measurement modality.
- `random_feasible`: uniformly random affordable action.
- `realised_oracle`: uses hidden outcomes to choose the action with greatest
  realised Brier improvement per cost; it is an unattainable upper bound and is
  excluded from superiority gates.

## Outcomes and uncertainty

Per case and strategy, report:

- candidate recall;
- final edge Brier score, PR-AUC and log loss;
- joint-hypothesis entropy reduction and entropy reduction per cost;
- cost spent, action count, abstention and selected measurement types; and
- regret relative to the realised oracle.

Use 5,000 paired case-bootstrap resamples for every proposed-versus-comparator
contrast. Report medians and 95% percentile intervals. Common random numbers
must be used across strategies.

## Claim gate

The bounded active-learning claim is supported only if all of the following
hold on the untouched locked test:

1. mean candidate recall is at least 0.80;
2. the proposed method is actionable in at least 90% of cases;
3. its paired Brier-score difference versus random has a strictly negative 95%
   interval, and its upper interval versus the strong legacy uncertainty policy
   is below the predeclared absolute noninferiority margin of +0.01;
4. its paired entropy-reduction-per-cost difference versus random has a
   strictly positive 95% interval, and the lower interval versus the strong
   legacy policy exceeds the predeclared absolute noninferiority margin of
   -0.01 nats per relative-cost unit;
5. its PR-AUC difference versus both random and the legacy policy is not
   materially harmful: the lower interval bound must exceed -0.01; and
6. every registered artifact regenerates byte-identically from the locked
   source and calibration-model hashes.

Mean-EIG and oracle comparisons are diagnostic rather than pass/fail gates. A
failure of any gate requires a negative or mixed claim decision; no threshold
may be changed after the locked test is inspected.

## Claim boundary

A passing result permits only this class of statement: under a controlled
synthetic, independently generated aquifer ensemble and declared relative-cost
model, scenario-robust Bayesian acquisition improved Brier score and topology
uncertainty per cost over random acquisition, and was noninferior in both
outcomes to the strong uncertainty-only chemistry policy. It does not establish
field effectiveness, universal optimality, or correctness of supplied
likelihoods and costs.

Temporal-series, three-dimensional and vadose-zone capabilities remain outside
scope by instruction.
