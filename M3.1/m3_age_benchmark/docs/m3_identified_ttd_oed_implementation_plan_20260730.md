# M3 identified-TTD and measurement-design implementation plan

**Status:** development integration implemented; confirmatory validation remains pending  
**Date:** 2026-07-30  
**Scope:** local transit-time-distribution evidence, held-out tracer prediction, graph compatibility, and measurement design  

## Implementation record (2026-07-30)

The development integration now includes:

- `hydrosheaf/nuclear/ttd_identified.py`: deterministic identified-set
  inference, sharp functional bounds, witness TTDs, rank/nullity diagnostics,
  censored constraints, explicit abstention, and held-out predictive bounds;
- `hydrosheaf/nuclear/joint_lpm.py::tracer_response_kernel`: the shared linear
  tracer forward operator on a supplied age grid;
- `hydrosheaf/nuclear/ttd_graph.py`: non-negative, mass-conserving graph
  compatibility and obstruction audits with tightening disabled;
- `hydrosheaf/nuclear/ttd_design.py`: a tracer/well/time adapter to the existing
  cost-aware information-gain engine, guarded by declared probability
  semantics;
- `configs/identified_ttd_protocol.yaml`: a hashable development protocol with
  locked age grid, functionals, correction modes, leakage controls, and claim
  boundary; and
- `scripts/run_m3_identified_ttd_benchmark.py`: local-only, tracer-withholding
  evaluation with native-unit coverage/width outputs and a hash-bearing
  manifest.

This is implementation evidence, not confirmatory scientific evidence. The
posterior generator, prospective synthetic-truth claim decision, full USGS
held-out run, and adverse-control graph experiment remain future work. No new
field-validation, true-TTD, graph-benefit, or optimal-design claim is authorized
until those gates pass.

## 1. Decision

HydroSheaf is not missing all of the required machinery. It is missing the integration layer that makes the machinery answer one coherent question:

> Which decision-relevant properties of a groundwater transit-time distribution are supported by the observed tracer panel, which remain unresolved, how does independently sourced topology change that support, and which feasible measurement is expected to reduce the remaining uncertainty?

The current implementation should not be replaced wholesale. The following assets should be retained:

- `hydrosheaf/nuclear/joint_lpm.py` already supplies tracer forward models, fitted LPM parameterizations, transit-time PDFs, and point age-fraction predictions.
- `hydrosheaf/nuclear/diagnostics.py` already diagnoses acceptable-model ambiguity and tracer-removal sensitivity.
- `hydrosheaf/nuclear/network_aging.py` already separates sampler convergence from tracer identifiability and performs posterior-predictive checks.
- `M3/m3_age_benchmark/scripts/run_m3_usgs_benchmark.py` already implements a reference-free scalar reportability guard.
- `hydrosheaf/calibration/bayesian_active_learning.py` already implements probability-bearing, cost-aware, scenario-sensitive expected-information-gain ranking, abstention, posterior updating, and batch redundancy control.
- `M3/m3_age_benchmark/scripts/run_m3_cross_validation_benchmark.py` already provides the leakage-control pattern for held-out tracer prediction.

The older `rank_next_measurements()` function must remain labeled as categorical topology triage. It is not the engine to extend for tracer/well/time design. The newer `bayesian_active_learning.py` engine is the reusable design engine, but it still needs a groundwater-tracer adapter and probability-bearing TTD hypotheses. Its controlled-synthetic M8 evidence does not by itself validate groundwater-age measurement design.

## 2. Non-negotiable scientific contract

### 2.1 Identified sets and posterior ensembles are different objects

HydroSheaf must never use the terms interchangeably.

- An **identified set** contains every TTD satisfying declared observation-error, forward-model, and physical constraints. Its reported limits are bounds, not credible intervals.
- A **posterior ensemble** contains probability-weighted TTD draws under a declared likelihood and prior. Its reported limits may be credible intervals, but prior sensitivity must be shown when the tracer panel is weak.
- AIC weights, uniform weights over optimizer candidates, and arbitrary weights over feasible LP vertices are not posterior probabilities. If used for sensitivity or design, their semantics must be named explicitly.

The identified set is the primary output. A posterior ensemble is optional and is needed only for Bayesian expected-information-gain design.

### 2.2 Primary mathematical object

Represent a local TTD on a fixed fine age grid as a non-negative probability-mass vector `w`:

```text
w_j >= 0
sum_j w_j = 1
```

For tracer `i`, the forward response is:

```text
y_hat_i = sum_j A[i, j] * w_j
```

`A[i, j]` is the tracer response to unit recharge mass at age-grid point `j`, including the locked input history, radioactive decay, sampling date, and only those correction factors supported for that observation. This fine-grid formulation avoids hiding a uniform-within-bin assumption.

For the first implementation, define the feasible set with auditable interval constraints:

```text
y_i - k*sigma_i <= A[i, :] @ w <= y_i + k*sigma_i
```

Use one-sided constraints for censored observations. Do not silently convert contaminated-mixture likelihoods into interval constraints. Mark such rows unsupported or evaluate them as separately named model-discrepancy scenarios.

Decision-relevant quantities are linear functionals `c @ w`, including:

- fraction younger than a locked threshold;
- fraction within each locked age band;
- fraction older than a locked threshold;
- mean residence time over the finite locked age domain; and
- CDF values at locked ages.

For each functional, solve two linear programs to obtain the minimum and maximum over the feasible set. Store the extremal TTDs as witnesses and retain solver feasibility certificates. Do not report a median age in version 1 because a median is not a linear functional; CDF bounds are the defensible precursor.

### 2.3 Reportability

The report must use three states:

- `IDENTIFIED`: every required decision functional has a feasible bound narrower than its predeclared threshold.
- `PARTIALLY_IDENTIFIED`: the feasible set exists and at least one functional is useful, but one or more required functionals remain too broad.
- `ABSTAIN`: the feasible set is empty under all supported scenarios, the observation panel is unsupported, or no required functional is informative.

The existing M3 scalar reportability guard remains available for legacy comparisons. It must not be used as the definition of set identification. A 100-year scalar likelihood-width cutoff must not be reused as a universal threshold. Functional-specific thresholds belong in the locked protocol.

### 2.4 Graph information

Local tracer inference must be completed and frozen before graph information is introduced. An edge may narrow a local identified set only when:

- the edge and its direction have provenance independent of every tracer being evaluated;
- the restriction/transport operator is declared before evaluation;
- the operator is mass conserving and non-negative;
- graph hyperparameters are chosen without the held-out tracer; and
- reversed, randomized, and identity-map controls are evaluated.

Without those conditions, graph information may be used only for compatibility or falsification, not for tightening age claims.

## 3. Target architecture

### 3.1 New module: `hydrosheaf/nuclear/ttd_identified.py`

Add immutable public data structures:

```python
@dataclass(frozen=True)
class AgeFunctional:
    name: str
    coefficients: np.ndarray
    maximum_reportable_width: float
    units: str

@dataclass(frozen=True)
class IdentifiedBound:
    name: str
    lower: float
    upper: float
    status: str
    lower_witness: np.ndarray
    upper_witness: np.ndarray
    solver_diagnostics: Mapping[str, Any]

@dataclass(frozen=True)
class TtdEvidenceReport:
    status: str
    age_grid_years: np.ndarray
    bounds: Mapping[str, IdentifiedBound]
    response_rank: int
    nullity: int
    feasible_scenarios: Sequence[str]
    excluded_scenarios: Mapping[str, str]
    assumptions: Mapping[str, Any]
    provenance: Mapping[str, Any]
```

Add public functions:

```python
build_tracer_response_matrix(...)
build_age_functionals(...)
solve_ttd_identified_set(...)
infer_ttd_evidence(...)
predict_tracer_bounds(...)
```

Use `scipy.optimize.linprog(method="highs")` for the first implementation. Verify every returned witness independently against non-negativity, mass balance, and tracer constraints. A solver success flag without a primal-feasibility audit is not sufficient.

### 3.2 Minimal refactor: `hydrosheaf/nuclear/joint_lpm.py`

Do not break `JointLpmFit`. Extract and expose only the shared forward-kernel operations needed by `ttd_identified.py`:

- a public age grid/integration-grid builder;
- a public tracer-response kernel on a supplied age grid; and
- a public observation-normalization path that preserves likelihood type, units, sigma, correction provenance, and input-history version.

Keep `age_fraction_predictions()` for point LPM outputs. Add a separate helper that applies the same locked age-band definitions to a probability-mass vector. Do not make point predictions appear set-valued.

### 3.3 New module: `hydrosheaf/nuclear/ttd_posterior.py`

This is phase 4, not phase 1. Implement an optional probability-bearing ensemble on a deliberately low-dimensional, predeclared basis. Requirements:

- simplex-valued age-bin masses;
- explicit likelihood and prior;
- prior-sensitivity scenarios;
- convergence checks separated from identifiability;
- posterior-predictive checks in native tracer units; and
- an export that clearly labels draws as posterior draws.

If posterior results change materially across reasonable priors, the identified-set output remains primary and the posterior must be labeled prior-sensitive.

### 3.4 New module: `hydrosheaf/nuclear/ttd_graph.py`

Add:

```python
@dataclass(frozen=True)
class MassTransportMap:
    edge_id: str
    matrix: np.ndarray
    provenance_tier: str
    source: str

audit_ttd_graph_compatibility(...)
condition_ttd_sets_on_graph(...)
```

Every transport matrix must be non-negative and mass conserving. Version 1 should implement graph compatibility and obstruction diagnostics only. Set tightening is version 2 and is disabled by default until the independent-provenance gate passes.

The identity restriction map is a graph-smoothing baseline, not evidence of a non-trivial sheaf contribution. A publishable sheaf claim requires a competence-matched graph baseline and directly tested non-identity transport maps.

### 3.5 New module: `hydrosheaf/nuclear/ttd_design.py`

Add a groundwater-specific adapter to the existing Bayesian acquisition engine:

```python
@dataclass(frozen=True)
class CandidateTracerMeasurement:
    option_id: str
    well_id: str
    tracer: str
    sample_year: float
    sigma: float
    cost: float
    feasible: bool
    input_history_id: str
    metadata: Mapping[str, Any]

build_measurement_options_from_ttd_ensemble(...)
rank_ttd_measurements(...)
```

For every posterior/design-ensemble member, forward-predict each candidate tracer/well/time action. Convert those predictions and declared observation errors into `PredictiveScenario` and `MeasurementOption` objects, then call `rank_measurement_options()` or `select_measurement_batch()`.

Pass decision-functional values, such as young fraction or threshold exceedance, as `decision_values`. This makes acquisition target a decision rather than entropy alone.

The adapter must refuse EIG ranking if the supplied ensemble has no defensible probability semantics. For a pure identified set, return `ABSTAIN_NO_PROBABILITY_MODEL` and offer a separately named worst-case bound-reduction analysis. Do not invent uniform probabilities over the feasible set.

## 4. Work packages

### WP0: protect the current evidence baseline

Before coding:

1. Commit or otherwise isolate the current intended M3 and nuclear changes. The present worktree contains uncommitted edits in `joint_lpm.py`, the M3 benchmark, tests, and other project areas.
2. Create a clean worktree or branch named `codex/m3-identified-ttd` from the approved baseline.
3. Record the commit, Python environment, dependency lock, input-history hashes, USGS source hashes, and the existing result-manifest hash.
4. Freeze `M3/m3_age_benchmark/configs/identified_ttd_protocol.yaml` before confirmatory evaluation.

Stop if the intended baseline is ambiguous. Do not overwrite or normalize unrelated user changes.

### WP1: deterministic identified-set core

Files:

- add `hydrosheaf/nuclear/ttd_identified.py`;
- minimally refactor `hydrosheaf/nuclear/joint_lpm.py`;
- export the new public API from `hydrosheaf/nuclear/__init__.py` and, lazily, from `hydrosheaf/__init__.py`;
- add `tests/test_ttd_identified.py`.

Required tests:

1. probability-mass, non-negativity, and finite-value invariants;
2. exact recovery of bounds for a small hand-solvable response matrix;
3. two distinct TTDs with the same tracer response remain unresolved;
4. an added informative tracer narrows or leaves unchanged every feasible bound;
5. duplicating a tracer does not create information when its dependence is declared;
6. affine unit conversion leaves bounds unchanged;
7. censored observations create the correct one-sided constraints;
8. infeasible observations cause explicit abstention;
9. lower and upper witnesses independently satisfy every constraint; and
10. grid refinement does not materially change locked synthetic bounds.

WP1 is complete only when the existing focused baseline and all new tests pass.

### WP2: M3 held-out tracer evaluation

Files:

- add `M3/m3_age_benchmark/scripts/run_m3_identified_ttd_benchmark.py`;
- add `M3/m3_age_benchmark/tests/test_m3_identified_ttd_benchmark.py`;
- add result, manifest, QA, table, and figure generators only after the protocol is locked.

For each eligible sample and held-out tracer:

1. remove the target tracer before model construction;
2. build the local identified set from the remaining tracers;
3. propagate the set through the held-out tracer forward model;
4. score whether the observed concentration falls within the prediction bounds;
5. record bound width in native units and standardized by observation sigma;
6. preserve all abstentions and failure reasons; and
7. prevent reference ages and reported age fractions from entering fitting, selection, thresholds, graph construction, or hyperparameter tuning.

Primary real-data metrics:

- reportable coverage;
- empirical held-out concentration coverage;
- median and upper-quantile prediction width;
- abstention rate and reason distribution;
- calibration stratified by tracer, age-information regime, and supported correction scenario; and
- comparison with best-fit scalar LPM and a deliberately non-informative full-simplex baseline.

USGS ages may be used only for post-fit parity/emulation context. M3 cannot validate true TTD recovery because the release does not contain true TTDs.

### WP3: graph compatibility and conditional tightening

Files:

- add `hydrosheaf/nuclear/ttd_graph.py`;
- add `tests/test_ttd_graph.py`;
- extend the M3 runner with a graph-disabled local baseline and separately named graph candidates.

Version 1 tests graph compatibility without narrowing local sets. Version 2 may intersect constraints across nodes only after the independent-provenance gate passes.

Required controls:

- no graph;
- identity-map graph;
- independently sourced candidate graph;
- reversed graph;
- degree-matched randomized graph; and
- edge-removal sensitivity.

A graph benefit requires both:

- non-inferior held-out concentration coverage under a prospectively frozen margin; and
- narrower decision-functional bounds or prediction intervals with paired uncertainty excluding no change.

If wrong or randomized graphs also improve, the graph contribution is not identified. If the candidate graph only raises incompatibility, report topology falsification rather than predictive benefit.

### WP4: posterior/design ensemble and OED adapter

Files:

- add `hydrosheaf/nuclear/ttd_posterior.py`;
- add `hydrosheaf/nuclear/ttd_design.py`;
- add `tests/test_ttd_posterior.py` and `tests/test_ttd_design.py`.

Reuse `hydrosheaf/calibration/bayesian_active_learning.py`; do not duplicate its EIG, cost, scenario, abstention, or batch-selection logic.

Required OED tests:

1. a candidate with identical predictions under all hypotheses has zero EIG;
2. EIG increases with signal-to-noise ratio;
3. rankings are invariant to units and hypothesis order;
4. high cost can demote a more informative action;
5. infeasible tracer/well/time actions are excluded and audited;
6. prior-sensitive recommendations are labeled and may trigger abstention;
7. batch selection rejects redundant time points or tracers;
8. the selected action improves realized decision loss in a locked independent synthetic test more than random acquisition; and
9. no field or optimality claim is emitted from a retrospective benchmark.

### WP5: prospective validation and claim decision

Use an independent synthetic-truth generator that does not call the HydroSheaf inverse code. Include TTDs outside the fitted LPM families, mixed-age distributions, input-history perturbations, correction errors, detection limits, wrong topology, and missing tracers.

Freeze before the confirmatory run:

- synthetic generator version and seeds;
- age grid and decision functionals;
- observation-error and model-discrepancy scenarios;
- graph candidates and adverse controls;
- acquisition costs and feasibility;
- baselines;
- primary metrics, uncertainty intervals, non-inferiority margins, and failure rules; and
- all file hashes and regeneration commands.

The confirmatory runner must write immutable manifests and refuse to overwrite a completed run unless an explicit development flag is supplied.

## 5. Proposed protocol gates

The exact values must be approved and frozen before the confirmatory run. A defensible starting proposal is:

- nominal 95% set-valued held-out prediction coverage with the lower 95% binomial confidence bound at or above 0.90 in controlled synthetic truth;
- no more than a 0.02 absolute coverage loss for graph-conditioned predictions relative to the local-only baseline;
- a paired median width reduction whose 95% interval excludes zero before claiming graph tightening;
- no claimed graph benefit if a reversed or degree-matched randomized graph passes the same rule;
- OED actionability only above the frozen minimum robust EIG and only when the recommendation is stable across declared prior/model-discrepancy scenarios;
- improvement over random acquisition and non-inferiority to a strong sensitivity/uncertainty baseline on untouched synthetic cases; and
- exact regeneration of every registered result artifact.

Development data may be used to set these thresholds, but confirmatory data must remain untouched until the protocol is locked.

## 6. Claim ladder

Code completion alone permits only:

> HydroSheaf implements set-valued TTD inference, explicit abstention, and a conditional adapter for measurement design.

Passing independent synthetic-truth validation may permit:

> Under the locked synthetic regimes, HydroSheaf's identified bounds achieved the declared coverage and width criteria, and the declared measurement-design policy improved the locked decision metric relative to its comparators.

Passing M3 held-out tracer evaluation may permit:

> On the public USGS release, set-valued fits were evaluated by leakage-guarded prediction of withheld tracer concentrations, with explicit coverage, width, and abstention reporting.

M3 still does not permit claims of true-age recovery, true-TTD recovery, verified field flow paths, universal graph benefit, or field-validated measurement optimization.

## 7. Definition of done

The integration is complete only when:

- local output is a set or a probability-bearing ensemble, not only a best scalar age;
- every reported functional has an auditable support interval and witness distributions;
- unsupported samples abstain with machine-readable reasons;
- held-out tracer concentrations are the real-data validation target;
- graph information enters only after local inference is frozen and only with independent provenance;
- graph candidates are tested against adverse controls;
- tracer/well/time actions are passed into the existing Bayesian acquisition engine through a dedicated adapter;
- EIG is never computed from invented probabilities;
- controlled synthetic truth, retrospective public-data prediction, and future prospective field evidence are reported as separate evidence tiers;
- focused tests, the full suite, manifests, and artifact-regeneration checks pass; and
- a claim-decision document records both passed and failed gates.

## 8. Agent execution order

Do not ask one agent to implement all work packages in one undifferentiated patch. Use the following sequence, reviewing and merging each package before starting the next:

1. WP0 and WP1 only.
2. Review mathematical constraints, API stability, witness feasibility, and tests.
3. WP2 only, followed by a development run and protocol revision.
4. Freeze the confirmatory protocol.
5. WP3 compatibility-only version, then decide whether graph tightening is scientifically eligible.
6. WP4 posterior/OED adapter.
7. WP5 untouched confirmatory run and claim decision.

### Copy-paste prompt for the first implementation agent

```text
Implement WP0 and WP1 from
M3/m3_age_benchmark/docs/m3_identified_ttd_oed_implementation_plan_20260730.md.

First inspect git status. The current repository may contain uncommitted user work in
joint_lpm.py, M3 benchmark scripts, and tests. Do not overwrite, reset, format, or
otherwise alter unrelated changes. If the approved baseline is not clean or cannot be
isolated safely, stop and report the exact overlapping files.

Create the deterministic identified-set core in a new
hydrosheaf/nuclear/ttd_identified.py module. Use a fine age-grid probability-mass
representation and scipy.optimize.linprog(method="highs"). Implement independently
verified primal-feasibility checks and return extremal witness TTDs for each functional.
Keep identified-set bounds distinct from posterior intervals. Do not implement graph
conditioning, posterior sampling, active learning, manuscript edits, or scientific claims
in this first package.

Make only the minimal joint_lpm refactor needed to expose a public tracer-response kernel
and observation normalization without breaking JointLpmFit. Add public exports and the
complete tests listed under WP1. Preserve all existing behavior.

Run at minimum:
  .\.venv\Scripts\python.exe -m pytest tests\test_joint_lpm.py tests\test_network_aging.py tests\test_bayesian_active_learning.py tests\test_m3_usgs_benchmark.py tests\test_ttd_identified.py -q

Then run the full suite if the focused set passes. Report files changed, exact tests,
solver/feasibility assumptions, remaining limitations, and any divergence from the plan.
Do not commit or push unless explicitly authorized.
```
