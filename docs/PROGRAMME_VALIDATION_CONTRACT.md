# HydroSheaf programme validation contract

**Status:** working contract; synthetic lock and adjudication implemented; protocol release pending
**Prepared:** 2026-08-01  
**Purpose:** turn the M2-M8 lessons into one auditable, end-to-end validation programme.

## 1. Frozen system-level question

Under independently generated groundwater systems with hydraulic, tracer, and
hydrochemical observations, can HydroSheaf:

1. report calibrated uncertainty over flow connectivity, groundwater age, and
   reaction family;
2. identify when those quantities are non-identifiable or mutually
   inconsistent; and
3. select the next feasible measurement that reduces joint decision risk per
   unit cost more safely than strong specialist baselines?

The primary claim is conditional on the declared data, generator, observation
model, costs, and discrepancy scenarios. It is not a claim of universal
groundwater-model superiority and it is not field validation.

## 2. Primary decision output

For each case and decision round, the system must return one of:

- `ACTION`: a named measurement, target, cost, expected utility, and the
  posterior state before and after the hypothetical measurement;
- `SET_REPORT`: a compatible set of flow, age, or reaction hypotheses when a
  unique answer is not supported; or
- `ABSTAIN`: no action is sufficiently informative, feasible, or robust across
  the declared model-discrepancy scenarios.

An apparently precise point estimate without a calibrated uncertainty or an
explicit identifiability status is not a successful system output.

## 3. Validation layers

| Layer | Question | Current state | What is needed for a gate |
|---|---|---|---|
| L1 component | Do topology, age, chemistry, kinetics, and active-learning components behave as specified? | Substantial locked M2-M8 evidence exists | Common registry, common metrics, common provenance check, and no component result promoted to a system claim |
| L2 integrated synthetic | Does the complete decision workflow recover withheld edge, age, and reaction truth under independent generators? | Executable active package. The independent M7 v3 generator supplies a two-layer MODFLOW/MODPATH field, fixed branch/merge observation mixtures, weighted age truth, nonlinear chemistry, tracers, and sealed truth; analytic-lattice and independent-mixing provide additional families | Larger locked ensembles, competence-matched specialist baselines, calibrated uncertainty, and a preregistered decision-utility rule |
| L3 external retrospective | Does the workflow transfer to data with genuinely independent reference information? | Partly available, but M3 reference ages and several field products are model-derived or lack process truth | Separate emulation from truth validation; add external datasets with independent age, path, or reaction reference information |
| L4 prospective field | Does a pre-registered HydroSheaf decision improve a real sampling campaign? | Not currently possible from the repository alone | New campaign with pre-issued predictions, repeated heads, age tracers, chemistry, screen metadata, and independent outcome measurements |

## 4. Required comparators

All methods must receive the same observations, missingness, uncertainty, cost
budget, and tuning allowance.

1. Hydraulic/head-only candidate graph.
2. Strong edge-local graph model.
3. Local tracer-age model without global sheaf refinement.
4. Conventional PHREEQC inverse or an explicitly defined thermodynamic inverse
   baseline.
5. HydroSheaf-Core and the full evidence-gated HydroSheaf workflow.
6. The existing M8 uncertainty-priority policy.
7. Random feasible acquisition and permuted-evidence controls.
8. An oracle diagnostic using withheld truth, excluded from superiority claims.

The comparator registry must record the input channels, parameter budget,
uncertainty output, convergence rule, abstention rule, and computational cost
for every method.

## 5. Model-discrepancy and model-averaging contract

The workflow must not silently assume that one interpretation is correct.
Each case must be evaluated under at least these declared scenario families:

- correct or nominal topology;
- false-positive and false-negative graph edges;
- scalar-age versus distributed-age mismatch;
- reaction-family misspecification or equivalence classes;
- shared tracer nuisance such as recharge temperature, C14 initial activity,
  or CFC degradation;
- observation noise, censoring, and structured missingness.

For each scenario, the system produces a predictive distribution and a
decision utility. Development-only cases may estimate model weights using
proper scoring, AICc/LOO, or another frozen rule. Locked-test truth must never
choose the weights. When posterior predictions disagree materially across
models, the result must be reported as `MODEL_DISAGREEMENT`, converted to a
compatible set, or abstained. Model averaging must not be used to hide
disagreement.

## 6. Primary metrics and gates

The primary metrics are:

- coverage and width of age and reaction compatible sets;
- Brier score, log loss, and calibration error for edge and family outputs;
- selective risk at each abstention rate;
- false-positive and false-negative rates for topology;
- entropy reduction and expected decision-risk reduction per unit cost;
- actionability and abstention frequency;
- recovery of withheld truth for age and reaction families; and
- reproducibility of all registered artifacts.

An integrated superiority claim requires a predeclared joint rule. At minimum,
the proposed workflow must improve the primary cost-adjusted decision metric
against random acquisition, remain non-inferior to the strongest specialist
baseline within a locked margin, and show no material degradation in
calibration, false commitment, or coverage. If the joint rule fails, the
allowed result is conditional, mixed, or negative.

## 7. Provenance and change control

Before locked testing, freeze:

- generator source and independent-import audit;
- generator families, development seeds, tuning seeds, and locked seeds;
- observation and missingness models;
- baseline implementations and tuning budgets;
- model-discrepancy scenarios and weights;
- primary metrics, non-inferiority margins, and stopping rule;
- source hashes, dependency lock, executable hashes, environment hash, and
  input manifest hash.

Every run receives a new immutable run ID. Failed runs remain recorded. No
threshold, comparator, seed, or outcome may be changed after inspecting locked
test results.

## 8. What can be done now

The following are executable from the current repository:

1. Build L1's common component registry from the locked M2-M8 claim decisions.
2. Use the independent v3 generator as the active L2 integrated truth source; retain v2 only as a historical extension baseline.
3. Wrap the existing topology, tracer-age, reaction, and active-learning
   outputs into the `ACTION` / `SET_REPORT` / `ABSTAIN` contract.
4. Add the model-discrepancy scenario interface, conservative scenario
   envelope, and model-disagreement flag.
5. Run development cases, freeze the baselines and seeds, then run a new
   locked integrated synthetic test.
6. Rewrite the M8 interpretation so its result is a component of this
   decision framework, not evidence of a universal sampling optimum.

### 8.1 Synthetic implementation checkpoint — 2026-08-01

The first executable development gate is now implemented in the authoritative
package:

- `hydrosheaf/validation/programme_contract.py` defines the allowed
  `ACTION`, `SET_REPORT`, and `ABSTAIN` outputs, identifiability states, stage
  records, and sealed-truth run contract.
- `hydrosheaf/validation/programme_gate.py` owns shared truth-blind row
  preparation, stage checks, topology scoring, and age interval scoring.
- `hydrosheaf/validation/discrepancy.py` defines scenario predictions and a
  conservative, explicitly non-calibrated discrepancy envelope that can emit
  `MODEL_DISAGREEMENT` and `ABSTAIN`.
- `scripts/run_programme_synthetic_gate.py` is a thin adapter around the
  independent M7 v3 generator and the public `fit_network_pipeline` API.
- `tests/test_programme_contract.py` and `tests/test_programme_gate.py` cover
  the contract and shared gate helpers.
- `tests/test_programme_discrepancy.py` covers weighting, envelope formation,
  disagreement detection, and invalid scenario records.

The development records are generated under `.codex_work` and are not treated
as locked evidence:

- `RUN-PROGRAMME-SYNTHETIC-DEV-20260801-01-SMOKE` passed the execution gate for
  one case.
- `RUN-PROGRAMME-SYNTHETIC-DEV-20260801-01` passed the execution gate for
  seeds 9101–9103 across hydraulic-only, local-age, full-sheaf, and
  age-permuted conditions.

The current development readout is diagnostic rather than a performance pass.
Across the three cases, mean topology candidate precision was 0.263, mean
candidate recall was 0.963, and mean selected F1 was 0.413. The local-age
condition had mean MAE 5.47 years with 66.6-year mean 95% interval width; the
full-sheaf condition had mean MAE 5.57 years with 66.5-year mean interval
width. The age-permuted negative control degraded to mean MAE 19.61 years and
0.528 interval coverage.

The public network-fit reaction outputs are now also recorded and scored as a
separate chemistry/reaction stage; this first runner does not yet treat that
component result as a complete reaction-family claim. Across the three cases,
the coarse reaction-family accuracy was 0.620 on 26 of 27 true-process edges;
one true edge had no fitted output. This is useful as a failure-facing
diagnostic, not evidence that the reaction module is field-ready.

An explicit discrepancy report is now wired into the runner for a declared
input-history stress pair: the nominal steady history and a tropical history.
The tropical history is extended to the present before inference so that the
comparison is not caused by a truncated lookup table. Across the three seeds,
36 node-level reports were written; 11 returned `MODEL_DISAGREEMENT` and
`ABSTAIN`. The reports use equal predeclared scenario weights and a
conservative scenario envelope marked `not_calibrated`. This is a useful
disagreement-safety layer, but it is not calibrated model averaging and it is
not an ensemble of independent generator families.

Two limitations are explicit in this first result. First, the current
candidate graph and sheaf refinement both use the same three-neighbour cap, so
the full sheaf retained the same 33 candidate edges in these cases; it has not
yet demonstrated topology discrimination. Second, the v2 generator declares
an analytic tritium response that is intentionally different from HydroSheaf's
input-history lookup. The age result therefore tests model discrepancy, not a
fair matched-forward-model age claim. The current discrepancy layer is limited
to the declared input-history pair with fixed weights; calibrated selective
risk for the age interval output is now available in a separate two-family
development/held-out runner, but full calibration across topology, chemistry,
and decisions, independent-generator model averaging, and joint
next-measurement selection remain unimplemented in this root-level programme
runner.

Accordingly, the current status is:

- software execution gate: **pass**;
- controlled input-history discrepancy reporting: **pass**;
- integrated scientific performance gate: **not yet locked or passed**;
- field validation: **deferred until field data and a prospective protocol are
  available**.

### 8.2 Independent-generator and calibration checkpoint — 2026-08-01

The next synthetic layer is now executable in the authoritative package:

- `scripts/independent_lattice_generator.py` supplies a second generator
  family with a branching/merging multilayer analytic geometry, an analytic
  head field, segment-wise chemistry, different tracer equations, and a
  distinct random stream. It imports neither HydroSheaf nor the
  M7/MODFLOW/MODPATH generator.
- `hydrosheaf/validation/calibration.py` fits an additive split-conformal age
  interval expansion on development truth only, applies the frozen calibrator
  to blind predictions, and reports selective-risk curves.
- `scripts/run_programme_generator_ensemble.py` runs both families through the
  same public inference pipeline, records family-specific provenance, and
  scores held-out cases after calibration without retuning.
- `tests/test_programme_calibration.py` and
  `tests/test_independent_lattice_generator.py` cover the new contracts.

The development/held-out execution record
`RUN-PROGRAMME-ENSEMBLE-CALIBRATION-DEV-20260801-01` passed its execution gate.
That record is retained as the first two-family historical checkpoint; the
current three-family result is recorded below in section 8.3.
The development fit used 48 full-sheaf age predictions from two M7 v2 cases
and two analytic-lattice cases. The frozen 95% calibrator required an additive
radius of 0.0 years because the development intervals already covered the
development truth at that rule; this is a valid calibration result, but it is
not evidence that the intervals are informative or calibrated in general.
Two additional cases were then scored without refitting. The held-out
full-sheaf age MAE was 5.54 years for the M7 v2 family and 3.36 years for the
analytic-lattice family; both had 1.00 interval coverage in this small locked
sample, with mean widths of 60.18 and 69.61 years respectively.

Selective risk decreased from MAE 5.14 when accepting all 24 held-out age
predictions to 1.18 when accepting the least-uncertain 25% (six predictions).
The corresponding interval coverage was 1.00 at every reported acceptance
level, with zero false commitments in this small sample. This is a diagnostic
curve, not a claim of calibrated selective superiority.

The failure-facing results remain important. The revised analytic-lattice
held-out full-sheaf topology precision was 0.23 with selected F1 0.32, and its
coarse reaction-family accuracy was 0.14; the M7 v2 held-out values were 0.27,
0.43, and 0.56 respectively. Across both families, 72 input-history discrepancy
reports were written, with 26 emitting `MODEL_DISAGREEMENT`/`ABSTAIN`. The
software can now expose cross-family transfer and calibrated selective-risk
failures, but the integrated scientific performance gate remains **not
passed**.

### 8.3 Adversarial generator-critic checkpoint — 2026-08-01

The generator is now reviewed by a separate adversarial critic in
`hydrosheaf/validation/generator_critic.py`. It audits actual source imports
with the Python AST and checks, independently of the inference run:

- hidden truth or route fields in observation rows;
- finite and physically non-negative tracer values;
- truth-edge referential integrity, cycles, age direction, and head direction;
- branching, merging, layer heterogeneity, geometry duplication, and dynamic
  tracer/chemistry range;
- duplicated measurement channels and absent missingness/censoring stress;
- exact same-seed reproducibility and different-seed variation; and
- reuse of an earlier generator core.

The ensemble writes `generator_critic.json` and separates the ordinary
execution gate from the generator-quality critic gate. The current
`RUN-PROGRAMME-ENSEMBLE-CALIBRATION-DEV-20260801-03` record covers three
families, four observation scenarios, 72 inference-condition records, and
zero recorded execution errors. Structured missingness, left censoring, and
combined stress are now actually passed through the full inference pipeline.
The missing-redox path was repaired so absent `NO3`/`Fe` evidence produces an
ambiguous state rather than a crash.

The adversarial critic reports zero blockers. The two primary independent
families (analytic lattice and independent mixing) have zero major findings.
The nine remaining major findings are confined to the deliberately retained
M7 v2 baseline: disjoint-chain topology, single-layer geometry, and reuse of
the v1 core. Those findings keep the global critic gate closed, but they are
now explicitly separated from the primary generator-release gate. Highly
correlated head channels are recorded as an information finding after the
declared covariance/discrepancy model is consumed and the secondary channel
is removed from independent evidence.

The independent mixing family is a held-out mechanism family with three
layers, branch/merge/shortcut connectivity, weighted endmember mixing, and a
two-component age-mixture forward model. Its full-sheaf locked diagnostic was
not uniformly strong: age MAE was about 11.9 years and reaction-family
accuracy was about 0.17. This is a useful negative result, not evidence of
general HydroSheaf superiority.

At that checkpoint, the next release requirements were fair specialist
baselines, topology/reaction probability calibration, development-fitted model
averaging, and a cost-adjusted next-measurement benchmark. The comparator and
policy slice is now implemented below, while the current claim remains
controlled-synthetic software validation, not field validation or universal
inference superiority.

### 8.4 Held-out comparator and policy implementation checkpoint — 2026-08-01

The next required software slice is now implemented and exercised in a scoped
analytic-family run:

- `hydrosheaf/validation/baselines.py` provides an auditable registry for the
  hydraulic-only directional-gradient and edge-local topology comparators,
  including fixed tuning, uncertainty, abstention, cost, and permutation
  metadata.
- `hydrosheaf/validation/controls.py` provides deterministic no-flow head,
  tracer, and chemistry permutation controls. Controls preserve row identity
  and channel marginals and are checked for truth blindness.
- `hydrosheaf/validation/metrics.py` now reports Brier score, log loss, ECE,
  and reliability-bin records for binary probabilities. These are diagnostic
  probability scores; they do not silently claim calibration where no
  development calibrator has been fitted.
- `hydrosheaf/validation/decision_utility.py` implements a one-step,
  cost-adjusted robust information-gain policy with explicit abstention and
  serialisable audit records.
- `scripts/run_programme_generator_ensemble.py` now records baseline outputs,
  probability diagnostics, no-flow/permutation conditions, policy decisions,
  and a truth-released synthetic target diagnostic. The policy itself is
  generated before truth scoring; prospective post-measurement utility remains
  deferred.

The scoped run
`RUN-PROGRAMME-ENSEMBLE-CALIBRATION-DEV-20260801-04-ANALYTIC-SMOKE` used the
analytic-lattice and independent-mixing families, completed 44 condition
records with zero recorded errors, and passed its scoped execution gate. The
no-flow control reduced selected topology F1 in both locked cases, while the
full integrated performance gate correctly remained deferred because the
scoped run omitted the MODFLOW/MODPATH family and no prospective measurement
outcome was simulated. This is implementation evidence, not a superiority
result.

The historical three-family smoke
`RUN-PROGRAMME-ENSEMBLE-CALIBRATION-DEV-20260801-04-FULL-SMOKE` confirmed that
the existing local executables run end to end: the external-solver execution
check passed, both MODFLOW/MODPATH development and locked cases recorded the
expected executable versions and hashes, and `errors.json` contained zero
records. Its overall gate remained closed only because the adversarial critic
flagged the retained M7 v2 generator's simple topology, single-layer geometry,
and v1-core reuse. This is a generator-quality finding, not an executable
configuration failure. The v2 result is retained as a failure-facing historical
baseline rather than an active required family.

Three guardrails were also made fail-closed: generator-critic findings now
close the execution gate, nested and explicit latent-age fields are rejected
from blind rows, and negative tracer/detection-limit parsing is checked
numerically. The full regression suite passes (`895 passed, 2 skipped`).

### 8.5 Observation-model and held-out-mechanism implementation checkpoint — 2026-08-01

The following software components are now authoritative for the synthetic
programme:

- `hydrosheaf/validation/observation_stress.py` defines reproducible blind
  complete, structured-missingness, left-censoring, and combined scenarios.
- `hydrosheaf/validation/head_evidence.py` applies declared GLS covariance or
  primary-channel-plus-discrepancy handling and prevents double counting.
- `scripts/independent_mixing_generator.py` supplies the held-out mechanism
  family without importing HydroSheaf or either earlier generator.
- `hydrosheaf/models/redox.py` and `hydrosheaf/models/reactions.py` now treat
  missing diagnostic ions conservatively instead of comparing `None` values.

The historical three-family run is reproducible under its recorded source
hashes and remains `PASS_WITH_CRITIC_FINDINGS` only because the M7 v2 baseline
is retained for comparison. It is not a locked scientific performance pass;
the newer scoped comparator/policy checkpoint is recorded in section 8.4.

### 8.6 Active multilayer MODFLOW/MODPATH generator checkpoint — 2026-08-01

The active MODFLOW/MODPATH family is now
`M7/m7_nonuniqueness_benchmark/scripts/independent_modflow_generator_v3.py`.
It is a separate source module: it imports neither HydroSheaf nor the v1/v2
MODFLOW generators. Its fixed synthetic contract includes:

- a two-layer heterogeneous MODFLOW 6 grid with layer-specific hydraulic
  properties and vertical conductance;
- four MODPATH 7 particle releases in both layers;
- six unique observed nodes formed by a predeclared particle-mixture mapping
  with two branch nodes and three merge nodes;
- explicit particle-release weights, component ages, and a weighted-mean scalar
  age target for the current scalar age scorer; and
- solver executable paths, hashes, versions, model dimensions, and topology
  provenance in the sealed generator record.

The adversarial generator tests in
`tests/test_independent_modflow_v3.py` verify source independence, two-layer
observations, duplicate-free coordinates, branch/merge truth, pathline-edge
support, age/head direction, mixture-weight normalization, and exact same-seed
reproducibility.

The immutable full lock
`RUN-PROGRAMME-ENSEMBLE-CALIBRATION-LOCK-20260801-02` completed all three
active families with two development seeds and one locked seed per family at
`age_samples=600`. Its 99 case-metric rows and 87 discrepancy reports were
written with `errors.json` empty; the MODFLOW/MODPATH executable gate passed;
all nine generator-critic records passed with zero major or blocker findings;
and both the execution and integrated performance gates passed. The lock
manifest records the run ID, source hashes, seed registries, sample count,
solver provenance, and per-artifact SHA-256 hashes; an independent audit found
all recorded hashes consistent. This remains a controlled-synthetic
software-validation result; it is not field validation or evidence of
universal HydroSheaf superiority.

### 8.7 Locked-result adjudication — 2026-08-01

The immutable lock was adjudicated by
`scripts/adjudicate_locked_ensemble.py`, which reads the lock without writing
into it and emits a separate JSON/Markdown record under
`.codex_work/programme-validation-adjudication/`.

- The execution verdict remained **PASS** and all input source and artifact
  hashes remained consistent.
- Conditional topology comparison was **PASS** for the analytic-lattice and
  independent-mixing cases because HydroSheaf tied the strongest registered
  baseline, but **FAIL** for MODFLOW v3 (selected F1 0.700 versus 0.737).
- Pooled calibrated age coverage was 0.966 against the declared 0.95 target;
  the independent-mixing family was below target at 0.909. No age specialist
  comparator is registered, so the comparative age claim is **ABSTAIN**.
- Reaction-family accuracy was failure-facing in all three locked cases
  (0.167, 0.167, and 0.143), and no reaction-family specialist comparator is
  registered.
- All locked complete-case next-measurement policies returned `ABSTAIN`; no
  prospective measurement outcome was simulated.

Accordingly, the system-level superiority claim is **ABSTAIN**. The current
baseline comparison is conditional rather than end-to-end because the runner
supplies the baselines with HydroSheaf's candidate universe and derived local
support. This is a candid synthetic adjudication, not field validation or
evidence of universal HydroSheaf superiority.

### 8.8 Specialist baseline extension — 2026-08-01

The baseline contract now covers the age and reaction output spaces without
pretending that they are edge-probability outputs:

- `hydrosheaf/validation/specialist_baselines.py` defines a fixed local
  tracer-decay age comparator and a fixed local concentration-difference
  reaction-family comparator.
- Both comparators consume only declared observed channels, record input,
  tuning, uncertainty, abstention, cost, and fingerprint metadata, and run
  through the same truth-blind row boundary as the topology baselines.
- The age comparator uses fixed modern activities and physical half-lives,
  reports an interval, and abstains when no age tracer is available.
- The reaction comparator emits a probability distribution over declared
  families and abstains when paired chemistry is missing or the maximum rule
  probability is below the fixed threshold.
- `run_programme_generator_ensemble.py` now records and truth-scores these
  specialist outputs separately from topology scores. The diagnostic smoke
  `RUN-PROGRAMME-SPECIALIST-BASELINES-SMOKE` completed 44 records with no
  execution errors and passed the execution gate.

This is progress toward WP3, not a completed fair-superiority comparator. The
age rule is not a full atmospheric-history inverse, the reaction rule is not a
PHREEQC inverse, and neither probability output is calibrated. A competence-
matched PHREEQC/thermodynamic inverse, development-only probability
calibration, and specialist-specific tuning/uncertainty audits remain required
before WP3 can pass.

### 8.9 Development-only model averaging checkpoint — 2026-08-01

The first statistically valid model-combination layer is now implemented in
`hydrosheaf/validation/model_averaging.py`:

- `DiscreteModelObservation` requires a complete model-by-target probability
  matrix and an explicit phase.
- `fit_discrete_model_weights` fits non-negative weights with a fixed floor
  using a case-blocked mean log score, so a case with more observed nodes
  cannot dominate the fit.
- `ModelWeightFit` records model IDs, prior and fitted weights, case IDs,
  iteration/convergence metadata, score rule, fit scope, and a development-data
  hash. Locked or field observations are rejected at the fit boundary.
- `apply_discrete_model_average` computes the predictive mixture while
  retaining pairwise total-variation disagreement. If the frozen threshold is
  exceeded, the mixture remains an audit object but the aggregate is explicitly
  non-reportable and returns `ABSTAIN` with `MODEL_DISAGREEMENT`.
- `score_locked_model_average` reports locked log score, Brier score, and
  disagreement rate without refitting.

The ensemble runner now applies this layer to the native HydroSheaf candidate
probabilities and the two topology baselines. The smoke
`RUN-PROGRAMME-MODEL-AVERAGING-SMOKE` fitted on 57 development targets and
scored 57 locked targets with no execution errors. This remains a conditional
diagnostic because every component shares HydroSheaf's native candidate
universe; it is not an independent topology-generation comparison, and the
weights are forecast-combination weights rather than model posterior
probabilities. The layer is not yet connected to the next-measurement utility,
and age/reaction discrepancy scenarios still require their own proper
predictive outputs and calibration.

### 8.10 Independent specialist candidate-universe checkpoint — 2026-08-01

The end-to-end comparator path now has an independent candidate generator in
`hydrosheaf/validation/specialist_candidate_generation.py`. It consumes only
blind site identifiers, geometry, and observed heads; it does not import
HydroSheaf inference, read `p_uv`/edge confidence, or receive generator truth.
The fixed k-nearest-neighbour graph orients clear head drops downstream and
retains both directions for tied or missing heads. Its input hash, algorithm
version, parameters, generated edges, and dropped-coordinate nodes are stored
in every case record.

The runner now executes the hydraulic-only and edge-local topology baselines,
the local tracer-age and reaction controls, and the stronger multi-tracer age
and stoichiometric reaction specialists on this independent universe. Age
remains a node-level comparator because it does not require candidate edges;
reaction is scored only for generated candidate edges, so omitted true edges
are counted rather than silently excluded. A separate candidate-universe
report records candidate recall and precision against truth only after the
blind inference boundary.

The `RUN-PROGRAMME-INDEPENDENT-SPECIALIST-SMOKE` run completed with no errors
and an execution gate of `true`. On its complete analytic-lattice and
independent-mixing cases, independent candidate recall was 0.92 and 0.93,
respectively. These are implementation diagnostics, not a superiority result:
the specialist outputs are still bounded diagnostic comparators, and the
candidate generator is a fixed geometry/head method rather than a full
MODFLOW/MODPATH or PHREEQC specialist implementation.

### 8.11 Stronger specialist baselines and calibration checkpoint — 2026-08-01

The age and reaction comparator layer now has explicit null controls and
stronger specialist diagnostics. `local_tracer_decay_age` and
`local_thermodynamic_reaction_rules` remain fixed controls. The new
`multitracer_atmospheric_history_age` evaluates a bounded age grid against
fixed atmospheric-history curves plus tritium/argon-39 decay, reports
posterior multimodality, and abstains when the evidence is unsupported. The
new `stoichiometric_reaction_inverse` projects observed chemistry onto a
fixed signed template library with non-negative extents, residual norms, and
an explicit zero-change/null reaction model. Neither is promoted to a
site-specific LPM/excess-air/degradation inverse or a coupled PHREEQC inverse.

The runner now fits frozen calibrators on development cases only and applies
them to locked cases without truth: case-blocked bias plus finite-sample
interval expansion for age, and case-blocked temperature scaling for reaction
probabilities. Calibrated predictions carry their fit scope, fingerprint,
calibration status, abstention reason, and uncertainty metadata. The revised
quick smoke completed with no errors and `execution_gate=true`.

The final three-family run
`RUN-PROGRAMME-ENSEMBLE-SPECIALIST-CALIBRATION-20260801-06` completed all 12
development and 6 locked cases, including the
MODFLOW/MODPATH executable path, and passed both `execution_gate` and
`integrated_performance_gate` with zero recorded runner errors. Its separate
adjudication also passed input-integrity checks, while retaining
`system_level_claim_status=ABSTAIN`. The locked specialist record reports
calibrated age/reaction outputs, but the reaction calibrator reached its
declared grid boundary and selected zero reaction edges at the fixed 0.60
threshold; this is a warning and selective-risk result, not evidence of
reaction-specialist success.

This is a calibration implementation checkpoint, not a performance pass. In
the small two-family smoke, the reaction temperature fit was deliberately
conservative and the 0.60 decision threshold produced zero selected reaction
edges. That result is retained as an auditable selective-risk failure signal;
it is not hidden by relabelling abstentions as correct predictions. Larger
generator/family holdouts and additional independent mechanisms are still
required before generalising coverage, calibration, or specialist ranking.

### 8.12 Bounded synthetic topology component claim — 2026-08-01

The full lock
`RUN-PROGRAMME-ENSEMBLE-BOUNDED-COMPONENT-FINAL-20260801-01` now separates
execution, component performance, integrated performance, and field-validation
status. The topology component is adjudicated on a truth-blind common
all-pairs selection universe (`common_all_pairs_v1`) against one fixed,
predeclared hydraulic-only directional-gradient comparator. Candidate-generation
recall is retained as a separate metric and is not folded into the selection
F1 comparison.

The six locked cases contain two cases from each required generator family. The
family-stratified HydroSheaf-minus-comparator selection-F1 margins are:

| Generator family | HydroSheaf F1 | Fixed comparator F1 | Margin |
|---|---:|---:|---:|
| analytic-lattice v1 | 0.558 | 0.329 | +0.229 |
| independent-mixing v1 | 0.585 | 0.406 | +0.180 |
| MODFLOW/MODPATH v3 | 0.700 | 0.636 | +0.064 |

The weakest case-level margin is +0.064 against the frozen non-inferiority
margin of -0.05. Execution, source hashes, artifact hashes, and recorded
errors passed. The independent candidate generator recovered 0.923, 0.929,
and 1.000 of reference edges in the three families, respectively; these are
candidate-generation diagnostics, not evidence that the specialist graph is
perfect.

Accordingly, the admissible bounded claim is:

> On the preregistered six-case controlled-synthetic lock, HydroSheaf was
> non-inferior to the fixed hydraulic-only directional-gradient comparator for
> topology-selection F1 on a common all-pairs evaluation universe, with two
> cases in each of the analytic-lattice, independent-mixing, and MODFLOW/MODPATH
> v3 generator families and a frozen per-case margin of -0.05.

This is a topology-selection component claim only. It does not claim age or
reaction superiority, integrated decision utility, field validity, or
universal HydroSheaf superiority. The synthetic integrated claim remains
**ABSTAIN**: all six locked policies abstained before a measurement, so no
realised post-measurement improvement was demonstrated, and the full model-
averaging fit did not converge. Field validation remains **DEFERRED**.

## 9. What cannot be completed without new data

The repository cannot create independent field truth. L4 requires a
prospective campaign in which HydroSheaf predictions are issued before new
measurements are taken. The minimum useful campaign should pair hydraulic
heads and pumping metadata with screen intervals, environmental age tracers,
major ions and isotopes, repeated sampling, and independent connectivity or
reaction evidence.

Until L4 exists, the strongest current positive claim is the bounded topology
component claim in section 8.12. The complete HydroSheaf workflow remains
unadjudicated for integrated decision performance; the next-measurement
simulator is implemented and records abstentions as unevaluated, but a positive
integrated claim requires selected actions with realised improvement under the
declared outcome scenarios.

## 10. Technical remediation work packages

The missing technical capabilities are converted into the following bounded
work packages. A work package is not complete because code exists; its gate
requires an independent test and an evidence-bearing result.

### WP1: integrated end-to-end benchmark

Create one orchestrator that takes a blind synthetic case through candidate
graph generation, tracer-age inference, chemistry/reaction inference,
uncertainty aggregation, discrepancy handling, and next-measurement selection.
The orchestrator must write stage status, predictions, compatible sets,
abstentions, action logs, and truth-blind run records before the truth files
are opened for scoring.

**Pass condition:** every required stage completes on every locked case; no
truth field appears in inference inputs; all methods use the same case and
observation budget; all registered outputs are finite and reproducible.

### WP2: independent generator ensemble

Use the historical v2 generator as a bounded diagnostic baseline, and use the
active v3 generator plus independent families that do not reuse its process
equations, geometry, or random-stream design. A held-out mechanism family and a
deliberately misspecified observation model are minimum stress tests. The v2
generator is independent of HydroSheaf, but v1 and v2 are not sufficient as two
independent generator families because v2 preserves the v1 core.

**Pass condition:** the inference code is not imported by either generator;
development, tuning, and locked cases are disjoint; results are reported by
generator family rather than pooled into one headline number.

**Current implementation:** execution-complete but not a superiority pass. The
M7 v3, analytic-lattice, and independent-mixing families satisfy the
independent-import, disjoint development/locked split, family-specific
reporting, held-out mechanism, and observation-stress requirements. The active
generator critic gate is clean. The historical v2 baseline remains available
for comparison and is intentionally not part of the required active-family
gate because its v1-core reuse and simple topology are known limitations.

### WP3: fair specialist baseline registry

Define a machine-readable baseline contract for hydraulic-only MODFLOW/MODPATH
information, edge-local graph inference, local tracer-age inference,
conventional PHREEQC or thermodynamic inversion, HydroSheaf-Core, and the full
HydroSheaf workflow. Record input channels, tuning data, parameter budget,
convergence rule, uncertainty output, abstention rule, runtime, and cost.

**Pass condition:** each comparator receives identical observations, missingness,
uncertainty, action budget, and development allowance. A baseline cannot be
weakened merely to make HydroSheaf look better.

**Current implementation:** partial but materially strengthened. The registry
contains the hydraulic-only directional-gradient and edge-local topology
comparators, local controls, an independent bounded exponential/gamma age
specialist with fixed tracer histories, and an independent constrained
reaction-family candidate generator with an explicit null model. All record
fixed tuning, uncertainty, abstention, cost, input channels, and fingerprints.
A full local atmospheric-history/excess-air/degradation inverse and a coupled
PHREEQC/thermodynamic inverse remain outside the present scope; therefore this
work package is not a complete specialist-superiority gate.

### WP4: calibrated uncertainty and selective risk

Define one output schema for probability, interval, compatible-set,
identifiability status, calibration status, and abstention reason. Add
reliability diagrams, Brier score, log loss, coverage, interval width,
selective-risk curves, false-commitment rate, and risk conditional on
`ACTION`, `SET_REPORT`, and `ABSTAIN`.

**Current implementation:** partial. Development-only age interval fitting,
case-blocked reaction temperature scaling, held-out interval/probability
scoring, uncertainty-ranked selective-risk curves, and diagnostic
topology/baseline Brier-log-loss-ECE records are implemented. Topology
probabilities are not yet development-calibrated, reaction selection remains
conservative in the quick diagnostic, and decision-conditional risk by all
three programme outputs remains to be added.

**Pass condition:** calibration is fitted only on development cases; locked
cases evaluate coverage and selective risk without retuning. A method that is
accurate only after silently excluding difficult cases fails the gate.

**Reaction RAPM addition:** the programme now includes a separate
development-only regularized adjusted reaction-family model. It uses a frozen
positive chemistry allowlist, case-blocked ridge/on-off-weight selection,
cross-fitted development logits for temperature calibration, an explicit
`mixing` class distinct from `none`, and the independent geometry/head
candidate-edge universe. Locked scoring includes decoy edges, unconditional
error, selective risk, false-commitment rate, class-wise coverage, and rank,
coherence, and parameter-stability diagnostics. The on/off term is an
adjusted evidence diagnostic, not a causal reaction effect; the model is not
a PHREEQC inverse and cannot establish unique stoichiometric identification.

The first full RAPM run is execution-complete but remains a bounded diagnostic:
the calibrated layer reduced multiclass log loss relative to the fixed
reaction comparators on the locked candidate universe, while accepted
coverage remained low and the predeclared reaction-superiority claim gate was
not passed. This result supports further development and honest abstention,
not a claim that HydroSheaf is already a strong reaction-inference engine.

### WP5: explicit model discrepancy and model averaging

Implement a common scenario interface for topology error, age-model mismatch,
reaction misspecification, shared tracer nuisance, and structured observation
error. Estimate model or scenario weights only on development data using a
predeclared proper-scoring rule. Propagate scenario-specific predictions into
the decision utility.

When model predictions materially disagree, emit `MODEL_DISAGREEMENT`; return a
compatible set or abstain. The averaging layer must expose disagreement and
cannot turn incompatible models into a falsely precise mean.

**Current implementation:** partial. The root-level runner provides a reusable
scenario record and conservative envelope for the declared nominal versus
tropical input-history pair, plus structured observation stress and
head-channel discrepancy handling. It deliberately uses fixed equal weights
for that age-history pair and labels the interval as non-calibrated. A separate
development-only, case-blocked log-score layer now combines the conditional
topology probability outputs and exposes model disagreement, but it does not
yet average calibrated age or reaction predictive distributions. The
three-family run has development-only HydroSheaf age calibration and
selective-risk scoring, and the specialist age/reaction layer now has its own
development-only calibrators. Independent candidate-universe comparison is
implemented as a separate diagnostic, but proper scenario-family weights,
topology predictive calibration, calibrated age/reaction distribution
averaging, and connection of the final mixture to decision utility are still
required for the full WP5 gate.
The new one-step decision utility is propagated from declared
nominal/discrepant scenarios, but its post-measurement outcome utility is
explicitly deferred.

**Pass condition:** the model-discrepancy layer improves or preserves locked
calibration and selective risk under stress scenarios. If it only improves a
point metric while worsening calibration or false commitment, the result is
negative or mixed.

### WP6: prospective next-measurement validation

Extend the new one-step policy from a declared synthetic edge-support model to
a simulated joint action space containing age, connectivity, and reaction
measurements. Each strategy must choose using observations available before
the hidden outcome is revealed. The simulator then reveals the outcome,
updates the posterior, records cost, and measures realised uncertainty or
decision-risk reduction.

**Current implementation:** partial. Truth-blind cost-adjusted robust action
selection and explicit abstention are implemented. The programme now simulates
declared modality-specific outcomes after the blind policy decision, updates a
binary posterior, and records expected post-measurement entropy reduction under
nominal and discrepant scenarios. Abstention is recorded as unevaluated, not as
zero risk or success. The current lock still has no selected actions, so the
integrated decision-performance gate remains abstained.

**Pass condition:** HydroSheaf improves the predeclared cost-adjusted primary
metric against random acquisition and is non-inferior to the strongest
specialist policy. The result remains controlled-synthetic until a field
campaign is run.

### WP7: clean reproducibility suite and stable release

Create one authoritative Track-A package for the programme-level benchmark:
canonical input manifests, dependency lock, generator and executable hashes,
immutable run IDs, complete failed-run logs, artifact hashes, environment
records, and a restoration script. Repair missing legacy fixtures and separate
environment records when different milestones were produced under different
software versions.

The final release check must include the full test suite, the integrated smoke
run, locked benchmark regeneration, source-hash verification, and byte-stable
artifact comparison. The current M7/M8 subpackage provenance is useful but does
not yet constitute a clean root-level programme release.

**Pass condition:** a clean checkout can regenerate the locked results without
undocumented workstation state, and the release snapshot identifies every run,
input, environment, and output hash.

### WP8: prospective field validation

Pre-register one field decision problem, issue HydroSheaf and comparator
predictions before sampling, record the proposed actions and costs, and collect
independent follow-up observations. The field protocol must include repeated
heads and sampling dates, pumping conditions, screen intervals, environmental
age tracers, chemistry, and independent path or reaction evidence.

**Pass condition:** the field outcome is scored against measurements that were
not used to construct the initial graph, tune the model, or select the action.
This work package is externally blocked until the data and campaign exist.

## 11. Recommended execution order

The safest order is:

1. WP7 minimum provenance and baseline registry;
2. WP1 integrated orchestrator using the active v3 generator;
3. WP4 common uncertainty and selective-risk outputs;
4. WP5 discrepancy and model-averaging layer;
5. WP2 second generator family and held-out stress cases;
6. WP6 joint next-measurement benchmark;
7. WP8 prospective field validation.

No universal-superiority manuscript claim should be drafted before WP1-WP5
have passed. If the joint gate fails, the publication should report the
conditional strengths and failure modes rather than expanding the model until
the claim passes.
