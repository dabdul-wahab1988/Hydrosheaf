# Age evidence and direct-adjacency benchmark protocol

**Protocol identifier:** `age-adjacency-v1`  
**Purpose:** distinguish what an age observation can identify from what must be
supplied by topology or transport information. This protocol is the evaluation
design for the available-data and controlled-synthetic panels.

## Why this protocol is needed

An age difference is naturally evidence about temporal order and, with a
transport model, about a travel-time increment.  It is not by itself evidence
that two observed nodes are consecutive in a graph.  If `u -> w -> v` and
`u -> v` have compatible endpoint ages, the same age ordering can support both
the two-hop path and the direct edge.  Direct adjacency is a **cover-relation**
estimand; reachability is a path-existence estimand.  They must therefore be
labelled and scored separately.

The existing direction-gate result remains valid for its stated estimand.  It
should be described as a test of downstream-older compatibility, not as a
test of direct-edge identification.  The benchmark below adds the missing
estimand and makes an explicit abstention path available when the metadata
cannot identify it.

## Estimands and labels

Every candidate ordered pair receives one evaluation relation label from the
truth graph.  Labels are joined only after prediction and are never model
features.

| Estimand | Positive relation | What age can contribute | Required caution |
| --- | --- | --- | --- |
| Direction | `downstream_older` | A one-sided probability that the proposed downstream node is older | It does not establish a path or an edge |
| Reachability | `transitive_reachable` or `direct_adjacent` | Compatibility of an accumulated age increment with a path | A multi-hop path can have the same endpoint ages as a direct edge |
| Direct adjacency | `direct_adjacent` only | A direct-segment travel-time likelihood, if edge-specific transport information is available | Without an independent segment/RTD model, the result is `ABSTAIN`, not a negative edge |

The disjoint directed relation labels are:

1. `direct_adjacent`: `(u, v)` is a truth edge;
2. `transitive_reachable`: `v` is reachable from `u` by a directed path of
   length at least two and `(u, v)` is not a truth edge;
3. `reverse_incompatible`: `u` is reachable from `v`;
4. `cross_path`: path metadata identify different generating trajectories;
5. `unrelated`: neither direction is reachable and the endpoints are on the
   same supplied trajectory (or trajectory identity is not available).

If the generator has connected-component or flow-path annotations, an
additional `cross_component` annotation may be retained, but it must not
replace the directed relation label.  A relation is missing/unknown whenever
the complete truth graph is unavailable.

For each case let `C` be the generated candidate set, `E` the complete direct
edge set, and `R` the complete directed reachability set.  Report both
candidate coverage and conditional detection:

\[
  \mathrm{candidate\ coverage}_E = |E \cap C|/|E|,
  \qquad
  \mathrm{conditional\ detection}_E = TP/|E \cap C|,
\]

and, when `E` is complete,

\[
  \mathrm{all\ truth\ detection}_E = TP/|E|.
\]

The same denominators are reported for `R`.  A metric calculated only on `C`
must be named `candidate-contained`; it must never be presented as recall
over the complete graph.

## Phased implementation and execution

### Phase 0 — freeze the audit boundary

Before changing a locked result, record a new run identifier, the protocol
hash, generator hash, inference-code hash, environment versions, random seeds,
and input hashes.  Existing locked outputs are immutable.  The current
supporting-validation directory is an input/reference artifact; a new
age-adjacency run writes to its own directory.

Freeze these decisions before inspecting locked-test outcomes:

- direct adjacency is the primary new estimand;
- direction and reachability are secondary, separately scored estimands;
- the candidate set is fixed before labels are joined;
- model selection, calibration, and threshold selection use development cases
  only;
- case-level paired resampling is used for uncertainty;
- missing edge-specific travel information produces `ABSTAIN`.

### Phase 1 — re-score the existing result table by relation

This phase is an evaluation-only decomposition of the existing predictions.
It does not refit the direction gate and does not alter its locked values.
Generate relation labels from the independent generator truth and calculate,
for each arm:

- direction compatibility metrics;
- reachability metrics where the complete closure is available;
- direct-edge metrics restricted to `direct_adjacent` versus the explicitly
  chosen negative relations;
- candidate coverage, candidate-contained detection, and all-truth detection;
- per-case results before any macro/micro aggregation.

The comparison must include the frozen baseline, the age condition, and an
age-permuted control with the same candidate rows.  A lower direct-adjacency
score in this phase is evidence that the direction gate is not a cover-relation
detector; it is not evidence that age observations are incorrect.

Applied to the current 12-case locked test table (825 candidate pairs; 108
complete direct truths; source-table SHA-256
`7dd10a2afe75a278e761efb108fa97358ba73187811d401c1c2f1ca9a3d5955e`), the
post-inference audit gives direction consistency 1.00, direct-adjacency
ROC-AUC 0.437 and PR-AUC 0.104, and false two-hop rejection 0.00.  Candidate
coverage is 0.954 for direct edges and 0.977 for all directed reachable pairs.
These figures are a diagnostic of the historical direction-compatible score,
not a new calibrated directness result.  Because the locked arm has no
separate reachability probability, its reachability ranking/calibration
metrics are intentionally undefined; only the truth-closure coverage is
reported until a dedicated reachability arm is supplied.

### Phase 2 — add an edge-specific direct-versus-indirect likelihood

For a candidate `u -> v`, retain the observed age difference

\[
  \Delta a = a_v-a_u,
  \qquad
  \sigma_\Delta^2 = \sigma_u^2+\sigma_v^2-2\,\mathrm{Cov}(a_u,a_v)
  +\sigma_\mathrm{process}^2,
\]

when endpoint age estimates are correlated.  A
direct hypothesis requires a segment-specific travel-time distribution
`T_direct`; for example, a path-length and velocity distribution with an RTD
dispersion.  An indirect hypothesis requires a specified multi-segment or
intermediate-node distribution `T_indirect`.  Score the two hypotheses using
the same observed `Delta a`, and report

\[
  \log BF_{direct:indirect} =
  \log p(\Delta a\mid T_\mathrm{direct})-
  \log p(\Delta a\mid T_\mathrm{indirect}).
\]

The implementation must expose the component terms, uncertainty, and a
status.  `ABSTAIN` is required when path length, velocity/RTD information, or
the intermediate-path alternative is absent.  Replacing absent information by
an arbitrary default and calling the result a directness probability is not a
valid identification strategy.  The likelihood is an edge score; it must not
be substituted for the one-sided direction gate without an explicit arm
definition.

The repository implementation is opt-in.  The pure functions
`hydrosheaf.validation.age_adjacency.compute_direction_evidence` and
`compute_age_adjacency_evidence` never infer distance, velocity, paths, or
truth labels.  `Config(sheaf_age_adjacency_enabled=True)` enables the optional
edge term in `refine_edges_with_sheaf`; an edge must provide an independent
`indirect_travel_years` (or equivalent) attribute before a direct-versus-
indirect Bayes factor can affect its score.  With missing information the
edge is annotated as `ambiguous` or `insufficient_information` and receives no
adjacency cost.  Optional `direct_travel_sigma_years`,
`indirect_travel_sigma_years`, and `age_covariance_years2` attributes preserve
RTD dispersion and correlated endpoint-error information.  The default is
disabled, so existing locked runs are unchanged.

### Phase 3 — independent controlled-synthetic validation

Use a generator that is independent of HydroSheaf and has truth for both the
direct graph and the complete directed closure.  Keep development and locked
test cases disjoint.  Vary the following prospectively:

- spatially variable velocity and edge lengths;
- heteroskedastic age error and correlated endpoint error;
- non-zero RTD dispersion, skew/tails, and mixing;
- censored or interval age observations;
- branching and converging paths, including deliberate two- and three-hop
  skips;
- candidate-radius and candidate-recall regimes.

The generator must not choose a favourable velocity or RTD after observing
locked predictions.  A replicate may be scored only if its metadata and truth
graph pass the provenance checks.  Results are model-conditioned synthetic
evidence and must remain labelled as such.

### Phase 4 — outside the current execution scope

This phase is retained as protocol history but is not an active work item. The
current package works with the supplied Ghana panels, Aiken model reference,
and controlled synthetic generator; no external-data acquisition step is part
of the present run.

## Required result metadata

Each prediction row must include `edge_id`, `u`, `v`, `case_id` (or a stable
`split` + `seed` pair), split, prediction arm, prediction score, and an
explicit status (`scored`, `ABSTAIN`, or `invalid`).  Each age record must
include the posterior/measurement mean and uncertainty, units, censoring
status, and—where applicable—the covariance with other endpoint ages.

Each directness record must additionally identify:

- edge/segment length and coordinate source;
- velocity distribution or independently estimated travel-time distribution;
- RTD family, scale, and dispersion parameters;
- process/mixing variance and any age-retardation assumptions;
- intermediate nodes or the path set used for the indirect hypothesis;
- generator version, case seed, graph hash, and code/environment hashes.

The complete truth artifact must state the number of direct edges and the
number of reachable ordered pairs, including truths absent from `C`.  A result
table containing only candidate rows cannot supply the all-truth denominator
by itself.

## Leakage and calibration controls

The following are hard gates, not optional reporting preferences:

- `is_true_edge`, relation labels, true process labels, graph closure, and
  future/intermediate truth are evaluation-only columns;
- the prediction function receives no truth-bearing column, including one
  under a renamed alias;
- the locked test is not used for feature selection, RTD/velocity selection,
  calibration, threshold selection, or stopping decisions;
- all arms use the same candidate rows and case splits;
- an age-permuted or otherwise adverse control uses the same calibration and
  scoring budget;
- any post-hoc change to the directness model starts a new run identifier and
  new locked split.

Calibration is fitted on development cases only.  Thresholds are locked
before test scoring.  Report threshold-free ranking metrics in addition to any
thresholded F1, because a threshold selected for direct edges is not
necessarily suitable for reachability or direction.

## Metrics and uncertainty

For every estimand and arm, report PR-AUC, ROC-AUC (when both classes are
present), Brier score, log loss, and a fixed-bin expected calibration error.
Report precision, recall, and F1 only with the locked threshold.  Add:

- relation-stratified confusion tables;
- the fraction and accuracy of `ABSTAIN` predictions;
- direct-edge retention and false two-/three-hop acceptance;
- candidate coverage and all-truth detection separately;
- per-case macro summaries and pooled summaries with their aggregation rule.

Uncertainty is estimated by resampling complete cases, not individual edges,
because candidates from one hydraulic/age realization are dependent.  Use a
fixed-seed paired bootstrap (or a predeclared exact paired permutation) for
arm differences and report 95% intervals.  When multiple estimands or arms
are treated as confirmatory, state the multiplicity family and use a
simultaneous interval or adjusted decision rule.  If a case has no positive
example for an estimand, record that metric as undefined for that case and
retain the case count; do not silently convert it to zero.

## Realistic transport controls

The directness model is most vulnerable when endpoint ages are almost
deterministic functions of the same hydraulic quantity used by the baseline.
The independent benchmark must therefore include controls that break this
redundancy without making the task artificial:

1. draw velocity and effective porosity independently by segment, with the
   distribution and coefficient of variation frozen in the protocol;
2. draw direct travel time from the segment RTD and indirect travel time from
   the convolution of its component segments;
3. include age observation error, mixing, and censoring after transport rather
   than treating the simulated age as exact;
4. retain intermediate nodes and score skipped pairs as reachable but not
   direct;
5. repeat at several uncertainty/dispersion levels and report the interaction
   rather than selecting the most favourable level;
6. keep a direction-only arm so any gain can be attributed to edge-specific
   timing rather than to simply increasing the age weight.

These controls test identifiability.  Increasing the age-gate weight, tuning a
cutoff on locked results, or adding a scalar age feature does not create the
missing cover-relation information.

## Runnable audit surface

The dependency-light evaluator
[`M7/m7_nonuniqueness_benchmark/scripts/age_adjacency_protocol.py`](../M7/m7_nonuniqueness_benchmark/scripts/age_adjacency_protocol.py)
checks a completed CSV and emits deterministic JSON.  It verifies required
candidate fields, finite probabilities, duplicate candidate pairs, and the
truth-bearing score-name guard; it computes ranking, calibration, threshold,
relation-count, and candidate-versus-all-truth summaries.  It is intentionally
descriptive and does not replace the preregistered case bootstrap.

Example:

```powershell
.venv\Scripts\python.exe M7\m7_nonuniqueness_benchmark\scripts\age_adjacency_protocol.py `
  --input M7\m7_nonuniqueness_benchmark\results\RUN-AGE-ADJACENCY-YYYYMMDD-01\edge_results.csv `
  --score-column probability_age_adjacency `
  --output M7\m7_nonuniqueness_benchmark\results\RUN-AGE-ADJACENCY-YYYYMMDD-01\audit.json
```

Add `--truth-count <complete-graph direct-edge count>` from the truth manifest
when it is available; otherwise interpret the all-truth recall as unavailable.
The audit must never be pointed at an existing locked directory with an
overwrite flag.

For the current locked M7 supporting-validation output, the separate
post-inference audit is run with:

```powershell
.venv\Scripts\python.exe M7\m7_nonuniqueness_benchmark\scripts\age_adjacency_locked_audit.py `
  --results-dir M7\m7_nonuniqueness_benchmark\results\supporting_validation `
  --output-dir M7\m7_nonuniqueness_benchmark\results\RUN-AGE-ADJACENCY-AUDIT-01
```

This command reads the locked table and its per-case truth/pathline files,
writes relation-labelled copies to a new directory, and evaluates the
historical age-compatibility score as a *diagnostic* against direct labels.
That score is deliberately not reported as a calibrated direct-adjacency
probability; a new run with independent segment/RTD hypotheses is required
for the Phase 2 directness arm.

## Run isolation and reproducibility

Use a new directory such as
`M7/m7_nonuniqueness_benchmark/results/RUN-AGE-ADJACENCY-YYYYMMDD-01/` for
every protocol execution.  Write at least:

- `protocol.lock.json` with this protocol identifier, settings, and hash;
- `manifest.json` with input/output SHA-256 hashes, seeds, versions, and
  generator/inference commits;
- candidate predictions and complete truth metadata in separate files;
- relation-stratified metrics and case-level metrics;
- the deterministic audit JSON.

Do not overwrite `results/supporting_validation`, `results/m7_3_locked`, or
any other locked artifact.  A rerun after viewing a locked result is a new
exploratory run unless its protocol, split, and decision gates were already
frozen.

## Evidence ceiling and interpretation

Success in Phase 2 or Phase 3 would support a bounded statement such as:
“under the specified independent transport generator, edge-specific timing
information improved discrimination of direct adjacency.”  It would not prove
that age alone identifies adjacency in aquifers.

Failure of the directness model under a generator with realistic independent
RTDs is also informative: it supports non-identifiability under those tested
conditions.  Failure when the required path/RTD metadata are absent is an
`ABSTAIN` condition, not evidence against the edge.  Field chemistry
hold-forward, component age accuracy, and controlled synthetic topology
performance must remain separate claims.  No phase of this protocol alone
establishes universal graph superiority, general field topology recovery, or
reaction-truth recovery.

## Executed implementation (2026-09-08)

The new controlled run is archived at
`.codex_work/runs/RUN-AGE-BF-CONTROLLED-20260908-03/`.  It contains 24
generated cases (12 development and 12 locked-test), 360 held-out candidate
pairs per arm, a complete truth graph, a sealed truth sidecar, the four
prespecified arms, and a passing two-tier package audit.  The complete
held-out truth has 192 direct edges and 864 reachable ordered pairs; the
candidate table contains the one-hop direct pairs and two-hop skips only, so
candidate-contained detection is not the same denominator as all-truth
detection.

For an endpoint pair, (a_u) and (a_v) are the reported ages (years) at
upstream node (u) and downstream node (v); (sigma_u) and (sigma_v)
are their one-standard-deviation uncertainties (years); and
(operatorname{Cov}(a_u,a_v)) is their error covariance (years²).  The
observed increment is (Delta a=a_v-a_u).  With process uncertainty
(sigma_p) (years), the propagated increment uncertainty is

\[
\sigma_\Delta = \sqrt{\sigma_u^2+\sigma_v^2-2\operatorname{Cov}(a_u,a_v)+\sigma_p^2}.
\]

For the controlled generator, the covariance is explicitly set to zero and
the transport hypotheses are generated independently of the *observed* ages
(they share only the latent segment travel time, as required by a physical
transport model).  (T_D) and (T_I) denote the supplied direct and indirect
travel-time means (years), and (s_D) and (s_I) their RTD/model standard
deviations (years).  The normal likelihood for hypothesis (H) is

\[
p(\Delta a\mid H)=\mathcal{N}\!\left(\Delta a;\,T_H,\sqrt{\sigma_\Delta^2+s_H^2}\right),
\]

where (mathcal{N}(x;\mu,s)) is the normal probability density at value
(x), with mean (mu) and standard deviation (s>0).  The reported
(log BF_{D:I}=\log p(\Delta a\mid T_D)-\log p(\Delta a\mid T_I)) is therefore
conditional on the two transport hypotheses; it is not an age-only adjacency
probability.

The locked result is consistent with the identifiability theory:

- order-only age compatibility had directness PR-AUC 0.4313, below the
  equal-prior baseline (0.5333 prevalence in this candidate design);
- the paired full Bayes-factor arm had PR-AUC 0.8847, with 70.6% of rows
  scorable and 29.4% explicitly abstaining;
- the overlapping-transport stratum abstained for every row, while the
  separable stratum reached PR-AUC 0.9987;
- full Bayes factors exceeded the permuted-transport control by +0.0420
  PR-AUC, but the paired 95% complete-case bootstrap interval was
  [-0.0004, +0.0839], so the global advantage is uncertain at this sample
  size; the contrast against order-only was +0.4513 [0.3669, 0.4990].

These values support only a conditional controlled-synthetic component claim:
the directness model can work when independent segment-level transport
hypotheses are sufficiently separated, and it must abstain when they overlap.
The Aiken panel is reported separately as model-conditioned transport and
age-concordance evidence.

The separated two-panel replay is archived at
`.codex_work/runs/RUN-AGE-ADJACENCY-TWO-TIER-20260908-03/`.  Its Aiken child
contains the same 61-row explicit crosswalk, 61 endpoint hypotheses, 5,261
pathline segments and 20 CFC screening intervals.  The Aiken writer hashes
each generated CSV (SHA-256 and byte count) but intentionally does not hash the
multi-gigabyte source archives or the manifest itself; its field-scoring flag
remains `false`.
