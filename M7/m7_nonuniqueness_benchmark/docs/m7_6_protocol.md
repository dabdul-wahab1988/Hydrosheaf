# M7.6 multi-stream integration and shared-nuisance benchmark — locked protocol

**Status:** declared, not yet executed. No M7.6 result exists at the time of writing.

## Relationship to M7.3

M7.3 is locked and its results stand unchanged. M7.6 does not revise, re-run or
reinterpret any M7.3 number. It adds two evidence streams that M7.3 did not
contain and one factor that M7.3 did not vary, and it uses fresh seeds so the
two benchmarks never share a case.

M7.3 established, in its own generator, that chemistry contributed strongly to
edge ranking (PR-AUC +0.447), hydraulics weakly but calibratedly (+0.0091), age
not at all (−0.0060, interval excluding zero), that a correct graph improved age
estimation while a reversed graph harmed it, and that permuted streams reduced
entropy while worsening every predictive metric.

Three limitations motivate M7.6:

1. **Only three streams.** Environmental isotopes were absent entirely, and the
   nuclear panel was a two-tracer surrogate (3H, 39Ar) that omits every tracer
   a practitioner normally holds — 14C, CFCs, SF6, 4He, 3H/3He.
2. **One inference target.** "Does stream X help?" was answered only for edge
   ranking. M7.3's own results show the answer is target-dependent: age added
   nothing to topology ranking while topology improved age. That dependence was
   incidental rather than designed.
3. **Negative transfer only by permutation.** Every harmful case in M7.3 was
   produced by permuting a stream, which is adversarial by construction. Real
   integration failures arise from streams that share a latent nuisance and are
   therefore wrong *together* and in a correlated way. M7.3 could not produce
   that failure mode.

## Thesis question

When does adding an evidence stream improve inference, add nothing, or actively
weaken it — and does the answer depend on which quantity is being inferred?

M7.6 additionally asks whether a shared, unmodelled nuisance affecting a family
of tracers produces detectable incompatibility, or silently produces confident
error.

## Independence and analysis split

The generator continues to import no HydroSheaf code and continues to execute
official MODFLOW 6 and MODPATH 7 binaries. The new isotope and nuclear
formulations are independent analytic constructions, not HydroSheaf's
input-history or correction routines.

- Development seeds: **5401–5406**.
- Locked test seeds: **5501–5512**.
- Non-claim-bearing smoke seeds: development **5951–5952**, test **5961–5963**.
- Seeds do not overlap M7.3 (5201–5206, 5301–5312) or its smoke seeds.
- All fusion coefficients, thresholds and stopping rules are fitted on
  development cases only.
- Test truth is joined only after truth-blind inference.
- No test result may change a feature, coefficient, condition, threshold or
  stopping rule.

New randomness is drawn from a generator stream independent of the one used by
M7.3, so that M7.3 cases remain bit-reproducible. This is a release criterion,
not an aspiration, and is asserted in the generator's tests.

## Evidence streams

| Code | Stream | Contents |
|---|---|---|
| `H` | hydraulic | heads, gradients |
| `C` | hydrochemistry | 12 major ions, pH, temperature |
| `N` | nuclear | 14C, CFC-11/12/113, SF6, 4He, 3H/3He, 3H, 39Ar |
| `E` | environmental isotopes | d18O, d2H, 87Sr/86Sr |

`d13C` is deliberately **not** assigned to a single stream. It is coupled by
construction to the carbonate reaction path (`C`) and to the 14C initial
activity (`N`). It is the mechanism by which two nominally independent streams
share a latent nuisance, and it is the object of Experiment 3.

## Inference targets

The factor M7.3 lacked. All stream combinations are evaluated against each:

- `T1` **topology**: edge-existence ranking. M7.3's target; retained for
  continuity.
- `T2` **age**: groundwater residence time at each node.
- `T3` **mixing**: recharge-source fraction, recoverable from conservative
  tracers.

## Predeclared predictions

Stated before execution. Their purpose is falsifiability: a result contradicting
a prediction is informative, and a pipeline that "improves" on a target where a
stream carries no information by construction indicates a defect, not a finding.

| Stream | T1 topology | T2 age | T3 mixing |
|---|---|---|---|
| `E` environmental isotopes | little | **none by construction** | strong |
| `N` nuclear | little | strong | little |
| `C` chemistry | strong | moderate | moderate |
| `H` hydraulic | small, calibrated | small | little |

The `E`/`T2` cell is the strongest prediction in this protocol. d18O and d2H are
generated at recharge and transported conservatively, with no age dependence
whatsoever. Any reported improvement in age inference from adding `E` alone is
therefore attributable to leakage, overfitting or a coding defect, and must be
investigated as such before any other M7.6 result is reported.

This is a **physically justified negative control**, and it is stronger than
M7.3's permutation controls: it destroys no structure and requires no
adversarial manipulation. The stream is genuinely informative — about something
else.

## Experiment 1: stream contribution by target

All 15 non-empty subsets of {H, C, N, E} are fitted on development cases and
evaluated on locked test cases, for each of T1, T2, T3.

Metrics follow M7.3 so results remain comparable:

- **T1**: PR-AUC, ROC-AUC, Brier, log loss, normalized edgewise Bernoulli
  entropy, fraction of edges in [0.1, 0.9], expected edge count, calibration
  error, overconfident error.
- **T2**: MAE, interval width, empirical coverage, ESS-based convergence rule.
- **T3**: mixing-fraction MAE, interval width, coverage.

Case-block bootstrap intervals, 10,000 replicates.

**M7.3's guardrail is retained and is binding.** An entropy or width reduction
counts as beneficial only if paired Brier score and log loss (T1) or MAE and
coverage (T2, T3) do not materially worsen. Confidently wrong output is negative
transfer, not resolution.

## Experiment 2: redundancy versus complementarity

Streams are evaluated for whether they constrain overlapping or disjoint
subspaces of each target. For each pair, report the incremental contribution of
the second stream given the first, in both orders. Asymmetry is expected and is
the signature of complementarity; near-zero increments in both directions
indicate redundancy.

Redundancy is reported as a distinct outcome from harm. M7.3 conflated them
under "no positive incremental value".

## Experiment 3: shared nuisance and correlated failure

The experiment M7.3 could not run, and the one connecting M7 to M3.

A per-case systematic is applied to a whole tracer family at once, at three
predeclared levels — `none`, `mild`, `severe`:

- `recharge_temperature_error`: perturbs CFC-11, CFC-12, CFC-113 and SF6
  together through their solubility conversion, i.e. common-mode within the CFC
  family.
- `c14_a0_error`: perturbs the 14C initial activity through the d13C
  relationship, coupling the nuclear stream to the chemistry stream.

The true perturbation is retained in ground truth and is never exposed to
inference.

Reported per level:

1. accuracy and calibration on T2 and T1;
2. whether entropy or interval width falls while error rises — the false
   confidence signature;
3. **the rate at which the local tracer constraint set becomes infeasible**;
4. whether the affected family is identifiable as the source, e.g. by
   leave-one-family-out reconciliation.

### Connection to M3

M3's national USGS benchmark found that 27.85% of eligible folds admit no
feasible transit-time distribution, that the failure is 83% pairwise, and that
it concentrates on tracers sharing a correction — CFC-CFC pairs at 19–33%, 3H
versus 3H/3He at 17.9%, against 0.8% for the independent 14C/3H pair. Seven
candidate explanations were tested and none accounted for it. M3 cannot resolve
this because it has no ground truth.

M7.6 has ground truth. If the `severe` level reproduces M3's signature — a
pairwise-dominated infeasibility concentrated within the perturbed family, at a
comparable rate — then the mechanism is demonstrated and M3's observation has a
validated explanation. If it does not reproduce it, shared correction error is
excluded as the explanation and M3's cause lies elsewhere.

**This is a diagnostic experiment, not a claim about the USGS data.** A positive
result licenses the statement that shared correction error *can* produce the
observed signature, never that it *did* in Ghana or in the USGS release.

## Experiment 4: adverse controls

Retained from M7.3 and extended:

1. `native`;
2. `permuted` — each stream permuted within case, singly and jointly;
3. `identity_relation` — the graph-smoothing baseline;
4. **`isotope_age_control`** — is `E` reported as improving `T2`? It must not
   be. A positive result here fails the run and blocks all other M7.6 reporting
   until resolved.

## Decision rules and guardrails

A stream is recorded as **improving** a target only if a predictive metric
improves with a bootstrap interval excluding zero, and no calibration metric
materially worsens.

A stream is recorded as **adding no value** if incremental effects are
indistinguishable from zero in both orders of Experiment 2.

A stream is recorded as **weakening** if any predictive or calibration metric
worsens with an interval excluding zero, *or* if uncertainty falls while error
rises.

No claim of universal superiority for any stream, framework or integration
strategy may be drawn from this benchmark. Synthetic truth is
model-conditioned. It is not evidence about any real aquifer, and specifically
not about the Northern Ghana workbook or the USGS release.

## Execution order

1. Commit this protocol and the executable runner **before** any 5501-series
   test case is generated. M7.3 set this precedent with `d336e87`; M7.6 must
   match it or it is not confirmatory.
2. Generate and analyse development cases (5401–5406). Fit all coefficients.
3. Freeze. Any change after this point invalidates the run.
4. Generate locked test cases (5501–5512), infer truth-blind, then join truth.
5. Write results, then the claim decision.
6. Only then update the manuscript.

Results written before step 4 completes are not confirmatory and must not be
described as such.

The pre-test implementation choices and M3-specific feasibility diagnostic are
frozen in [docs/m7_6_execution_amendment_20260731.md](m7_6_execution_amendment_20260731.md).
