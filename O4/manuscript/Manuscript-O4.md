# Fit quality is not identifiability: a harmonized robustness, integration-value and calibration audit of groundwater inference under data limitation

Dickson Abdul-Wahab

HydroSheaf PhD Chapter 1, Objective 4. Sibling article to `O3` (Objective 3)
and companion to the HydroSheaf framework paper (Objective 2, same target
journal). Submitted to Computers & Geosciences.

## Abstract

Groundwater inverse frameworks are usually judged by whether they converge,
fit the data, and report a plausible uncertainty interval. Whether those
internally-generated signals track anything externally verifiable once
evidence is limited, an evidence stream is misspecified, or the calibrating
model's form does not match the truth-generating process is rarely asked of
the same framework twice. This paper harmonizes three already-locked
HydroSheaf benchmarks -- a five-level evidence-tier field-transfer audit
across 160 Northern Ghana wells (`M6`), a conditional evidence-integration
test against fresh independent MODFLOW 6/MODPATH 7 truth (`M7`), and three
prospectively locked controlled-synthetic calibration protocols against
known parameters (`M8`) -- under one question: does an internal confidence
signal track an external truth signal under stress? It does not, in any of
the three layers. Stripping measurement tiers leaves `M6`'s mean
mechanism-resolution score flat (68.4-71.0) while the fraction of wells
correctly flagged non-identifiable rises from 0.0% to 60.0% [95% Wilson CI
52.3, 67.3], with a 51.3 percentage-point cliff at the single step that
removes strontium/silica (0.6% [0.1, 3.5] to 51.9% [44.2, 59.5]).
Permuting an evidence stream in `M7` reduces posterior-edge entropy by more
than the native condition while simultaneously degrading ranking skill
(joint misspecification: entropy -0.071 nats, PR-AUC -0.139). The same 50 d
observation design that gives `M8`'s matched analytical model near-nominal
95% coverage (0.832-0.836) collapses to 0.02-0.004 once truth is generated
by an independently implemented numerical solver instead, even as the
dispersivity point estimate itself improves under that same independent
truth. Negative and structural
controls confirm each pattern is a property of the stress test, not an
artefact of the reporting convention: an evidence-gate-off ablation returns
0% non-identifiable at every `M6` tier, and every `M7` native condition beats
its permuted-map control. The contribution is this cross-component
demonstration -- unavailable from any one underlying benchmark alone -- that
numerical health is not evidence of identifiability, integration value, or
calibration validity, offered as a companion to the HydroSheaf framework
article rather than a fourth validation of any one component.

**Keywords:** identifiability; uncertainty calibration; robustness; data-limited aquifers; equifinality; MODFLOW; MODPATH; PEST; optimal experimental design; reproducible research software

## 1. Introduction

A groundwater inverse framework typically reports three kinds of signal
before a reader ever sees an external comparison: whether the optimizer
converged, how well the model fits the observations, and how confident its
own uncertainty interval claims to be. These signals are computed entirely
from the model and the data used to fit it. They are also, by construction,
unable to detect a class of failure that has been documented across
environmental modelling for decades: a model can fit available observations
well without its structure or parameters being uniquely determined by those
observations [@beven2001glue; @beven2006manifesto], and predictions built on
an unidentified model still carry uncertainty that a converged fit does not
reveal on its own [@neuman2003bma]. The Prediction in Ungauged Basins
programme made the same argument for data-scarce catchments specifically:
inference under data limitation requires reporting the boundary of what the
evidence supports, not only a single calibrated best estimate
[@sivapalan2003pub; @hrachowitz2013pub]. Whether a framework's own
internally-generated confidence signals actually respect that boundary, or
instead stay flat, healthy, or even improve while the boundary is being
crossed, is an empirical question that has to be asked of a specific
framework under a specific stress, not assumed from the general argument
alone.

HydroSheaf is a claim-gated evidence-integration framework whose three
inference layers -- residence time, connectivity and reaction -- were
benchmarked against external references for detection versus correctness in
a companion article (Objective 3 of this thesis chapter, `O3`)
[@hydrosheaf_o3_inprep]. That comparison left a different question open:
under stress -- not simply "is the estimate right," but "does the framework
know when it might be wrong" -- do its internally-generated signals of
confidence behave honestly? Answering that question required three further,
physically distinct stress designs, each already built, locked and reported
independently, on its own timetable, against its own kind of ground truth.
`M6` [@hydrosheaf_m6_inprep] field-transfers the reaction-inference workflow
to three real Ghanaian chemistry datasets under a five-level measurement-tier
ladder and asks whether the framework degrades gracefully as evidence is
removed. `M7` [@hydrosheaf_m7_inprep] tests conditional evidence integration
-- combining or deliberately misspecifying hydraulic, tracer and chemical
evidence streams -- against fresh, independently generated MODFLOW 6/MODPATH
7 truth, and asks when integration helps, does nothing, or actively
misleads. `M8` [@hydrosheaf_m8_inprep] runs three prospectively locked
controlled-synthetic calibration protocols with a production PEST-style
Gauss-Levenberg-Marquardt engine against known parameters, and asks whether
optimiser convergence and nominal uncertainty intervals track parameter
recovery -- the same dissociation between a well-posed calibration and a
well-identified one that motivates identifiability and null-space analysis
in the wider groundwater calibration literature [@hilltiedeman2007effective;
@dohertyhunt2009identifiability].

Each of these three benchmarks already reports its own internal-versus-
external divergence as a standalone finding: `M6`'s own headline is that
"good numerical fit and stable support persist even where identifiability
has collapsed" (`M6/docs/m6_results_summary.md`); `M7`'s own decision table
records that lower uncertainty does not always mean better integration; `M8`
titles its own contribution "calibration is not identification." No single
one of the three states this as a cross-component pattern, because none was
designed to be read against the other two, and each uses an unrelated
metric, an unrelated notion of ground truth, and an unrelated stress
mechanism -- a real field-measurement ladder for `M6`, an independent
numerical generator for `M7`, and a known synthetic parameter vector for
`M8`. That separation is a strength for each benchmark's own internal
validity. It is also the reason no reader can currently tell whether "fit
quality is not identifiability" is a `M6`-specific finding about mineral
saturation indices, an `M7`-specific finding about entropy heuristics, or
something that recurs whenever a HydroSheaf-style framework is placed under
any of these three kinds of stress.

This paper harmonizes the three already-locked benchmarks under one
taxonomy and asks a single question: across data limitation, evidence
misspecification, and a model-form shift away from the calibrating model,
does an internally-generated confidence signal track an externally-
verifiable truth signal? The answer reported here is no, consistently, but
not identically. In `M6`, an evidence-tier ablation strips real field
measurements from 160 Northern Ghana wells and shows the internal signal
(mean mechanism-resolution score) essentially flat while the external
signal (verified identifiability class) collapses at one specific,
identifiable step. In `M7`, adverse controls that deliberately misspecify
an evidence stream reduce the internal signal (posterior-edge entropy) by
more, not less, than the honest native condition, while the external signal
(ranking skill against independent truth) moves in the opposite direction.
In `M8`, a shift from the calibrating model's own form to an independently
generated numerical truth shows the internal signal (nominal interval
coverage) collapsing for a parameter whose external signal (point-estimate
recovery) simultaneously improves. Three further objectives follow: to
report the negative and structural controls that bound each pattern, so it
is shown to be a property of the stress test and not a reporting artefact;
to compare the three benchmarks' design scale, since a benchmark with 8,500
calibrations and a benchmark with 800 well-tier rows are not making
comparably powered claims by construction; and to state plainly that the
three layers' field-transfer scope is uneven, `M8` in particular having
none.

The paper does not introduce new software architecture, does not revisit
HydroSheaf's integrated field demonstration, and does not claim a fourth,
independent validation of any component beyond what each component's own
result package already supports. It is a benchmark-synthesis article: its
contribution is the cross-component demonstration itself, the taxonomy that
makes the comparison possible, and the explicit statement -- available from
no single underlying paper -- that an internal confidence signal not
tracking an external truth signal is not a quirk of one benchmark's design
but a recurring property of this framework under three unrelated kinds of
stress.

## 2. Data

This paper introduces no new dataset. It reports a comparison built on
already-existing datasets and already-generated synthetic designs, each
already used as the stress mechanism or the truth reference for exactly one
of the three component benchmarks. This section describes each in enough
detail for a reader unfamiliar with `M6`, `M7` or `M8` to judge what kind of
evidence underlies every number reported later; the full inventory is given
as a table in Section 4.5.

### 2.1 Robustness: Northern Ghana, Talensi and Lower Anayari field chemistry

`M6` is stressed against a five-level measurement-tier ladder built from the
canonical Northern Ghana workbook: 160 boreholes across four administrative
regions in the Voltaian and crystalline basement aquifer systems typical of
semi-arid West Africa [@banoengyakubo2011ghana; @macdonald2012africa],
sampled once each, providing major-ion chemistry, nitrate,
fluoride, strontium, silica, stable water isotopes, and Hydrosheaf-computed
PHREEQC saturation indices [@parkhurst2013phreeqc]. The tier ladder removes
measurement classes one at a time from a full panel (T4) down to major ions
alone (T0), and every well is scored at every tier by the same evidence-gated
mechanism-resolution classifier. Two further field datasets, Talensi
District (63 samples) [@chegbeleh2020talensi] and Lower Anayari catchment (41
samples) [@zakaria2021anayari], provide an external sparse-transfer check at
fixed native tiers (Talensi at Tier 1, Lower Anayari at Tier 2); both are
independently collected field surveys with no reconstructed attribute. A
synthetic extended model with known dominant-process truth, built
independently of the field data, provides the one ground-truth check
available for this layer (Section 3.2).

### 2.2 Identifiability: an independent MODFLOW 6/MODPATH 7 generator

`M7` is stressed against fresh cases from an independent MODFLOW 6/MODPATH 7
generator that imports no HydroSheaf code [@langevin2017modflow6;
@pollock2016modpath]. `M7.3`, the core locked benchmark, uses six development
cases and twelve untouched locked-test cases, each carrying hydraulic, tracer
(tritium, argon-39) and chemical (reaction-family) evidence that can be
combined, withheld, or deliberately permuted. `M7.4` and `M7.5` use a
separate, code-independent scalar-section generator over 64 held-out and 128
locked-test cases respectively to compare an affine-sheaf solver against
edge-local and identity-Laplacian graph baselines under four scenario strata
(native, identity-limit, heterogeneous-affine, incompatible-cycles,
noisy-missing). `M7.6` uses a further independent generator, informed only by
declared synthetic nuisance parameters, to test one specific mechanism
hypothesis about a shared-nuisance explanation for the age layer's own
infeasibility pattern. None of these four generators shares code with the
HydroSheaf inference it evaluates; the Northern Ghana workbook is audited
separately, at component-diagnostic level only, and is not a source of
topology, age, or connectivity ground truth for this layer (Section 4.5).

### 2.3 Calibration: three controlled-synthetic protocols with known truth

`M8` is stressed against parameters that are known exactly by construction.
A matched-model transport experiment calibrates the production
`TransportCalibrationAdapter` in one-dimensional advection-dispersion mode
against synthetic dispersivity (2.0 m) and first-order decay (0.005 d^-1)
truth, across 16 fixed four-observation designs and a further optimal-design
experiment (4,000 and 4,500 calibrations respectively). A kinetic experiment
drives the production PHREEQC kinetic adapter against known calcite rate
constant and reactive surface area. An independent-model experiment
regenerates truth with an implicit finite-volume/upwind numerical solver
sharing no code with the calibration forward model, gated on a locked
grid-convergence check, and recalibrates the same production adapter against
it (1,000 locked-test calibrations). A separate frontier active-learning
experiment uses a fourth, independent MODFLOW 6/MODPATH 7 heterogeneous-
aquifer generator with nonlinear synthetic geochemistry, over 24 untouched
locked-test cases. No field dataset enters `M8` anywhere; this is a
deliberate design choice reported directly in Section 4.5 and Section 5.4,
not a gap discovered after the fact.

### 2.4 What this comparison does and does not add to the underlying data

No file listed in the dataset inventory (Section 4.5) was re-read,
re-fitted, or re-simulated to produce this paper. Every number reported in
Sections 4 and 5 traces to a result file that `M6`, `M7` or `M8` had already
written and locked before this comparison was constructed (Methods, Section
3.1); the only new computation performed here is the harmonization
arithmetic described in Section 3 -- selecting rows, relabelling into a
common internal-signal/external-signal schema, and passing through
already-computed paired bootstrap contrasts unchanged -- applied to those
existing files.

## 3. Methods

### 3.1 Computation and figure authority; no re-simulation; staleness check

This paper performs no PHREEQC, MODFLOW/MODPATH, PEST-GLM, or Bayesian
active-learning re-run. Every number below was read from a result file that
`M6`, `M7` or `M8` had already written and locked under each component's own
predeclared design. Because three commits made after those results were
locked modified shared `hydrosheaf` modules (`inference/network_fit.py`,
`sheaf/topology_refine.py`, `models/reactions.py`, `calibration/*`, the
`validation/` package), residual staleness risk was checked directly, not
assumed away: the per-script import graph of every `M6`/`M7`/`M8` script
producing a number cited here was traced against the literal diff of every
changed shared module. The changed `network_fit`/`topology_refine` branches
are gated behind an opt-in flag left at its default in every calling script;
the changed reaction-dictionary logic is imported only by `M5`'s
fully-specified synthetic benchmark (where it is inert) and by two `M6`
scripts explicitly excluded from `M6`'s own submission pipeline; and three
headline numbers (`M6`'s tier-ablation table, `M7.3`'s native/adverse-control
contrasts, `M7.4`'s sheaf-versus-graph contrasts) were independently
recomputed from raw locked CSVs and matched their components' own prose
summaries exactly. One number chain was found genuinely stale -- `M7`'s
public-pipeline system-acceptance run, which imports a module changed after
that run was locked -- and is excluded from every claim in this paper. Full
detail is in `O4/DECISIONS.md`, Decision D3. A Python harmonization layer
reads already-locked files and performs only disclosed, auditable arithmetic
-- selecting a row, relabelling a column, or passing through an
already-computed paired bootstrap contrast -- and writes the results as tidy
CSV exports. R consumes those exports and owns every figure; no figure
recomputes a reported statistic.

### 3.2 Each component's own stress design

**Robustness (`M6`).** Every one of 160 Northern Ghana wells is classified at
five nested measurement tiers (T4 full panel down to T0 majors-only) by the
same evidence-gated mechanism-resolution classifier, yielding a resolution
class (identifiable, partially identifiable, non-identifiable), a
mechanism-resolution score (MRS, a continuous 0-100 fit-quality index), and a
bootstrap support-stability estimate at every tier. A negative control
disables the evidence gate entirely and re-scores every tier, isolating
whether an apparent identifiability collapse reflects the classifier's
conservative prior or the underlying evidence. An independent synthetic
extended model, built with declared true dominant processes and no access to
the field classifier, provides a second check: whether the direction of an
identifiability change (for example, adding strontium/silica) matches a
known truth, via structural rank and exact-mineral recovery, rather than only
trusting the field-side classification.

**Identifiability (`M7`).** `M7.3` scores directed-edge classification
(precision-recall against the independent generator's true edges) and
groundwater-age accuracy under seven evidence conditions: hydraulics alone,
age alone, chemistry alone, each pairwise combination, the full three-stream
combination, and three adverse controls in which one stream (age, hydraulics,
or all three jointly) is deliberately permuted before being combined with the
others. Every condition reports both a posterior-uncertainty summary (mean
edge entropy) and a predictive-skill summary (PR-AUC, Brier score, log-loss)
against the same held-out true edges, with 10,000 case-block paired bootstrap
resamples per contrast. `M7.4` and `M7.5` extend this design to isolate what
an affine-sheaf restriction-map solver contributes beyond an edge-local
weighted graph and an identity-coupled graph Laplacian, under equal inputs,
optimisation routine, and calibration form, across four scenario strata
designed to stress heterogeneity, incompatible cycles, and noisy or missing
observations.

**Calibration (`M8`).** The production `PESTGLM` engine -- a
Gauss-Levenberg-Marquardt optimiser with log10-parameter bounds, parallel
Jacobian evaluation, and covariance by matrix inverse with an
SVD-pseudoinverse fallback -- calibrates against synthetic transport and
kinetic truth under three linked designs. A fixed 16-design sweep varies
observation placement (early/late/clustered/log-spread sampling) at 250
paired replicates per design. An optimal-design experiment scores eight
candidate additional observations by parameter-marginal-variance,
D-/A-/E-optimality, and a balanced criterion in the whitened log-parameter
Fisher information matrix, then verifies each candidate by actually adding it
and recalibrating, under both the matched analytical forward model and an
independently generated numerical truth accepted only after passing a
grid-convergence gate. A kinetic experiment tests whether one or four
residence-time observations can separate a rate constant from a reactive
surface area, with structural support for non-identifiability required to
show near-collinear sensitivities, numerical rank one, and invariant
predictions along a constant-product ridge. A separate frontier
active-learning experiment scores a prospectively specified Bayesian
value-of-information policy against random acquisition and a strong legacy
heuristic under a frozen superiority/noninferiority gate.

### 3.3 A common stress/signal taxonomy

Every retained experiment cited in this paper was classified by three axes:
the stress mechanism (data limitation, evidence misspecification, or a
model-form shift away from the calibrating model), the internal signal
reported (a quantity computed without reference to ground truth: fit
quality, posterior-entropy reduction, nominal interval coverage, or
optimiser convergence), and the external signal reported (a quantity
verified against ground truth: true identifiability class, predictive skill
against an independent generator, realised parameter-recovery error, or
structural information rank). This taxonomy (Table 1) is applied without
altering any component's own claim; it adds a shared vocabulary across `M6`,
`M7` and `M8`, not a re-interpretation of any single result.

### 3.4 What counts as divergence, and its limits

The three benchmarks do not share a statistic, any more than `O3`'s three
components did. `M6`'s internal signal (mean MRS, a 0-100 index) and external
signal (percent of wells verified non-identifiable) are read on the same
five-tier x-axis but are not the same kind of quantity. `M7`'s internal
signal (a change in posterior edge entropy, in nats) and external signal (a
change in PR-AUC, a probability-scale ranking statistic) are two different
units entirely; the divergence reported is a sign-and-magnitude comparison of
two changes, not a single combined statistic. `M8`'s internal signal (a
coverage fraction) and external signal (a log-scale point error) are read
for the same parameter under the same design, which makes their divergence
the most directly interpretable of the three, but the comparison is still
between two different statistics, not one ratio. This is stated as a
disclosed methodological choice, following the precedent `O3` set for its
own capture/correctness comparison: the claim is that an internal signal and
an external signal move in different directions or at different rates under
stress, not that they are the same measurement examined twice.

### 3.5 Reproducibility

The full harmonization layer -- five Python scripts and six R figure scripts
-- is archived under `O4/analysis/`, runs end to end from the repository
root with fixed, already-locked inputs, and writes nothing back to any `M6`,
`M7` or `M8` result file.

## 4. Results

### 4.1 Applying the common taxonomy

Table 1 classifies twelve retained experiments across `M6`, `M7` and `M8` by
stress mechanism and by the internal/external signal pair each reports.
Every experiment carries a negative or structural control except `M6`'s
edge-set perturbation and external-transfer rows, which are themselves used
as diagnostics rather than as claims requiring their own control (Section
4.3). No experiment's classification here contradicts the claim already
recorded for it in its own source component's `DECISIONS.md`.

![Figure 1. Three independently locked HydroSheaf stress tests and the common internal-versus-external signal taxonomy this paper adds. `M6` strips measurement tiers from real field chemistry; `M7` combines or deliberately misspecifies evidence streams against an independent generator; `M8` shifts from a matched calibration model to an independently generated numerical truth. Each pipeline was locked independently, on its own timetable, before this comparison was constructed.](artifacts/figures/FIG-1.png)

**Table 1.** Common stress/signal taxonomy applied to twelve retained `M6`, `M7` and `M8` experiments (abridged; full table in Supplementary Table S1).

| Component | Experiment | Stress axis | Internal signal | External signal | Control |
|---|---|---|---|---|---|
| M6 | Tier ablation T4→T0 | data limitation | mean MRS | verified identifiability class | evidence-gate-off ablation |
| M6 | Edge-set perturbation | model-form (connectivity) | per-edge MRS | network process-composition shift | — |
| M6 | External sparse transfer | data limitation + domain shift | mean MRS | % non-identifiable | — |
| M6 | Synthetic validation, known truth | model-form | exact-mineral F1 | true dominant-process recovery | independent extended model |
| M7 | Native evidence-panel integration | evidence combination | entropy reduction | PR-AUC/Brier/log-loss | permuted-stream controls |
| M7 | Adverse controls (permuted age/hydraulics/joint) | evidence misspecification | entropy reduction | PR-AUC/Brier/log-loss | self (each IS the control) |
| M7 | M7.4 sheaf vs graph | model-form (restriction maps) | log-loss/ECE | PR-AUC/selected-F1 | permuted-map control |
| M7 | M7.5 robust hybrid | evidence misspecification | Brier/log-loss | PR-AUC | permuted-map control |
| M8 | Fixed 16-design sweep | data limitation | success rate; coverage | recovery error | — |
| M8 | Optimal design, independent truth | model-form | coverage | recovery error | no-new-measurement; random |
| M8 | Kinetic structural confound | data limitation (measurement type) | convergence; objective value | Fisher-information rank | off-ridge product-doubling check |
| M8 | Frontier active learning | evidence-value/decision design | entropy reduction per cost | Brier/PR-AUC | random acquisition; realised oracle |

### 4.2 The central pattern: internal signals do not track external signals under stress

Figure 2 and Table 2 report the paper's central result, one panel per
component on its own native scale. In `M6`, mean MRS is essentially flat
across all five evidence tiers (71.0, 70.7, 70.7, 68.6, 68.4), a range of
2.6 points, while the percent of wells independently verified non-
identifiable rises from 0.0% [95% Wilson CI -0.0, 2.3] to 60.0% [52.3, 67.3],
a range of 60 points whose endpoint intervals do not overlap, with the
single largest step (T3→T2, removing strontium/silica: 0.6% [0.1, 3.5] to
51.9% [44.2, 59.5], also non-overlapping) accounting for 51.3 of those 60
points on its own. In `M7`, adding chemistry to the hydraulics-plus-
chemistry-absent baseline reduces mean edge entropy by 0.0827 nats (95% CI
-0.0985, -0.0653) while improving PR-AUC by 0.447 (0.357, 0.540) -- entropy
and skill move together here, the one condition in this comparison where
they do. Every other native or adverse condition breaks that alignment:
adding age reduces entropy by 0.0006 nats (-0.0009, -0.0003) while *reducing*
PR-AUC by 0.0060 (-0.0122, -0.0011); permuting age reduces entropy by 0.0207
nats (-0.0236, -0.0175) -- more than three times the entropy reduction from
honestly adding age -- while reducing PR-AUC by 0.0754 (-0.1353, -0.0148),
more than twelve times the honest condition's own loss. In `M8`, the design
with the lowest realised dispersivity error (`c90_s45`, median absolute
log10 error 0.0172) and the design with the highest (`c120_s5`, 0.3818) both
report success rate 1.0 and linearised 95% coverage within nine points of
each other (0.804 versus 0.840); success rate and coverage, the two signals
available before a reader ever sees the recovery error, cannot distinguish a
21-fold difference in realised accuracy. A second, sharper `M8` instance
appears under a model-form shift rather than a design change: adding one 50 d
observation to the same three-point late-time design gives the matched
analytical calibration model near-nominal coverage for both parameters
(dispersivity 0.832, decay 0.836), but recalibrating the identical production
adapter against an independently implemented numerical solver, with the same
50 d observation, collapses coverage to 0.02 for dispersivity and 0.004 for
decay -- a signal that stayed inside nine points of nominal under the
matched model gives almost no warning of a model-form mismatch it cannot see
from inside the calibration it is reporting on.

![Figure 2. Central result: paired internal confidence signal (grey/dark) versus external truth signal (red), one panel per component, on each component's own native scale. Panels are not numerically comparable to one another. M6: mean MRS versus percent of wells verified non-identifiable, five evidence tiers. M7: entropy change versus PR-AUC change, native evidence increments and one adverse control. M8: linearised 95% coverage versus median absolute log10 recovery error, lowest- and highest-error fixed designs for each parameter.](artifacts/figures/FIG-2.png)

**Table 2.** Headline internal-signal versus external-signal pairs, selected conditions, with 95% intervals where the source component reports them.

| Component | Condition | Internal signal | Value [95% CI] | External signal | Value [95% CI] |
|---|---|---|---|---|---|
| M6 | Tier T4 → T0 | Mean MRS (range across tiers) | 71.0 → 68.4 (Δ2.6) | % wells non-identifiable | 0.0% → 60.0% (Δ60) |
| M6 | Tier T3 → T2 (the cliff) | Mean MRS | 70.7 → 70.7 (Δ0.0) | % wells non-identifiable | 0.6% [0.1, 3.5] → 51.9% [44.2, 59.5] (Δ51.3) |
| M7 | Native: age added to HC | Mean edge-entropy change | −0.0006 [−0.0009, −0.0003] | PR-AUC change | −0.0060 [−0.0122, −0.0011] |
| M7 | Adverse control: age permuted | Mean edge-entropy change | −0.0207 [−0.0236, −0.0175] | PR-AUC change | −0.0754 [−0.1353, −0.0148] |
| M7 | Native: chemistry added to HA | Mean edge-entropy change | −0.0827 [−0.0985, −0.0653] | PR-AUC change | +0.4471 [+0.3575, +0.5401] |
| M8 | 50 d obs., dispersivity: matched model → independent truth (same design) | Coverage | 0.832 → 0.02 | Median \|log10 error\| (independent truth, no obs. → 50 d obs.) | 0.826 → 0.167 |
| M8 | 50 d obs., decay: matched model → independent truth (same design) | Coverage | 0.836 → 0.004 | Median \|log10 error\| (independent truth, no obs. → 50 d obs.) | 0.137 → 0.154 |

Notes. M6 percentages carry 95% Wilson score intervals computed from the
underlying well counts (`O4/analysis/python/derive_robustness_gradient.py`;
n=160 per tier), added after internal review found bare point estimates on
proportions with a known denominator (`manuscript/review/
O4_manuscript_review.md`); MRS is a continuous index with no equivalent
closed-form interval reported by `M6` and is shown without one. M6 values
recomputed directly from `results/m6_tier_ablation.csv` (800
well-tier rows); M7 values pass through `results/m7_3_locked/
evidence_case_bootstrap_contrasts.csv` unchanged; M8 values pass through
`manuscript/artifacts/m8_transport_oed_summary_reviewed.csv` and
`results/RUN-M8-INDEPENDENT-20260728-01/strategy_summary.csv` unchanged. The
two M8 rows report two different comparisons deliberately placed side by
side, not one: the coverage column contrasts the matched analytical model
against an independently generated numerical truth calibrated with the
identical 50 d design (Section 4.2); the error column contrasts, within the
independent-truth arm only, no added observation against the same 50 d
addition (Section 4.3). Read together, they show a signal (coverage) that
looks nominal under the matched model and fails under a model-form shift,
next to a point estimate that the same added observation moves in opposite
directions for the two parameters. Full three-way values (matched;
independent-truth baseline; independent-truth with the 50 d observation) are
in Supplementary Table S1's source export,
`calibration_model_form_shift.csv`.

### 4.3 The pattern does not run in only one direction

The `M8` model-form shift is the clearest case that this is a divergence
pattern, not simply "internal signals are always too optimistic." Under the
independently generated numerical truth, the same 50 d observation that
improves dispersivity recovery (median absolute log10 error 0.8262 to 0.1674,
paired difference −0.6690, 95% CI −0.7374 to −0.5757) simultaneously
*degrades* decay recovery (0.1367 to 0.1541, paired difference +0.0210,
+0.0092 to +0.0276) -- and coverage collapses for both parameters regardless
of which direction the point estimate moved (dispersivity 0.788 to 0.02;
decay 0.64 to 0.004). A reader who used coverage alone to decide whether to
trust the 50 d recommendation would have no way to tell, from the internal
signal, that one parameter's point estimate had just become more trustworthy
while the other's had become less so. `M7.4`'s sheaf-versus-weighted-graph
contrast shows a related, milder asymmetry: pooled across all four scenario
strata, PR-AUC is statistically indistinguishable between the two methods
(+0.0097, 95% CI −0.0054 to +0.0248), but log-loss is significantly worse for
the sheaf (+0.0117, 95% CI +0.0008 to +0.0231) -- the external ranking signal
reports a tie while the internal calibration signal reports a measurable
cost, in the same comparison, for the same sixty-four cases.

### 4.4 Negative and structural controls bound each pattern

Table 3 reports the controls that rule out the simplest alternative
explanation for each finding -- that the divergence is an artefact of the
classifier, the entropy heuristic, or the coverage formula, rather than a
real property of the stress test. `M6`'s evidence-gate-off ablation returns
0% non-identifiable at every one of the five tiers (compare Table 2's
0.0-60.0% range with the gate on), confirming the tier "collapse" is the
framework's conservative prior being triggered by genuinely reduced evidence,
not a classifier artefact; the same ablation's structural-rank diagnostic
moves in the theoretically expected direction (rank 8 at tiers T0-T2, rank 11
at tiers T3-T4; silicate coherence falls from 0.707 to 0.500 as strontium/
silica are added), an independent, code-separate confirmation that the T3→T2
step is where real information is lost, not merely where the classifier's
threshold happens to sit. `M7`'s three adverse controls are each beaten by
every native evidence condition on PR-AUC: permuted age (external signal
−0.0754) is worse than native age (−0.0060), permuted hydraulics reduces
PR-AUC by 0.0686 against a native hydraulics gain of +0.0091, and joint
misspecification reduces PR-AUC by 0.139, the largest degradation of any
condition tested, while producing the largest entropy reduction (−0.0706
nats) of any condition tested -- the sharpest instance of the central pattern
in this paper. `M8`'s kinetic experiment confirms its rank-one finding is
structural rather than a dead experiment: predictions along the constant-
product ridge differ by at most 1.33×10⁻⁶ declared observation standard
deviations, while doubling the product changes predictions by 6.42 standard
deviations, and an independent surface-area observation restores numerical
rank two (condition number 45.9, parameter correlation −0.957).

**Table 3.** Negative- and structural-control results bounding each layer's central finding.

| Component | Finding being bounded | Control | Result |
|---|---|---|---|
| M6 | Tier collapse (0.0%→60.0% non-identifiable) | Evidence-gate-off ablation | 0% non-identifiable at every tier |
| M6 | T3→T2 cliff is where information is lost | Independent structural-rank/silicate-coherence diagnostic | Rank 8→11, silicate coherence 0.707→0.500, exactly at the T3/T2 boundary |
| M7 | Entropy reduction ≠ predictive gain | Every native condition vs. its adverse-control counterpart | Native beats permuted on PR-AUC in all three matched pairs |
| M7 | Joint misspecification is the sharpest instance | Joint-permutation adverse control | Largest entropy reduction (−0.0706) and largest PR-AUC loss (−0.139) of any condition |
| M8 | Kinetic rank deficiency is structural | Off-ridge product-doubling; independent surface-area observation | Off-ridge response 6.42 SD; rank restored to 2 with the added observation |
| M8 | Frontier policy's entropy-reduction gain is a real decision benefit | Random-acquisition baseline; realised-oracle ceiling | Beats random on Brier (−0.0248) and entropy-reduction-per-cost (+0.0358); below oracle as expected |

### 4.5 Benchmark scale and field-/archive-transfer scope

Table 4 reports the three benchmarks' design scale and transfer scope on
their own terms. `M8` is the largest by calibration count: 8,500 matched-model
transport calibrations, 1,000 independent-model calibrations, and 24 frontier
active-learning cases, none against field data. `M7` is the largest by
resampling intensity: 10,000 case-block bootstrap replicates per contrast
across `M7.3`, `M7.4` and `M7.5`, against fresh independent MODFLOW 6/MODPATH
7 cases throughout. `M6` is the smallest by design count but the only layer
with a genuine field-transfer footprint: 160 Northern Ghana wells at five
tiers, plus 63 Talensi and 41 Lower Anayari samples as external sparse-
transfer checks. This asymmetry -- `M6` transfers to three real datasets,
`M7` audits the Northern Ghana workbook for component-diagnostic readiness
only, and `M8` has no field-transfer component at all -- is reported directly
rather than smoothed into a single "field-tested" label for the framework as
a whole.

![Figure 6. Benchmark scale and field-/archive-transfer scope across the three components. Bar presence indicates the design/replicate scale in Table 4; shading indicates whether any field or archive dataset enters that specific benchmark row.](artifacts/figures/FIG-6.png)

**Table 4.** Benchmark scale and field-/archive-transfer scope (abridged; full table in Supplementary Table S4).

| Component | Primary scale | Replicate/bootstrap count | Field/archive transfer |
|---|---|---|---|
| M6 robustness | 160 wells × 5 tiers = 800 rows | not applicable (deterministic, seed 1234) | Northern Ghana (160), Talensi (63), Lower Anayari (41) |
| M7.3 identifiability | 6 dev + 12 locked-test cases | 10,000 case-block bootstrap replicates | Northern Ghana audited for readiness only |
| M7.4/M7.5 identifiability | 64 + (64 dev + 128 locked-test) cases | 10,000 paired bootstrap replicates | none (controlled synthetic) |
| M8 matched-model | 8,500 calibrations | 250 replicates/design; 2,000-resample bootstrap | none |
| M8 independent-model | 1,000 locked-test calibrations | 80 dev replicates; 5,000 paired bootstrap | none |
| M8 frontier active learning | 24 locked-test cases | 5,000 paired case-bootstrap | none (controlled synthetic) |

## 5. Discussion

### 5.1 Principal finding

Across three physically unrelated stress designs -- stripping real field
measurements, misspecifying synthetic evidence streams, and shifting to an
independently generated numerical truth -- an internally-generated
confidence signal failed to track an externally-verifiable truth signal in
every layer, though not in the same way twice. `M6`'s signal stayed flat
while the truth collapsed; `M7`'s signal improved (entropy fell) exactly
where the truth got worse, most sharply under deliberate misspecification;
`M8`'s signal collapsed for a parameter whose truth-tracking point estimate
simultaneously improved. A framework that reported only convergence, fit
quality, entropy, or nominal coverage, without independently verifying what
those signals do under stress, would have no internal way to notice any of
these three patterns. This paper's contribution is not that internal signals
can be wrong in general -- that is well established [@beven2001glue;
@neuman2003bma] -- but that this specific framework's specific
internally-generated signals were checked against ground truth under three
unrelated stress mechanisms and found not to track it in all three,
including one case (`M8`'s decay parameter) where the divergence runs in the
direction opposite to the framework's own better-known result for the
companion parameter in the same experiment.

### 5.2 Relation to the equifinality, identifiability and optimal-design literature

The pattern reported here instantiates, with three quantified cases, a
general argument that a converged, well-fitting model is not thereby an
identified one [@beven2001glue; @beven2006manifesto]. `M8`'s condition-
number and coverage diagnostics operationalise the same distinction the
groundwater calibration literature draws between a well-posed inverse
problem and a well-identified one, where sensitivity and null-space analysis
are used specifically to separate parameters and combinations of parameters
the data can and cannot resolve, independently of whether the optimiser
converges [@hilltiedeman2007effective; @dohertyhunt2009identifiability];
this paper's contribution is to show the same dissociation directly, via
coverage collapse under a model-form shift, rather than only via a static
sensitivity or null-space diagnostic computed at one fixed model. It also
intersects a narrower, more recent literature on the field-transferability
of optimal experimental designs: designs chosen by a linearised Fisher-information
criterion under one forward-model form are not guaranteed to remain optimal,
or even beneficial, under a different model form
[@sreekanth2017monitoring], a difficulty `M8`'s independent-numerical-truth
result reproduces directly and quantitatively rather than only in principle.
`M7`'s finding that reduced uncertainty can accompany reduced accuracy is a
specific, adversarially controlled instance of the broader caution that
model-averaging and multi-source evidence combination reduce reported
uncertainty without a guarantee that the combined estimate is more accurate
[@neuman2003bma]; this paper's contribution is to show the divergence occurs
even in the well-behaved native evidence-addition case (age added to
hydraulics-plus-chemistry), not only under an adversarial permutation, and to
quantify how much larger the divergence becomes once misspecification is
deliberately introduced.

### 5.3 Relation to M2.3 and O3

The companion framework article (Objective 2) [@hydrosheaf_m2_3_inprep]
reports HydroSheaf's claim-gated architecture and an integrated field
demonstration; it does not
attempt this paper's cross-component stress comparison. `O3` (Objective 3)
harmonizes a disjoint evidence triad (`M3`, `M4`, `M5`) on a disjoint
question -- whether detection exceeds correctness under independent
evaluation -- and reports no finding about internal-versus-external signal
divergence under stress. This paper's question could not have been asked of
`O3`'s evidence, because none of `M3`, `M4` or `M5` was designed as a stress
test in the sense used here: each compares one estimate against one external
reference under a fixed evidence condition, not the same estimator's own
confidence signal across a deliberately varied stress level. Together, the
two O-series papers cover what their six underlying components already
established independently -- `O3` the detection/correctness axis, this
paper the confidence/truth axis -- without either restating the other's
figures or claiming the other's result as its own.

### 5.4 Practical implication: report the divergence check, not only the signal

A calibration or field-transfer report that states convergence, fit quality,
or a nominal coverage figure without independently checking whether that
signal tracks a verifiable truth under the specific stress the application
will face is, on this paper's evidence, not distinguishable from one that
checked and found the signal untrustworthy. The three checks reported here
are inexpensive relative to the calibrations or field campaigns they audit:
`M6`'s gate-off ablation and structural-rank cross-check reuse the same 160
wells; `M7`'s adverse controls reuse the same generator and cases as the
native conditions; `M8`'s independent-numerical-truth arm reuses the same
production calibration engine against a solver that took a few hundred lines
to implement independently. The practical recommendation is not a new
diagnostic but a reporting discipline already implicit in each component's
own design: state the stress condition under which a confidence signal was
verified, and do not report a coverage figure, an entropy reduction, or a
fit-quality index as evidence of identifiability without that verification
attached.

### 5.5 Limitations

Four limitations bound this comparison. First, as with any comparison of a
model to a reference, agreement should be read as consistency under stated
conditions, not proof of a correct underlying process
[@oreskes1994validation]; that general caution applies to every synthetic
ground-truth check reported here, not only to the field-side `M6` result.
Second, the internal-versus-external divergence reported for each layer uses
that layer's own native metric pair; Section 3.4's disclosed analogy across
MRS-versus-identifiability-class, entropy-versus-PR-AUC, and
coverage-versus-recovery-error is a methodological choice, not a proof that
the three axes measure one underlying quantity, and Figure 2's three panels
should be read as three parallel single-component results on adjacent axes,
not one merged statistic. Third, `M6` and `M7`'s result files were locked
before the 2026-08-01/02 `hydrosheaf/validation` consolidation; the
per-script import-graph check reported in Section 3.1 and `DECISIONS.md` D3
found the code paths these results depend on unaffected, but a change in a
shared dependency that neither this nor any other check can fully rule out
remains a residual risk. Fourth, and most importantly, `M6` is the only layer
with genuine field-transfer evidence; `M7`'s Ghana-workbook audit is
component-diagnostic, not topology or age validation, and `M8` has none. No
result in this paper is field validation of robustness, non-uniqueness, or
calibration accuracy; the divergence pattern reported is demonstrated under
controlled synthetic or field-transfer-diagnostic conditions specific to
each layer's own design, not a general claim about groundwater inverse
frameworks beyond HydroSheaf.

## 6. Conclusion

Three already-locked HydroSheaf stress tests -- a real-data evidence-tier
ablation (`M6`), an adversarially controlled evidence-integration test
against independent synthetic truth (`M7`), and three prospectively locked
controlled-synthetic calibration protocols (`M8`) -- were harmonized under
one taxonomy without re-running any of the underlying computation. In all
three, an internally-generated confidence signal failed to track an
externally-verifiable truth signal under stress: flat fit quality against a
60-point identifiability collapse in `M6`; entropy reduction exceeding, in
adverse conditions, three times its honest-condition value while predictive
skill degraded twelve times as much in `M7`; and interval coverage
collapsing for a parameter whose point-estimate accuracy simultaneously
improved in `M8`. Negative and structural controls in every layer confirm
the pattern is a property of the stress test, not an artefact of the
classifier, the entropy heuristic, or the coverage formula. None of this is
a fourth, independent validation of any component, and none of it is a
second HydroSheaf framework contribution; both remain the domain of the
components' own result packages and of the companion framework article.
What this comparison adds is the demonstration, unavailable from any one
underlying benchmark, that numerical health is not evidence of
identifiability, integration value, or calibration validity in this
framework, and a reporting discipline -- verify the signal against ground
truth under the specific stress an application will face -- that each
component already practises but that no prior paper stated as a shared
requirement across all three.

## Code and data availability

The harmonization scripts, R figure scripts, tidy CSV exports, and the
component result files this comparison reads (`M6/m6_field_transfer_benchmark/`,
`M7/m7_nonuniqueness_benchmark/`, `M8/m8_calibration_benchmark/`) are archived
at `https://doi.org/10.5281/zenodo.PLACEHOLDER-DOI` and developed at
`https://github.com/PLACEHOLDER-ORG/hydrosheaf`. No PHREEQC, MODFLOW/MODPATH,
PEST-GLM, or Bayesian active-learning re-run is required to reproduce any
number in this paper; only the already-published result files listed in
Table 4 and Python/R are required (`O4/README.md`).

## Author contributions

To be completed with named CRediT roles before submission: conceptualisation,
methodology, software, formal analysis, data curation, writing -- original
draft, writing -- review and editing.

## Ethics, competing interests and funding

This comparison used only already-locked result files and involved no new
human or animal subject procedures, and no new field sampling. The authors
declare no competing interests. Funding is to be declared in the final
submission metadata.
