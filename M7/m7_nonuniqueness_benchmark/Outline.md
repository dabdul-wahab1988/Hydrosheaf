# Q1 Journal Manuscript Outline

Development status: drafted 24 July 2026. Word budget fixed at **6,500 words of
main text (Introduction through Conclusion)** on explicit user instruction,
consistent with the M5 revised budget. Grounded entirely in the locked M7.3
package (`docs/m7_3_protocol.md`, `docs/m7_3_results.md`,
`docs/figures_and_tables.md`, `tables/publication/*.csv`,
`figures/publication/*`) and its retained M7.2 supporting-validation audit
trail (`docs/supporting_validation_*.md`). No new computation is introduced by
this outline; it organises already-executed, already-locked results into a
manuscript.

## Recommended Journal Positioning

Primary target: **Water Resources Research** (AGU/Wiley), Research Article.

Rationale: M7 is fundamentally an identifiability/equifinality and
multi-evidence-integration benchmark for groundwater inverse problems — the
same class of contribution as M5, which already targets WRR at this word
budget. AGU's 25-publication-unit allowance (1 PU = 500 words or one display
item) accommodates a 6,500-word main text plus six main figures and one main
table (13 PU) with substantial headroom, so the paper can carry four linked
experiments and their adverse controls without an excess-length fee. WRR
readers and reviewers are also the correct audience for a paper whose central
finding is conditional, partly negative, and organised around explicit
decision rules rather than a single positive claim — that framing is
consistent with WRR's history of publishing rigorous non-uniqueness and
equifinality studies (e.g. Beven and Freer, 2001; Beven, 2006) and joint/data
fusion evaluations (e.g. Linde et al., 2006).

Strong alternative: **Journal of Hydrology** (Elsevier), Research Paper (word
allowance up to 10,000 words; explicitly welcomes negative/conditional
transferability results).

Shorter-form fallback if the target is later constrained to a Nature Portfolio
Article format (~5,000 words): **Communications Earth & Environment**. A
condensed section allocation is archived at the end of this outline for that
scenario; it is not the working budget while the target is WRR.

Stretch target after an additional independent field/basin validation:
**Nature Water**.

## Working Title

**Evidence integration reduces groundwater interpretive non-uniqueness only conditionally: an independent benchmark with adverse controls**

Alternative titles:

1. **When more evidence helps and when it deceives: a truth-blind benchmark of hydraulic, age and chemical evidence integration**
2. **Conditional value of evidence integration in groundwater inference: complementarity, redundancy and negative transfer**
3. **False confidence is a measurable failure mode of multi-evidence groundwater inference**

## Role Within the Hydrosheaf PhD

M7 is not a new inference method and must not be written as a second topology,
age, or reaction paper. It is the cross-layer integration audit of the
Hydrosheaf inference chain:

1. **M4 — Where can groundwater move?** Candidate directed edges and their evidential support.
2. **M3 — How long has groundwater moved?** Residence-time distributions and age-ordering constraints.
3. **M5 — What happened chemically along a supported edge?** Reaction classes, extents, alternatives, and identifiability.
4. **M6 — How stable is the resulting interpretation under real, data-limited field transfer?**
5. **M7 — Does combining these independent evidence streams actually reduce interpretive non-uniqueness, and when does it not?**

M3, M4 and M5 each validate one evidential layer in isolation, largely against
internally generated or PHREEQC-based synthetic truth. M6 tests whether the
resulting workflow transfers to real, uneven Ghanaian field data. M7 asks a
different, prior question that the other papers assume rather than test:
*is combining these layers actually beneficial, and can it be told apart from
merely looking more confident?* It answers this with a truth generator that
imports no Hydrosheaf code (official MODFLOW 6/MODPATH 7 plus an independent
nonlinear chemistry/tracer model), truth-blind inference, and predeclared
adverse controls that can expose false confidence.

**Claim boundary.** M7 validates whether integrating hydraulic, age and
chemical evidence changes topological and age non-uniqueness in a defensible,
calibration-checked direction. It does not re-validate the internal accuracy
of the M3 age model, the M4 topology sampler, or the M5 reaction inversion in
isolation — those are each paper's own domain. It does not use Ghanaian field
data as truth for topology, age or reaction mechanism.

## Central Scientific Claim

Reduced posterior uncertainty is not, by itself, evidence that integrating an
additional evidence stream improved a groundwater interpretation. Evidence
streams can be complementary (jointly informative beyond any single stream),
redundant (numerically inert once other streams are present), or actively
harmful when misspecified — and misspecification can look like an
improvement, because it lowers apparent uncertainty while degrading
predictive skill and calibration. Correct flow-topology information reduces
groundwater-age uncertainty, particularly when the tracer panel is weak, but
age evidence alone does not improve topology ranking in the reverse
direction, and an incorrect topology assumption is not a merely-worse
alternative estimate — it can destabilise the inference outright. Reaction
non-uniqueness is only partly resolved by adding chemical indicators:
denitrification, sulfate reduction and silicate weathering respond to an
enhanced ion panel, while carbonate weathering and precipitation remain
non-identifiable regardless of panel richness.

This claim makes M7 essential to Hydrosheaf rather than a redundant fourth
validation study. It supplies the framework's central integration guardrail:
a decision rule that separates genuine complementarity from redundancy and
negative transfer, so that combining evidence streams is a defensible choice
rather than an assumed one.

## Research Questions

1. Does combining hydraulic, age and chemical evidence reduce topological interpretive non-uniqueness, and does every pairing behave the same way?
2. How do correct, partial and reversed flow-topology assumptions change groundwater-age accuracy, interval width and coverage under an informative versus a weak tracer regime?
3. Does expanding the hydrochemical indicator panel resolve reaction-mechanism non-uniqueness uniformly across processes, or only for some?
4. When an evidence stream is misspecified, does the resulting reduction in posterior uncertainty still indicate improved inference, or can it manufacture false confidence?
5. What can and cannot be defensibly claimed about a real, data-limited aquifer system (Northern Ghana) given the evidence it actually contains?

## Testable Hypotheses

- **H1.** Combining independent evidence streams reduces topological non-uniqueness only when the added stream is complementary (jointly informative beyond the others already present); it is not a monotonic function of the number of streams combined.
- **H2.** Correct flow-topology information reduces groundwater-age error and interval width relative to no topology assumption, with a larger benefit under a weak (single-tracer) regime than an information-rich (paired-tracer) regime.
- **H3.** An incorrect (reversed) topology assumption degrades age accuracy or coverage under an informative tracer regime, and produces importance-sampling degeneracy under a weak regime, rather than yielding a competing but merely less accurate estimate.
- **H4.** Predeclared adverse permutation controls that preserve each evidence stream's marginal distribution while destroying its case-specific meaning reduce posterior entropy while simultaneously degrading discrimination (PR-AUC) and calibration (Brier score, log loss, overconfident error).
- **H5.** Expanding the hydrochemical indicator panel from a core to an enhanced set improves reaction-family recovery only for a subset of redox and silicate processes; carbonate weathering and precipitation remain non-identifiable under both panels.
- **H6.** Northern Ghana field data support component-level chemical and isotopic diagnostics but not field validation of topology, groundwater age or a unique reaction mechanism, given the absence of environmental age tracers, screen intervals, repeated head series and independently verified connectivity.

---

## Article Architecture and Word Budget

The target is **6,500 words from Introduction through Conclusion**, excluding
the Abstract, references, figure legends, code/data availability statements,
and Supplementary Information. Word counts are prose only and exclude
headings, equations, and citation markers, consistent with
`manuscript/methods/methods_manifest.json` policy.

| Section | Target words |
|---|---:|
| 1. Introduction | 1,050 |
| 2. Methods | 1,550 |
| 3. Results | 2,100 |
| 4. Discussion | 1,550 |
| 5. Conclusion | 250 |
| **Total** | **6,500** |

Results receives the largest share because four linked experiments (evidence
integration, topology-conditioned age, reaction non-uniqueness, Ghana scope
audit) each require a fully reported quantitative paragraph with uncertainty.
Methods and Discussion are weighted second so the adverse-control logic,
decision rules, and calibration guardrail can be explained without
compression. Introduction and Conclusion answer five research questions and
one decision framework respectively and do not require additional length.

Abstract: 180–200 words, unstructured, no references, word-capped
independently of the main-text budget. Write only after all analyses are
locked (they already are; the Abstract may be drafted alongside the other
sections).

### Abstract sequence (180–200 words, excluded from the 6,500-word total)

1. State the problem: combining independent groundwater evidence streams is
   widely assumed to reduce interpretive uncertainty, but reduced uncertainty
   alone does not prove improved inference.
2. State the methodological advance: a truth-blind benchmark using an
   independent MODFLOW 6/MODPATH 7 generator, four linked experiments, and
   predeclared adverse (permutation) controls that can distinguish
   complementarity, redundancy, and negative transfer.
3. State the benchmark scale: development/locked-test case counts, particle
   and bootstrap counts, and the Northern Ghana data-scope audit.
4. Report the principal quantitative results: the chemistry/hydraulics
   contribution to topology ranking, the null incremental-age result, the
   topology-conditioned age contrasts, the carbonate non-identifiability
   result, and the adverse-control entropy/PR-AUC divergence.
5. State the significance: a decision rule that lets practitioners tell
   genuine integration benefit apart from false confidence.

---

# 1. Introduction — 1,050 words

## 1.1 Groundwater interpretation from incomplete evidence — 150 words

Establish that hydraulic heads, environmental tracers, and hydrochemistry are
routinely combined to interpret groundwater flow paths, residence times and
reaction pathways, because no single evidence type is individually sufficient
for water-resource, contaminant-attribution or aquifer-vulnerability
decisions. Frame this as the general motivation for evidence integration in
subsurface hydrology.

## 1.2 Current integration practice assumes more evidence helps — 190 words

Summarise the standard justification for combining evidence: joint and
hydrogeophysical inversion is reported to reduce model ambiguity (Linde et
al., 2006); multi-tracer environmental dating combines several
radioisotopes to constrain age distributions more tightly than any single
tracer (Visser et al., 2013); Bayesian model averaging treats conceptual and
parametric uncertainty jointly across candidate models (Neuman, 2003).
Existing groundwater equifinality work (Beven and Freer, 2001; Beven, 2006)
establishes that a single model can fit well without being correct, but
typically evaluates one evidence type or one model family at a time rather
than testing whether adding a stream is itself beneficial, redundant, or
harmful.

## 1.3 Problem statement: uncertainty reduction is not proof of improvement — 200 words

State the central problem: an added evidence stream can appear to resolve
ambiguity by lowering posterior entropy while being numerically redundant, or
by being actively misspecified in a way that manufactures false confidence
rather than genuine resolution. Environmental tracer dating is known to carry
systematic limitations that are not always visible from fit quality alone
(McCallum et al., 2015). No widely used groundwater benchmark evaluates
evidence integration with an external, code-independent truth generator, a
locked development/test split, and predeclared adverse controls that can
separate a genuinely complementary stream from one that is redundant or
harmful. Absent such controls, apparently more confident integrated
interpretations may simply be more confidently wrong.

## 1.4 Novelty: a truth-blind conditional-integration benchmark with adverse controls — 180 words

Introduce the shift from assuming integration helps to testing when it does.
Define the principal innovations: (1) an external synthetic-truth generator
built on official MODFLOW 6/MODPATH 7 simulation (Langevin et al., 2017;
Pollock, 2016) plus an independently coded nonlinear chemistry/tracer model
that imports no Hydrosheaf code; (2) predeclared permutation-based adverse
controls that preserve each stream's marginal distribution while destroying
its case-specific meaning; (3) four linked experiments spanning topology
ranking, topology-conditioned age inference, reaction non-uniqueness, and a
real-data scope audit; (4) explicit, predeclared decision rules that
classify an evidence stream as complementary, redundant, or negative
transfer using paired discrimination and calibration metrics rather than
uncertainty alone (Brier, 1950; Davis and Goadrich, 2006).

## 1.5 Aim and objectives — 180 words

State the aim:

> To determine, using an independent and truth-blind synthetic benchmark,
> when combining hydraulic, age and hydrochemical evidence reduces
> interpretive non-uniqueness in groundwater topology, age and reaction
> inference, when it is redundant, and when a misspecified stream produces
> false confidence.

State four objectives, one per linked experiment:

1. Quantify the incremental topology-ranking contribution of hydraulic, age
   and chemical evidence, alone, in pairs and combined, under native and
   adverse (permuted) conditions.
2. Quantify how correct, partial and reversed topology assumptions change
   groundwater-age accuracy, interval width and coverage under an
   informative and a weak tracer regime.
3. Quantify reaction-family recovery under core and enhanced hydrochemical
   indicator panels, with carbonate processes reported separately.
4. Audit a real, data-limited aquifer system (Northern Ghana) to separate
   supportable component diagnostics from non-identifiable field claims.

## 1.6 Significance — 150 words

State the scientific significance: the study replaces an assumed
integration benefit with a measurable, falsifiable decision rule. State the
methodological significance: predeclared adverse controls make false
confidence a detectable, reportable outcome rather than an invisible failure
mode. State the practical significance: the results identify which evidence
combinations are worth collecting and which additional streams risk
degrading, not improving, an interpretation. State the Hydrosheaf
significance: M7 supplies the cross-layer guardrail that the M3–M6 papers
individually assume, closing the inference chain with a defensible answer to
whether integration itself is beneficial.

---

# 2. Methods — 1,550 words

## 2.1 Independent synthetic-truth generator and blinding — 250 words

Describe the external generator: official MODFLOW 6 and MODPATH 7 executables
(Langevin et al., 2017; Pollock, 2016) driven through FloPy, combined with an
independently coded nonlinear hydrochemistry and tracer model that imports no
Hydrosheaf code. State the development/locked-test split (six development
cases, twelve locked test cases), that development cases alone fit all
evidence-fusion coefficients, and that test truth is joined only after
truth-blind inference with no downstream tuning permitted. Report candidate
edge recall on the locked test set. Reference Supplementary Methods for full
generator and blinding-contract detail.

## 2.2 Evidence streams and topology-ranking fusion — 230 words

Define the three evidence streams: hydraulic log odds (H), uncertainty-aware
age compatibility (A), and a negative constrained-chemistry log objective
(C). Describe the seven unweighted logistic fusion models fitted on
development cases (H, A, C and all four combinations) and state why fitting
is unweighted (posterior entropy and calibration would be distorted by
artificial class balancing). Reference Supplementary Methods for the fusion
model specification.

## 2.3 Adverse controls: permutation design — 180 words

Define the four locked test conditions: native, age-permuted,
hydraulic-permuted, and jointly misspecified. Explain that within-case
permutation preserves each stream's marginal evidence distribution while
destroying its edge-specific meaning, making it a generic negative control
for misspecification rather than a model of one particular field failure
mode.

## 2.4 Topology-conditioned groundwater-age inference — 280 words

Describe the local ³H/³⁹Ar likelihood and lognormal age prior, the two
predeclared tracer regimes (informative: ³H and ³⁹Ar; weak: ³H only), and the
four topology conditions (no assumption, 50% partial true graph, complete
true graph, completely reversed graph). Describe the importance-sampling
reweighting by a soft downstream-older topology potential and state that its
effective sample size and log mean weight are reported so an incompatible
graph cannot silently appear precise. State particle count per case/regime.
Reference Supplementary Methods for the full importance-sampling formulation
and convergence rules.

## 2.5 Reaction-family bootstrap under evidence panels — 220 words

Describe scoring reaction mechanisms only on externally true flow edges after
inference, the lognormal analytical-noise bootstrap, and the two compared
hydrochemical panels (core; enhanced with additional indicator ions). State
that both panels retain PHREEQC-derived thermodynamic direction bounds
(Parkhurst and Appelo, 2013) and that carbonate weathering and precipitation
are reported separately rather than folded into an aggregate accuracy figure.

## 2.6 Northern Ghana data-scope audit — 170 words

Describe the audit as descriptive rather than a truth benchmark: the
workbook is checked for measured chemistry, environmental age tracers, screen
intervals, coordinates, hydraulic observations and independent connectivity
labels, and the permitted claim is restricted to which integrated
interpretations are supportable under available data.

## 2.7 Statistical reporting and decision rules — 220 words

State the reported metrics (PR-AUC, ROC-AUC, Brier score, log loss,
normalised edgewise entropy, calibration error, and, for age, MAE, bias, and
95% interval coverage), noting that PR-AUC is prioritised for a sparse
true-edge class (Davis and Goadrich, 2006) and Brier score/log loss jointly
characterise calibration (Brier, 1950). State the case-block bootstrap
(10,000 replicates) used for all paired contrasts and the three predeclared
decision rules distinguishing complementarity, redundancy, and negative
transfer. Reference Supplementary Methods for exact metric definitions and
the full decision-rule specification.

---

# 3. Results — 2,100 words

## 3.1 Benchmark scale and candidate recall — 150 words → **Table 1**

Report development/test case counts, particle and bootstrap counts, and
locked-test candidate recall.

## 3.2 Evidence-panel performance on native locked-test data — 350 words → **Figure 2a**, **Table 2**

Report PR-AUC, ROC-AUC, Brier score, log loss, entropy and calibration error
for H, A, C and all combinations. State that chemistry is the strongest
individual stream, hydraulics weaker, age alone below a useful ranking level,
and that hydraulics-plus-chemistry is the strongest pair.

## 3.3 Case-block evidence contrasts and adverse controls — 400 words → **Figure 2b–d**, **Table 3**

Report the paired case-block contrasts for native age, chemistry and
hydraulic addition, and for the three adverse (permuted) conditions. State
the central negative result (age does not improve topology ranking beyond
hydraulics and chemistry) and the central guardrail result (permuted streams
reduce entropy while worsening PR-AUC, Brier score and log loss).

## 3.4 Topology-conditioned groundwater-age inference — 350 words → **Figure 3**, **Table 4**

Report age MAE, coverage and interval-width contrasts for complete-true,
partial-true and reversed topology relative to no topology, under both
tracer regimes, including the importance-sampling ESS failure pattern for
the reversed/weak-tracer combination.

## 3.5 Reaction-family recovery under core and enhanced panels — 350 words → **Figure 4**, **Table 5**

Report modal-family accuracy, true-family probability, support entropy and
effective supported families per process and panel, with carbonate weathering
and precipitation reported separately.

## 3.6 Northern Ghana data-scope audit — 300 words → **Figure 5**, **Table 6**

Report which evidence types are available, absent, or masked in the
Northern Ghana workbook, and which integrated interpretations this evidence
does and does not support.

## 3.7 Summary decision table — 200 words

Synthesise the four experiments into the complementary/redundant/negative-
transfer classification for each contrast tested, without importing
Discussion-level explanation.

---

# 4. Discussion — 1,550 words

## 4.1 Evidence integration is conditional, not additive — 350 words

Answer RQ1/H1: chemistry and hydraulics are complementary topology evidence
in this generator; age is redundant in the reverse direction once hydraulics
and chemistry are present. Integration benefit therefore depends on which
streams are already available, not on the number of streams combined.

## 4.2 Topology assumptions change age uncertainty in both directions — 300 words

Answer RQ2–RQ3/H2–H3: correct topology conditionally reduces age
uncertainty, most when the tracer panel is weak; an incorrect topology
assumption is not a merely less-accurate alternative but can destabilise
inference outright, evidenced by importance-sampling degeneracy.

## 4.3 Entropy reduction without calibration is false confidence — 300 words

Answer RQ4/H4: the adverse controls show that reduced posterior uncertainty
can coexist with worse discrimination and calibration. Connect this to the
broader literature on overconfidence under model misspecification and argue
that uncertainty reduction alone is an unsafe stopping rule for evidence
integration.

## 4.4 Carbonate reactions remain non-identifiable regardless of panel richness — 250 words

Answer RQ3/H5: an enhanced indicator panel improves several redox/silicate
processes but not carbonate weathering or precipitation, illustrating that
low support entropy (confidence) can coexist with a wrong modal family.

## 4.5 What Ghana field data can and cannot support — 200 words

Answer RQ5/H6: state the defensible claim boundary for the Northern Ghana
data explicitly, consistent with the companion field-transfer analysis, and
reiterate that synthetic truth throughout the benchmark is model-conditioned,
not field truth.

## 4.6 Limitations — 150 words

State the single-generator/geometry scope, the model-conditioned nature of
all synthetic truth, the insensitivity of the predeclared univariate conflict
flag (reported as a failed heuristic rather than silently replaced), and the
descriptive-only status of the Ghana audit.

---

# 5. Conclusion — 250 words

Answer each of the five research questions directly in one to two sentences
each, then close with the decision-relevant statement that integration
benefit in groundwater inference must be demonstrated per evidence pairing
under predeclared adverse controls, not assumed from the number of streams
combined.

---

# Display-Item Plan

All items below are already generated, locked artifacts
(`figures/publication/`, `figures/supporting_validation/`,
`tables/publication/`); this outline registers them rather than requesting
new computation.

## Main figures

1. **Figure 1 — Benchmark architecture and claim boundary.** Independent
   synthetic-truth branch and Ghana claim-boundary branch.
2. **Figure 2 — Evidence integration is conditional.** Native panel
   performance; case-block incremental contrasts; adverse-control certainty/
   discrimination split.
3. **Figure 3 — Correct topology improves age inference.** Age MAE and
   interval-width changes; importance-sampling ESS; coverage changes.
4. **Figure 4 — Reaction recovery remains process-dependent.** Core versus
   enhanced panel modal accuracy, true-family probability, support entropy,
   effective support size.
5. **Figure 5 — Ghana data support component diagnostics, not complete
   field truth.** Evidence-to-claim map; M6 evidence-tier ablation link;
   truth-free seasonal hold-forward.

## Main tables

1. Locked benchmark design and computational scale.
2. Native locked-test performance for all seven evidence panels.
3. Case-block evidence contrasts, including adverse controls.
4. Topology-conditioned age contrasts.
5. Reaction-family recovery and non-uniqueness.
6. Ghana data scope and claim boundary.

## Supplementary figures

- Figure S1: locked synthetic model domain (representative realisation
  4101; synthetic model-space coordinates; not a geographic map).

## Supplementary tables

- Table S1: every evidence panel under every native/adverse condition.
- Table S2: all case-block evidence contrasts.
- Table S3: per-case topology-to-age sensitivity and importance-sampling ESS.
- Table S4: edgewise reaction support and entropy.
- Table S5: the predeclared conflict diagnostic, including its zero
  sensitivity.

## Supplementary Methods

Full generator and blinding contract; fusion-model specification and
coefficients; permutation-control construction; importance-sampling
formulation, topology-potential definition, and MCMC/ESS convergence rules;
PHREEQC constraint construction; reaction-dictionary and bootstrap detail;
full metric definitions (PR-AUC, Brier score, log loss, entropy, calibration
error, ESS); the three decision rules in full; and the Ghana audit criteria.
Target length: technically complete, not word-capped; expected to run
approximately 2,800–3,400 words of prose (excluding equations, headings, and
citations), matching the depth of four linked experiments.

---

# Claim Guardrails

Use:

> M7 evaluates whether combining independent hydraulic, age and hydrochemical
> evidence reduces groundwater interpretive non-uniqueness, using an
> independent synthetic-truth generator and predeclared adverse controls.

Do not use:

> M7 proves that integrating more evidence always improves groundwater
> interpretation, or that the Northern Ghana results validate field topology,
> age, or reaction truth.

Do not equate:

- reduced posterior entropy with improved inference;
- a synthetic MODFLOW 6/MODPATH 7 truth with field truth;
- Ghana chemistry/isotope evidence with residence-time, connectivity, or
  reaction-mechanism validation;
- low reaction-support entropy with correct mechanistic recovery.

---

# Data & Analysis Status — LOCKED (committed `d336e87`, `771388a`, `8ca2036`)

The M7.3 package is built and run at `M7/m7_nonuniqueness_benchmark/`. All
four experiments are executed and locked; no new computation is required to
write this manuscript.

## Implemented evidence

- Independent MODFLOW 6/MODPATH 7 generator (`independent_modflow_generator.py`)
  importing no Hydrosheaf code; six development cases (5201–5206), twelve
  locked test cases (5301–5312); candidate recall 1.000.
- Seven-panel logistic evidence-fusion benchmark with four adverse
  (permutation) conditions; 10,000-replicate case-block bootstrap.
- Topology-conditioned Bayesian age inference across four topology
  conditions and two tracer regimes; 50,000 particles per case/regime.
- Reaction-family bootstrap under core/enhanced panels; 64 bootstraps per
  case; PHREEQC direction bounds.
- Northern Ghana descriptive data-scope audit, cross-linked to the M6
  evidence-tier ablation.
- Six main figures, one supplementary figure, six main tables and five
  supplementary tables generated as vector PDF, 600-dpi PNG and
  LZW-compressed 300-dpi TIFF (figures) or CSV/Markdown (tables).

## Principal locked results

See `docs/m7_3_results.md` for the full decision table. Headline results:

- Chemistry is the strongest individual topology-ranking stream (PR-AUC
  0.459 native); age alone is below a useful ranking level (0.111).
- Age adds no positive incremental topology-ranking value beyond hydraulics
  and chemistry (HAC − HC PR-AUC = −0.0060, 95% CI [−0.0122, −0.0011]).
- Correct topology reduces age MAE by 0.062 y (³H+³⁹Ar) and 0.164 y
  (³H only) relative to no topology assumption; reversed topology worsens
  MAE by 0.282 y with two tracers and causes importance-sampling ESS
  failure in 8/12 tritium-only cases.
- Permuted (misspecified) evidence reduces entropy while worsening PR-AUC,
  Brier score, log loss and overconfident error in every adverse condition
  tested.
- Carbonate weathering and precipitation recovery remained 0/36 edges under
  both core and enhanced chemistry panels.
- The Northern Ghana workbook supports component chemical/isotopic
  diagnostics but not field validation of age, exact topology, or reaction
  truth.

---

# Archived 5,000-word section allocation (Nature Portfolio Article fallback)

Retain only if the target journal is later constrained to ~5,000 words.

| Section | Target words |
|---|---:|
| 1. Introduction | 800 |
| 2. Methods | 1,150 |
| 3. Results | 1,650 |
| 4. Discussion | 1,150 |
| 5. Conclusion | 250 |
| **Total** | **5,000** |

Do not use this allocation while the working target is Water Resources
Research or Journal of Hydrology.
