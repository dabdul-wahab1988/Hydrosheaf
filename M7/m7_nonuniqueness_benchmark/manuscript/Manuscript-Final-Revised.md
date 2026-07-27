# Evidence integration reduces groundwater interpretive non-uniqueness only conditionally in a single-generator benchmark with adverse controls

**Abstract**

Combining hydraulic, age, and hydrochemical evidence is widely assumed to
reduce groundwater interpretive non-uniqueness, but lower uncertainty does
not prove improvement. We asked, within one independent MODFLOW 6/MODPATH 7
generator sharing no code with the evaluated inference method, when evidence
integration is complementary, redundant, or produces false confidence. Four
linked experiments used predeclared permutation controls, six development and
twelve locked test cases, 50,000 age particles per case/regime, 64 reaction
bootstraps per case, and a descriptive Ghana audit. Chemistry was robustly
complementary topology-ranking evidence. Age and hydraulics each showed a
directional PR-AUC change that did not survive exact permutation testing with
multiplicity correction, although their entropy changes, and hydraulics'
calibration changes, did. Correct topology reduced age error by 0.062-0.164
years depending on tracer regime; reversed topology degraded accuracy or
destabilised inference. Carbonate reactions remained unrecovered under both
chemistry panels. Every adverse control lowered apparent uncertainty while
worsening discrimination and calibration, robustly to correction. The
Northern Ghana audit supported component diagnostics but not field truth for
topology, age, or unique reactions. Specific to this tested generator, an
integration benefit must be demonstrated per pairing against calibration and
adverse-control checks, not assumed.

## Introduction

### Groundwater interpretation from incomplete evidence

Hydraulic heads, environmental tracer concentrations, and hydrochemistry are
routinely combined to interpret where groundwater moves, how long it has
travelled, and which chemical reactions have shaped its composition along
the way. No single evidence type is individually sufficient for the
decisions that depend on these interpretations, including
contaminant-source attribution, aquifer-vulnerability assessment, and the
design of additional monitoring. A hydraulic head survey constrains
plausible flow directions but says little about residence time or reaction
history; an isotopic or radiometric tracer panel constrains age but is
silent on connectivity; and a hydrochemical survey constrains reaction
plausibility but cannot, on its own, establish which flow path produced the
observed composition. Practitioners therefore integrate hydraulic,
isotopic, and chemical evidence as a matter of routine practice, on the
general expectation that combining independent observations narrows the
range of interpretations consistent with the data.

### Current integration practice assumes more evidence helps

This expectation has methodological support from several directions,
although each addresses a narrower problem than combining hydraulic, age,
and chemical evidence for topology inference. Joint inversion of two
continuous geophysical fields measuring the same subsurface volume,
crosshole electrical resistance and ground-penetrating-radar traveltime,
has been reported to reduce model ambiguity relative to inverting either
field alone [@linde2006jointinversion]; whether that conclusion transfers
from jointly inverting two continuous physical fields to fusing discrete,
heterogeneous evidence streams for a classification task is not established
by that literature and is one of the questions this study addresses.
Multi-tracer environmental dating combines several radioisotopes, such as
tritium, krypton-85, and argon-39, to constrain groundwater age
distributions more tightly than any individual tracer permits, particularly
where the tracers cover complementary timescales [@visser2013multitracer].
Bayesian model averaging treats conceptual and parametric uncertainty
jointly across candidate models rather than committing to a single
best-fitting one [@neuman2003bma]. Equifinality research has separately
established that a single hydrological or geochemical model can fit
observed data well without being the correct model [@beven2001glue;
@beven2006manifesto]. Separately, posterior inference over directed graph
structure, the problem underlying this study's topology-conditioning
experiment, has been developed extensively in the Bayesian network and
causal-discovery literature using Markov-chain and generative-flow
approaches to sample plausible graphs under uncertainty; this study does
not attempt full posterior structure learning and instead uses a
deliberately simplified importance-sampling reweighting of a fixed particle
set, a design choice that trades posterior completeness for a fast,
auditable diagnostic (effective sample size and log mean weight) capable of
flagging an incompatible topology assumption, which is the property the
four linked experiments need rather than a complete graph posterior. Each
of the four literatures above evaluates one evidence type, one model
family, or one narrower structure-inference problem in isolation; none
directly tests whether adding a further evidence stream to an existing
interpretation is itself beneficial, redundant, or actively harmful to that
interpretation.

### Problem statement: uncertainty reduction is not proof of improvement

An added evidence stream can appear to resolve interpretive ambiguity
simply by lowering the posterior uncertainty attached to an interpretation,
without that reduction reflecting genuine new information about the
system. This can occur for at least two reasons: the added stream may be
numerically redundant with evidence already present, contributing no
information that was not already available; or it may be misspecified
relative to the true system, in a way that manufactures apparent confidence
rather than genuine resolution. Even well-established environmental
tracers carry systematic limitations that are not always visible from
goodness of fit alone [@mccallum2015limitations], and an inference pipeline
that reports posterior uncertainty without a paired predictive-accuracy or
calibration check has no way to distinguish a genuinely informative
addition from one of these failure modes. To our knowledge, a targeted
literature search did not
identify a groundwater benchmark that combines all three of an external,
code-independent truth generator, a locked development/test split, and
predeclared adverse controls capable of separating a genuinely
complementary evidence stream from one that is redundant or actively
harmful; this combination, rather than any one element in isolation, is
the specific gap this study addresses. In the absence of such controls, an
integrated interpretation that looks more confident may simply be more
confidently wrong, and a practitioner following standard reporting practice
would have no way to tell the difference.

### Novelty: a truth-blind conditional-integration benchmark with adverse controls

This study reframes evidence integration from an assumed benefit to a
testable, per-pairing property. The central methodological contribution is
a pair of predeclared devices that make false confidence detectable rather
than merely possible in principle: permutation-based adverse controls that
preserve each evidence stream's marginal distribution within a case while
destroying its case-specific meaning, and an explicit three-way decision
rule that classifies every tested evidence addition as complementary,
redundant, or negative transfer using paired discrimination and calibration
metrics [@davisgoadrich2006prcurves; @brier1950verification] rather than
posterior uncertainty alone, so that an entropy reduction can never, by
itself, count as evidence of an improved interpretation. This contribution
is supported by two further design choices rather than standing alongside
them as independently novel: an external synthetic-truth generator built on
official MODFLOW 6 groundwater-flow simulation and MODPATH 7 particle
tracking [@langevin2017modflow6; @pollock2016modpath], combined with an
independently coded nonlinear chemistry and tracer model, so that the
generator shares no code with the inference method it evaluates, and four
linked experiments spanning topology-ranking evidence integration,
topology-conditioned groundwater-age inference, reaction-mechanism
non-uniqueness, and a real-data scope audit, so that the conditional
nature of integration benefit could be examined across several evidence
layers rather than only one. Synthetic benchmarks built on official flow
simulators with independently coded observation models are established
practice in groundwater methods research; what this study adds to that
practice is the adverse-control and decision-rule layer described above.

### Aim and objectives

The aim of this study was to determine, using an independent and
truth-blind synthetic benchmark, when combining hydraulic, age, and
hydrochemical evidence reduces interpretive non-uniqueness in groundwater
topology, age, and reaction inference; when it is redundant; and when a
misspecified stream produces false confidence. Four objectives followed
directly from this aim. The first was to quantify the incremental
topology-ranking contribution of hydraulic, age, and chemical evidence,
alone, in pairs, and combined, under native and adverse conditions. The
second was to quantify how correct, partial, and reversed topology
assumptions change groundwater-age accuracy, interval width, and coverage
under an informative and a weak tracer regime. The third was to quantify
reaction-family recovery under core and enhanced hydrochemical indicator
panels, reporting carbonate processes separately rather than folding them
into an aggregate accuracy figure. The fourth was to audit a real,
data-limited aquifer system in Northern Ghana, where groundwater is a
primary potable resource drawn from a geologically unevenly distributed
resource base [@macdonald2012africa] and where, as the audit itself
documents, the available monitoring evidence is far from complete, to
separate supportable component diagnostics from non-identifiable field
claims.

### Significance

This study replaces an assumed integration benefit with a measurable,
falsifiable decision rule, directly addressing whether combining evidence
streams should be treated as an assumed good practice or as a claim that
requires its own supporting evidence. Predeclared adverse controls make
false confidence a detectable and reportable outcome rather than an
invisible failure mode of routine evidence fusion, a property of practical
relevance wherever an integrated interpretation informs a monitoring or
investigation decision, although translating the present benchmark's
metrics into a specific field decision protocol remains future work rather
than something demonstrated here. The results, obtained from one synthetic
generator and reported with that scope in mind, identify which evidence
combinations were worth collecting in this benchmark and which additional
streams risked degrading, rather than improving, an existing
interpretation. Taken together, the four linked experiments supply a
candidate cross-layer integration guardrail that groundwater studies
combining hydraulic, age, and chemical evidence can test against their own
data before treating a lower reported uncertainty as evidence of a better
interpretation.

## Methods

### Independent synthetic-truth generator and blinding

Interpretive non-uniqueness was benchmarked against an external synthetic
truth generator that shared no code with the inference method under test.
Groundwater flow and pathlines were simulated with the official MODFLOW 6
groundwater-flow model [@langevin2017modflow6] and the MODPATH 7
particle-tracking model [@pollock2016modpath], executed through FloPy.
Nonlinear hydrochemistry and environmental-tracer concentrations were
produced by an independently coded model that received no candidate-edge,
process, or age information from the inference side. This separation is
summarised in the benchmark architecture (Figure 1). Six development
cases fitted every evidence-fusion coefficient and classification threshold
used later; twelve locked test cases were withheld until after truth-blind
inference, and no locked-test outcome was permitted to alter a feature,
coefficient, condition, threshold, or stopping rule. The protocol, feature
definitions, and decision rules were frozen in version control (commit
d336e87) before any locked test case was generated, and this commit is the
verifiable record of the predeclaration referred to throughout this
section. Candidate edge generation retained every true test edge, giving a
locked-test candidate recall of 1.000. The full generator geometry,
chemistry archetypes, computational environment, and blinding contract are
reported in Supplementary Methods.

### Evidence streams and topology-ranking fusion

Three evidence streams were derived for every candidate directed edge: a
hydraulic log-odds term (H); an uncertainty-aware negative age-cost term (A);
and a negative constrained-chemistry log-objective term (C). Seven evidence
panels were formed from these streams and their combinations (H, A, C, HA,
HC, AC, HAC). For each panel, a logistic regression was fitted only on
development-case edges, using the panel's raw features standardised to
development-sample mean and standard deviation, with no class weighting.
Unweighted fitting was a deliberate choice: class-balanced weighting would
have distorted the posterior entropy and Brier/log-loss calibration values
reported later. For a panel with standardised features
$\tilde{x}_{e,1},\dots,\tilde{x}_{e,J}$ on edge $e$, fitted intercept
$\beta_0$, and fitted coefficients $\beta_1,\dots,\beta_J$, the fused
edge-truth probability was

<!-- EQ:EQ-1 -->
$$
p_e = \frac{1}{1 + \exp\!\left(-\beta_0 - \sum_{j=1}^{J} \beta_j\, \tilde{x}_{e,j}\right)}.
$$

Fitted coefficients from the six development cases were applied unchanged
to every locked test case.

### Adverse control construction

Four locked test conditions were scored for every panel: native evidence;
age-permuted; hydraulic-permuted; and joint misspecification. Permutation
was applied within each independent case, so every case's marginal
distribution of hydraulic or age evidence was preserved while its
edge-specific correspondence to that case's true topology was destroyed.
This made the adverse conditions a generic negative control for
misspecification rather than a simulation of one specific field failure
mode. A stream was accepted as complementary only when its native addition
improved discrimination or calibration with a 95% case-block bootstrap
interval excluding zero, or reduced entropy without a paired worsening of
Brier score and log loss; an addition that changed neither materially was
treated as redundant; and an addition that worsened skill or calibration, or
lowered uncertainty while raising overconfident error, was treated as
negative transfer. The full decision rule is stated in Supplementary
Methods.

### Topology-conditioned groundwater-age inference

Node-level groundwater age was inferred independently of any assumed
topology, from a local likelihood combining tritium and, in an informative
tracer regime, argon-39 against a lognormal age prior; a weak tracer regime
used tritium alone [@visser2013multitracer; @mccallum2015limitations]. Fifty
thousand common local age-posterior particles were drawn per case and
tracer regime. These particles were then reweighted under four topology
conditions: no connectivity assumption; a five-of-nine partial true graph
(55.6% of true edges); the
complete true graph; and the completely reversed graph. Reweighting applied
a soft downstream-older importance potential to every assumed edge, so that
an incompatible topology assumption degraded the effective sample size
rather than silently producing an artificially narrow posterior. Effective
sample size and log mean importance weight were reported alongside age
mean absolute error (MAE), bias, 95% interval coverage, and interval width,
so that an unstable reweighting could be identified rather than
misinterpreted as a precise result. The exact likelihood, prior, and
importance-weight specification are given in Supplementary Methods.

### Reaction-family bootstrap under evidence panels

Reaction-mechanism recovery was scored only on candidate edges confirmed as
true flow paths by the external generator, so that mechanism inference was
not credited or penalised using topologically incorrect edges. Chemistry
was perturbed by 3% lognormal analytical noise across 64 bootstrap
replicates per case, under two hydrochemical indicator panels: a core panel
and an enhanced panel with additional indicator ions. Both panels retained
thermodynamic reaction-direction bounds derived from PHREEQC speciation and
saturation calculations [@parkhurst2013phreeqc]. Modal-family accuracy, the
probability assigned to the true reaction family, family-support entropy,
and the effective number of supported families were computed per edge and
aggregated by process. The reaction solver was evaluated inference, not part
of the synthetic generator, preserving generator-inference separation.
Carbonate weathering and carbonate
precipitation reported separately from the aggregate rather than folded
into it.

### Northern Ghana data-scope audit

The Northern Ghana workbook was audited descriptively rather than treated
as a truth benchmark. The audit recorded whether measured hydrochemistry,
stable water isotopes, environmental age tracers, screen intervals,
repeated head observations, site coordinates, and independent connectivity
labels were present, absent, or masked, and mapped each available evidence
type to a defensible use. The audit's permitted objective was to determine
which integrated interpretations are supportable under the available
Ghanaian evidence and which remain non-identifiable, not to validate field
topology, age, or reaction truth.

### Statistical reporting and decision rules

Topology-ranking discrimination was reported as precision-recall area under
the curve (PR-AUC) alongside receiver-operating-characteristic area under
the curve (ROC-AUC). PR-AUC was prioritised as the primary discrimination
metric because the locked-test candidate edge set is class-imbalanced and
because Davis and Goadrich (2006) argued that PR-AUC is more informative
than ROC-AUC under imbalance [@davisgoadrich2006prcurves]; more recent work
has questioned how generally that preference holds, so PR-AUC priority is
reported here as a stated design choice rather than an uncontested
property, and ROC-AUC is retained throughout as a check. Calibration was
reported using Brier score [@brier1950verification], log loss, and a
ten-bin expected calibration error, so that any reduction in posterior
uncertainty (normalised edgewise Bernoulli entropy) could be checked
against predictive accuracy rather than accepted on its own. All paired
contrasts used a case-block bootstrap with 10,000 replicates, resampling
whole independent MODFLOW cases and reporting 95% percentile confidence
intervals; because this yields twenty-four such intervals across the six
predeclared contrasts and four metrics reported in Table 3, an exact
permutation test and a Benjamini-Hochberg false-discovery-rate correction
across that full family of twenty-four tests were performed as a
post-hoc robustness check and are reported in Supplementary Methods and
Supplementary Table S6. This design followed the general expectation that
combining independent evidence sources can reduce model ambiguity
[@linde2006jointinversion; @neuman2003bma], but tested that expectation as a
falsifiable, per-pairing decision rule rather than assuming it. Complete
metric definitions, the importance-sampling stability rule, and the full
complementarity/redundancy/negative-transfer decision rule are reported in
Supplementary Methods.

### Estimator clarification and reviewer-requested sensitivities

PR-AUC was computed on all candidate edges pooled within each test case and
then contrasted by resampling the 12 independent cases as blocks. This preserves
the edge-ranking estimand while preventing edges within a generated aquifer
from being treated as independent replicates. The reported interval is the
percentile interval of the paired case-block contrast; direct pooled-edge
contrasts are supplied only as descriptive checks. Posterior graph sampling by
generative flow networks provides a relevant modern alternative
[@deleu2022daggflownet], but this benchmark intentionally evaluates a fixed,
auditable candidate-edge scoring task rather than full graph-posterior learning.

The effective-sample-size threshold was varied from 200 to 400 and 1,000. All
conditions remained stable in every case except reversed topology: under
tritium only, the stable-case fraction changed from 0.583 to 0.333 and 0.083;
under two tracers it changed from 1.000 to 1.000 and 0.833. Thus, the qualitative
diagnosis of reversed-topology incompatibility is stable and becomes more
conservative as the threshold rises.

Logistic L2 regularisation was varied over C = 0.1, 0.3, 1, 3, and 10. Native
HAC PR-AUC changed narrowly from 0.477 to 0.486; Brier score improved from
0.096 to 0.083 and log loss from 0.305 to 0.247. Standardised H/A/C
coefficients are archived for every C in
`results/m7_3_locked/review_sensitivity.csv`, so the stability claim is not
based on performance values alone.

The probability-span conflict threshold was varied from 0.10 to 0.50. At 0.10
it flagged 39.8% of native edges and 59.5%-61.2% under adverse conditions, but
its error rate remained condition-dependent; thresholds of 0.30 or higher
flagged no edges. This confirms that the originally predeclared univariate
diagnostic is not a reliable standalone safeguard. The 3% reaction-noise level
and five-year topology order scale remain predeclared fixed settings, not
validated universal constants, and are retained as explicit single-generator
limitations.

A reviewer-requested case-block bootstrap placed the enhanced-minus-core modal
family accuracy difference at 0.0556 (95% CI 0.0278-0.0833; 12 cases). Under
tritium-only evidence, reversed topology also produced a falsely narrow mean
interval relative to complete topology (-3.394 years; 95% CI -5.777 to
-1.069), while failing the ESS diagnostic in most cases. Values in prose are
rounded only after analysis; CSV outputs retain full precision.

Recent hydrochemical metamodelling has used machine-learning ensembles to
estimate groundwater-age distributions [@tschritter2023age]. That work predicts
ages from chemistry, whereas the present benchmark tests whether adding an
evidence stream improves or harms a separately scored topology, age, or
reaction inference under adverse controls.

## Results

### Benchmark scale and candidate recall

The locked benchmark comprised six development cases (seeds 5201-5206),
used only to fit evidence-fusion coefficients and any classification
threshold, and twelve independent locked test cases (seeds 5301-5312),
generated by the official MODFLOW 6/MODPATH 7 simulator combined with an
independently coded nonlinear chemistry/tracer model that shared no code
with the inference method it evaluated. The benchmark architecture is
summarised in Figure 1.

![Figure 1. Benchmark architecture and claim boundary. The independent synthetic-truth branch separates official MODFLOW 6/MODPATH 7 flow and pathline generation, nonlinear chemistry and tracer generation, truth-blind inference, and locked adverse-control scoring from the Northern Ghana branch, which distinguishes supportable component diagnostics from field interpretations that remain non-identifiable under the available evidence.](figures/publication/figure1_benchmark_and_claim_design.png)

Topology-conditioned age inference used 50,000 importance-sampling
particles per case and tracer regime; reaction-family inference used 64
chemistry-perturbation bootstrap replicates per case; and every paired
contrast used a 10,000-replicate case-block bootstrap resampling whole
independent cases, as summarised in Table 1.

**Table 1.** Locked benchmark design and computational scale, including development and locked-test case counts, age-importance particle count, reaction-bootstrap count, case-block bootstrap count, and locked-test candidate recall.

| Design item | Value | Scope |
| --- | --- | --- |
| Development cases | 6 | 5201–5206 |
| Locked test cases | 12 | 5301–5312 |
| Age particles | 50000 | per case/regime |
| Reaction bootstraps | 64 | per test case |
| Case-block bootstraps | 10000 | independent case resampling |
| Candidate recall | 1.000 | locked test |

Candidate edge generation retained every true test edge, giving a
locked-test candidate recall of 1.000, so none of the incremental-value
contrasts reported below could be confounded by a missing true edge in the
candidate set.

### Evidence-panel performance on native locked-test data

On native locked-test data, chemistry was the strongest individual
topology-ranking evidence stream (PR-AUC 0.459, ROC-AUC 0.892, Brier score
0.094, log loss 0.295, mean edge entropy 0.455), hydraulics was
substantially weaker (PR-AUC 0.176, ROC-AUC 0.650, Brier score 0.108), and
age alone was below a useful ranking level (PR-AUC 0.111, ROC-AUC 0.452),
only marginally distinguishable from the base rate implied by chance
ranking (see Table 2 below). Among the three evidence pairs,
hydraulics-plus-chemistry (HC) was the strongest (PR-AUC 0.485, ROC-AUC
0.909, Brier score 0.088, log loss 0.268, mean edge entropy 0.422),
exceeding both age-plus-chemistry (AC; PR-AUC 0.471, Brier score 0.089) and
hydraulics-plus-age (HA; PR-AUC 0.111, effectively unchanged from age alone
and far below either single-stream chemistry result). The fully integrated
seven-feature panel (HAC) achieved PR-AUC 0.480, ROC-AUC 0.908, Brier score
0.088, log loss 0.268, mean edge entropy 0.421, and expected calibration
error 0.060, a profile marginally below HC on PR-AUC rather than above it,
despite combining strictly more evidence, as shown in Figure 2.

**Table 2.** Native locked-test evidence-panel performance for all seven hydraulic (H), age (A) and chemistry (C) panel combinations, reporting PR-AUC, ROC-AUC, Brier score, log loss, mean edge entropy and expected calibration error.

| Panel | PR-AUC | ROC-AUC | Brier | Log loss | Edge entropy | ECE |
| --- | --- | --- | --- | --- | --- | --- |
| H | 0.176 | 0.650 | 0.108 | 0.347 | 0.503 | 0.006 |
| A | 0.111 | 0.452 | 0.109 | 0.356 | 0.524 | 0.012 |
| C | 0.459 | 0.892 | 0.094 | 0.295 | 0.455 | 0.078 |
| HA | 0.111 | 0.453 | 0.108 | 0.347 | 0.504 | 0.005 |
| HC | 0.485 | 0.909 | 0.088 | 0.268 | 0.422 | 0.060 |
| AC | 0.471 | 0.904 | 0.089 | 0.274 | 0.436 | 0.069 |
| HAC | 0.480 | 0.908 | 0.088 | 0.268 | 0.421 | 0.060 |

![Figure 2. Evidence integration is conditional. Panel a compares native hydraulic (H), age (A), chemistry (C), pairwise and fully integrated topology-ranking performance. Panel b shows case-block incremental PR-AUC contrasts with 95% bootstrap confidence-interval error bars. Panels c and d separate the misspecification stress test into certainty and discrimination shifts, with error bars representing 95% bootstrap confidence intervals throughout, showing that permuted evidence can reduce entropy while degrading PR-AUC.](figures/publication/figure2_evidence_integration.png)

The full native and adverse condition grid for all seven panels, extending
this comparison beyond the headline HAC/HC/HA/AC contrasts, is reported in
Supplementary Table S1.

### Case-block evidence contrasts and adverse controls

Paired case-block contrasts isolated the incremental value of each stream
within the fully integrated panel across the twelve locked test cases
(Table 3). Adding age to hydraulics-plus-chemistry produced a
small decrease in PR-AUC (HAC minus HC = -0.006) whose bootstrap interval
excluded zero (95% CI [-0.0122, -0.0011]), with log loss and Brier score
changes both including zero, alongside an entropy reduction (-0.0006) that
did not offset the ranking loss. An exact permutation test on the same
twelve paired differences, corrected for the full family of contrasts
tested (Supplementary Table S6), found that this PR-AUC contrast did not
survive multiplicity correction (adjusted p = 0.070), whereas the
accompanying entropy reduction did (adjusted p = 0.004); age's incremental
effect on topology ranking is therefore better described as a robust
entropy reduction unaccompanied by a statistically robust change in
discrimination, precisely the pattern the adverse-control framework below
treats as insufficient evidence of improvement. Adding chemistry to
hydraulics-plus-age produced a large positive PR-AUC gain (HAC minus HA =
+0.447, 95% CI [0.357, 0.540]) that survived multiplicity correction
(adjusted p < 0.01), with Brier score and log loss both improving
(-0.0196, 95% CI [-0.0213, -0.0176]; and -0.0791, 95% CI [-0.0850,
-0.0720], respectively) and entropy falling by 0.0827, all likewise robust
to correction. Adding hydraulics to age-plus-chemistry produced a smaller
positive PR-AUC gain (HAC minus AC = +0.0091, 95% CI [0.001, 0.020]) that,
like the age contrast, did not survive multiplicity correction (adjusted
p = 0.070), while its Brier score, log loss, and entropy improvements all
did (adjusted p < 0.01 each); hydraulics' incremental value is accordingly
supported by calibration rather than by a statistically robust ranking
change at this sample size.

**Table 3.** Case-block evidence contrasts, including native incremental-evidence additions and adverse (permuted) controls, each reported as a paired 10,000-replicate case-block bootstrap mean difference with 95% confidence interval.

| Contrast | PR-AUC Δ [95% CI] | Brier Δ [95% CI] | Log-loss Δ [95% CI] | Entropy Δ [95% CI] |
| --- | --- | --- | --- | --- |
| Native age added | -0.006 [-0.012, -0.001] | 0.0001 [-0.0000, 0.0001] | 0.0000 [-0.0001, 0.0002] | -0.0006 [-0.0009, -0.0003] |
| Native chemistry added | 0.447 [0.357, 0.540] | -0.0196 [-0.0213, -0.0176] | -0.0791 [-0.0850, -0.0720] | -0.0827 [-0.0985, -0.0653] |
| Native hydraulics added | 0.009 [0.001, 0.020] | -0.0010 [-0.0012, -0.0008] | -0.0059 [-0.0070, -0.0050] | -0.0146 [-0.0185, -0.0118] |
| Permuted age added | -0.075 [-0.135, -0.015] | 0.0034 [0.0012, 0.0055] | 0.0105 [0.0030, 0.0176] | -0.0207 [-0.0236, -0.0175] |
| Permuted hydraulics added | -0.069 [-0.112, -0.027] | 0.0110 [0.0055, 0.0164] | 0.0745 [0.0404, 0.1091] | -0.0482 [-0.0593, -0.0372] |
| Age + hydraulics misspecified | -0.139 [-0.204, -0.074] | 0.0106 [0.0055, 0.0163] | 0.0730 [0.0401, 0.1073] | -0.0706 [-0.0827, -0.0586] |

Every predeclared adverse control moved in the same qualitative direction:
entropy fell while discrimination and calibration worsened, and every one
of these degradations survived multiplicity correction (Supplementary
Table S6, adjusted p < 0.05 throughout). Permuting age within each case
reduced mean edge entropy by 0.0207 (95% CI [-0.0236, -0.0175]) while
reducing PR-AUC by 0.075 (95% CI [-0.135, -0.015]) and increasing Brier
score by 0.0034. Permuting hydraulics reduced entropy by 0.0482 while
reducing PR-AUC by 0.069 (95% CI [-0.112, -0.027]) and increasing log loss
by 0.0745. Joint misspecification of both streams reduced entropy the most
of any condition (-0.0706, 95% CI [-0.0827, -0.0586]) while reducing
PR-AUC by 0.139 (95% CI [-0.204, -0.074]) and raising overconfident error
by 0.0387. The predeclared univariate probability-span conflict diagnostic
(threshold 0.50) flagged no edges as conflicting under any of the four
conditions, including joint misspecification, despite the clear
degradation that the case-level discrimination and calibration metrics
detected in every adverse condition; this comparison is reported in full
in Supplementary Table S5.

### Topology-conditioned groundwater-age inference

Correct and partial topology assumptions narrowed age-posterior intervals
without a coverage penalty in both tracer regimes (Figure 3,
Table 4). Under the informative (³H + ³⁹Ar) regime, complete-true
topology reduced age mean absolute error (MAE) by 0.062 years relative to no
topology assumption (95% CI [-0.071, -0.053]), narrowed 95% intervals by
0.252 years (95% CI [-0.281, -0.222]), and left coverage effectively
unchanged (+0.014, 95% CI [0.000, 0.035]); the five-of-nine partial-true graph
(55.6%) gave a
smaller benefit (MAE -0.025 years, interval width -0.145 years). Under the
weak (³H-only) regime, the same complete-true contrast reduced MAE by 0.164
years (95% CI [-0.184, -0.145]) and narrowed intervals by 0.912 years (95%
CI [-0.964, -0.863]), a substantially larger absolute benefit than under
the informative regime, consistent with topology contributing information
that a single-tracer panel cannot supply on its own.

![Figure 3. Correct topology improves age inference. Panels a and b show changes in age MAE and 95% interval width for correct, partial and reversed topology assumptions, with 95% bootstrap confidence-interval error bars on each contrast. Panel c reports importance-sampling effective-sample-size fractions (point values, no interval), including the reversed-graph failures. Panel d shows coverage changes with confidence-interval error bars. Reversed tritium-only accuracy contrasts are not interpreted because 8 of 12 cases failed the predeclared effective-sample-size rule.](figures/publication/figure3_topology_conditions_age.png)

**Table 4.** Topology-conditioned age contrasts for complete-true, partial-true and reversed topology relative to no topology assumption, reported for both the informative (³H+³⁹Ar) and weak (³H-only) tracer regimes.

| Tracer regime | Comparison | Metric | Difference [95% CI] | Cases |
| --- | --- | --- | --- | --- |
| ³H + ³⁹Ar | Complete true − none | Age MAE (years) | -0.062 [-0.071, -0.053] | 12 |
| ³H + ³⁹Ar | Complete true − none | 95% coverage | 0.014 [0.000, 0.035] | 12 |
| ³H + ³⁹Ar | Complete true − none | 95% interval width (years) | -0.252 [-0.281, -0.222] | 12 |
| ³H + ³⁹Ar | Partial true − none | Age MAE (years) | -0.025 [-0.040, -0.009] | 12 |
| ³H + ³⁹Ar | Partial true − none | 95% coverage | -0.007 [-0.028, 0.014] | 12 |
| ³H + ³⁹Ar | Partial true − none | 95% interval width (years) | -0.145 [-0.167, -0.123] | 12 |
| ³H + ³⁹Ar | Reversed − complete true | Age MAE (years) | 0.282 [0.016, 0.545] | 12 |
| ³H + ³⁹Ar | Reversed − complete true | 95% coverage | -0.028 [-0.049, -0.007] | 12 |
| ³H + ³⁹Ar | Reversed − complete true | 95% interval width (years) | 0.296 [0.102, 0.469] | 12 |
| ³H only | Complete true − none | Age MAE (years) | -0.164 [-0.184, -0.145] | 12 |
| ³H only | Complete true − none | 95% coverage | 0.014 [0.000, 0.035] | 12 |
| ³H only | Complete true − none | 95% interval width (years) | -0.912 [-0.964, -0.863] | 12 |
| ³H only | Partial true − none | Age MAE (years) | -0.057 [-0.074, -0.040] | 12 |
| ³H only | Partial true − none | 95% coverage | 0.000 [-0.021, 0.021] | 12 |
| ³H only | Partial true − none | 95% interval width (years) | -0.456 [-0.489, -0.422] | 12 |
| ³H only | Reversed − complete true | Age MAE (years) | 0.024 [-0.889, 0.888] | 12 |
| ³H only | Reversed − complete true | 95% coverage | -0.056 [-0.090, -0.021] | 12 |
| ³H only | Reversed − complete true | 95% interval width (years) | -3.394 [-5.777, -1.069] | 12 |

Reversed topology behaved asymmetrically to this benefit rather than as a
merely weaker version of it. Under the informative regime, reversed
topology increased MAE by 0.282 years relative to complete-true topology
(95% CI [0.016, 0.545]), widened 95% intervals by 0.296 years, and reduced
coverage by 0.028 (95% CI [-0.049, -0.007]). Under the weak regime, the
reversed-topology importance-sampling effective sample size fell below the
predeclared stability rule (ESS < 400 of 50,000 particles) in 8 of the 12
locked test cases; the nominal reversed-minus-correct MAE contrast under
this regime carried a 95% interval of [-0.889, 0.888], spanning almost two
years in each direction, and is reported for completeness in Table 4 but is
not interpreted as a stable estimate given this diagnostic failure.
Per-case topology-to-age sensitivity and effective-sample-size values
underlying this comparison are reported in full in Supplementary Table S3.

### Reaction-family recovery under core and enhanced panels

Reaction-family recovery was strongly process-dependent under both the core
and enhanced hydrochemical panels (Figure 4, Table 5).
Denitrification and sulfate reduction were recovered reliably in both
panels (modal accuracy 1.000 for both processes under both panels; true-
family probability 0.994 and 0.997 respectively under the core panel, and
0.969 and 0.986 under the enhanced panel, each with low support entropy).
Silicate weathering was also recovered well under both panels (accuracy
1.000, mean true-family probability rising from 0.864 under the core panel
to 0.951 under the enhanced panel, with support entropy falling from 0.349
to 0.161). Iron reduction was mixed and improved materially with the
enhanced panel, from 0.000 modal accuracy and 0.341 support entropy under
the core panel to 0.500 accuracy and 0.469 true-family probability under the
enhanced panel. In contrast, carbonate weathering and carbonate
precipitation were recovered in neither panel: modal accuracy remained
0.000 for both processes under both the core and enhanced panels, across
the 36 scored edges per panel (24 carbonate weathering and 12 carbonate
precipitation edges).

![Figure 4. Reaction recovery remains process-dependent. The 2x2 panels compare core and enhanced evidence panels for modal-family recovery, probability assigned to the true family, support entropy and effective support size. Carbonate weathering and precipitation remain unrecovered despite low support entropy, demonstrating that numerical stability can be confidently wrong.](figures/publication/figure4_reaction_nonuniqueness.png)

**Table 5.** Reaction-family recovery and non-uniqueness by process and hydrochemical panel (core; enhanced), reporting edge counts, modal-family accuracy, true-family probability, support entropy and effective supported families.

| Panel | Process | Edges | Modal accuracy | True-family probability | Support entropy | Effective families |
| --- | --- | --- | --- | --- | --- | --- |
| Core | Carbonate Precipitation | 12 | 0.000 | 0.000 | 0.034 | 1.042 |
| Core | Carbonate Weathering | 24 | 0.000 | 0.000 | 0.116 | 1.138 |
| Core | Denitrification | 24 | 1.000 | 0.994 | 0.019 | 1.023 |
| Core | Iron Reduction | 12 | 0.000 | 0.000 | 0.341 | 1.449 |
| Core | Silicate Weathering | 24 | 1.000 | 0.864 | 0.349 | 1.466 |
| Core | Sulfate Reduction | 12 | 1.000 | 0.997 | 0.013 | 1.014 |
| Enhanced | Carbonate Precipitation | 12 | 0.000 | 0.000 | 0.000 | 1.000 |
| Enhanced | Carbonate Weathering | 24 | 0.000 | 0.000 | 0.090 | 1.110 |
| Enhanced | Denitrification | 24 | 1.000 | 0.969 | 0.062 | 1.080 |
| Enhanced | Iron Reduction | 12 | 0.500 | 0.469 | 0.526 | 1.715 |
| Enhanced | Silicate Weathering | 24 | 1.000 | 0.951 | 0.161 | 1.199 |
| Enhanced | Sulfate Reduction | 12 | 1.000 | 0.986 | 0.058 | 1.064 |

Overall modal-family accuracy rose from 0.556 under the core panel to 0.611
under the enhanced panel across all 108 scored edges. A post-hoc case-block
bootstrap estimated the paired increase as 0.0556 (95% CI 0.0278-0.0833;
12 independent cases). Because this interval was requested after review rather
than predeclared, it is labelled post hoc and does not alter the locked primary
decision rule.
Several carbonate rows showed low family-support entropy (as low as 0.000
for enhanced-panel carbonate precipitation) despite zero accuracy,
indicating that the reaction solver was frequently confident in an
incorrect family rather than merely undecided between candidate families.
Edgewise detail underlying this pattern is reported in Supplementary Table
S4.

### Northern Ghana data-scope audit

The canonical Northern Ghana workbook audit found major hydrochemistry and
stable water isotopes available, alongside a single static water-level
measurement per well, but no independent aquifer-type classification, no
processed graph edges, no environmental age-tracer panel, no screen-interval
data, no repeated (time-varying) head series, no intra-season sampling-date
field, and masked rather than exact site coordinates (Figure 5,
Table 6).

![Figure 5. Ghana data support component diagnostics, not complete field truth. Panels a and b map observed evidence to defensible claims. Panel c shows a companion study's predeclared five-tier field-transfer ablation of reaction identifiability as diagnostic chemistry is withdrawn; panel d shows a truth-free supporting seasonal hold-forward test. Synthetic panels elsewhere use 12 locked cases; these Ghana panels use 160 wet-to-dry well pairs and 140 complete pairs, respectively.](figures/publication/figure5_ghana_supportability_boundary.png)

**Table 6.** Northern Ghana data scope and claim boundary, listing evidence availability status and the corresponding defensible use for each evidence type.

| Evidence | Status | Defensible use |
| --- | --- | --- |
| Major hydrochemistry | Available | Component inference and QC |
| Stable water isotopes | Available | Recharge/source evidence |
| Environmental age tracers | Absent | Residence time non-identifiable |
| Screen intervals | Absent | Vertical connectivity non-identifiable |
| Repeated head series | Absent | Dynamic head validation unavailable |
| Coordinates | Masked | No site-scale connectivity truth |
| Processed graph edges | Absent | Edge sets are self-generated, not supplied |
| Intra-season sampling dates | Absent | Batch order is disclosed and arbitrary |
| Independent aquifer-type classification | Absent | Stratified reporting uses region instead |
| Independent reaction truth | Absent | Unique mechanisms unvalidated |

Consequently, the audit supported component-level chemical and isotopic
diagnostics, within-campaign seasonal chemistry hold-forward under a
disclosed arbitrary well-revelation order, reaction-family plausibility and
equivalence-class reporting, sensitivity to alternative assumed edge sets,
and explicit non-identifiability mapping.
Residence-time inversion, exact directed flow paths, screen-resolved
vertical connectivity, unique reaction mechanisms, and a fully observed
field digital twin were each classified as non-identifiable or unsupported
from this dataset, and no result from Experiments 1-3 was used, or could be
used, to validate any of these field-scale claims.

### Summary decision table

Table 7 synthesises the classification of every tested contrast
across the four experiments under the predeclared decision rule, refined
by the multiplicity-corrected re-analysis reported above and in
Supplementary Table S6. Chemistry was classified as complementary
topology-ranking evidence on every metric; hydraulics was complementary on
calibration and entropy but not, after correction, on PR-AUC alone; and age
was likewise supported only by a robust entropy reduction rather than by a
statistically robust ranking change. Correct and partial topology
assumptions were classified as complementary age evidence, conditional on
tracer regime and strongest under a weak tracer panel, while reversed
topology was classified as negative transfer under the informative regime
and as numerically unstable, rather than interpretably worse, under the
weak regime. Every one of the three adverse (permutation) conditions was
classified as negative transfer under the predeclared decision rule and
survived multiplicity correction on every reported metric. The enhanced
chemistry panel was classified as a partial, process-dependent improvement
for reaction-family recovery, benefiting denitrification, sulfate
reduction, silicate weathering, and iron reduction, while carbonate
weathering and precipitation remained non-identifiable under either panel.
The Northern Ghana audit supported component diagnostics and non-
identifiability mapping only, with no extension to field topology, age, or
reaction-truth validation.

**Table 7.** Summary decision table synthesising the complementary, redundant, or negative-transfer classification of every tested contrast across the four experiments, refined by the multiplicity-corrected re-analysis in Supplementary Table S6.

| Experiment | Contrast or condition | Classification | Basis |
| --- | --- | --- | --- |
| Evidence integration | Chemistry added to hydraulics-plus-age | Complementary | PR-AUC and calibration gains survive multiplicity correction |
| Evidence integration | Hydraulics added to age-plus-chemistry | Complementary (calibration); PR-AUC not distinguishable from zero after correction | Calibration and entropy gains survive correction; PR-AUC gain does not |
| Evidence integration | Age added to hydraulics-plus-chemistry | Redundant (entropy); PR-AUC not distinguishable from zero after correction | Entropy reduction survives correction; PR-AUC decrease does not |
| Evidence integration | Age permuted | Negative transfer | PR-AUC/Brier/log-loss degradation and entropy reduction all survive correction |
| Evidence integration | Hydraulics permuted | Negative transfer | PR-AUC/log-loss degradation and entropy reduction all survive correction |
| Evidence integration | Age and hydraulics jointly misspecified | Negative transfer | PR-AUC/Brier/log-loss degradation and entropy reduction all survive correction |
| Topology-conditioned age | Complete-true topology vs none (informative regime) | Complementary | MAE and interval-width improvement with unchanged coverage |
| Topology-conditioned age | Complete-true topology vs none (weak regime) | Complementary (larger effect) | MAE and interval-width improvement several times larger than informative regime |
| Topology-conditioned age | Partial-true topology vs none | Complementary (smaller effect) | Smaller but directionally consistent MAE and interval-width improvement |
| Topology-conditioned age | Reversed topology vs complete-true (informative regime) | Negative transfer | MAE increase and coverage decrease |
| Topology-conditioned age | Reversed topology vs complete-true (weak regime) | Numerically unstable (not interpreted) | Importance-sampling ESS failure in 8 of 12 cases |
| Reaction-family recovery | Enhanced vs core panel (denitrification, sulfate reduction, silicate weathering, iron reduction) | Partial improvement (descriptive) | Accuracy and true-family probability improve; no case-block interval predeclared |
| Reaction-family recovery | Enhanced vs core panel (carbonate weathering and precipitation) | Non-identifiable under either panel | Zero modal accuracy under both panels despite low support entropy |
| Northern Ghana audit | Component chemical and isotopic diagnostics | Supportable | Major hydrochemistry and stable isotopes available |
| Northern Ghana audit | Field topology, age, and reaction-mechanism validation | Not supportable | Environmental age tracers, screen intervals, repeated head series, and independent connectivity truth all absent |

## Discussion

### Evidence integration is conditional, not additive

The first research question asked whether combining hydraulic, age, and
chemical evidence reduces topological non-uniqueness in this benchmark, and
whether every pairing behaves alike. It does not. Chemistry and hydraulics
behaved as complementary topology-ranking evidence in this generator:
adding chemistry to hydraulics-plus-age produced the largest PR-AUC gain
observed in the study (+0.447), and adding hydraulics to age-plus-chemistry
produced a smaller gain (+0.009) whose calibration and entropy components
survived multiplicity correction even though its PR-AUC component alone did
not (Supplementary Table S6). This complementarity is broadly consistent
with the general expectation that independent data types can reduce model
ambiguity when combined, an expectation drawn mainly from joint inversion
of continuous geophysical fields [@linde2006jointinversion] and shown here
to extend, at least in this generator, to discrete evidence-panel fusion
for a classification task. Age behaved differently: adding age to
hydraulics-plus-chemistry produced a small decrease in PR-AUC whose
bootstrap interval excluded zero but whose exact permutation test did not
survive correction for the twenty-four contrasts tested (adjusted
p = 0.070), while the accompanying entropy reduction did survive
(adjusted p = 0.004). The most defensible statement is therefore that, in
this generator, age's incremental effect on topology ranking is supported
as an entropy reduction but not as a statistically robust change in
discrimination, which is itself the pattern the adverse-control framework
below treats as insufficient evidence of genuine improvement. This
combination of outcomes is inconsistent with a model in which integration
benefit scales with the number of evidence streams combined; it instead
depends on which information a candidate stream contributes beyond what is
already available, consistent with model-averaging accounts of evidence
combination that weight sources by their distinct information content
rather than their count [@neuman2003bma].

A rival explanation deserves consideration before this null result is
attributed to age evidence in general: the age stream entering the fusion
model here is a single engineered feature, an uncertainty-aware negative
age-cost term, and its failure to add discriminative value could reflect a
limitation of this specific feature construction rather than of age
information as such. This concern is not hypothetical; a related result
from the same programme's earlier confirmatory audit found an age-based
compatibility gate that correctly screened false candidate edges without
improving the primary ranking metric, consistent with a feature-design
bottleneck rather than an evidentiary one. This study's design does not
distinguish between these two explanations, and doing so would require
testing an alternative age-feature construction, which is left as a
direction for the feature-engineering choice to be examined rather than
treated as settled. In this generator, where chemistry already carried
most of the edge-discriminating signal (native PR-AUC 0.459 for chemistry
alone versus 0.111 for age alone), the age evidence's clearest value lay in
constraining groundwater age itself rather than in ranking topology, and
extending an existing hydraulic-and-chemistry survey with an age-dating
campaign should not, on this evidence alone, be assumed to sharpen
topological interpretation without testing that assumption on the
practitioner's own data.

### Topology assumptions change age uncertainty in both directions

The second and third research questions concerned how topology assumptions
affect age inference in this generator. Correct and partial topology
assumptions narrowed age intervals without a coverage penalty, and the
benefit was several times larger under a weak, tritium-only tracer regime
(MAE improvement 0.164 years, interval-width reduction 0.912 years) than
under an informative, two-tracer regime (MAE improvement 0.062 years,
interval-width reduction 0.252 years). This asymmetry is consistent with
topology contributing information that a weak tracer panel cannot supply on
its own, while a strong, two-isotope tracer panel already resolves most of
what topology could add, leaving comparatively little additional
uncertainty for a correct topology assumption to remove. Reversed topology,
however, was not a merely less-accurate alternative estimate: under the
informative regime it degraded both accuracy and coverage outright (MAE
+0.282 years, coverage -0.028), and under the weak regime it caused
importance-sampling degeneracy, with the effective-sample-size stability
rule failing in 8 of 12 cases, a diagnostic signature of model-data
incompatibility rather than of a competing but weaker inference. This
distinction matters for practice: a wrong flow-topology assumption does not
merely add noise to an age estimate in this generator, it can silently
destabilise the inference in a way that becomes visible only if
effective-sample-size diagnostics are checked alongside the point estimate,
and a practitioner who reports only a mean and an interval width could
easily miss that the underlying reweighting had already broken down.

### Entropy reduction without calibration is a detectable form of false confidence

The fourth research question asked whether a reduction in posterior
uncertainty from a misspecified evidence stream still indicates improved
inference. It does not, and this pattern is, in one sense, an expected
consequence of permuting an informative feature under a proper scoring
rule: destroying a feature's case-specific correspondence to the truth
while preserving its marginal distribution should be expected to degrade
discrimination on general statistical grounds. What this benchmark adds
beyond that expectation is a quantitative account of how large the
resulting entropy-versus-skill divergence is in this generator, and a
demonstration that a plausible, predeclared univariate heuristic (the
probability-span conflict diagnostic) failed to catch the divergence in
every one of the four tested conditions even though case-level
discrimination and calibration metrics caught it clearly and every one of
those degradations survived multiplicity correction. Permuted age reduced
entropy by 0.0207 while reducing PR-AUC by 0.075; permuted hydraulics
reduced entropy by 0.0482 while increasing log loss by 0.0745; and joint
misspecification reduced entropy the most of any condition (-0.0706) while
reducing PR-AUC by 0.139 and raising overconfident error by 0.0387. An
interpretation can become measurably more confident and measurably more
wrong at the same time, and the size of the entropy reduction gave no
indication of which direction the accompanying skill change would take.
Reporting the conflict diagnostic's failure transparently, rather than
replacing it with a better-tuned threshold after seeing the locked-test
outcome, was necessary to avoid overstating this guardrail's reliability.
This generalises a concern raised for tracer-based age inference
specifically, that fit quality alone can obscure systematic limitations
[@mccallum2015limitations], to the broader practice of multi-evidence
integration: entropy reduction is an unsafe stopping rule for judging
whether an added evidence stream has improved, rather than merely appeared
to improve, an interpretation, and should be reported alongside
discrimination and calibration checks
[@davisgoadrich2006prcurves; @brier1950verification], a discipline this
manuscript's own revised age and hydraulics results (above) were held to
once multiplicity correction was applied.

### Carbonate reactions remain non-identifiable regardless of panel richness

The third research question also concerned whether an enhanced
hydrochemical indicator panel resolves reaction-mechanism non-uniqueness
uniformly across processes. It does not in this generator. Denitrification,
sulfate reduction, and silicate weathering were recovered reliably under
both panels, and iron reduction improved materially with the enhanced
panel (modal accuracy rising from 0.000 to 0.500), but carbonate weathering
and precipitation were recovered in neither, across all 36 scored
carbonate edges under each panel. The coexistence of zero accuracy with low
family-support entropy for several carbonate rows is a direct illustration
of the equifinality concern raised for mechanistic model fitting generally
[@beven2001glue; @beven2006manifesto]: numerical stability and a
confident-looking result are not evidence of a correct mechanism, and
adding chemical indicators does not automatically break a structural
non-identifiability that a richer ion panel does not actually constrain.
The practical implication is process-specific rather than general, and
specific to the reaction dictionary and indicator panels tested here: an
enhanced indicator panel of the kind tested is a defensible investment for
redox and silicate process resolution in this system, but carbonate
reaction attribution appears to require a different kind of evidence, such
as isotopic or thermodynamic constraints not represented in either panel
tested here, rather than simply more of the same class of ion measurement.

### What Ghana field data can and cannot support

The fifth research question asked what the Northern Ghana data can
defensibly support. The available major hydrochemistry and stable isotopes
support component-level diagnostics, within-campaign chemistry
hold-forward, and explicit non-identifiability mapping. The absence of an
environmental age-tracer panel, screen intervals, a repeated head series,
and independently verified connectivity means the dataset cannot support
field validation of residence time, exact directed topology, or a unique
reaction mechanism, regardless of how the synthetic-benchmark results in
Experiments 1 through 3 are interpreted. This evidence boundary may be
common to other data-limited groundwater monitoring programmes that were
not designed with residence-time or connectivity validation in mind, but
that is an inference from this single dataset's structure rather than a
claim verified against other datasets, and the audit's explicit mapping
from available evidence to defensible claims is offered as a template for
stating such a boundary rather than as evidence that the boundary is
universal. Synthetic truth throughout the first four experiments is
model-conditioned by the independent MODFLOW 6/MODPATH 7 generator; it is
not, and should not be read as, evidence about a real aquifer, and none of
the quantitative integration results reported above should be transferred
numerically to a field setting without independent verification that the
generator's chemistry archetypes and topology are representative of that
setting.

### Limitations

Five limitations bound these conclusions. First, all quantitative
integration results derive from one generator geometry, one set of
chemistry archetypes, and a fixed twelve-case locked test set; the specific
magnitude of each contrast, though not necessarily its qualitative
direction, may differ in other hydrogeological settings or at a different
sample size, and every practical-implication statement above is scoped to
the tested conditions rather than to groundwater evidence integration in
general. Second, all synthetic truth is model-conditioned, so the
benchmark tests consistency with a known generative process rather than
field reality, and the strength of a conclusion about, for example,
carbonate non-identifiability is bounded by how representative the
generator's carbonate chemistry is of real carbonate systems. Third, the
predeclared univariate conflict diagnostic was retained and reported as
insensitive rather than replaced after seeing locked-test outcomes, which
is a conservative but incomplete safeguard: a better single-threshold
heuristic may exist and was not searched for post hoc, so the diagnostic's
failure should be read as evidence against relying on any single univariate
threshold rather than as evidence against this specific one. Fourth, the
Ghana audit is descriptive rather than a truth benchmark, and its
conclusions are limited to which interpretations the available evidence can
support, not to the correctness of any specific interpretation drawn from
it; extending the audit to a multi-year or multi-basin field dataset with
independent age or connectivity truth would be required before any of the
synthetic-benchmark findings could be treated as field-validated. Fifth,
the importance-sampling effective-sample-size stability rule (400 of
50,000 particles) was adopted by analogy with a convention used for a
related but distinct multi-chain Markov chain Monte Carlo diagnostic
rather than derived independently for single-particle-set importance
sampling, and the constants governing the adverse-control and
decision-rule design (the conflict-diagnostic threshold, the topology
potential's order scale, the reaction-noise level) were fixed by
predeclaration rather than validated by a sensitivity analysis across
alternative values; whether the qualitative pattern of results is stable
under alternative reasonable choices for these constants is left for
future work.

## Conclusion

This truth-blind benchmark tested, in one independent synthetic generator,
whether combining hydraulic, age, and chemical evidence reduces groundwater
interpretive non-uniqueness. Combining evidence streams reduced topological
non-uniqueness only conditionally: chemistry was robustly complementary,
while age and hydraulics each added an entropy reduction (and, for
hydraulics, a calibration gain) that survived multiplicity correction
without a topology-ranking change that did. Topology assumptions changed
groundwater-age uncertainty in both directions: correct and partial
topology conditionally reduced age error and interval width, most strongly
under a weak tracer regime, while reversed topology degraded accuracy
under an informative regime and destabilised inference outright under a
weak one. An enhanced hydrochemical indicator panel improved several redox
and silicate reaction processes but left carbonate weathering and
precipitation non-identifiable regardless of panel richness. Every
predeclared adverse control demonstrated, robustly to multiplicity
correction, that a reduction in posterior uncertainty can coexist with
worse discrimination and calibration, establishing false confidence as a
measurable, reportable failure mode rather than an assumption to guard
against informally. The Northern Ghana audit showed that this real,
data-limited aquifer's available evidence can support component chemical
and isotopic diagnostics without supporting field validation of residence
time, exact topology, or reaction mechanism. Taken together, and bounded by
the single generator geometry and locked test-case count on which they
rest, these results indicate that the benefit of integrating independent
groundwater evidence streams must be demonstrated per evidence pairing,
against paired discrimination and calibration checks corrected for the
number of contrasts tested, rather than assumed from the number of streams
combined or from a reduction in reported uncertainty alone.

## Open Research

**Author contributions.** [To be completed by the submitting authors using
a CRediT-style statement, for example: conceptualisation, methodology,
software, formal analysis, investigation, data curation, writing (original
draft), writing (review and editing), and visualisation, attributed to
each named author individually rather than as "all authors contributed
equally."]

**Competing interests.** [To be completed by the submitting authors. State
any financial or non-financial competing interests, or state that none are
declared.]

**Data availability.** The locked benchmark result tables (evidence panel
and case-block contrast summaries, topology-conditioned age sensitivity
tables, reaction-family support tables, the Northern Ghana data-scope
audit record, and the multiplicity-correction robustness check introduced
in this revision) and the per-case simulation provenance records
(executable identity, checksums, and version strings for every generated
case) are held in the authors' version-controlled project repository at
`https://github.com/dabdul-wahab1988/Hydrosheaf`. Prior to submission, these result
tables and provenance records should be deposited in a persistent,
publicly accessible repository (for example Zenodo or a comparable
AGU-recommended repository) and the resulting dataset DOI substituted for
this placeholder. The Northern Ghana workbook audited in Experiment 4
(`data/FieldData/NorthenGhana/NorthernGhana.xlsx`) is this project's own
field dataset and is held in the same repository; it should be deposited
alongside the result tables above and cited by the same dataset DOI. An
earlier revision instead audited a different, antecedent study's own
derived reprocessing of the same boreholes
(`Aquifers_Dataset_Mendeley.xlsx`, documented in `data/FieldData/
NorthenGhana/SI.pdf` as the supplementary information for "Graph-inverted
Ghanaian aquifers under aridification"); that workbook is not this
project's data, is not redistributed by this manuscript, and has been
removed from the analysis (DECISIONS.md).

**Code availability.** The independent synthetic-truth generator, the
evidence-fusion and topology-conditioning analysis scripts, the reaction
bootstrap procedure, and the multiplicity-correction script used for the
robustness check reported in this revision are held in the same
version-controlled repository referenced above (protocol-freeze commit
d336e87). This repository should be archived at a persistent, versioned
release (for example via a Zenodo-linked GitHub release) and the resulting
software DOI substituted for this placeholder before submission. Software
versions used to generate the reported results were MODFLOW 6 (build
6.7.0), MODPATH 7 (version 7.2.001), FloPy 3.10.0, PHREEQC via the
IPhreeqc interface, scikit-learn 1.9.0, numpy 2.4.6, scipy 1.17.1, and
pandas 2.3.3, run under Python 3.14.6, as detailed in Supplementary
Methods.
