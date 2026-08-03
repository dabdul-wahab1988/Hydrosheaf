# Conditional evidence integration and the incremental contribution of sheaf structure in controlled-synthetic groundwater benchmarks

**Abstract**

**Purpose.** The problem is that groundwater evidence integration is often
assumed to reduce interpretive non-uniqueness, although lower uncertainty does
not prove better inference, and graph-based connectivity prediction does not
isolate a sheaf layer's contribution. We tested both propositions.

**Design and methods.** Seven linked audits were conducted. The process-based
integration and public-pipeline audits used an independent MODFLOW 6/MODPATH 7
flow, tracer and nonlinear-chemistry generator; the competence-matched
representation benchmark and follow-up estimator diagnostic used a separate
scalar affine graph generator. Locked tests comprised 12, 6, 64 and 128 cases,
respectively. A descriptive Northern Ghana audit defined the field-data claim
boundary. The full representation-benchmark and estimator-diagnostic contrast
families were evaluated with
10,000-resample max-z simultaneous intervals.

**Primary representation result.** The competence-matched representation
benchmark passed identity-limit nesting
and the native affine sheaf outperformed the identity graph on PR-AUC (+0.0854,
simultaneous 95% CI [0.0539, 0.1169]) and the permuted-map control (+0.0909
[0.0571, 0.1246]). Planted-cycle conflict localisation also improved over the
edge-local graph (+0.0689 [0.0318, 0.1061]). However, overall PR-AUC versus
edge-local was +0.0097 [-0.0149, 0.0343], and neither Brier score nor log loss
improved. The representation benchmark failed the prespecified complete gate.

**Follow-up result.** In the estimator diagnostic, development selected a local
weight of 1.0, so
the estimator was local-first/global-fallback. Against edge-local, differences
were +0.0200 [-0.0055, 0.0454] for PR-AUC, -0.00151 [-0.00691, 0.00389] for
Brier score and +0.00333 [-0.0102, 0.0168] for log loss. The estimator diagnostic failed the
prespecified complete gate. Native maps nevertheless beat the permuted-map
control on all three outcomes, and conflict localisation survived
correction.

**Scope and significance.** In the process-based integration benchmark,
chemistry improved the topology-ranking
outcomes, correct topology reduced age MAE by 0.062--0.164 years, adverse
controls reduced uncertainty while worsening skill, and carbonate reactions
were not recovered under either tested indicator panel. These controlled-
synthetic benchmarks support a conditional representation claim: affine maps
encode non-identity relations and global compatibility supplies a conflict
diagnostic and missing-endpoint fallback. They do not establish general
predictive superiority, field validity, or performance for temporal,
three-dimensional, vadose-zone, vector-stalk or active-learning capabilities.

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
linked audits need rather than a complete graph posterior. Each
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
addition from one of these failure modes. In a targeted, non-systematic search
of publisher pages, Crossref-indexed records, Google Scholar-indexed records
and a preprint index, we did not identify a groundwater benchmark combining
all three of an external, code-independent truth generator, a locked
development/test split, and predeclared adverse controls capable of separating
a complementary evidence stream from one that is redundant or actively
harmful. This is a bounded account of the recorded search, not a systematic-
review claim of exhaustive absence. The combination, rather than any one
element in isolation, is the specific gap this study addresses. In the absence
of such controls, an
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

A second unresolved issue concerns representation rather than evidence
count. A cellular sheaf augments a graph by attaching local state spaces and
restriction maps to its incidences, so that a global section must satisfy
relations that can differ from edge to edge [@robinson2017sheaves;
@hansen2019spectral]. In the scalar identity limit, a sheaf Laplacian reduces
to an ordinary weighted graph Laplacian; non-identity scalar gains or affine
offsets are therefore the additional structure that requires direct testing
[@hansen2019learning]. This nesting makes a graph-only ablation necessary:
success of a full Hydrosheaf pipeline cannot by itself show that the sheaf
layer, rather than its features, optimiser, calibration, or graph topology,
caused the result. A recent preprint finding competitive identity-sheaf
baselines in unrelated graph-learning tasks reinforces the need for such a
competence-matched comparison, but does not supply groundwater evidence
[@caralt2026necessity].

Recent non-sheaf alternatives clarify the comparison class. Hydrochemical
metamodels have predicted groundwater-age distributions using symbolic and
gradient-boosted regression, although their targets were lumped-parameter-
model age products rather than known age truth [@tschritter2023age]. Statistical
hydrochemistry has also been used to infer potential aquifer connectivity and
source mixing [@zhao2022connectivity]. Graph neural networks have represented
monitoring-network dependencies for contaminant-transport emulation and source
attribution [@pang2024agnn], and multiple spatial dependencies for groundwater-
level forecasting [@wu2025groundwatergraph]. These methods solve prediction,
surrogate-modelling or association tasks; they do not provide the present
ablation of identity, affine and permuted restriction maps under a common
scoring pipeline. They nevertheless show why the sheaf comparison must be
against competent structured and edge-local alternatives, not an intentionally
weak graph baseline.

### Aim and objectives

The aim of this study was to determine, using two independent controlled-
synthetic generator systems and a truth-blind inference design, when combining hydraulic, age, and
hydrochemical evidence reduces interpretive non-uniqueness in groundwater
topology, age, and reaction inference; when it is redundant; when a
misspecified stream produces false confidence; and what the sheaf layer adds
beyond an ordinary weighted graph. Five objectives followed
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
claims. The fifth was to isolate restriction-map value by comparing the
native affine sheaf with an identity-restriction graph Laplacian, a stronger
edge-local weighted-graph comparator, and a restriction-map permutation
control under identical observations, candidate edges, solver budget,
calibration model, development split, and locked-test cases.
The same objective included a second, fresh-seed diagnostic to determine
whether any loss against the edge-local graph arose from replacement of direct
endpoint evidence, candidate self-influence in the global section, or
representation-specific calibration. This diagnostic compared original and
leave-one-edge-out global residuals, local/global hybrids, separately
cross-fitted calibration, a shared-calibrator analysis, and a frozen
permuted-map control.

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
than something demonstrated here. The results identify which evidence
combinations improved the prespecified synthetic-benchmark outcomes and which
additional streams risked degrading, rather than improving, an existing
interpretation. They do not establish a value-of-information threshold or a
measurement-cost recommendation. Taken together, the seven linked audits---
(1) multi-evidence integration, (2) topology-conditioned age, (3) reaction-
family recovery, (4) Northern Ghana evidence supportability, (5) public-
pipeline execution, (6) competence-matched sheaf-versus-graph representation,
and (7) local-first/global-fallback estimator diagnosis---supply a
candidate cross-layer integration guardrail that groundwater studies
combining hydraulic, age, and chemical evidence can test against their own
data before treating a lower reported uncertainty as evidence of a better
interpretation. They also replace the broad statement that "the sheaf helps"
with two falsifiable propositions: non-identity restriction maps should
outperform their identity limit where edge relations are heterogeneous, and
any claimed benefit should survive comparison with a strong graph method that
receives the same local evidence. The estimator follow-up adds a third: a
robust or hybrid correction must improve both discrimination and calibration
under fresh locked seeds before it can replace the bounded diagnostic claim.

**Table 1.** Design and claim map for the seven audits, distinguishing generator, sample, claim-bearing comparator, outcomes, lock, multiplicity family and permitted claim.

| Audit | Generator or data | Development / test cases | Claim-bearing comparator | Primary outcomes and gate | Lock and multiplicity | Permitted claim |
| --- | --- | --- | --- | --- | --- | --- |
| Multi-evidence integration | Process-based MODFLOW 6/MODPATH 7 plus independent tracer/chemistry | 6 / 12 | HAC versus nested panels and within-case permutations | PR-AUC, Brier, log loss, entropy | Commit d336e87; 24-test BH family | Conditional incremental value of each stream |
| Topology-conditioned age | Same process-based generator | 6 / 12 | Complete, partial or reversed topology versus none | MAE, interval width, coverage, ESS | Commit d336e87; case-block intervals | Model-conditioned topology effect |
| Reaction-family recovery | Same process-based generator; six planted archetypes | 6 / 12 | Enhanced versus core indicator panel | Modal accuracy and family support | Commit d336e87; descriptive/post-hoc interval | Discrimination among planted archetypes only |
| Ghana scope audit | Northern Ghana workbook | Not applicable / 160 wells; 140 complete seasonal pairs | Observed support versus claim requirements | Availability and truth-free hold-forward | Data-scope audit; no truth-bearing multiplicity family | Component diagnostics and non-identifiability map |
| Public-pipeline execution | Fresh process-based-generator seeds | 0 / 6 | Full sheaf versus hydraulic, local-age and age-permuted arms | Execution plus PR-AUC, Brier, log loss and selected F1 | Protocol lock; predeclared paired contrasts | System execution; no general incremental claim |
| Sheaf-versus-weighted-graph representation | Independent scalar affine graph generator | 32 / 64 | Affine sheaf versus edge-local, identity and permuted-map arms | Complete discrimination/calibration gate; conflict localisation | Protocol/source locks; all 120 contrasts in one max-z family | Identity nesting, map semantics and conditional conflict localisation |
| Local-first/global-fallback estimator diagnosis | Same scalar generator; fresh seeds | 64 / 128 | Local-first/global-fallback versus edge-local; original, LOO and permuted controls | Complete PR-AUC/Brier/log-loss gate | Two-stage locks; all 560 contrasts in one max-z family | Estimator diagnosis and missing-endpoint fallback; no superiority claim |

## Methods

### Independent synthetic-truth generator and blinding

The process-based evidence-integration, topology-conditioned-age and reaction-
family audits, and the later
public-pipeline audit, were benchmarked against an external synthetic-truth
generator that shared no code with the inference method under test.
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
locked-test candidate recall of 1.000 for the process-based integration
benchmark. The public-pipeline audit used
six additional fresh cases from the same generator family and had candidate
recall 0.9815. The full generator geometry,
chemistry archetypes, computational environment, and blinding contract are
reported in Supplementary Methods.

### Evidence streams and topology-ranking fusion

Three evidence streams were derived for every candidate directed edge: a
hydraulic log-odds term (H); an uncertainty-aware negative age-cost term (A);
and a negative constrained-chemistry log-objective term (C). Seven evidence
panels were formed from these streams and their combinations (H, A, C, HA,
HC, AC, HAC). For each panel, a logistic regression was fitted only on
development-case edges, using the panel's raw features standardised to
development-sample mean and standard deviation, with no class weighting. True
edges represented 54 of 414 development candidates (prevalence 0.1304) and 108
of 827 locked-test candidates (0.1306). Unweighted fitting was the prespecified
design choice because the estimand was the probability of edge truth under
this generated candidate distribution. Class weighting would target a
different effective prevalence; weighted estimators could still be assessed
and recalibrated on untouched data, but that alternative was not part of the
locked protocol. For a panel with standardised features
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

### Auxiliary shared-nuisance mechanism diagnostic

After the primary process-based and representation locks, a separate
controlled-synthetic diagnostic was run for the unresolved M3 tracer-
infeasibility mechanism. The fresh generator added environmental isotopes,
the expanded nuclear panel and paired recharge-temperature and 14C-initial-
activity nuisance levels. Development seeds 5401--5406 fitted the 15-panel,
three-target fusion models; locked seeds 5501--5512 were scored truth-blind
with 10,000 case-block bootstrap replicates. A local non-negative
convex-mixture transit-time test was then applied to the nuclear panel,
including CFC-11-containing pairs, the CFC-12 specificity control and the
3H--3H/3He pair. Stable-isotope evidence was masked from age by the declared
target-conditional design and was required to produce an exact zero age
increment. Full implementation details, provenance and the claim decision are
given in Supplementary Methods and Supplementary Table S14.

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
aggregated by process. The reaction solver was the inference method under
evaluation and was not part of the synthetic generator, preserving generator-
inference separation. Generation and scoring nevertheless used the same six-
family reaction vocabulary. The estimand was therefore discrimination among
those six planted archetypes under the specified stoichiometry, mineral
assemblage, noise model and two indicator panels; out-of-dictionary reactions
were not tested. Carbonate weathering and carbonate precipitation were
reported separately from the aggregate rather than folded into it.

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

### Generator construct validity and non-transfer boundary

Two controlled-synthetic generator systems answer different questions. The
process-based MODFLOW 6/MODPATH 7 generator represents two-dimensional saturated flow,
pathline-derived age, nonlinear tracer responses and six planted reaction
archetypes. It can test evidence addition, topology conditioning and recovery
of those planted reactions, but it cannot validate a field aquifer. The scalar
affine graph-case generator used for the representation benchmark and estimator
diagnostic instead plants identity, heterogeneous-
affine, incompatible-cycle and noisy/missing-observation relations directly.
It can isolate restriction-map information and estimator behaviour without
hydrogeological simulator confounding, but it does not reproduce the process-based
flow, age or chemistry process. A result from either generator is not treated
as replication in the other. The Ghana workbook provides a descriptive scope
audit, not connectivity, age or reaction truth. Accordingly, no result is
transferred across generator systems or to field validity without an
independent claim-bearing replication.

### Public-pipeline acceptance and competence-matched sheaf ablation

System-level execution and representation-level contribution were evaluated
separately. The public-pipeline acceptance run passed six independently
generated cases through candidate generation, hydraulic evidence, local-age
inference, directed-section solving, and calibrated edge scoring without
private benchmark shortcuts. Full-sheaf predictions were compared with
hydraulic-only, local-age, and within-case age-permutation conditions using
paired case-block bootstrap intervals. Execution success was necessary but
not sufficient: an incremental full-sheaf claim required a predeclared gain
over the simpler comparators.

The subsequent competence-matched representation benchmark isolated the
mathematical representation. The
design treated sheaf restriction maps and
global sections according to their information-integration and spectral
formulations [@robinson2017sheaves; @hansen2019spectral], including the known
identity weighted-graph special case [@hansen2019learning]. An
independent generator importing no Hydrosheaf code produced scalar graph
cases in four equal scenarios: an identity limit, heterogeneous affine
relations, incompatible cycles, and noisy or missing observations. Thirty-two
development cases (seeds 7401-7432) fixed one common calibration model and all
thresholds; 64 locked cases (seeds 7501-7564; 16 per scenario) were evaluated
once after the protocol and source hashes were frozen. Every model received
the same nodes, candidate edges, edge weights, local observations, features,
regularisation and three-solve reweighting budget. The identity graph and
affine sheaf used the same directed least-squares implementation, thereby
isolating the restriction maps rather than the optimiser. For node states
$x_v$ and directed edge $e:u\rightarrow v$, the affine section estimate
minimised

<!-- EQ:EQ-7 -->
$$
\sum_e w_e\left(\alpha_e x_u+b_e-x_v\right)^2
+\sum_{v\in O}w_v^{\mathrm{obs}}(x_v-y_v)^2
+\varepsilon\sum_v x_v^2.
$$

The identity graph fixed $\alpha_e=1$ and $b_e=0$. The native affine sheaf
used the generated relation maps, while a negative control permuted those maps
within each case. A stronger edge-local weighted graph received the same
common prior and local affine residual but did not solve for a globally
compatible section. The generated edge prevalence was one-third in development
and locked test. The primary contrasts were affine sheaf minus each graph
comparator for PR-AUC, Brier score, log loss and selected F1; secondary
endpoints were expected calibration error, false-confidence rate,
abstention accuracy and planted-conflict localisation PR-AUC. Ten-thousand-
replicate case-block bootstrap intervals were reported. In the post-review
audit, all 120 published representation-benchmark scenario-comparator-metric
contrasts formed one
family and received two-sided max-z simultaneous 95% intervals from the same
shared case-block resamples. The identity-limit
equivalence gate required prediction differences below $10^{-12}$, and a
general incremental-sheaf claim required both discrimination and calibration
improvement over the strong edge-local graph, not merely over the deliberately
restricted identity baseline.

The prospectively locked estimator diagnostic then tested whether the
representation-benchmark loss
arose from discarded endpoint evidence, candidate self-influence, or
calibration. The unchanged independent generator produced 64 development
cases (seeds 8401--8464) and 128 fresh locked cases (8501--8628), balanced
across the same four scenarios. The protocol was hashed before implementation;
after development-only selection, the runner, fitted settings and development
manifest were frozen before the single permitted test execution. Every arm
received the prior logit, one residual scalar and a missing-endpoint indicator.
Original and leave-one-edge-out (LOO) global residuals were evaluated alone
and in hybrids with the local residual. For LOO scoring, each candidate was
evaluated against a section fitted without that candidate before robust weight
updating, preventing an edge from making itself appear compatible.

Hybrid weights and logistic regularisation were selected jointly from fixed
grids by eight-fold seed-group cross-validation minimising mean case log loss.
The generated edge prevalence was one-third in both development and locked
test. Arm-specific cross-fitted calibration was primary; a common shared calibrator
was secondary. The overall claim required favourable 95% intervals for
PR-AUC, Brier score and log loss against edge-local, identity-limit
non-degradation, and favourable native-versus-permuted PR-AUC and Brier
intervals. Failure of any component prohibited general superiority. Complete
LOO updates, grids, thresholds and secondary metrics are reported in
Supplementary Methods. The post-review multiplicity audit placed all 560
published estimator-diagnostic scenario-comparator-metric contrasts in one
family and used
10,000 shared case-block resamples to form two-sided max-z simultaneous 95%
intervals. Scenario statements were retained as supported only when those
simultaneous intervals excluded zero in the favourable direction.

### Statistical reporting and decision rules

Topology-ranking discrimination was reported as precision-recall area under
the curve (PR-AUC) alongside receiver-operating-characteristic area under
the curve (ROC-AUC). PR-AUC was prioritised as the primary discrimination
metric because the locked-test candidate edge set is class-imbalanced and
because precision-recall curves can give a more informative view under class
imbalance [@davisgoadrich2006prcurves]; more recent work
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
predeclared contrasts and four metrics reported in Supplementary Table S2, an exact
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

### Post-review precision planning and practical-magnitude audit

No minimum practically important differences or prospective power analysis
were prespecified before the locked runs. To avoid retroactively labelling a
post-review choice as confirmatory, the revision used transparent planning
margins only for interpretation and future replication: PR-AUC 0.02, Brier
score 0.01, log loss 0.02, age MAE 0.25 years, interval width 0.50 years,
coverage 0.05 and modal reaction-family accuracy 0.10. These are not field-
validated minimum important differences. Twenty-thousand empirical planning
simulations resampled development-case effects for the 12-, 64- and 128-case
designs. The process-based integration and representation-benchmark development
models were evaluated on their fitting cases, so those planning probabilities
may be optimistic; the estimator diagnostic used its
eightfold out-of-fold development predictions. Because development topology-
age and reaction summaries were not archived, their locked-test vectors were
used only to plan future replication, never to justify the completed tests.
Realised effects were also expressed relative to the comparator mean and
checked against the post-review margins. Full methods and results are reported
in Supplementary Tables S10--S12.

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
validated universal constants, and are retained as explicit process-based-generator
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

### Benchmark scale and audit design

The locked process-based integration benchmark comprised six development cases
(seeds 5201-5206),
used only to fit evidence-fusion coefficients and any classification
threshold, and twelve independent locked test cases (seeds 5301-5312),
generated by the official MODFLOW 6/MODPATH 7 simulator combined with an
independently coded nonlinear chemistry/tracer model that shared no code
with the inference method it evaluated. The benchmark architecture is
summarised in Figure 1.

![](figures/publication/figure1_benchmark_and_claim_design.png)

**Figure 1.** Benchmark architecture and claim boundary. The independent synthetic-truth branch separates official MODFLOW 6/MODPATH 7 flow and pathline generation, nonlinear chemistry and tracer generation, truth-blind inference and locked adverse-control scoring from the Northern Ghana branch, which distinguishes supportable component diagnostics from field interpretations that remain non-identifiable under the available evidence.

Topology-conditioned age inference used 50,000 importance-sampling
particles per case and tracer regime; reaction-family inference used 64
chemistry-perturbation bootstrap replicates per case; and every paired
contrast used a 10,000-replicate case-block bootstrap resampling whole
independent cases. The seven-audit design, generator systems, test sizes and
claim gates are mapped in Table 1.

Candidate edge generation retained every true test edge, giving a
locked-test candidate recall of 1.000, so none of the incremental-value
contrasts reported below could be confounded by a missing true edge in the
candidate set.

The post-review simulation audit showed why near-zero contrasts require
caution. Under development-derived effects, the 12-case chemistry-
addition design excluded zero and the post-review planning margin in 100% of
20,000 simulations for PR-AUC, Brier score and log loss, although development
resubstitution may be optimistic. For the competence-matched affine-sheaf-versus-edge-local
contrast, the corresponding probabilities of a favourable interval excluding
zero were 0.193, 0.0138 and 0.0001, and probabilities of clearing the planning
margin were 0.00015, 0 and 0. For the local-first/global-fallback estimator versus edge-
local, out-of-fold development estimates gave probabilities 0.548, 0.123 and
0.0022 of excluding zero, and 0.0006, 0 and 0 of clearing the margin. Thus, the
later designs were adequate for capability diagnosis but were not prospectively
powered for small general-superiority effects. Topology-age and reaction
planning based on locked vectors is reported only for future replication
(Supplementary Table S12).

### Auxiliary shared-nuisance mechanism diagnostic

After the primary integration and representation results had been frozen, an
auxiliary controlled-synthetic diagnostic was run to investigate the unresolved
M3 pattern of pairwise tracer infeasibility. It used fresh development seeds
5401--5406 and locked test seeds 5501--5512, the independent v2 MODFLOW 6/
MODPATH 7 generator, three predeclared nuisance levels, and a local
convex-mixture transit-time feasibility test over the nuclear tracer panel.
This diagnostic was not a field-validation experiment and did not use the
USGS observations as truth.

The binding conservative-isotope age control passed exactly: adding the E
stream to N for T2 changed MAE by 0.0000 years (95% CI [0.0000, 0.0000]).
Severe shared nuisance increased full-panel synthetic infeasibility by 0.2882
(95% CI [0.2118, 0.3646]) relative to none. However, the redox-stratified
CFC-11-containing pair contrast (+0.7188 [0.6667, 0.7500]) did not pass the
predeclared CFC-12 specificity control (+0.7396 [0.7240, 0.7500]). The
selective CFC-11 mechanism was therefore not supported under this generator.
The complete diagnostic tables and claim decision are archived in
`results/RUN-M7-6-M3-MECHANISM-20260731-01/`; these results show that the
tested nuisance can increase synthetic infeasibility but do not identify the
cause of the USGS pattern.

### Evidence-panel performance on native locked-test data

On native locked-test data, chemistry was the strongest individual
topology-ranking evidence stream (PR-AUC 0.459, ROC-AUC 0.892, Brier score
0.094, log loss 0.295, mean edge entropy 0.455), hydraulics was
substantially weaker (PR-AUC 0.176, ROC-AUC 0.650, Brier score 0.108), and
age alone was below a useful ranking level (PR-AUC 0.111, ROC-AUC 0.452),
only marginally distinguishable from the base rate implied by chance
ranking (Supplementary Table S1). Among the three evidence pairs,
hydraulics-plus-chemistry (HC) was the strongest (PR-AUC 0.485, ROC-AUC
0.909, Brier score 0.088, log loss 0.268, mean edge entropy 0.422),
exceeding both age-plus-chemistry (AC; PR-AUC 0.471, Brier score 0.089) and
hydraulics-plus-age (HA; PR-AUC 0.111, effectively unchanged from age alone
and far below either single-stream chemistry result). The fully integrated
seven-feature panel (HAC) achieved PR-AUC 0.480, ROC-AUC 0.908, Brier score
0.088, log loss 0.268, mean edge entropy 0.421, and expected calibration
error 0.060, a profile marginally below HC on PR-AUC rather than above it,
despite combining strictly more evidence, as shown in Figure 2.

![](figures/publication/figure2_evidence_integration.png)

**Figure 2.** Evidence integration is conditional. Panel a compares native hydraulic (H), age (A), chemistry (C), pairwise and fully integrated topology-ranking performance. Panel b shows case-block incremental PR-AUC contrasts with 95% bootstrap intervals. Panels c and d show that permuted evidence can reduce entropy while degrading discrimination and calibration.

The full native and adverse condition grid for all seven panels, extending
this comparison beyond the headline HAC/HC/HA/AC contrasts, is reported in
Supplementary Table S1.

### Case-block evidence contrasts and adverse controls

Paired case-block contrasts isolated the incremental value of each stream
within the fully integrated panel across the twelve locked test cases
(Supplementary Tables S2 and S6). Adding age to hydraulics-plus-chemistry produced a
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
without a coverage penalty in both tracer regimes (Figure 3 and
Supplementary Table S3). Under the informative (³H + ³⁹Ar) regime, complete-true
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

Relative to the no-topology baselines, the informative-regime MAE and interval-
width changes were 2.24% of 2.764 years and 2.02% of 12.503 years; neither
cleared the post-review planning margins of 0.25 and 0.50 years. In the weak
regime, the MAE change was 3.46% of 4.750 years and did not clear 0.25 years,
whereas the 0.912-year interval-width reduction was 3.55% of 25.701 years and
cleared the post-review 0.50-year planning margin. These margins were not
prespecified or field-validated and are used only to distinguish statistical
precision from practical scale (Supplementary Table S12).

![](figures/publication/figure3_topology_conditions_age.png)

**Figure 3.** Correct topology improves model-conditioned age inference. Panels a and b show MAE and interval-width changes, panel c reports effective-sample-size fractions and panel d shows coverage changes. Reversed tritium-only accuracy contrasts are not interpreted because 8 of 12 cases failed the stability rule.

Reversed topology behaved asymmetrically to this benefit rather than as a
merely weaker version of it. Under the informative regime, reversed
topology increased MAE by 0.282 years relative to complete-true topology
(95% CI [0.016, 0.545]), widened 95% intervals by 0.296 years, and reduced
coverage by 0.028 (95% CI [-0.049, -0.007]). Under the weak regime, the
reversed-topology importance-sampling effective sample size fell below the
predeclared stability rule (ESS < 400 of 50,000 particles) in 8 of the 12
locked test cases; the nominal reversed-minus-correct MAE contrast under
this regime carried a 95% interval of [-0.889, 0.888], spanning almost two
years in each direction, and is reported for completeness in Supplementary
Table S3 but is
not interpreted as a stable estimate given this diagnostic failure.
Per-case topology-to-age sensitivity and effective-sample-size values
underlying this comparison are reported in full in Supplementary Table S3.

### Recovery among six planted reaction archetypes

Discrimination among the six planted reaction archetypes was strongly process-
dependent under both the core and enhanced hydrochemical panels
(Figure 4 and Supplementary Table S4).
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

![](figures/publication/figure4_reaction_nonuniqueness.png)

**Figure 4.** Recovery among the six planted reaction archetypes is process-dependent. The core and enhanced panels are compared for modal accuracy, probability assigned to the planted family, support entropy and effective support. Carbonate reactions were not recovered under either tested panel.

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
field, and masked rather than exact site coordinates (Figure 5 and the
detailed scope matrix in the Supplement).

![](figures/publication/figure5_ghana_supportability_boundary_m7_only.png)

**Figure 5.** Northern Ghana evidence and claim boundary using evidence from this study only. Panels a and b map observed evidence to defensible claims; panel c reports 160 wells, 140 complete seasonal pairs and 320 seasonal observations; panel d reports the truth-free seasonal hold-forward comparison. No external field-transfer experiment is included.

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

Table 2 synthesises the classification of every tested contrast
across the four original experiments under the predeclared decision rule,
refined by the multiplicity-corrected re-analysis reported above and in
Supplementary Table S6. Chemistry was classified as complementary
topology-ranking evidence on every metric; hydraulics was complementary on
calibration and entropy but not, after correction, on PR-AUC alone; and age
was likewise supported only by a robust entropy reduction rather than by a
statistically robust ranking change. Correct and partial topology assumptions
were classified as complementary age evidence, conditional on tracer regime
and strongest under a weak tracer panel, while reversed topology was
classified as negative transfer under the informative regime and as
numerically unstable, rather than interpretably worse, under the weak regime.
Every adverse permutation condition was classified as negative transfer and
survived multiplicity correction on every reported metric. The enhanced
chemistry panel was a partial, process-dependent improvement, while carbonate
weathering and precipitation were not recovered under either tested indicator
panel. The Northern Ghana
audit supported component diagnostics and non-identifiability mapping only.

**Table 2.** Primary process-based integration decision table. Detailed metrics and complete machine-readable contrasts are reported in Supplementary Tables S1--S6.

| Audit | Supported result | Unsupported or adverse result | Claim boundary |
| --- | --- | --- | --- |
| Evidence integration | Chemistry improved PR-AUC and calibration; hydraulics improved calibration and entropy | Age did not retain a ranking change after correction; every permutation control degraded skill despite lower entropy | Incremental value is stream- and outcome-specific |
| Topology-conditioned age | Complete and partial true topology reduced MAE and interval width | Reversed topology degraded the informative regime and was unstable in 8 of 12 weak-regime cases | Model-conditioned topology effect only |
| Reaction-family recovery | Recovery was strong for selected planted redox and silicate archetypes | Carbonate reactions were not recovered under either tested panel | Six planted archetypes, one noise model and two panels |
| Northern Ghana scope | Chemistry/isotopes support component diagnostics and truth-free hold-forward | Exact topology, residence time and unique mechanisms lack independent truth | Scope audit, not field validation |

### Public-pipeline acceptance and incremental sheaf contribution

The strict public-pipeline run completed every predeclared stage on all six
cases, with mean candidate recall 0.9815, and therefore passed its execution
gate. It did not pass the stronger incremental full-sheaf claim gate
(Supplementary Table S8). Full-sheaf PR-AUC was 0.3075, compared with 0.3272
for hydraulic-only evidence (difference -0.0197, 95% CI -0.0355 to -0.0039),
0.2488 for local-age evidence (difference 0.0586, 95% CI 0.0386-0.0777), and
0.3211 after within-case age permutation (difference -0.0136, 95% CI
-0.0622-0.0347). Thus, the public workflow executed reproducibly and improved
on one weak component, but did not establish incremental value over every
simple or adverse comparator.

Selection was not based on a scalar probability threshold: the public pipeline
retained all 198 generated candidates in every arm. Consequently, each arm had
the same conditional confusion counts (53 true positives, 145 false positives,
0 false negatives) and the same conditional selected F1 (0.4223), despite
different probabilities and log losses. Minimum retained probabilities were
0.1921 for full sheaf, 0.0257 for hydraulic only, 0.3170 for local age and
0.0222 for age-permuted. Candidate generation missed one of 54 truth edges,
making end-to-end counts 53 true positives, 145 false positives and 1 false
negative, end-to-end F1 0.4206 and candidate recall 0.9815. The system audit
therefore tested scoring only among recovered candidates; it did not establish
complete candidate generation (Supplementary Table S13).

The separate 64-case representation benchmark passed all execution and identity-limit
equivalence checks (Figure 6, Table 3). In all 16 identity-
limit cases, affine-sheaf and identity-graph raw residuals and calibrated
predictions were exactly equal. Across all scenarios, the affine sheaf
improved PR-AUC over the identity graph Laplacian by 0.0854 (simultaneous 95%
CI [0.0539, 0.1169]), reduced Brier score by 0.0193 [-0.0261, -0.0126], and
reduced log loss by 0.0472 [-0.0638, -0.0306]. Native affine maps also improved
PR-AUC over the permuted-map sheaf by 0.0909 [0.0571, 0.1246], with favourable
Brier and log-loss intervals. These intervals control the family-wise error
rate across all 120 published representation-benchmark contrasts
(Supplementary Table S10).

![](figures/publication/figure6_sheaf_vs_weighted_graph.png)

**Figure 6.** Incremental contribution of affine sheaf structure. Overall and selected scenario contrasts compare edge-local, identity, native affine and permuted-map arms. Error bars are simultaneous 95% intervals controlling all 120 published representation-benchmark contrasts as one family. Intended printed width is 7.08 in; minimum label size is 8 pt.

**Table 3.** Locked-test performance of the competence-matched edge-local graph, identity graph Laplacian, affine sheaf and permuted-map control across 64 cases. Family-wise contrasts are in Supplementary Table S10.

| Model              |   Cases |   PR-AUC |   Brier |   Log loss |   Selected F1 |
|:-------------------|--------:|---------:|--------:|-----------:|--------------:|
| Affine sheaf       |      64 |   0.6762 |  0.1777 |     0.5306 |        0.6235 |
| Permuted-map sheaf |      64 |   0.5853 |  0.1992 |     0.5852 |        0.5506 |
| Identity Laplacian |      64 |   0.5908 |  0.197  |     0.5778 |        0.5647 |
| Edge-local graph   |      64 |   0.6665 |  0.1772 |     0.5189 |        0.6219 |

The sheaf did not, however, outperform the stronger edge-local weighted graph
overall. Its PR-AUC difference was +0.0097 (simultaneous 95% CI [-0.0149,
0.0343]), its Brier difference was +0.00053 [-0.00584, 0.00691], and its log-
loss difference was +0.0117 [-0.00665, 0.0301]. The mean PR-AUC gain was 1.45%
of the edge-local mean and did not clear the post-review 0.02 planning margin;
neither calibration metric cleared its margin. In planted incompatible-cycle
cases, PR-AUC improved by +0.0483 [0.00954, 0.0871] and conflict-localisation
PR-AUC by +0.0689 [0.0318, 0.1061]. Conversely, the heterogeneous-affine
scenario remained worse than edge-local on PR-AUC (-0.0443 [-0.0788,
-0.00991]), Brier score (+0.0216 [0.0119, 0.0313]) and log loss (+0.0755
[0.0474, 0.1037]). No noisy/missing representation-benchmark contrast against
edge-local survived the full-family correction. The representation benchmark
therefore failed the prespecified complete
gate. The supported contribution is exact identity nesting, information in
native versus permuted maps, and conditional planted-cycle localisation---not
overall predictive superiority.

The separately locked estimator diagnostic passed its execution/provenance gate on
128 fresh cases (Figure 7, Table 4). The exact locked runner
was restored and archived under SHA-256
`a0ef13bde5391af62698927211cb4e701123affebb108d331795ce8596e2e191`, matching
the confirmatory lock; stored locked-test files were hash-verified and were not
regenerated. Development assigned both selected arms a local weight of 1.0,
retaining the local residual wherever both endpoints were observed and using
the global residual only when an endpoint was missing. The selected method is
therefore described as local-first/global-fallback, not as a general local/
global blend. Under primary separately cross-fitted calibration, its overall
difference from edge-local was +0.0200 for PR-AUC (simultaneous 95% CI
[-0.00550, 0.0454]), -0.00151 for Brier score [-0.00691, 0.00389], and
+0.00333 for log loss [-0.0102, 0.0168]. The PR-AUC mean was 3.09% of the
edge-local mean but did not exceed the post-review 0.02 planning margin. The
estimator diagnostic failed the prespecified complete gate.

![](figures/publication/figure7_robust_hybrid_sheaf.png)

**Figure 7.** Local-first/global-fallback estimator diagnostic. The selected nominal hybrid had local weight 1.0 and is therefore local-first/global-fallback. Error bars are simultaneous 95% intervals controlling all 560 published estimator-diagnostic contrasts as one family. Intended printed width is 7.08 in; minimum label size is 8 pt.

**Table 4.** Locked-test performance of seven estimator-diagnostic arms across 128 fresh cases. The selected local-first arms had local weight 1.0; family-wise contrasts are in Supplementary Table S11.

| Model                |   Cases |   PR-AUC |   Brier |   Log loss |   Selected F1 |
|:---------------------|--------:|---------:|--------:|-----------:|--------------:|
| Edge-local           |     128 |   0.6459 |  0.1821 |     0.5306 |        0.6106 |
| Identity             |     128 |   0.5801 |  0.1966 |     0.5762 |        0.5648 |
| Original global      |     128 |   0.6638 |  0.1796 |     0.5312 |        0.6252 |
| Original local-first |     128 |   0.6666 |  0.179  |     0.5297 |        0.6324 |
| Robust global        |     128 |   0.659  |  0.1825 |     0.5385 |        0.6209 |
| Robust local-first   |     128 |   0.6658 |  0.1806 |     0.534  |        0.6263 |
| Permuted             |     128 |   0.6218 |  0.1929 |     0.5651 |        0.5832 |

Mechanism contrasts separated this outcome. The LOO robust global estimator
worsened Brier score by +0.00290 (simultaneous 95% CI [0.000951, 0.00485]) and
log loss by +0.00734 [0.00249, 0.0122] relative to the original global
estimator. Other small original-versus-local-first and robust-versus-original
contrasts did not survive the 560-contrast family correction. Against edge-
local, the local-first/global-fallback estimator retained a conflict-
localisation gain of +0.0770 [0.0389, 0.1151], but its scenario-specific
ranking gains did not survive. Native maps did beat the frozen permuted-map
control on PR-AUC (+0.0441 [0.0192, 0.0690]), Brier score (-0.01230 [-0.0175,
-0.00712]) and log loss (-0.03114 [-0.0440, -0.0183]). Complete unadjusted and
simultaneous intervals are reported in Supplementary Tables S9 and S11. The
result retains a conditional map-semantic, conflict-diagnostic and missing-
endpoint-fallback claim, while rejecting general superiority and the LOO
robust estimator as a replacement for edge-local scoring.

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

### What the sheaf layer contributes beyond an ordinary weighted graph

The direct answer is narrower and more useful than a claim of universal
superiority. An ordinary weighted graph represents which nodes are connected
and how strongly, while the tested affine sheaf additionally represents how a
state at one endpoint is translated to the other through edge-specific gains
and offsets, and evaluates whether those local translations admit a coherent
global section [@robinson2017sheaves; @hansen2019spectral]. Exact equality in
the identity-limit cases confirms the expected nesting: when every restriction
map is identity, the sheaf introduced no hidden algorithmic advantage over the
graph Laplacian [@hansen2019learning]. The positive native-versus-identity and
native-versus-permuted contrasts therefore identify non-identity relation
structure, rather than simply graph connectivity or solver choice, as the
source of the gain.

That contribution had two family-wise supported empirical forms. First, native
affine maps improved overall prediction and calibration over both the identity-
restriction Laplacian and the permuted-map control. Second, the global-section
residual localised planted cycle incompatibilities more accurately than the
edge-local comparator. This provides a defensible role for the sheaf layer as a
relation and consistency model: it can express edge relations that an identity
graph cannot, and identify where locally plausible relations fail to compose
globally.

The strong edge-local comparator prevents this result from being inflated.
The competence-matched representation benchmark did not improve overall PR-AUC
or calibration against edge-local after the complete 120-contrast family
correction. The fresh-seed estimator diagnostic
clarified rather than erased that result. Development selection placed both
candidate combinations at local weight 1.0, so global compatibility was used
only when an endpoint residual was unavailable. The selected method was
therefore local-first/global-fallback, not evidence that a general local/global
blend was superior. Its mean overall PR-AUC gain over edge-local was small and
did not survive the complete 560-contrast correction.

Stronger false-edge downweighting did not provide the anticipated remedy.
Leave-one-edge-out scoring prevented candidate self-influence but worsened
Brier score and log loss for the robust global estimator after family-wise
correction. Separately cross-fitted calibration also did not rescue the three-
outcome gate. The persistent heterogeneous-affine penalty therefore cannot be
attributed solely to the original shared calibration or to false edges pulling
their own fitted section. In this dense candidate generator, aggressive LOO
downweighting removed useful relational support as well as contamination.

The answer to the research question is consequently conditional but specific.
Compared with an identity graph, the sheaf contributes edge-specific gains and
offsets plus a test of whether those relations compose globally. Compared with
a strong edge-local graph, the supported contribution is a network-level
conflict signal plus the capability for global fallback when endpoint evidence
is missing; scenario-specific ranking gains did not survive full-family
correction. Native maps remained superior to permuted maps, but the method was
not a uniformly better probability estimator and was worse in the
representation benchmark's
heterogeneous-affine stratum. This distinction is consistent with the need for
strong identity and graph baselines in sheaf learning [@caralt2026necessity]
and with groundwater graph methods that learn spatial or transport dependencies
for different predictive tasks [@pang2024agnn; @wu2025groundwatergraph]. The
recent sheaf source is a preprint, those graph studies do not test the present
estimand, and the present result remains limited to controlled scalar cases.

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

### Carbonate reactions were not recovered under either tested indicator panel

The third research question also concerned whether an enhanced
hydrochemical indicator panel resolves reaction-mechanism non-uniqueness
uniformly across processes. It does not in this generator. Denitrification,
sulfate reduction, and silicate weathering were recovered reliably under
both panels, and iron reduction improved materially with the enhanced
panel (modal accuracy rising from 0.000 to 0.500), but carbonate weathering
and precipitation were recovered in neither, across all 36 scored
carbonate edges under each panel. This result is restricted to discrimination
among the six planted archetypes, the tested stoichiometry and mineral
assemblage, the 3% noise model, and the core and enhanced panels. The
coexistence of zero accuracy with low
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

### The auxiliary nuisance experiment narrowed, but did not solve, the M3 mechanism question

The auxiliary controlled-synthetic diagnostic supplied the ground truth that
the national M3 benchmark lacks, but its result was deliberately narrower than
an explanation of the USGS data. Severe shared nuisance increased the rate of
infeasible full nuclear panels, showing that a declared correlated error can
create additional age-constraint failure under the tested generator. The
redox-stratified CFC-11 contrast nevertheless failed the CFC-12 specificity
control: the CFC-12-containing pair contrast was at least as large. The
selective CFC-11 degradation mechanism was therefore not supported as a
controlled-synthetic explanation under this implementation.

This negative result is informative but not exhaustive. The diagnostic used an
analytic nuisance formulation and synthetic redox histories, not the upstream
USGS DGMETA correction process or co-located field redox measurements. It
therefore narrows the tested mechanism without resolving the M3 cause. The
diagnostic remains an auxiliary controlled-synthetic audit and provides no
field age, topology or correction validation.

### Limitations

Eight limitations bound these conclusions. First, the process-based quantitative
integration results derive from one MODFLOW/MODPATH generator geometry, one
set of chemistry archetypes and 12 locked cases, while the representation-
benchmark and estimator-diagnostic results derive from a separate scalar affine
generator with 64 and 128 locked cases. The second system is not a cross-
generator replication
of the first, and every practical statement is scoped to its tested generator.
Second, all synthetic truth is model-conditioned, so the
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
future work. Sixth, the representation benchmark and estimator diagnostic used
static scalar stalks
and controlled affine relations. It did not test temporal, three-dimensional,
vadose-zone, or field-validated vector-stalk behaviour, and the edge-local
comparator was deliberately strong but not an exhaustive survey of graph-
learning methods. The hybrid and LOO grids were finite, the LOO influence rule
was one specific robust estimator, and active learning was not evaluated in
these ablations. The representation claim must therefore remain conditional
on the tested relation family and cannot be extrapolated to general
superiority of sheaf models in groundwater science.

Seventh, prospective minimum important differences and sample-size
justifications were not specified before the locked tests. The revision's
planning margins, relative effects and empirical power simulations are
explicitly post-review diagnostics, not field decision thresholds and not a
retroactive basis for claiming adequate confirmatory power. The next claim-
bearing step is replication with a different generator family or a field
dataset having independently measured connectivity. No temporal, three-
dimensional, vadose-zone, vector-stalk or active-learning performance is
inferred from these results.

Eighth, the auxiliary M7.6 mechanism diagnostic used synthetic redox histories
and analytic shared-nuisance perturbations rather than the USGS DGMETA
correction process, co-located field redox observations or laboratory-specific
measurement covariance. Its negative CFC-11 specificity result therefore
narrows the tested synthetic explanation without identifying a field cause;
the accompanying increase in synthetic infeasibility is a mechanism-screening
diagnostic, not field validation.

## Conclusion

Across seven linked audits, evidence integration reduced groundwater
interpretive non-uniqueness only conditionally. In the process-based
MODFLOW/MODPATH
benchmark, chemistry improved the prespecified topology-ranking outcomes;
correct topology improved age inference by small relative amounts; reversed
topology degraded or destabilised it; carbonate reactions were not recovered
under either tested indicator panel; and adverse controls showed that lower
uncertainty can coexist with worse discrimination or calibration. The Northern
Ghana data supported component diagnostics and explicit non-identifiability
mapping, not field validation of exact topology, age or reaction mechanisms.

The auxiliary shared-nuisance diagnostic showed that the declared synthetic
nuisance increased full nuclear-panel infeasibility, but its redox-stratified
CFC-11 result failed the CFC-12 specificity control. It therefore did not
support a selective CFC-11 mechanism and did not resolve the corresponding M3
or USGS cause; this remains controlled-synthetic diagnostic evidence only.

The representation audits answer what the sheaf layer contributes beyond an
ordinary weighted graph. The competence-matched representation benchmark
demonstrated exact identity-limit nesting,
information in native affine maps relative to identity and permuted-map
controls, and family-wise supported localisation of planted cycle conflicts.
It did not outperform the strong edge-local graph overall and failed the
prespecified complete gate. In the estimator diagnostic, development selected
local weight 1.0, so
the selected method was local-first/global-fallback. Its overall ranking and
calibration contrasts against edge-local did not survive the complete 560-
contrast correction, and it also failed the prespecified complete gate. Native
maps still beat the permuted control and global conflict localisation survived;
stronger LOO robustification worsened global-estimator calibration.

The experiments support use of the sheaf layer as a prespecified model of non-
identity relations and as a global-compatibility diagnostic under the tested
scalar scenarios, with global fallback when endpoint evidence is missing.
They do not validate HydroSheaf as a whole or establish general predictive
superiority over weighted graphs. The next claim-bearing step is independent
replication under a different generator family or a field dataset with
independently measured connectivity; temporal, three-dimensional, vadose-zone,
vector-stalk and active-learning capabilities remain outside the evidence
evaluated here.

## Declarations and open research

**Author contributions.** Verified contributor-role assignments were not
present in the technical source package reviewed for this revision. A final
CRediT statement assigning conceptualisation, methodology, software, formal
analysis, investigation, data curation, writing and visualisation to each named
author is required from the authors before submission. No roles have been
inferred from commit history or licence metadata.

**Funding.** A verified funding statement was not present in the technical
source package. The authors must provide the funder names and award identifiers,
or confirm that the work received no specific funding, before submission.

**Competing interests.** A verified competing-interest declaration was not
present in the technical source package. Each author must confirm any financial
or non-financial interests, or confirm that none exist, before submission.

**Data availability.** The locked process-based integration tables, topology-
conditioned-age and reaction-family outputs, Northern Ghana data-scope audit,
public-pipeline records, immutable representation-benchmark and estimator-
diagnostic outputs, and the post-review multiplicity and precision
audit are held in the project repository at
`https://github.com/dabdul-wahab1988/Hydrosheaf`. The process-based integration
benchmark is identified by protocol-freeze commit `d336e87`. The later
immutable records are identified within the
current technical package by run IDs `RUN-M7-SHEAF-VS-GRAPH-20260729-01` and
`RUN-M7-ROBUST-HYBRID-SHEAF-20260729-01`, whose manifests contain per-file
SHA-256 hashes. No dataset DOI has yet been minted; consequently this package
is technically auditable locally but is not submission-ready as a persistent
data release. The project-owned Northern Ghana workbook
(`data/FieldData/NorthenGhana/NorthernGhana.xlsx`) is included in that pending
release scope. The antecedent `Aquifers_Dataset_Mendeley.xlsx` is not used or
redistributed.

**Code availability.** The process-based integration code and its predeclaration
are identified by commit `d336e87`. The representation-benchmark and estimator-
diagnostic manifests record repository revision
`53beb46034d5230c1a061341a5cf2175d9af858e`, but that historical commit does not
contain the later protocol, runner, test and immutable-result files and is not
claimed as their release identifier. In the current technical package, the
exact estimator-diagnostic runner is archived by content hash
`a0ef13bde5391af62698927211cb4e701123affebb108d331795ce8596e2e191` and matches
the confirmatory lock. A new versioned repository commit and release containing
the complete representation-benchmark and estimator-diagnostic package,
followed by a persistent software DOI, remain
required before submission; neither has been fabricated in this revision. The
repository software licence is MIT. Reported runtime versions are MODFLOW 6
6.7.0, MODPATH 7 7.2.001, FloPy 3.10.0, scikit-learn 1.9.0, numpy 2.4.6, scipy
1.17.1 and pandas 2.3.3 under Python 3.14.6; PHREEQC was accessed through the
IPhreeqc interface.
