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
