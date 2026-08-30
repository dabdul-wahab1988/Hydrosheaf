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
