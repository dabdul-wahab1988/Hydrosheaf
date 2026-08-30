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
summarised in the benchmark architecture ([[FIGREF:FIG-1]]). Six development
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
