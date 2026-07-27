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
summarised in the benchmark architecture ([[FIGREF:FIG-1]]). Six development
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
