# Supplementary Methods

## Independent generator and blinding contract

The external MODFLOW 6/MODPATH 7 generator used for the process-based
integration benchmark was inherited unchanged from an earlier independent
supporting-validation experiment and imported no Hydrosheaf code. Each
independent case was a heterogeneous, one-layer, ten-by-twenty-cell MODFLOW 6
aquifer containing a low-conductivity lens, forced by three forward MODPATH 7
particles. Four milestones were recorded per particle path, yielding twelve
observed nodes and nine true directed edges with exact model-conditioned
pathline ages per case. Every case exercised six reaction processes:
carbonate weathering, carbonate precipitation, silicate weathering,
denitrification, sulfate reduction, and iron reduction, generated in a
reducing sequence that followed nitrate, sulfate, then ferric-iron
electron-acceptor progression. Laboratory- and process-noise, non-linearity,
and conservative-ion perturbations were included so that exact recovery by
Hydrosheaf's linear reaction dictionary was not guaranteed by construction.

Tracer boundary forcing was treated as known, as it would be supplied by an
independent recharge-history measurement in a field application, but the
generator computed tritium and argon-39 responses through its own nonlinear
model rather than Hydrosheaf's. A development-calibrated argon-39 reference
scale and Student-t discrepancy likelihood, frozen before locked testing, were
used identically in the topology-conditioned-age model described below.

Development seeds were 5201 to 5206 and locked test seeds were 5301 to 5312
(twelve cases). Non-claim-bearing smoke seeds (development 5901 to 5902 and
test 5911 to 5913) were used only for pipeline smoke-testing and were
excluded from every reported result. Development cases alone fitted the
evidence-fusion coefficients described in the next section and any
classification threshold used at test time. Locked-test truth (true edges,
true ages, true reaction processes) was withheld from every inference step
and joined only by the benchmark scorer after inference was complete. No
locked-test outcome was permitted to change a feature definition, a fitted
coefficient, a topology or tracer condition, a decision threshold, or a
stopping rule; this restriction is the blinding contract referred to
throughout the main text.

## Evidence-fusion model specification

Three edge-level evidence features were computed for every candidate edge in
every case: a hydraulic log-odds term (`hydraulic_logit`), an
uncertainty-aware negative age-cost term (`negative_age_cost`), and a
negative constrained-chemistry log-objective term
(`negative_chemistry_log_objective`). Seven evidence panels were defined by
feature subset: H (hydraulic only), A (age only), C (chemistry only), and the
four combinations HA, HC, AC, and HAC.

For each panel, a logistic-regression estimator was fitted on development
cases only, using an L2 penalty of unit strength, no class weighting, and up
to 2,000 optimisation iterations. Each panel's raw feature matrix
$\mathbf{x}_e = (x_{e,1},\dots,x_{e,J})$ for candidate edge $e$ was first
standardised using the development-sample mean $\mu_j$ and standard deviation
$\sigma_j$ of feature $j$ (with any near-zero $\sigma_j$ floored at 1 to avoid
division instability):

$$
\tilde{x}_{e,j} = \frac{x_{e,j} - \mu_j}{\sigma_j}.
$$

The fitted intercept $\beta_0$ and coefficients $\beta_1,\dots,\beta_J$ were
then applied, unchanged, to every locked test edge, giving the fused
edge-truth probability reported as Equation 1 in the main text. True edges
represented 54 of 414 development candidates (prevalence 0.1304) and 108 of
827 locked-test candidates (0.1306). Unweighted fitting was the prespecified
choice because the probability estimand used the generated candidate
distribution. Class weighting would target a different effective prevalence;
a weighted model could still be evaluated and recalibrated on untouched data,
but that alternative was outside the locked protocol.

## Auxiliary shared-nuisance and M3-mechanism diagnostic

The auxiliary diagnostic used a fresh version of the independent generator and
was kept separate from the locked process-based integration numbers above.
Development seeds were 5401--5406 and locked test seeds were 5501--5512. Each
seed was evaluated at three paired nuisance levels: none (0.0), mild (0.5) and
severe (1.0), multiplying the generator's declared recharge-temperature and
14C-initial-activity error scales. All 15 non-empty subsets of hydraulic (H),
hydrochemical (C), nuclear (N) and environmental-isotope (E) evidence were
evaluated against topology (T1), age (T2) and recharge-source fraction (T3).

The T2 environmental-isotope mask was empty because the generated d18O and
d2H values were fixed at recharge and transported conservatively. The binding
control compared N with N+E and required an exact zero MAE contrast before the
mechanism result could be reported. The T3 synthetic source fraction was the
recharge x-position divided by the maximum recharge x-position in the case;
this quantity is a generator-defined target, not a field endmember estimate.

The M3 diagnostic used a non-negative convex mixture over 240 ages from 0.25
to 100 years. The independently generated 3H, 39Ar, 14C, CFC-11, CFC-12,
CFC-113, SF6, 4He and 3H/3He response curves were required to lie within
declared measurement intervals at $k=1.96$. Feasibility was tested for the
full nuclear panel and for every tracer pair, with CFC-11-containing pairs,
CFC-12-containing pairs and the 3H--3H/3He pair retained as the predeclared
mechanism contrasts. Redox strata were assigned from the withheld generator
process history only after truth-blind predictions had been written.

All locked-test inference outputs were written before the truth files were
opened. The claim-bearing output was the predeclared specificity decision; the
diagnostic was not a calibration of the USGS DGMETA correction, a field-age
validation, or a test of the USGS cause.

## Adverse control construction

Four locked test conditions were applied to every case and evaluated for
every panel: `native` (no transformation); `age_permuted` (the age feature
permuted); `hydraulic_permuted` (the hydraulic feature permuted); and
`joint_misspecified` (both features permuted). Permutation was performed
independently within each case using a case-specific pseudorandom generator
seeded from an affine function of the case's MODFLOW seed and a
condition-specific salt (17 for the age feature, 29 for the hydraulic
feature), so that the same case always received the same permutation and
different cases received independent permutations. Permutation reassigned
each case's own set of feature values across its own candidate edges,
preserving the marginal distribution of that feature within the case while
destroying its edge-specific correspondence to the true topology. No values
were drawn from outside the case, and no permutation used locked-test truth
beyond the feature values being permuted.

A predeclared univariate diagnostic was also computed for every condition:
for each candidate edge, the span between the maximum and minimum of the
three univariate panel probabilities (H, A, C) was compared against a fixed
threshold of 0.50; edges exceeding the threshold were flagged as
"conflicting". This diagnostic was locked before locked-test scoring and was
not retuned afterwards, even though it proved insensitive to the
misspecification conditions (Table S5).

## Topology-conditioned age inference

Local, topology-free groundwater age was inferred at every node from an
exact grid quadrature over an 800-point age grid spanning 0.25 to 200 years.
A lognormal age prior with median 20 years and log-scale parameter 1.2 was
combined with a Student-t discrepancy likelihood (four degrees of freedom)
between each node's observed tritium concentration and its decay-corrected
expectation; in the informative tracer regime, an equivalent Student-t term
for argon-39 was added. Writing $a$ for candidate age, $T(a)$ and
$Ar(a)$ for the generator's tritium and argon-39 decay expectations at age
$a$, and $t_{4}(\cdot;\sigma)$ for the Student-t(4) log-density with scale
$\sigma$, the local age log-posterior at node $n$ was

<!-- EQ:EQ-2 -->
$$
\log \pi_n(a) \;=\; \log\mathrm{LogNormal}(a; \mu=\log 20,\ \sigma=1.2)
\;+\; t_{4}\!\big(T_n - T(a);\ \sigma_{T,n}\big)
\;+\; \mathbb{1}_{\text{informative}}\, t_{4}\!\big(Ar_n - Ar(a);\ 1.8\big),
$$

with the tritium discrepancy scale $\sigma_{T,n}$ taken as the larger of 0.10
and 12% of the node's observed tritium concentration, and the resulting
density normalised over the 800-point grid. Fifty thousand particles per
case and tracer regime were drawn independently for each node by categorical
sampling from its own normalised local posterior; because the same particle
draws were reused across topology conditions, the four conditions differed
only in how the shared particle set was reweighted, not in the underlying
local age evidence.

Four topology conditions were then applied to the same particle set: no
connectivity assumption; a fixed 50% subset of the true generator edges
(every second true edge in generator order); the complete set of true
generator edges; and the same edge set with every edge's direction reversed.
For an assumed directed edge from node $u$ to node $v$ and particle draws
$a_u, a_v$, a soft downstream-older potential was applied with a five-year
order scale $s = 5$, incrementing the particle's log importance weight by

<!-- EQ:EQ-3 -->
$$
\Delta \log w \;=\; -\log\!\Big(1 + \exp\!\big(-\tfrac{a_v - a_u}{s}\big)\Big),
$$

summed over every edge assumed under that condition. Normalised importance
weights $w_i \propto \exp(\log w_i)$, with effective sample size

$$
\mathrm{ESS} = \left(\sum_i w_i^2\right)^{-1},
$$

were used to compute weighted posterior means, 95% weighted quantile
intervals, and normalised marginal age entropy per node. This effective
sample size is the classical importance-weighting diagnostic
[@kish1965survey]; a rule of $\mathrm{ESS}\ge 400$ out of 50,000 particles
was adopted as the stability threshold, chosen by analogy with the
per-chain effective-sample-size convention used for rank-normalised Markov
chain Monte Carlo diagnostics [@vehtari2021rhat], which addresses a related
but distinct diagnostic problem (bulk and tail effective sample size across
multiple chains, not single-particle-set importance-sampling effective
sample size) and is cited here only for that analogy, not as a direct
endorsement of this threshold for importance sampling. No independent
derivation of the 400 threshold for this specific importance-sampling
context was performed, and this is recorded as a limitation in the main
text. Cases and regimes failing this rule were reported explicitly rather
than interpreted as precise; the reversed-graph, tritium-only combination
failed in 8 of 12 locked test cases. Age mean
absolute error, bias, and 95% coverage were computed against each node's
true model-conditioned age only after reweighting, and paired case-block
bootstrap contrasts (10,000 replicates, resampling whole cases) were computed
between complete-true and no-topology, partial-true and no-topology, and
reversed and complete-true, separately for each tracer regime.

## Reaction-family bootstrap and thermodynamic constraints

Reaction-mechanism recovery was scored using Hydrosheaf's own sparse
linear reaction-inversion solver, applied only to candidate edges that the
external generator confirmed as true flow paths, so that mechanism scoring
never rewarded or penalised a topologically incorrect edge. For each locked
test case, PHREEQC speciation and saturation-index calculations were run
once on the unperturbed observations to derive thermodynamic reaction-
direction bounds [@parkhurst2013phreeqc]. Sixty-four bootstrap replicates
were then generated per case by perturbing every measured ion concentration
with independent lognormal multiplicative noise (multiplicative standard
deviation 3%), using a bootstrap- and case-specific pseudorandom generator so
that replicates were reproducible. For each replicate, the reaction solver
was fitted twice: once restricted to a core indicator-ion panel (calcium,
magnesium, sodium, potassium, bicarbonate, chloride, sulfate, and nitrate)
and once using an enhanced panel adding fluoride, iron, phosphate, and
silica, both under the same PHREEQC-derived direction bounds. For each fitted
edge, the reaction label with the largest absolute fitted extent was taken
as that replicate's dominant reaction; dominant reactions were mapped to one
of six process families (carbonate weathering, carbonate precipitation,
silicate/exchange, denitrification, sulfate reduction, and iron reduction)
for scoring.

Generation and scoring therefore shared the same six-family vocabulary. Code
independence prevented implementation leakage, but this design did not test
out-of-dictionary chemistry, alternative mineral assemblages or arbitrary
reaction combinations. All reaction claims are restricted to discrimination
among the six planted archetypes under the specified stoichiometry, 3% noise
model and two indicator panels.

Across the 64 replicates for a given edge, the empirical frequency of each
predicted family defined that edge's bootstrap family-support distribution
$\{\hat{p}_{e,f}\}_f$. The modal family was the most frequent label; the
probability assigned to the true family was $\hat{p}_{e,f^\ast}$ for the
edge's true process family $f^\ast$. Family-support entropy and the
effective number of supported families were

<!-- EQ:EQ-6 -->
$$
H_e = -\sum_{f} \hat{p}_{e,f} \log \hat{p}_{e,f}, \qquad
N^{\text{eff}}_e = \exp(H_e),
$$

with the sum taken over families with non-zero support. These edge-level
quantities were aggregated to per-process, per-panel means (Supplementary
Table S4). Carbonate weathering and carbonate precipitation were always reported
as their own rows rather than combined into an overall accuracy figure, so
that a high aggregate accuracy driven by other processes could not conceal
their non-recovery.

## Statistical metric definitions

For a candidate edge with fused probability $p_e \in (0,1)$ (numerically
bounded to $[10^{-8}, 1-10^{-8}]$) and binary truth label $y_e$, normalised
edgewise entropy was

<!-- EQ:EQ-4 -->
$$
H(p_e) = -\frac{p_e \log p_e + (1-p_e)\log(1-p_e)}{\log 2} \in [0,1].
$$

Expected calibration error (ECE) partitioned $[0,1]$ into ten equal-width
bins and summed, over non-empty bins, the edge-count-weighted absolute
difference between the bin's mean predicted probability and its mean
observed truth:

<!-- EQ:EQ-5 -->
$$
\mathrm{ECE} = \sum_{k=1}^{10} \frac{n_k}{n}\left| \bar{p}_k - \bar{y}_k \right|,
$$

where $n_k$ is the number of edges whose probability falls in bin $k$, $n$
is the total edge count, and $\bar{p}_k$, $\bar{y}_k$ are the bin's mean
predicted probability and mean observed truth. Discrimination was reported
as precision-recall area under the curve (PR-AUC) and receiver-operating-
characteristic area under the curve (ROC-AUC); PR-AUC was prioritised because
the candidate edge set is class-imbalanced (most candidate edges are not true
edges), and precision-recall curves can give a more informative view under
class imbalance [@davisgoadrich2006prcurves]. ROC-AUC was retained rather than
treating one scalar AUC summary as unconditionally superior. Calibration was
additionally reported as Brier score [@brier1950verification] and binomial
log loss. Overconfident error was the fraction of edges with a probability
at or beyond 0.1/0.9 whose 0.5-threshold classification disagreed with
truth. All aggregate statistics were also computed per case, and every
paired contrast between two evidence panels or two topology conditions used
a case-block bootstrap: for each of 10,000 replicates, cases were resampled
with replacement, and the mean paired difference across the resampled cases
was recorded; the reported 95% confidence interval was the 2.5th and 97.5th
percentile of that replicate distribution.

## Decision rules for complementarity, redundancy and negative transfer

Three decision rules were fixed before locked-test scoring and applied
identically to every contrast reported in the main text and in Tables 3 and
S2. An evidence addition was classified as complementary only if the
paired case-block bootstrap difference in PR-AUC or log loss for the fuller
panel, relative to the reduced panel, had a 95% confidence interval
excluding zero in the improving direction, or if entropy fell without a
paired worsening of Brier score and log loss. An evidence addition was
classified as redundant if it changed neither discrimination nor calibrated
uncertainty materially. An evidence addition was classified as negative
transfer if discrimination or calibration worsened, or if entropy fell while
overconfident error rose; this last pattern was the operational definition
of false confidence used throughout the study. A reduction in entropy alone,
without a paired calibration and discrimination check, was never sufficient
by itself to classify an addition as beneficial.

## Northern Ghana audit criteria

The canonical Northern Ghana workbook (`data/FieldData/NorthenGhana/
NorthernGhana.xlsx`, Dry/Wet sheets) was audited programmatically against a
fixed checklist rather than fitted or scored against any assumed truth. An
earlier revision instead audited a different, antecedent study's own derived
reprocessing of the same boreholes (`Aquifers_Dataset_Mendeley.xlsx`,
documented in `data/FieldData/NorthenGhana/SI.pdf` as the supplementary
information for a separate publication); that workbook is not this
project's own field data and has been removed from the audit (DECISIONS.md).
The audit checked for: major-ion hydrochemistry (calcium, magnesium, sodium,
potassium, bicarbonate, chloride, sulfate, and nitrate); stable water
isotopes (δ¹⁸O, δ²H); any environmental age-tracer column (tritium, ³H,
CFCs, SF₆, argon-39, or ¹⁴C, matched by column-name token); screen-interval
columns; a single static water-level measurement per well versus a
repeated, time-varying head series; masked or exact site coordinates; an
intra-season sampling-date field; an independent aquifer-type
classification; and processed graph edges. None of the last three exist in
the canonical workbook. Based on this checklist, the audit recorded which
interpretations were supportable: data-readiness and measurement-value
auditing; within-campaign seasonal chemistry hold-forward under a fixed,
disclosed, arbitrary well-revelation order (the workbook records one
dry-season and one wet-season observation per well with no intra-season
date field, so a real chronological issue sequence cannot be reconstructed);
reaction-family plausibility and equivalence-class reporting; sensitivity to
alternative assumed edge sets; and explicit non-identifiability
classification. It also recorded which interpretations were not supportable
from this workbook: field residence-time inversion; exact directed field
flow paths; screen-resolved vertical connectivity; unique field reaction
mechanisms; and a fully observed field digital twin. The audit's permitted
objective, fixed before the workbook was inspected for this study, was to
apply the framework and its component diagnostics to Ghanaian aquifer
datasets to determine which integrated interpretations are supportable under
available data and which remain non-identifiable, not to validate field
topology, age, or reaction mechanisms.

## Multiplicity-correction robustness check

The six predeclared case-block contrasts reported in Supplementary Table S2 are each
evaluated on four metrics (PR-AUC, Brier score, log loss, and mean edge
entropy), giving twenty-four confidence intervals reported without
correction for multiple comparisons in the primary analysis. As a
post-hoc robustness check performed in response to pre-submission review,
an exact two-sided sign-flip permutation test was computed for each of
these twenty-four contrasts directly from the per-case paired differences
(twelve independent cases per contrast, giving $2^{12}=4096$ equally likely
sign assignments under the null hypothesis of no systematic difference),
and a Benjamini-Hochberg false-discovery-rate correction was applied across
the full family of twenty-four tests at $q=0.05$. Full results are reported
in Supplementary Table S6.

Two results changed classification under this stricter check relative to
the bootstrap-interval screening used in the primary analysis. The native
incremental-age PR-AUC contrast (mean difference -0.0060) had an exact
permutation $p=0.063$, Benjamini-Hochberg-adjusted $p=0.070$, and did not
survive correction at $q=0.05$, whereas its accompanying entropy contrast
(mean difference -0.0006) had an exact permutation $p=0.002$,
adjusted $p=0.004$, and did survive correction. The native
incremental-hydraulics PR-AUC contrast (mean difference +0.0091) likewise
had an exact permutation $p=0.064$, adjusted $p=0.070$, and did not survive
correction, whereas its Brier score, log loss, and entropy contrasts all
survived correction with adjusted $p<0.01$. Every contrast in the native
incremental-chemistry row and every contrast in the three adverse
(permuted and jointly misspecified) rows survived correction with adjusted
$p<0.05$. These results are reported in the main text with the
corresponding qualification: PR-AUC evidence for the age and hydraulics
increments is directional but not distinguishable from zero once
multiplicity is corrected at this sample size, while the entropy, and for
hydraulics also the calibration, evidence for both increments is robust to
correction.

## Strict public-pipeline system acceptance

The system-acceptance audit used the released Hydrosheaf interfaces rather
than private benchmark-scoring utilities. For each of six locked cases, candidate
edges were generated from the public graph representation and passed through
the hydraulic evidence, local groundwater-age inference, directed affine
section solver, residual extraction, iterative reweighting, and development-
calibrated edge classifier. The four frozen conditions were full sheaf,
hydraulic only, local age only, and age permuted within case. Candidate recall,
PR-AUC, Brier score, log loss and selected F1 were recorded per case, and
10,000 case-block bootstrap replicates were used for the three full-sheaf
contrasts. An execution gate required all stages, invariants, provenance and
artifact hashes to pass. A separate scientific gate required full-sheaf
incremental value over the simple and adverse comparators. Passing the former
did not imply passing the latter; Supplementary Table S8 records both the
successful execution and failed overall incremental claim.

The selected-edge output retained all 198 generated candidates in every arm;
there was no scalar probability threshold. Each arm therefore had 53 true
positives, 145 false positives and no false negatives conditional on candidate
generation, giving conditional selected F1 0.4223. One of 54 truth edges was
missed before scoring, so the end-to-end false-negative count was one, F1 was
0.4206 and candidate recall was 0.9815. Minimum retained probabilities differed
by arm despite identical selected sets. Supplementary Table S13 gives the
complete selection-rule and confusion-count audit.

## Competence-matched sheaf-versus-graph benchmark

The representation-benchmark protocol, runner, independent generator and native directed-section
implementation were hashed before the first claim-bearing run. The generator
created 96 cases in total: 32 development and 64 locked test cases, with the
test set balanced across identity-limit, heterogeneous-affine,
incompatible-cycle, and noisy/missing-observation scenarios. Truth labels were
unavailable to prediction functions and were joined only during scoring.
Temporal, three-dimensional and vadose-zone capabilities were outside this
protocol by design because corresponding data were unavailable.

All four models used the same scalar candidate graph, observations, edge
weights, common-prior feature, three iterative reweighting solves, logistic
regularisation ($C=1$), development-fitted global calibration and frozen
selection threshold. The identity graph used the main-text affine objective
with $(\alpha_e,b_e)=(1,0)$. The affine sheaf used the generated edge maps. The
permutation control reassigned the same native maps among edges within each
case, retaining the marginal map distribution while destroying map-edge
correspondence. The edge-local weighted graph received both the common-prior
feature and the same local affine residual available to the sheaf calibrator,
but no global-section residual; it was therefore a stronger comparator than an
identity graph and was designated as the claim-bearing baseline.
This comparator discipline also addressed recent, non-peer-reviewed evidence
that identity-sheaf baselines can remain competitive in unrelated graph-
learning benchmarks [@caralt2026necessity]; that preprint motivated baseline
strength but was not treated as groundwater validation.

The exact identity-limit equivalence test compared raw residuals and calibrated
predictions with tolerance $10^{-12}$. Case-level metrics comprised PR-AUC,
ROC-AUC, Brier score, log loss, ten-bin expected calibration error, selected
F1, false-confidence rate, abstention coverage and abstention accuracy. For
the planted incompatible-cycle scenario, edge-level conflict localisation was
scored by PR-AUC against the planted incompatible-edge indicator. All
contrasts used the sign convention affine sheaf minus comparator and 10,000
case-block bootstrap resamples. The 120-row published contrast matrix was also
treated as a single family using shared resamples and a two-sided max-z
simultaneous 95% interval. Supplementary Tables S7 and S10 report the
unadjusted and family-wise results, respectively.

## Robust and hybrid sheaf estimator diagnostic

The estimator diagnostic was specified as a separate two-stage experiment
rather than a post-hoc alteration of the representation benchmark. The first
lock bound the protocol, unchanged independent generator, original
representation-benchmark runner, production directed-section solver, four
scenario definitions, development seeds 8401--8464, test seeds 8501--8628,
three robust iterations, candidate grids and 10,000 bootstrap replicates
before the estimator-diagnostic runner was implemented. Only the 64 development
cases were then generated. After model and threshold selection, the estimator-
diagnostic runner,
development manifest and complete serialised settings were SHA-256 bound in a
second confirmatory lock. The exact runner has been restored and archived under
SHA-256 `a0ef13bde5391af62698927211cb4e701123affebb108d331795ce8596e2e191`,
which matches that lock. The stored locked-test artifacts were hash-verified
without rerunning the confirmatory test. The 128 test cases had been generated
and scored in one execution. Development, locked-test, representation-
benchmark and smoke seeds were
disjoint. Test truth was present only in scoring records and was not accessed
by residual construction, standardisation, calibration or threshold fitting.

Seven arms were evaluated. `edge_local` used the observed endpoint affine
residual. `identity_graph` and `original_affine_global` used the original
three-iteration global-section residual with identity and native affine maps,
respectively. `original_hybrid` blended log-transformed local and original-
global residuals. `robust_affine_global` and `robust_hybrid` substituted a
leave-one-edge-out (LOO) global residual. `robust_hybrid_permuted` applied the
native robust-hybrid settings and calibrator unchanged to deterministically
permuted maps. Every arm received exactly three fitted inputs: the common
prior logit, one scalar residual and the local-missing indicator. Hybrids
therefore did not receive extra fitted coefficients. For a locally observed
edge, the hybrid residual was the convex combination of the log-transformed
local and global residuals; where either endpoint was missing, only the global
residual was used.

The LOO procedure removed the candidate being scored, solved the directed
section using all remaining current-weight edges and fixed observation/ridge
terms, and evaluated the held-out affine discrepancy. At each of three fixed
iterations, the median held-out residual supplied the scale, and each edge
weight was updated to its prior probability multiplied by
$\exp[-0.5(r_e/\mathrm{median}(r))^2]$, clipped to $[10^{-4},1]$. Because the
candidate was absent from the section used to calculate $r_e$, it could not
reduce its own diagnostic residual by pulling the fitted section towards its
restriction map. The same algorithm and permutation seed rule were used for
the adverse map control.

For the primary calibration regime, local weights in
$\{0,0.25,0.50,0.75,1\}$ and logistic regularisation values
$C\in\{0.1,1,10\}$ were selected using eight deterministic seed-group folds,
minimising mean per-case log loss. Ties were resolved towards a local weight
of 0.50 and then the smaller value. Selection thresholds maximising F1 were
derived from out-of-fold development probabilities and resolved nearest 0.50.
Each selected calibrator was then refitted to all development cases and
applied unchanged to the locked test. For the secondary regime, each arm was
standardised on development data, all native arm blocks were stacked, and one
common logistic model and out-of-fold threshold were fitted. The permuted
control was not separately recalibrated in either regime.

Per-case PR-AUC, Brier score and log loss were primary; ROC-AUC, expected
calibration error, selected F1, false-confidence rate, abstention
coverage/accuracy and conflict-localisation PR-AUC were secondary. The primary
claim required the selected local-first/global-fallback estimator minus edge-
local PR-AUC to have an interval above
zero and Brier/log-loss intervals below zero, identity-limit non-degradation
within the frozen margins, and favourable native-versus-permuted PR-AUC and
Brier intervals. Failure of any component prohibited overall superiority.
Mechanism contrasts were interpreted only in their prespecified directions:
original hybrid versus original global tested restoration of local evidence;
robust global versus original global tested candidate self-influence; and
robust hybrid versus original hybrid tested the increment due to LOO
robustification. Development selected local weight 1.0 for both nominal hybrid
arms, so they are described in explanatory text as local-first/global-fallback
estimators. The 560-row published contrast matrix was treated as one family
using 10,000 shared resamples and a two-sided max-z simultaneous 95% interval.
The complete unadjusted and family-wise matrices are Supplementary Tables S9
and S11.

## Generator construct validity and non-transfer boundary

The process-based MODFLOW 6/MODPATH 7 generator and the scalar affine generator
used for the representation benchmark and estimator diagnostic were separate
construct-validity tests. The former represented
two-dimensional saturated flow, pathline-derived age, nonlinear tracer
responses and six planted reaction archetypes. The latter directly planted
identity, heterogeneous-affine, incompatible-cycle and noisy/missing relations
without simulating groundwater flow, age or chemistry. Results from either
generator were not treated as replication in the other.

Neither generator supplied field truth. The Northern Ghana workbook was used
only for data-support and truth-free hold-forward audits, not for validating
connectivity, groundwater age or reaction mechanisms. The claim-bearing next
step is an independent generator family or a field dataset with separately
measured connectivity.

## Post-review precision and practical-magnitude audit

No prospective minimum important differences or power analysis were specified
before the locked tests. The revision therefore labels the following as post-
review planning margins, not field-validated decision thresholds: PR-AUC 0.02,
Brier score 0.01, log loss 0.02, age MAE 0.25 years, interval width 0.50 years,
coverage 0.05 and modal reaction-family accuracy 0.10. Twenty-thousand
empirical simulations resampled development-case effects at the planned test
sizes. The process-based integration and representation-benchmark models were
evaluated on their development fitting cases, so their planning estimates may
be optimistic. The estimator diagnostic used eightfold out-of-fold
development predictions. Development topology-age and reaction summaries had
not been archived; locked vectors were used only for future replication
planning and were never used to justify the completed tests. Realised effects
were also divided by the comparator mean and checked against the post-review
margins. Supplementary Table S12 reports all assumptions and probabilities.

## Computational environment

Groundwater flow was simulated with MODFLOW 6 (build 6.7.0, dated
05 February 2026) and particle tracking with MODPATH 7 (version 7.2.001),
both executed through FloPy 3.10.0; executable identity was confirmed by
SHA-256 checksum and recorded per case in the run provenance. PHREEQC
speciation and saturation-index calculations were performed through the
IPhreeqc interface. The evidence-fusion logistic-regression models were
fitted with scikit-learn 1.9.0; numpy 2.4.6, scipy 1.17.1, and pandas 2.3.3
were used for array operations, statistical functions, and tabular data
handling; the pipeline was run under Python 3.14.6. The locked result
tables, case-level provenance records, and analysis code are archived
alongside the manuscript's data availability statement.

# Supplementary Figures and Tables

**Conditional evidence integration and the incremental contribution of sheaf structure in controlled-synthetic groundwater benchmarks**

The complete authoritative supplementary tables are the version-controlled CSV files in `tables/publication/`. The compact views below preserve the claim-bearing comparisons while keeping the Word and PDF tables legible; omitted columns and per-case or per-edge rows remain available without loss in the cited CSV. Tables S1-S6 derive from locked process-based integration outputs, Table S7 from the prospectively locked representation benchmark, Table S8 from the strict public-pipeline acceptance run, Table S9 from the prospectively locked estimator diagnostic, and Tables S10-S13 from the post-review family-wise, precision and selection audit. Table S14 is the auxiliary M7.6 controlled-synthetic M3-mechanism diagnostic.

## Figure S1

![](figures/supporting_validation/figure_s1_model_domain_map.png)

**Figure S1.** Locked synthetic model domain (representative realization 4101). The MODFLOW/MODPATH truth network, observation nodes and hydraulic heads are shown in synthetic model-space coordinates. This is not a geographic map.

## Table S1

**Table S1.** Evidence-panel performance in every locked process-based integration condition. Compact view of `tables/publication/tableS1_all_evidence_conditions.csv`.

| Condition          | Panel   |   PR-AUC |   Brier |   Log loss |   Mean entropy |
|:-------------------|:--------|---------:|--------:|-----------:|---------------:|
| native             | H       |   0.1765 |  0.1076 |     0.3467 |         0.5029 |
| native             | A       |   0.1112 |  0.1087 |     0.3558 |         0.524  |
| native             | C       |   0.4587 |  0.0944 |     0.2951 |         0.4551 |
| native             | HA      |   0.1114 |  0.1076 |     0.3467 |         0.5042 |
| native             | HC      |   0.4846 |  0.088  |     0.2676 |         0.4221 |
| native             | AC      |   0.471  |  0.089  |     0.2736 |         0.4361 |
| native             | HAC     |   0.4805 |  0.088  |     0.2677 |         0.4215 |
| age_permuted       | H       |   0.1765 |  0.1076 |     0.3467 |         0.5029 |
| age_permuted       | A       |   0.1397 |  0.1179 |     0.4539 |         0.524  |
| age_permuted       | C       |   0.4587 |  0.0944 |     0.2951 |         0.4551 |
| age_permuted       | HA      |   0.1825 |  0.1096 |     0.3588 |         0.4726 |
| age_permuted       | HC      |   0.4846 |  0.088  |     0.2676 |         0.4221 |
| age_permuted       | AC      |   0.3521 |  0.1016 |     0.3588 |         0.4236 |
| age_permuted       | HAC     |   0.4165 |  0.0914 |     0.2781 |         0.4014 |
| hydraulic_permuted | H       |   0.1302 |  0.1191 |     0.467  |         0.5029 |
| hydraulic_permuted | A       |   0.1112 |  0.1087 |     0.3558 |         0.524  |
| hydraulic_permuted | C       |   0.4587 |  0.0944 |     0.2951 |         0.4551 |
| hydraulic_permuted | HA      |   0.1173 |  0.1155 |     0.4358 |         0.4713 |
| hydraulic_permuted | HC      |   0.3803 |  0.1032 |     0.3717 |         0.4054 |
| hydraulic_permuted | AC      |   0.471  |  0.089  |     0.2736 |         0.4361 |
| hydraulic_permuted | HAC     |   0.3971 |  0.1001 |     0.3482 |         0.3879 |
| joint_misspecified | H       |   0.1302 |  0.1191 |     0.467  |         0.5029 |
| joint_misspecified | A       |   0.1397 |  0.1179 |     0.4539 |         0.524  |
| joint_misspecified | C       |   0.4587 |  0.0944 |     0.2951 |         0.4551 |
| joint_misspecified | HA      |   0.1371 |  0.1201 |     0.463  |         0.4729 |
| joint_misspecified | HC      |   0.3803 |  0.1032 |     0.3717 |         0.4054 |
| joint_misspecified | AC      |   0.3521 |  0.1016 |     0.3588 |         0.4236 |
| joint_misspecified | HAC     |   0.3305 |  0.105  |     0.3682 |         0.3845 |

## Table S2

**Table S2.** Predeclared evidence contrasts for discrimination and calibration (10,000 case-block bootstrap replicates). Complete five-metric record: `tables/publication/tableS2_case_block_contrasts.csv`.

| Contrast                      | Metric   |   Difference |   CI low |   CI high |   Cases |
|:------------------------------|:---------|-------------:|---------:|----------:|--------:|
| native incremental age        | PR-AUC   |      -0.006  |  -0.0122 |   -0.0011 |      12 |
| native incremental age        | Brier    |       0.0001 |  -0      |    0.0001 |      12 |
| native incremental age        | Log loss |       0      |  -0.0001 |    0.0002 |      12 |
| native incremental chemistry  | PR-AUC   |       0.4471 |   0.3575 |    0.5401 |      12 |
| native incremental chemistry  | Brier    |      -0.0196 |  -0.0213 |   -0.0176 |      12 |
| native incremental chemistry  | Log loss |      -0.0791 |  -0.085  |   -0.072  |      12 |
| native incremental hydraulics | PR-AUC   |       0.0091 |   0.001  |    0.0198 |      12 |
| native incremental hydraulics | Brier    |      -0.001  |  -0.0012 |   -0.0008 |      12 |
| native incremental hydraulics | Log loss |      -0.0059 |  -0.007  |   -0.005  |      12 |
| permuted age increment        | PR-AUC   |      -0.0754 |  -0.1353 |   -0.0148 |      12 |
| permuted age increment        | Brier    |       0.0034 |   0.0012 |    0.0055 |      12 |
| permuted age increment        | Log loss |       0.0105 |   0.003  |    0.0176 |      12 |
| permuted hydraulic increment  | PR-AUC   |      -0.0686 |  -0.112  |   -0.0271 |      12 |
| permuted hydraulic increment  | Brier    |       0.011  |   0.0055 |    0.0164 |      12 |
| permuted hydraulic increment  | Log loss |       0.0745 |   0.0404 |    0.1091 |      12 |
| joint misspecification        | PR-AUC   |      -0.139  |  -0.2041 |   -0.0742 |      12 |
| joint misspecification        | Brier    |       0.0106 |   0.0055 |    0.0163 |      12 |
| joint misspecification        | Log loss |       0.073  |   0.0401 |    0.1073 |      12 |

## Table S3

**Table S3.** Topology-to-age sensitivity aggregated over twelve locked cases. Values are means except median ESS; the complete CSV also reports bias, entropy, order violations and importance weights. Complete 96-row case audit: `tables/publication/tableS3_topology_age_sensitivity.csv`.

| Regime      | Topology   |   n |   MAE (yr) |   Coverage |   Width (yr) |       ESS |   Stable |
|:------------|:-----------|----:|-----------:|-----------:|-------------:|----------:|---------:|
| Informative | None       |  12 |     2.7644 |     0.9306 |      12.5035 | 50000     |   1      |
| Informative | Partial    |  12 |     2.7399 |     0.9236 |      12.3584 | 48866.4   |   1      |
| Informative | Complete   |  12 |     2.7025 |     0.9444 |      12.2511 | 48301.1   |   1      |
| Informative | Reversed   |  12 |     2.9846 |     0.9167 |      12.5469 |  2558     |   1      |
| 3H only     | None       |  12 |     4.7501 |     0.9583 |      25.7014 | 50000     |   1      |
| 3H only     | Partial    |  12 |     4.6927 |     0.9583 |      25.2452 | 47183.5   |   1      |
| 3H only     | Complete   |  12 |     4.586  |     0.9722 |      24.7896 | 45117.1   |   1      |
| 3H only     | Reversed   |  12 |     4.6101 |     0.9167 |      21.396  |   227.926 |   0.3333 |

## Table S4

**Table S4.** Reaction-family non-uniqueness aggregated by benchmark tier and true process (C, core; E, enhanced). Complete 216-row edgewise audit: `tables/publication/tableS4_reaction_edge_nonuniqueness.csv`.

| Tier   | Process            |   n |   Accuracy |   True prob. |   Entropy |   Eff. families |
|:-------|:-------------------|----:|-----------:|-------------:|----------:|----------------:|
| C      | Carbonate weather. |  24 |        0   |       0      |    0.1155 |          1.1378 |
| E      | Carbonate weather. |  24 |        0   |       0      |    0.0902 |          1.1103 |
| C      | Carbonate precip.  |  12 |        0   |       0      |    0.0338 |          1.0417 |
| E      | Carbonate precip.  |  12 |        0   |       0      |    0      |          1      |
| C      | Silicate weather.  |  24 |        1   |       0.8639 |    0.349  |          1.4663 |
| E      | Silicate weather.  |  24 |        1   |       0.9505 |    0.1613 |          1.1987 |
| C      | Denitrification    |  24 |        1   |       0.9941 |    0.0191 |          1.0226 |
| E      | Denitrification    |  24 |        1   |       0.9694 |    0.0623 |          1.0805 |
| C      | Sulfate reduction  |  12 |        1   |       0.9974 |    0.0134 |          1.014  |
| E      | Sulfate reduction  |  12 |        1   |       0.9857 |    0.0577 |          1.0636 |
| C      | Iron reduction     |  12 |        0   |       0      |    0.3413 |          1.4493 |
| E      | Iron reduction     |  12 |        0.5 |       0.4688 |    0.5256 |          1.7147 |

## Table S5

**Table S5.** Conflict diagnostics under native and adverse evidence conditions. The complete CSV additionally reports error rates within flagged and concordant subsets: `tables/publication/tableS5_conflict_diagnostics.csv`.

| Condition      |   N |   Flagged |   Span |   Error |   Overconf. |
|:---------------|----:|----------:|-------:|--------:|------------:|
| Native         | 827 |         0 | 0.092  |  0.1306 |      0      |
| Age perm.      | 827 |         0 | 0.1193 |  0.1306 |      0.0036 |
| Hyd. perm.     | 827 |         0 | 0.1205 |  0.1306 |      0.0351 |
| Joint misspec. | 827 |         0 | 0.123  |  0.1306 |      0.0387 |

## Table S6

**Table S6.** Exact paired-permutation tests with Benjamini-Hochberg correction for the predeclared process-based integration contrast family. Complete record: `tables/publication/tableS6_multiplicity_correction.csv`.

| Contrast       | Metric   |   n |   Difference |   Exact p |   BH p |
|:---------------|:---------|----:|-------------:|----------:|-------:|
| Native +age    | PR-AUC   |  12 |      -0.006  |    0.0625 | 0.0703 |
| Native +age    | Brier    |  12 |       0.0001 |    0.1875 | 0.1957 |
| Native +age    | Log loss |  12 |       0      |    0.8491 | 0.8491 |
| Native +age    | Entropy  |  12 |      -0.0006 |    0.002  | 0.0043 |
| Native +chem.  | PR-AUC   |  12 |       0.4471 |    0.0005 | 0.0012 |
| Native +chem.  | Brier    |  12 |      -0.0196 |    0.0005 | 0.0012 |
| Native +chem.  | Log loss |  12 |      -0.0791 |    0.0005 | 0.0012 |
| Native +chem.  | Entropy  |  12 |      -0.0827 |    0.0005 | 0.0012 |
| Native +hyd.   | PR-AUC   |  12 |       0.0091 |    0.0645 | 0.0703 |
| Native +hyd.   | Brier    |  12 |      -0.001  |    0.0005 | 0.0012 |
| Native +hyd.   | Log loss |  12 |      -0.0059 |    0.0005 | 0.0012 |
| Native +hyd.   | Entropy  |  12 |      -0.0146 |    0.0005 | 0.0012 |
| Permuted +age  | PR-AUC   |  12 |      -0.0754 |    0.041  | 0.0492 |
| Permuted +age  | Brier    |  12 |       0.0034 |    0.0181 | 0.0241 |
| Permuted +age  | Log loss |  12 |       0.0105 |    0.0264 | 0.0333 |
| Permuted +age  | Entropy  |  12 |      -0.0207 |    0.0005 | 0.0012 |
| Permuted +hyd. | PR-AUC   |  12 |      -0.0686 |    0.0151 | 0.0214 |
| Permuted +hyd. | Brier    |  12 |       0.011  |    0.0049 | 0.0073 |
| Permuted +hyd. | Log loss |  12 |       0.0745 |    0.0029 | 0.005  |
| Permuted +hyd. | Entropy  |  12 |      -0.0482 |    0.0005 | 0.0012 |
| Joint misspec. | PR-AUC   |  12 |      -0.139  |    0.0044 | 0.007  |
| Joint misspec. | Brier    |  12 |       0.0106 |    0.0024 | 0.0045 |
| Joint misspec. | Log loss |  12 |       0.073  |    0.0024 | 0.0045 |
| Joint misspec. | Entropy  |  12 |      -0.0706 |    0.0005 | 0.0012 |

## Table S7

**Table S7.** Claim-bearing representation-benchmark sheaf-versus-graph contrasts, including the identity limit and incompatible-cycle diagnostics (10,000 case-block bootstrap replicates). Complete 120-row contrast matrix: `tables/publication/tableS7_sheaf_vs_graph_contrasts.csv`.

| Scenario       | Comparator   | Metric          |   Cases |   Difference |   CI low |   CI high |
|:---------------|:-------------|:----------------|--------:|-------------:|---------:|----------:|
| All            | Edge-local   | PR-AUC          |      64 |       0.0097 |  -0.0054 |    0.0248 |
| All            | Edge-local   | Brier           |      64 |       0.0005 |  -0.0033 |    0.0044 |
| All            | Edge-local   | Log loss        |      64 |       0.0117 |   0.0008 |    0.0232 |
| All            | Edge-local   | Selected F1     |      64 |       0.0016 |  -0.0156 |    0.0183 |
| All            | Edge-local   | Conflict PR-AUC |      16 |       0.0689 |   0.0466 |    0.0914 |
| All            | Identity     | PR-AUC          |      64 |       0.0854 |   0.0666 |    0.105  |
| All            | Identity     | Brier           |      64 |      -0.0193 |  -0.0235 |   -0.0152 |
| All            | Identity     | Log loss        |      64 |      -0.0472 |  -0.0573 |   -0.0372 |
| All            | Identity     | Selected F1     |      64 |       0.0588 |   0.0418 |    0.0758 |
| All            | Identity     | Conflict PR-AUC |      16 |       0.1098 |   0.0912 |    0.1277 |
| All            | Permuted-map | PR-AUC          |      64 |       0.0909 |   0.0705 |    0.1117 |
| All            | Permuted-map | Brier           |      64 |      -0.0215 |  -0.0263 |   -0.0169 |
| All            | Permuted-map | Log loss        |      64 |      -0.0546 |  -0.066  |   -0.0433 |
| All            | Permuted-map | Selected F1     |      64 |       0.0729 |   0.0491 |    0.0964 |
| All            | Permuted-map | Conflict PR-AUC |      16 |       0.052  |   0.0323 |    0.07   |
| Identity limit | Identity     | PR-AUC          |      16 |       0      |   0      |    0      |
| Identity limit | Identity     | Brier           |      16 |       0      |   0      |    0      |
| Identity limit | Identity     | Log loss        |      16 |       0      |   0      |    0      |
| Identity limit | Identity     | Selected F1     |      16 |       0      |   0      |    0      |
| Incompatible   | Edge-local   | PR-AUC          |      16 |       0.0483 |   0.0258 |    0.0712 |
| Incompatible   | Edge-local   | Conflict PR-AUC |      16 |       0.0689 |   0.0467 |    0.0912 |

## Table S8

**Table S8.** Strict public-pipeline system acceptance. The execution criterion passed, whereas a general full-sheaf incremental-performance claim did not.

### Condition summary

| Condition      |   Cases |   Recall |   PR-AUC |   Brier |   Log loss |   Selected F1 |
|:---------------|--------:|---------:|---------:|--------:|-----------:|--------------:|
| age permuted   |       6 |   0.9815 |   0.3211 |  0.2378 |     0.6953 |        0.4222 |
| full sheaf     |       6 |   0.9815 |   0.3075 |  0.2171 |     0.6243 |        0.4222 |
| hydraulic only |       6 |   0.9815 |   0.3272 |  0.6068 |     7.764  |        0.4222 |
| local age      |       6 |   0.9815 |   0.2488 |  0.2256 |     0.6433 |        0.4222 |

### Paired case-block contrasts

| Left       | Comparator     | Metric      |   Cases |   Difference |   CI low |   CI high |
|:-----------|:---------------|:------------|--------:|-------------:|---------:|----------:|
| full sheaf | hydraulic only | PR-AUC      |       6 |      -0.0197 |  -0.0355 |   -0.0039 |
| full sheaf | hydraulic only | Brier       |       6 |      -0.3897 |  -0.4019 |   -0.3801 |
| full sheaf | hydraulic only | Selected F1 |       6 |       0      |   0      |    0      |
| full sheaf | local age      | PR-AUC      |       6 |       0.0586 |   0.0386 |    0.0777 |
| full sheaf | local age      | Brier       |       6 |      -0.0085 |  -0.0098 |   -0.0072 |
| full sheaf | local age      | Selected F1 |       6 |       0      |   0      |    0      |
| full sheaf | age permuted   | PR-AUC      |       6 |      -0.0136 |  -0.0622 |    0.0347 |
| full sheaf | age permuted   | Brier       |       6 |      -0.0207 |  -0.0429 |    0.0014 |
| full sheaf | age permuted   | Selected F1 |       6 |       0      |   0      |    0      |

Complete record: `tables/publication/tableS8_public_pipeline_acceptance.csv`.

## Table S9

**Table S9.** Claim-bearing estimator-diagnostic robust/hybrid contrasts under the primary separately cross-fitted calibration regime. Differences are left minus comparator; lower Brier score and log loss are favourable. Estimator abbreviations are EL, edge-local; OG, original affine-global; OH, original hybrid; RG, robust affine-global; RH, robust hybrid; and PRH, permuted robust hybrid. PRAUC denotes precision-recall area under the curve. The full 560-row matrix, including the shared-calibrator regime and secondary metrics, is `tables/publication/tableS9_robust_hybrid_contrasts.csv`.

| Scenario         | Comparison   | Metric   |   n | Difference [95% CI]        |
|:-----------------|:-------------|:---------|----:|:---------------------------|
| All              | RH vs EL     | PRAUC    | 128 | +0.0200 [+0.0073, +0.0324] |
| All              | RH vs EL     | Brier    | 128 | -0.0015 [-0.0042, +0.0011] |
| All              | RH vs EL     | LogLoss  | 128 | +0.0033 [-0.0034, +0.0101] |
| All              | RH vs PRH    | PRAUC    | 128 | +0.0441 [+0.0321, +0.0568] |
| All              | RH vs PRH    | Brier    | 128 | -0.0123 [-0.0149, -0.0098] |
| All              | RH vs PRH    | LogLoss  | 128 | -0.0311 [-0.0376, -0.0250] |
| All              | OH vs OG     | PRAUC    | 128 | +0.0029 [+0.0010, +0.0047] |
| All              | OH vs OG     | Brier    | 128 | -0.0006 [-0.0011, -0.0002] |
| All              | OH vs OG     | LogLoss  | 128 | -0.0015 [-0.0025, -0.0005] |
| All              | RG vs OG     | PRAUC    | 128 | -0.0048 [-0.0097, +0.0005] |
| All              | RG vs OG     | Brier    | 128 | +0.0029 [+0.0019, +0.0038] |
| All              | RG vs OG     | LogLoss  | 128 | +0.0073 [+0.0049, +0.0097] |
| All              | RH vs OH     | PRAUC    | 128 | -0.0008 [-0.0056, +0.0044] |
| All              | RH vs OH     | Brier    | 128 | +0.0017 [+0.0008, +0.0025] |
| All              | RH vs OH     | LogLoss  | 128 | +0.0043 [+0.0021, +0.0065] |
| Identity limit   | RH vs EL     | PRAUC    |  32 | +0.0126 [-0.0086, +0.0338] |
| Identity limit   | RH vs EL     | Brier    |  32 | +0.0004 [-0.0042, +0.0049] |
| Identity limit   | RH vs EL     | LogLoss  |  32 | +0.0014 [-0.0087, +0.0113] |
| Heterog. affine  | RH vs EL     | PRAUC    |  32 | -0.0100 [-0.0299, +0.0121] |
| Heterog. affine  | RH vs EL     | Brier    |  32 | +0.0040 [-0.0011, +0.0087] |
| Heterog. affine  | RH vs EL     | LogLoss  |  32 | +0.0251 [+0.0121, +0.0373] |
| Incompat. cycles | RH vs EL     | PRAUC    |  32 | +0.0437 [+0.0125, +0.0751] |
| Incompat. cycles | RH vs EL     | Brier    |  32 | -0.0033 [-0.0088, +0.0020] |
| Incompat. cycles | RH vs EL     | LogLoss  |  32 | -0.0056 [-0.0194, +0.0073] |
| Noisy/missing    | RH vs EL     | PRAUC    |  32 | +0.0335 [+0.0110, +0.0559] |
| Noisy/missing    | RH vs EL     | Brier    |  32 | -0.0072 [-0.0127, -0.0017] |
| Noisy/missing    | RH vs EL     | LogLoss  |  32 | -0.0076 [-0.0222, +0.0061] |

## Table S10

**Table S10.** Selected representation-benchmark contrasts with simultaneous 95% intervals controlling all 120 published contrasts as one family. Complete record: `tables/publication/tableS10_m7_4_multiplicity_adjusted.csv`.

| Scenario            | Comparison                            | Metric                       |   n | Difference [simultaneous 95% CI]   | FWER support   |
|:--------------------|:--------------------------------------|:-----------------------------|----:|:-----------------------------------|:---------------|
| all                 | affine sheaf vs weighted graph        | pr_auc                       |  64 | +0.0097 [-0.0149, +0.0343]         | No             |
| all                 | affine sheaf vs weighted graph        | brier                        |  64 | +0.0005 [-0.0058, +0.0069]         | No             |
| all                 | affine sheaf vs weighted graph        | log_loss                     |  64 | +0.0117 [-0.0066, +0.0301]         | No             |
| all                 | affine sheaf vs graph laplacian       | pr_auc                       |  64 | +0.0854 [+0.0539, +0.1169]         | Yes            |
| all                 | affine sheaf vs graph laplacian       | brier                        |  64 | -0.0193 [-0.0261, -0.0126]         | Yes            |
| all                 | affine sheaf vs graph laplacian       | log_loss                     |  64 | -0.0472 [-0.0638, -0.0306]         | Yes            |
| all                 | affine sheaf vs affine sheaf permuted | pr_auc                       |  64 | +0.0909 [+0.0571, +0.1246]         | Yes            |
| all                 | affine sheaf vs affine sheaf permuted | brier                        |  64 | -0.0215 [-0.0292, -0.0139]         | Yes            |
| all                 | affine sheaf vs affine sheaf permuted | log_loss                     |  64 | -0.0546 [-0.0734, -0.0357]         | Yes            |
| incompatible_cycles | affine sheaf vs weighted graph        | pr_auc                       |  16 | +0.0483 [+0.0095, +0.0871]         | Yes            |
| incompatible_cycles | affine sheaf vs weighted graph        | conflict_localisation_pr_auc |  16 | +0.0689 [+0.0318, +0.1061]         | Yes            |

## Table S11

**Table S11.** Selected estimator-diagnostic contrasts with simultaneous 95% intervals controlling all 560 published contrasts as one family. The selected robust-hybrid arm had local weight 1.0 and is the local-first/global-fallback estimator in the main text. Complete record: `tables/publication/tableS11_m7_5_multiplicity_adjusted.csv`.

| Scenario   | Comparison                                     | Metric                       |   n | Difference [simultaneous 95% CI]   | FWER support   |
|:-----------|:-----------------------------------------------|:-----------------------------|----:|:-----------------------------------|:---------------|
| all        | robust hybrid vs edge local                    | pr_auc                       | 128 | +0.0200 [-0.0055, +0.0454]         | No             |
| all        | robust hybrid vs edge local                    | brier                        | 128 | -0.0015 [-0.0069, +0.0039]         | No             |
| all        | robust hybrid vs edge local                    | log_loss                     | 128 | +0.0033 [-0.0102, +0.0168]         | No             |
| all        | robust hybrid vs edge local                    | conflict_localisation_pr_auc |  32 | +0.0770 [+0.0389, +0.1151]         | Yes            |
| all        | robust hybrid vs robust hybrid permuted        | pr_auc                       | 128 | +0.0441 [+0.0192, +0.0690]         | Yes            |
| all        | robust hybrid vs robust hybrid permuted        | brier                        | 128 | -0.0123 [-0.0175, -0.0071]         | Yes            |
| all        | robust hybrid vs robust hybrid permuted        | log_loss                     | 128 | -0.0311 [-0.0440, -0.0183]         | Yes            |
| all        | robust hybrid vs robust hybrid permuted        | conflict_localisation_pr_auc |  32 | +0.0505 [+0.0154, +0.0855]         | Yes            |
| all        | robust affine global vs original affine global | brier                        | 128 | +0.0029 [+0.0010, +0.0049]         | No             |
| all        | robust affine global vs original affine global | log_loss                     | 128 | +0.0073 [+0.0025, +0.0122]         | No             |

## Table S12

**Table S12.** Post-review empirical precision and future-replication planning audit (20,000 simulations). Margins were not prespecified or field validated; POST_TEST rows are not evidence for completed tests. Complete record: `tables/publication/tableS12_precision_and_power.csv`.

| Design                   | Metric         |   Source n |   Planned n |   Margin |   P(CI favourable) |   P(CI clears margin) | Status               |
|:-------------------------|:---------------|-----------:|------------:|---------:|-------------------:|----------------------:|:---------------------|
| Evidence panel           | pr_auc         |          6 |          12 |     0.02 |             1      |                1      | Development planning |
| Evidence panel           | brier          |          6 |          12 |     0.01 |             1      |                1      | Development planning |
| Evidence panel           | log_loss       |          6 |          12 |     0.02 |             1      |                1      | Development planning |
| Representation           | pr_auc         |         32 |          64 |     0.02 |             0.1928 |                0.0001 | Development planning |
| Representation           | brier          |         32 |          64 |     0.01 |             0.0138 |                0      | Development planning |
| Representation           | log_loss       |         32 |          64 |     0.02 |             0.0001 |                0      | Development planning |
| Estimator                | pr_auc         |         64 |         128 |     0.02 |             0.5475 |                0.0006 | Development planning |
| Estimator                | brier          |         64 |         128 |     0.01 |             0.1228 |                0      | Development planning |
| Estimator                | log_loss       |         64 |         128 |     0.02 |             0.0022 |                0      | Development planning |
| Topology-age: two-tracer | age MAE        |         12 |          12 |     0.25 |             1      |                0      | Replication planning |
| Topology-age: two-tracer | interval width |         12 |          12 |     0.5  |             1      |                0      | Replication planning |
| Topology-age: tritium    | age MAE        |         12 |          12 |     0.25 |             1      |                0      | Replication planning |
| Topology-age: tritium    | interval width |         12 |          12 |     0.5  |             1      |                1      | Replication planning |
| Reaction recovery        | modal accuracy |         12 |          12 |     0.1  |             0.9268 |                0.0004 | Replication planning |

## Table S13

**Table S13.** Public-pipeline selection and confusion-count audit. All generated candidates were retained; no scalar probability threshold was applied. Complete record: `tables/publication/tableS13_public_pipeline_selection.csv`.

| Condition      |   Candidates |   Selected |   Min probability |   TP |   FP |   Conditional FN |   End-to-end FN |   End-to-end F1 |   Candidate recall |
|:---------------|-------------:|-----------:|------------------:|-----:|-----:|-----------------:|----------------:|----------------:|-------------------:|
| age_permuted   |          198 |        198 |            0.0222 |   53 |  145 |                0 |               1 |          0.4206 |             0.9815 |
| full_sheaf     |          198 |        198 |            0.1921 |   53 |  145 |                0 |               1 |          0.4206 |             0.9815 |
| hydraulic_only |          198 |        198 |            0.0257 |   53 |  145 |                0 |               1 |          0.4206 |             0.9815 |
| local_age      |          198 |        198 |            0.317  |   53 |  145 |                0 |               1 |          0.4206 |             0.9815 |

## Table S14

**Table S14.** Auxiliary M7.6 controlled-synthetic M3-mechanism diagnostic. The complete result record is `tables/publication/tableS14_m7_6_m3_mechanism.csv`; this table is not field validation and does not identify the USGS cause.

| Contrast                                            |   Difference |   CI low |   CI high |   Cases |   Bootstrap | Decision                                           |
|:----------------------------------------------------|-------------:|---------:|----------:|--------:|------------:|:---------------------------------------------------|
| severe minus none full infeasibility                |       0.2882 |   0.2118 |    0.3646 |      12 |       10000 | shared nuisance increased full-panel infeasibility |
| reducing minus nonreducing cfc11 pair rate severe   |       0.7188 |   0.6667 |    0.75   |      12 |       10000 | CFC-11 contrast positive                           |
| reducing minus nonreducing cfc12 pair rate severe   |       0.7396 |   0.724  |    0.75   |      12 |       10000 | CFC-12 specificity control failed                  |
| reducing minus nonreducing tritium pair rate severe |       0.3229 |   0.1354 |    0.4896 |      12 |       10000 | positive association not selective to CFC family   |
| E added to N T2 MAE none                            |       0      |   0      |    0      |      12 |       10000 | binding isotope-age control passed                 |
