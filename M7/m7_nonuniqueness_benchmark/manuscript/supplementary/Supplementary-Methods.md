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
