# Supplementary Methods

## S1. Derivation of the coordinate-descent update

The estimator solved by the main text (Equation 2) minimised a bound-
constrained, elastic-net-penalised weighted least-squares objective combining
the lasso penalty [@tibshirani1996lasso] with the ridge-type term of the
elastic net [@zou2005elasticnet]. Writing
$\mathbf{u} = \mathbf{W}^{1/2}\mathbf{r}$ and $\mathbf{A} =
\mathbf{W}^{1/2}\mathbf{S}^{\mathsf T}$, the objective became

<!-- EQ:EQ-S1 -->
$$
J(\mathbf{z}) = \frac{1}{2}\|\mathbf{u} - \mathbf{A}\mathbf{z}\|_2^2
+ \lambda_1 \sum_j p_j |z_j| + \frac{\lambda_2}{2}\|\mathbf{z}\|_2^2 .
$$

$J$ is separable in each coordinate $z_j$ conditional on the current values
of the other coefficients, which permits coordinate descent. Holding all
coordinates except $z_j$ fixed, define the partial residual $\mathbf{u}_{-j}
= \mathbf{u} - \sum_{k \neq j} z_k \mathbf{a}_k$, where $\mathbf{a}_k$ is the
$k$-th column of $\mathbf{A}$. The one-dimensional subproblem in $z_j$ is

<!-- EQ:EQ-S2 -->
$$
\min_{z_j} \; \tfrac{1}{2}\|\mathbf{u}_{-j} - z_j \mathbf{a}_j\|_2^2
+ \lambda_1 p_j |z_j| + \tfrac{\lambda_2}{2} z_j^2 .
$$

Expanding the quadratic term and collecting coefficients of $z_j$ gives a
one-dimensional lasso-type problem with slope $c_j = \mathbf{a}_j^{\mathsf
T}\mathbf{u}_{-j}$ and curvature $\|\mathbf{a}_j\|_2^2 + \lambda_2$. Its
closed-form minimiser is the soft-thresholding operator

<!-- EQ:EQ-S3 -->
$$
z_j \leftarrow \frac{
  \mathrm{soft}\!\left(c_j,\; \lambda_1 p_j\right)
}{
  \|\mathbf{a}_j\|_2^2 + \lambda_2
}, \qquad
\mathrm{soft}(c, t) = \mathrm{sign}(c)\max(|c| - t, 0),
$$

applied after each coordinate update, followed by projection onto the
reaction-specific thermodynamic bounds $[z_j^{\min}, z_j^{\max}]$. Bound
projection was applied after the soft-thresholding step rather than folded
into it, because the elastic-net shrinkage and the thermodynamic feasibility
constraint address different sources of ill-posedness (statistical
regularisation and physical plausibility respectively) and combining them
inside a single closed form would obscure which constraint was active for a
given coefficient. Coordinates were updated in a fixed cyclic order over a
maximum of 250 iterations, with convergence declared when the maximum
absolute change in $\mathbf{z}$ across a full cycle fell below a
predeclared tolerance on the same scale as the locked objective convention
described in the main text. Non-convergent scenarios were flagged rather
than silently accepted, and their convergence status was carried into the
reported runtime and reliability metrics rather than discarded.

Coefficients were initialised at zero for the first point on the
regularisation grid and warm-started from the previous grid point's solution
for subsequent, less-penalised or more-penalised points, which reduced the
number of coordinate cycles needed to reach convergence when scanning a
grid of $(\lambda_1, \lambda_2)$ pairs and made the resulting regularisation
paths (Supplementary Figures S6, S7) smoother and less sensitive to solver
initialisation than independent cold starts at every grid point would have
been. The cyclic coordinate order was fixed rather than randomised across
scenarios so that any residual order-dependence in the converged solution
was itself a diagnosable and reproducible property of the reaction
dictionary rather than a source of run-to-run variability. Where a scenario
failed to converge within the iteration budget, its diagnostics were still
computed on the best iterate reached, and the non-convergence flag was
propagated to every downstream metric (phase recovery, extent error,
held-out-ion prediction, and Mechanism Resolution Score inputs) that
depended on that fit, so that a reader could distinguish a genuine recovery
failure from a solver failure to reach its own stopping criterion.

## S2. Bound handling and sign conventions

Each reaction $j$ was assigned an admissible extent interval
$[z_j^{\min}, z_j^{\max}]$ derived from its thermodynamic feasibility under
the local water chemistry: a mineral undersaturated with respect to the
initial water was permitted to dissolve but not precipitate, a
supersaturated mineral was permitted to precipitate but not dissolve, and
minerals close to saturation equilibrium were bounded near zero net
transfer. Redox and denitrification vectors were bounded to respect electron
and mass balance and to prevent the solver from selecting a thermodynamically
infeasible direction purely because it improved the residual. Ion-exchange
reactions were represented as signed pairs so that the two directions of a
given exchange process occupied opposite signs of the same coefficient
rather than being encoded as two independent non-negative variables; this
choice matters for identifiability, because two oppositely signed exchange
reactions are, by construction, exactly collinear and must be reported as a
single equivalence class rather than as two separately resolvable phases.
The sign convention fixed dissolution and oxidative release as positive and
precipitation, sorption, or reductive uptake as negative, applied uniformly
across the mineral, exchange, and redox vectors so that a reported negative
extent always had the same qualitative meaning regardless of reaction family.
Bound violations at the unconstrained coordinate-descent solution were
resolved by projection rather than by re-solving a constrained quadratic
program at every iteration, which kept the per-scenario solve time low
enough for the full factorial design while preserving exact feasibility of
the reported solution.

The six compared methods used bounds inconsistently by design, because the
ablation compared the effect of thermodynamic screening rather than
assuming it. Bounded least squares and the Hydrosheaf guarded and
Hydrosheaf-Core models applied the full thermodynamic bound set described
above; plain lasso and elastic net were fitted without saturation-index
bounds so that their unconstrained equifinality could be measured directly;
and the thermodynamic elastic net applied bounds derived only from
saturation index sign, without the additional geologic penalty priors used
by the guarded model. The conventional PHREEQC inverse baseline used its own
native uncertainty-bound mechanism (Supplementary Methods S6) rather than the
coordinate-descent bound-projection scheme, so its feasibility criterion is
not numerically identical to the guarded models' feasibility criterion; both
were reported on their own terms rather than forced onto a shared numeric
scale, and the comparison in the main text is limited to phase and
equivalence-class recovery, not to bound tightness. Saturation-index
thresholds separating "dissolving", "near equilibrium", and "precipitating"
states were fixed before any scenario was generated and were identical
across all four hydrochemical archetypes, so that archetype-specific
differences in recovery could not be attributed to an archetype-specific
threshold choice.

## S3. Rank, null-space, singular-value, condition-number, mutual-coherence, and reaction-equivalence-class diagnostics

For each ion panel, the corresponding columns of the design matrix
$\mathbf{A}$ were extracted and its singular value decomposition
$\mathbf{A} = \mathbf{U}\boldsymbol{\Sigma}\mathbf{V}^{\mathsf T}$ was
computed. Numerical rank was defined as the number of singular values
exceeding a fixed relative tolerance, and nullity as the number of reaction
columns minus the rank. A non-zero null space indicates the existence of at
least one non-trivial reaction combination $\mathbf{v}$ such that
$\mathbf{A}\mathbf{v} = \mathbf{0}$: adding any multiple of $\mathbf{v}$ to a
feasible reaction-extent vector leaves the fitted concentration change
unchanged, so recovering the true $\mathbf{v}$-component of the extent
vector from concentration data alone is impossible in principle, independent
of solver quality. The condition number, the ratio of the largest to
smallest non-zero singular value, summarised how close near-null directions
came to true collinearity, and mutual coherence, the largest normalised
inner product between any two reaction columns, summarised pairwise
near-collinearity that can degrade sparse recovery even when the matrix is
formally full rank. Exact reaction-equivalence classes were constructed by
identifying reaction subsets whose columns were linearly dependent to within
numerical tolerance and grouping them so that any feasible reallocation of
extent within the subset left the fitted concentrations unchanged. Approximate
equivalence classes extended this construction to reaction pairs or subsets
whose coherence exceeded a predeclared threshold, capturing near-collinear
cases in which the data technically permit discrimination but only at noise
levels well below those realised in practice. Every reported recovery metric
was computed against both the single true reaction and its equivalence
class, and the class definition itself, including its null-space basis, was
exported so that the boundary between formally provable non-identifiability
and merely practical ambiguity remained inspectable rather than asserted.

The full 11-ion, 16-reaction dictionary has rank 8 and nullity 4 under the
full ion panel, so at least four independent reaction directions can be
reallocated among themselves without changing any predicted concentration,
regardless of noise level, solver choice, or regularisation strength; this
property is a structural fact about the dictionary and the ion panel, not an
estimation error, and it sets a hard ceiling on exact-phase recoverability
that no comparator method can exceed under that panel. The dominant
contributor to this nullity is a pair of oppositely signed ion-exchange
reactions that are exact negatives of one another in the design matrix,
together with a smaller number of near-collinear mineral pairs whose
coherence exceeds the approximate-class threshold under representative noise
levels. Reduced ion panels (Supplementary Methods S5) generally increase
both nullity and coherence relative to the full panel, because removing a
diagnostic ion removes constraints on the reaction space rather than adding
any; the panel-by-panel rank, nullity, condition number, and coherence
values are reported in full for every combination of archetype and panel
(Supplementary Table S2, Supplementary Figure S2) so that the identifiability
consequence of a given monitoring decision can be read directly from the
diagnostic table rather than inferred from downstream recovery metrics
alone.

## S4. Stability-selection and Mechanism Resolution Score algorithms

Support stability was assessed by bootstrap resampling of the observed ion
vector under its assumed noise model, refitting the guarded model at each
resample, and recording the frequency with which each reaction was selected
across resamples and across a predeclared grid of penalty values, following
the stability-selection principle that a variable's importance should be
judged by how consistently it is selected across a range of regularisation
strengths and data perturbations rather than at a single, potentially
favourable, operating point [@meinshausen2010stability]. Reactions selected
with high and stable frequency across both axes were treated as robustly
supported; reactions whose selection frequency was sensitive to penalty
choice or resample were treated as candidates belonging to an ambiguous
equivalence class even if they appeared in the single best-fit solution.

The Mechanism Resolution Score (MRS) combined six predeclared diagnostics
per scenario: panel-specific rank deficiency, nullity, mutual coherence,
bootstrap support stability, held-out-ion prediction error, and
thermodynamic-feasibility violation rate. Each diagnostic was rescaled onto
a common range using statistics computed only on the training archetypes
(carbonate, crystalline, evaporitic), and the rescaled diagnostics were
combined by a calibrated classifier trained to predict the ordinal
resolution class (non-identifiable, partially identifiable,
equivalence-class identifiable, identifiable) defined from the known
ground-truth recovery outcome on those same training archetypes. Calibration
was assessed by reliability curves comparing predicted and observed class
frequencies, and discrimination was assessed by ordinal rank-correlation
between the predicted score and the realised recovery outcome. The mixed
archetype was withheld entirely from calibration and used only once, as a
transfer test, to obtain an untouched estimate of four-class classification
accuracy; this test was run a single time per model configuration to avoid
implicit tuning against the nominally held-out set. Threshold values
separating the four resolution classes were fixed from the training
archetypes before the mixed-archetype evaluation and were not adjusted
afterward.

Bootstrap resamples were drawn by perturbing the observed ion vector under
its assumed heteroscedastic noise model and refitting at each of a
predeclared grid of penalty values, so that the resulting stability path
recorded, for every reaction, its selection frequency jointly across the
resampling and regularisation axes rather than at one arbitrarily chosen
operating point; a reaction with high, flat selection frequency across the
grid was treated as robustly identified, while a reaction whose selection
frequency changed sharply with a small change in penalty was treated as
sitting on an unstable decision boundary regardless of its status in the
single reported best-fit solution. The four ordinal resolution classes used
by the Mechanism Resolution Score were defined operationally from the
training-archetype ground truth: "identifiable" scenarios were those in
which the exact true reaction set was recovered with high bootstrap support
and low held-out error; "equivalence-class identifiable" scenarios were
those in which the true reaction was recovered only as part of a multi-member
equivalence class but that class itself was recovered reliably;
"partially identifiable" scenarios recovered a strict subset of the true
active reactions with unstable support for the remainder; and
"non-identifiable" scenarios showed high held-out-ion error or unstable
support across the majority of the true active set. The classifier combining
the six rescaled diagnostics into a single ordinal prediction was selected
and its hyperparameters were fixed using only the three training archetypes,
with the mixed archetype excluded from every stage of that selection,
including any implicit model-family comparison, so that the reported
four-class transfer accuracy reflects genuine extrapolation to an
untouched hydrochemical setting rather than in-sample fit.

## S5. Held-out-ion falsification and next-best-measurement algorithm

Held-out-ion falsification withheld one ion at a time from the fitting
objective by removing its row from $\mathbf{W}$ and $\mathbf{r}$, refit the
reaction-extent vector on the remaining ions, and then used the fitted
stoichiometry to predict the withheld ion's concentration change. The
prediction error on the withheld ion, rather than the residual on the ions
used for fitting, was treated as the primary falsification signal, because a
reaction combination that fits the retained ions well but predicts a
withheld ion poorly demonstrates that the retained ions under-determine the
true mechanism. This differs from a conventional missing-ion sensitivity
test, in which an omitted ion's weight is simply set to zero and the
omitted ion is never scored, a procedure that can appear favourable purely
because the omitted ion was removed from the objective rather than because
the fitted reaction actually predicts it.

Next-best-measurement selection used the same held-out-ion mechanism in
reverse. For a given scenario and current ion panel, each candidate
additional ion was hidden from the current fit, its expected value under
the current best-supported equivalence class was compared against its
expected value under competing class members, and candidates were ranked by
the expected reduction in class ambiguity that observing them would produce.
The realised value of a candidate measurement was then computed by actually
adding it to the panel, refitting, and recording the observed reduction in
equivalence-class ambiguity and held-out-ion prediction error. Comparing
expected and realised value tested whether the ranking procedure was itself
informative rather than only self-consistent, and comparing the
identifiability-aware next-best selection against fixed and randomly chosen
panels of the same size tested whether the procedure outperformed
convention-based monitoring design. The data-tier experiment (core,
core-plus-lite, and enhanced diagnostic panels) applied the same
held-out/realised-value logic to controlled synthetic SiO2, strontium, water
isotope, bromide, dissolved oxygen, dissolved organic carbon, and
sulphate/nitrate isotope channels, each generated from the known scenario
truth with an added measurement-noise model, and reported an entropy-based
evidence-lifted resolution index that scored how far conditioning on the
additional diagnostic separated the evidence supporting each member of an
equivalence class, without asserting that the added diagnostics created new
mass-balance uniqueness.

The evidence-lifted resolution index (ELRI) was constructed as an entropy
statistic over the posterior-like evidence weights assigned to the members
of an equivalence class after conditioning on the available diagnostics: a
value of zero corresponds to a uniform distribution of evidence across class
members, meaning the additional diagnostics provided no discriminating
information, and a value approaching one corresponds to one class member
carrying essentially all of the evidence weight. ELRI was reported
separately from the exact- and equivalence-class recovery metrics because it
answers a different question: recovery metrics ask whether the correct
class was found at all, while ELRI asks, conditional on having found the
correct class, how far the available evidence narrows the choice within it.
Reporting both together prevented an improvement in ELRI from being
mistaken for an improvement in mass-balance identifiability itself, since
the underlying null space of the reaction dictionary was unchanged by adding
diagnostic ions; only the practical evidence available to prefer one member
of a class over another changed. The fluoride channel, evaluated within the
Northern Ghana transfer analysis, obtained the highest retrospective
measurement-value score of the tested ions, reflecting its diagnostic
specificity for fluorite/apatite dissolution reactions that are otherwise
weakly constrained by the major-ion panel alone; this result is reported as
a retrospective, dataset-specific finding rather than a general claim that
fluoride is the most informative measurement in every hydrochemical setting.

## S6. PHREEQC synthetic-data generator and species extraction

The synthetic benchmark used PHREEQC 3.7.3 [@parkhurst2013phreeqc] to
generate scenario ground truth by explicit forward simulation rather than by
a mass-balance proxy. For each of the four hydrochemical archetypes
(carbonate, crystalline/silicate, evaporitic, and mixed), 60 replicate
scenarios were generated by sampling initial water chemistry (pH,
temperature, ionic strength, and initial saturation state around an
archetype-specific centre), a reaction set of one to five phases drawn from
the 16-reaction dictionary, a signed extent for each active phase, and,
where applicable, a mixing or evaporation transport confounder. Each
scenario was executed as a PHREEQC batch-reaction run from the sampled
initial solution through the sampled reaction set to a simulated
"downstream" equilibrium composition, and the resulting major-ion
concentrations, pH, and saturation indices were extracted programmatically
from the PHREEQC output using a fixed species-mapping table so that solver
output could be joined unambiguously to the corresponding synthetic ion
vector used by the inverse-modelling pipeline. Analytical noise was then
added to the extracted concentrations at each of the three noise levels
(0%, 3%, 8%) using a heteroscedastic noise model in which the standard
deviation scaled with the concentration magnitude subject to a floor, so
that trace species were not assigned unrealistically small absolute
uncertainties. Random seeds were derived deterministically from the
scenario index and noise level so that every noise realisation was
reproducible from the archived seed table alone. Numerical validity was
checked before any inverse fit was attempted: the PHREEQC-simulated
concentration change was compared against the concentration change implied
by the sampled stoichiometry and extents, and scenarios whose generation
residual exceeded a strict tolerance were flagged and excluded from the
ground-truth set, giving a maximum PHREEQC-to-stoichiometric generation
residual across all retained scenarios that was several orders of magnitude
smaller than the smallest analytical-noise level applied. All PHREEQC input
decks, raw output listings, and the extracted scenario-level ground-truth
table were preserved so that any single scenario's provenance could be
re-examined independently of the aggregated results.

Each PHREEQC input deck comprised a `SOLUTION` block encoding the sampled
initial water chemistry, one or more `EQUILIBRIUM_PHASES` or `REACTION`
blocks encoding the sampled active phases and their signed extents, and, for
mixing-confounded scenarios, a `MIX` block combining the reacted solution
with a second end-member solution at a sampled mixing fraction. Custom phase
definitions for reactions without a native PHREEQC database entry (for
example, the signed exchange and denitrification vectors) were declared
explicitly in the input deck with fixed equilibrium constants consistent
with the reaction dictionary used by the inverse-modelling side of the
pipeline, so that the same stoichiometry governed both forward generation
and inverse recovery. A single fixed thermodynamic database was used across
all scenarios and archetypes to avoid conflating database-driven differences
in equilibrium constants with archetype-driven differences in recoverability;
database sensitivity was treated as a separate, secondary experiment
(Supplementary Figure S10) rather than folded into the primary factorial
design. Species concentrations, pH, and saturation indices were parsed from
the PHREEQC output listing using a fixed field-mapping table keyed on
PHREEQC's internal species and phase identifiers, and every parsed value was
cross-checked against the corresponding line in the raw output listing
during generation so that a parsing error would surface as a scenario-level
generation-residual failure rather than propagate silently into the
downstream ground-truth table.

## S7. Northern Ghana preprocessing, unit conversion, and 2025 sampling-date verification

The Northern Ghana dataset comprised 320 wet/dry hydrochemical records from
160 boreholes across four administrative regions, with major-ion chemistry,
water stable isotopes, strontium, and silica available from the raw
workbook (`data/FieldData/NorthenGhana/NorthernGhana.xlsx`, `Dry`/`Wet`
sheets). This raw workbook carries no independent aquifer-type,
geology-group, or lithology classification and no precomputed saturation
indices, so stratified reporting uses region and district rather than
aquifer or lithology, and saturation indices for calcite, dolomite, gypsum,
halite, and fluorite were computed independently by live PHREEQC 3.8.6 from
each borehole's own wet-season pH, temperature, and major-ion panel (no Al
or PO4 was measured, so silicate and apatite saturation states were not
computable). All concentration fields were converted to a single consistent
unit system (millimoles per litre for the species entering the reaction
dictionary) before any fitting, and charge-balance error was computed for
every record from the converted major ions. Records were paired by
borehole, retaining boreholes with a wet- and dry-season record and at
least six commonly measured ions; every one of the 160 boreholes met this
criterion, so no separate quantitative-subset filter was applied. Exact
per-record sampling dates are not present in the raw workbook; the
dry-season campaign (5-24 March 2025) and wet-season campaign (5-24 August
2025) covering these same 160 boreholes are documented in the accompanying
`SI.pdf` (Supplementary Table 19), which was consulted for this date range
only, not for any borehole-level attribute. Graph edges connecting
boreholes were generated by the Hydrosheaf topology layer from spatial and
hydraulic proximity rather than observed directly in the field, and
heuristic hydraulic residence times attached to those edges were derived
from a hydraulic heuristic rather than measured with an age tracer; both are
reported as processed analytical outputs and are explicitly excluded from
any claim of independently observed flow-path or groundwater-age ground
truth. Region, district, and season labels used for stratified reporting
were taken directly from the raw workbook without modification.

Two provenance issues identified during the pre-manuscript data-readiness
audit were carried forward explicitly into the reporting boundary rather
than resolved by assumption. First, all 1,019 candidate Hydrosheaf-generated
edges over the corrected workbook carried `phreeqc_ok = False` and
`rt_validated = False` flags and no populated edge-confidence field at the
time of the field-transfer analysis, so no edge-level topology or
residence-time weighting was applied to the chemistry-only results reported
here; where a figure or table shows edge-stratified results, the
stratification is by region, district, or season metadata, not by
topology confidence. Second, hydraulic residence-time estimates attached to
those edges spanned a median of approximately 1.74 million days and a
maximum of approximately 3.14 billion days, values consistent with an
unvalidated hydraulic heuristic rather than a measured groundwater age, and
were excluded from every age-conditioned analysis. The external field-
transfer extension repeated the ELRI audit described in Supplementary
Methods S5 on three independently sourced chemistry datasets covering the
same general region (`NorthernGhana.xlsx`, 160 edges; Talensi, 85 edges;
Lower Anayari, 49 edges), constructed using wet-dry or nearest-neighbour
pairing rules rather than an independently observed flow path, and is
reported strictly as a plausibility and transfer audit rather than as
additional reaction-truth validation. Talensi and Lower Anayari measure
complete pH, temperature and major-ion panels, so saturation indices were
computed independently for them by live PHREEQC exactly as for
`NorthernGhana.xlsx`, rather than left uncomputed on the assumption that
saturation indices required strontium or silica (neither of which this
computation uses); Talensi measures no fluoride, so its fluorite index is
not computable and is treated as absent rather than estimated
(`DECISIONS.md`).

## S8. Complete reproducibility protocol

The full analysis, table, and figure workflow was executed from a single
entry-point script that ran the live-PHREEQC factorial benchmark, the
factorial inversion across all methods, noise levels, and ion panels, the
Mechanism Resolution Score calibration and transfer test, the held-out-ion
and next-best-measurement analyses, the bootstrap stability analysis, and
the Northern Ghana field-transfer analysis in a fixed order with archived
per-stage logs. The complete set of generated results tables was exported to
a single evidence database together with a table catalogue, and a
run-specific artefact manifest recorded every figure, table, and generating
script produced by that run. Software versions, the environment lock file,
and the PHREEQC executable and database version were recorded for every run.
Two reduced-cost reproduction modes were retained alongside the full run: a
display-only mode that regenerated tables and figures from already-generated
result files without repeating the factorial inversion, and a
model-comparison mode that reused existing PHREEQC ground truth while
regenerating only the sparse-inversion outputs, so that a reviewer without
access to a licensed PHREEQC installation could still reproduce the
downstream statistical and reporting pipeline from the archived ground-truth
tables. Code, the environment lock file, generated inputs and outputs,
figure-generation scripts, and a model-card-style statement of claim
boundaries were deposited together in a permanent repository at
submission.

A run ledger recorded, for every executed run, its status, start and end
timestamps, the corresponding source-code revision, a hash of the resolved
configuration, a hash of the input manifest, and a hash of the resolved
software environment, so that any reported figure or table could be traced
back to the specific run that produced it rather than to the pipeline in the
abstract. The complete set of generated result tables was consolidated into
a single evidence database with an accompanying table catalogue describing
every table's columns, units, and generating stage, and a run-specific
artefact manifest linked every main and supplementary figure and table to
its underlying result table and the plotting or table-generation script that
produced it. This manifest, rather than narrative cross-referencing, is the
authoritative map from a displayed number in the manuscript to the
computation that generated it, and it was regenerated automatically at the
end of every full run so that it could not drift out of step with the
displays it describes. The claim-boundary statement accompanying the
deposited materials restates, in machine-readable form, the distinctions
enforced throughout this Supplementary Methods document: which results are
synthetic-benchmark evidence with known ground truth, which are
retrospective measurement-value simulations on controlled synthetic
diagnostics, and which are chemistry-only field-transfer results without
independent validation of connectivity, groundwater age, or reaction truth.
