# Methods

## Conceptual model and mass-balance formulation

The basic inferential unit was a directed Hydrosheaf edge linking an
upstream and a downstream groundwater composition. Where available, each
edge carried topology confidence from the graph/sheaf layer and a
residence-time estimate from the temporal layer; where these upstream
quantities were unavailable, the analysis was labelled explicitly as
chemistry-only rather than fully integrated Hydrosheaf inference. The
transport-corrected concentration residual was defined as

<!-- EQ:EQ-1 -->
$$
\mathbf{r} = \mathbf{x}_v - \mathbf{x}_{\mathrm{transport}},
$$

where $\mathbf{x}_v$ was the observed downstream ion-concentration vector and
$\mathbf{x}_{\mathrm{transport}}$ was the concentration expected from
conservative transport alone. Reaction extents $\mathbf{z}$ were estimated by
solving a bounded, penalised sparse regression

<!-- EQ:EQ-2 -->
$$
\hat{\mathbf{z}}
=
\arg\min_{\mathbf{z}}
\left[
\frac{1}{2}
\|\mathbf{W}^{1/2}(\mathbf{r}-\mathbf{S}^{\mathsf T}\mathbf{z})\|_2^2
+\lambda_1\sum_j p_j|z_j|
+\frac{\lambda_2}{2}\|\mathbf{z}\|_2^2
\right],
$$

subject to thermodynamic lower and upper bounds on each $z_j$, where
$\mathbf{S}$ was the signed stoichiometric reaction matrix, $\mathbf{W}$ was
a diagonal ion-specific analytical-uncertainty weight built from the assumed
measurement precision of each species, $p_j$ was a per-reaction penalty
scale that down-weighted geologically implausible reactions relative to the
locally expected assemblage, and $\lambda_1$, $\lambda_2$ were the lasso and
ridge regularisation weights of an elastic-net-type penalty
[@tibshirani1996lasso; @zou2005elasticnet]. Dissolution and precipitation
corresponded to opposite signs of $z_j$ under a single locked sign
convention applied consistently across every reaction in the dictionary, so
that a negative extent always denoted removal of the associated ions from
solution regardless of whether the underlying process was mineral
precipitation, sorption, or reductive uptake. The objective, penalty, and
bound scaling were fixed to one internally consistent convention before any
scenario was fitted, and the coordinate-descent convergence tolerance was
defined on that same scale so that stopping behaviour did not depend on an
arbitrary unit choice. Full symbol definitions, units, and the
coordinate-descent solution to Equation 2 are given in Supplementary
Methods S1-S2.

## Reaction dictionary, solver, and identifiability diagnostics

The stoichiometric reaction matrix comprised 16 signed reactions over an
11-ion panel, spanning mineral dissolution/precipitation, cation exchange,
nitrate-source, and denitrification vectors. Coefficients were estimated by
coordinate descent with soft thresholding, a fixed iteration budget, and
convergence and bound-feasibility checks (Supplementary Methods S1-S2).
Because several reaction directions are stoichiometrically collinear or
exactly opposed (for example, competing ion-exchange directions), the
dictionary was not assumed to be identifiable a priori. Rank, nullity,
singular values, condition number, and mutual coherence of the
panel-specific design matrix were computed for every ion panel to diagnose
identifiability directly, following the general principle that residual fit
alone cannot distinguish structurally equivalent solutions
[@christophersen1992emma; @beven2001glue]. Reactions or reaction combinations
that were linear-algebraically indistinguishable under a given ion panel
were grouped into exact reaction-equivalence classes; near-collinear
combinations were grouped into approximate classes using a coherence
threshold (Supplementary Methods S3). This distinction separated two
different failure modes that a single residual-based metric conflates: a
solver failing to find any adequate reaction combination, and a solver
correctly finding an adequate combination that the observed ion panel simply
cannot distinguish from an equally consistent alternative. Recovery was
therefore always reported at two levels: exact-phase recovery against the
single true reaction, and equivalence-class recovery against the class
containing the true reaction, so that an apparent failure at the exact-phase
level could be re-evaluated as an information-limited but structurally
correct result at the equivalence-class level.

## PHREEQC synthetic benchmark and comparators

Ground truth was generated with a live PHREEQC 3.7.3 factorial benchmark
[@parkhurst2013phreeqc] of 240 scenarios (60 replicates across four
hydrochemical archetypes: carbonate, crystalline/silicate, evaporitic, and a
mixed archetype combining carbonate and silicate weathering with cation
exchange in a single water), each with one to five simultaneous known
phases, dissolution or precipitation direction, and a reaction extent within
a realistic groundwater range. Scenarios were run at three analytical-noise
levels (0%, 3%, 8% of concentration) and under five ion panels ranging from
the full 11-ion set to a reduced 5-ion core panel, factorially combined with
mixing/evaporation transport confounders, giving 21,600 factorial inverse
fits. Seven methods were compared: bounded least squares, lasso, elastic
net, a thermodynamic elastic net, the Hydrosheaf Guarded model
(thermodynamically bounded, penalty-scaled elastic net), a Hydrosheaf-Core
evidence-gated variant that adds sparse geochemical evidence gates (Gibbs,
exchange, redox, charge-balance, and saturation-index checks) without
requiring rare tracers, and a conventional PHREEQC `INVERSE_MODELING`
baseline [@parkhurst2013phreeqc] run under strict 5% and relaxed 20%
analytical-uncertainty bounds. Hydrosheaf Guarded used the identical
coordinate-descent solve as the thermodynamic elastic net (same fit method
and penalty weights) and differed from it only in the support threshold
applied to convert continuous reaction extents into a selected phase set;
consequently, continuous quantities such as extent RMSE, reconstruction
RMSE, and runtime are identical between the two methods by construction, and
only phase-support and equivalence-class metrics, which depend on the
threshold, differ between them. The strict bound followed conventional
inverse-modelling practice; the relaxed bound represented the wider
tolerance sometimes adopted when a strict solution is infeasible, so both
feasibility and recovery accuracy could be reported as a function of the
declared uncertainty. Thermodynamic bounds, geologic penalty priors,
indicator-ion gates, and the L1/L2 penalty weights were each ablated
individually to isolate their marginal contribution to phase recovery,
false-discovery rate, and extent error. Regularisation was selected by a
predeclared grid with a one-standard-error rule, which selected the most
parsimonious penalty whose cross-validated performance was within one
standard error of the minimum-residual penalty, rather than by minimum
residual alone, because a minimum-residual choice systematically favours the
member of an equivalence class with the smallest incidental numerical
advantage rather than the member best supported by independent evidence
(Supplementary Methods S4). The conventional PHREEQC inverse baseline was
run without a sparsity penalty, so it reported every phase combination
consistent with the declared uncertainty bounds rather than a single
preferred solution; its full feasible-model count and its first minimal
model were both retained so that solver multiplicity itself could be
compared against the guarded models' equivalence-class reporting.

## Evaluation metrics and the Mechanism Resolution Score

Primary metrics were phase-support precision, recall, F1, and false-discovery
rate; equivalence-class precision, recall, and coverage; signed
reaction-direction accuracy; extent RMSE, MAE, and interval coverage;
concentration-reconstruction RMSE; held-out-ion RMSE and interval coverage;
thermodynamic-violation rate; bootstrap support stability across resamples
and penalty values [@meinshausen2010stability]; and convergence/runtime.
Held-out-ion prediction withheld a predeclared ion from the fitting objective
and tested whether the inferred reaction predicted its concentration, in
place of the weaker convention of merely zero-weighting an omitted ion.
Whether residual error predicts support accuracy was tested with regression,
calibration curves, and rank correlation. A calibrated Mechanism Resolution
Score (MRS), combining rank, nullity, coherence, support stability, held-out
error, and thermodynamic feasibility, was fitted on three training
archetypes (carbonate, crystalline, evaporitic) and evaluated, untouched,
on the mixed archetype, reporting four ordered resolution classes:
non-identifiable, partially identifiable, equivalence-class identifiable,
and identifiable (Supplementary Methods S4). Held-out-archetype testing was
used specifically because within-archetype cross-validation would let the
score learn archetype-specific shortcuts rather than a transferable relation
between structural diagnostics and phase-recovery reliability. Next-best
measurement was evaluated retrospectively by masking observed variables,
ranking the expected ambiguity reduction from each candidate ion, and
comparing the realised reduction in held-out prediction error against fixed
and random measurement panels, following the general principle that a
monitoring design should be scored by its measured value of information
rather than by convention or availability alone
[@sreekanth2017monitoring] (Supplementary Methods S5). An analogous
data-tier experiment compared a core major-ion panel against panels enriched
with controlled synthetic SiO2, Sr, and water-isotope diagnostics, and
further with Br, dissolved oxygen, dissolved organic carbon, and
sulphate/nitrate-isotope diagnostics, reporting an entropy-based
evidence-lifted resolution index that quantifies conditional evidence
separation within an equivalence class rather than new mass-balance
uniqueness (Supplementary Methods S5).

## Northern Ghana chemistry-only field demonstration

The framework was applied to 320 wet/dry hydrochemical records from 160
Northern Ghana boreholes sampled in March and August 2025, spanning four
administrative regions in terrain dominated by crystalline basement and
Voltaian sedimentary geology in which silicate weathering and cation
exchange are established regional controls [@anku2008ghana;
@banoengyakubo2011ghana]. No aquifer-type or lithology classification is
available, so results are reported by region. Records were paired by
borehole, retaining boreholes with both a wet-season and a dry-season record
and at least six ions measured in common between the two records; every one
of the 160 boreholes met this criterion, yielding 160 wet-to-dry pairs.
Support stability and held-out-ion performance were stratified by region
and season, and reaction-equivalence classes were reported wherever the
measured ions did not distinguish individual phases (Supplementary Methods
S7). Saturation indices for calcite, dolomite, gypsum, halite, and fluorite
were computed by live PHREEQC from each borehole's wet-season chemistry (no
Al/PO4, so silicate/apatite SI were not computable), giving the
thermodynamic/SI evidence gates real field-specific constraints. The
identifiability-aware output was compared directly against a conventional
single best-fit interpretation of the same field pairs so that the practical
difference between the two reporting styles could be shown on the same
data. Generated graph edges and heuristic hydraulic residence times were not
treated as independent ground truth; this component is reported as a
chemistry-only transfer demonstration, not as validation of connectivity,
groundwater age, or reaction mechanism.

## Reproducibility

Code, the environment lock file, the PHREEQC database, generated scenario
inputs and outputs, random seeds, and figure-generation scripts are archived
at `https://github.com/dabdul-wahab1988/Hydrosheaf` and will be deposited,
with data and claim limitations, in a permanent repository before submission
(Supplementary Methods S8). Software versions, fixed random seeds, and
predeclared hyperparameter grids were recorded for every stage of the
pipeline, and scenario-level PHREEQC input and output files were preserved
unmodified so any ground-truth value could be traced to the simulation that
produced it. The distinction between synthetic benchmark
evidence, retrospective measurement-value simulation, and the Northern Ghana
field transfer was preserved throughout so no synthetic result could be
misread as field validation.
