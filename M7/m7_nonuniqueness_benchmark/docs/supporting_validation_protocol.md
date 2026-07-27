# Supporting strong-integration validation — Locked protocol

## Question and estimands

M7.2 tests whether HydroSheaf integration recovers externally generated flow
edges and reaction families, whether Bayesian ages are calibrated and
incrementally useful, whether PHREEQC constraints materially alter the inverse
problem, and whether topology uncertainty is numerically converged.

The primary topology estimands on the locked cases are PR-AUC, ROC-AUC, Brier
score, precision, recall, F1, and MCC for:

- hydraulics plus constrained chemistry;
- hydraulics plus constrained chemistry plus age; and
- an adverse control in which age evidence is permuted within each aquifer.

Positive incremental age value is the paired locked-case difference between
the full and no-age models. Both F1 at the development-frozen threshold and
threshold-free PR-AUC are reported with a 10,000-replicate case-block
bootstrap. A positive F1 alone does not override a non-positive PR-AUC result.

Reaction-family recovery, Bayesian-age MAE/bias/coverage, topology-posterior
diagnostics, and PHREEQC constraint activity are secondary estimands.

## Independence and blinding

`independent_modflow_generator.py` imports no HydroSheaf module. It calls the
official MODFLOW 6 and MODPATH 7 executables through FloPy, then applies
standalone nonlinear reaction and tracer equations. The generator and
inference use different implementations and include model discrepancy.

HydroSheaf receives only the observation rows. Hidden age, edge, process, and
pathline fields are forbidden by the inference contract. Truth is joined only
by the benchmark scorer. The candidate graph must contain all true edges for
recovery to be estimable; candidate recall is reported rather than assumed.

Development seeds are 2101–2106. They alone fit the standardized,
class-balanced logistic fusion coefficients and all classification thresholds.
Locked test seeds are 3101–3112. No locked-test outcome may alter features,
coefficients, thresholds, priors, or stopping rules.

## External synthetic experiment

Each case is a heterogeneous one-layer 10 x 20 MODFLOW 6 aquifer with a
low-conductivity lens and three forward MODPATH 7 particles. Four milestones
per path yield 12 observed nodes, nine true directed edges, and exact
model-conditioned pathline ages.

Every case exercises carbonate weathering, carbonate precipitation, silicate
weathering, denitrification, sulfate reduction, and iron reduction. The
reducing sequence follows nitrate, sulfate, then ferric-iron electron-acceptor
progression. Laboratory/process noise, nonlinearity, and conservative-ion
perturbations prevent exact recovery by HydroSheaf's linear dictionary.

The tracer boundary forcing is treated as known, as it would be supplied by an
independent recharge-history measurement. The generator nevertheless uses
different nonlinear 3H/39Ar responses. A development-calibrated 39Ar reference
scale and Student-t discrepancy likelihood are frozen before locked testing.

## Bayesian age audit

The blind edge topology is not supplied to the age model. Node ages are
estimated independently by exact grid quadrature from 3H and 39Ar, which
separates dating performance from topology truth. Four independent posterior
draw chains must have:

- maximum rank-normalized split R-hat no greater than 1.01;
- minimum bulk and tail ESS at least 400; and
- zero divergences.

Sampler convergence and tracer identifiability are reported separately. Age
MAE, mean bias, and 95% interval coverage are scored against MODPATH only after
inference.

## Topology posterior audit

The posterior uses four dispersed feasible chains, a hard DAG constraint, a
minimum of nine edges, and maximum out-degree three. Its transition mixture
combines exact edge full-conditionals with symmetric parent swaps, general
swaps, two-edge flips, and single-edge flips. Sixteen transitions are made
between retained draws.

Convergence requires, for graph size and every edge indicator:

- rank-normalized split R-hat no greater than 1.01; and
- bulk and tail ESS at least 400.

No posterior topology claim is made for a case that fails these rules.

## Active geochemical constraints

HydroSheaf fits every candidate edge twice, with and without PHREEQC-derived
thermodynamic direction bounds. Bounds are constructed in the exact
sample/layer-specific reaction-label order. Precipitation-only phases can take
negative extents. A constraint is considered materially active when it changes
an edge objective by more than `1e-6`; success fractions, bound hits, and
objective changes are retained.

The reaction audit reports constrained and unconstrained dominant-family
accuracy for all six generator processes. PHREEQC activity is not, by itself,
evidence that family classification improved.

## Field prequential audit

The field component uses complete dry/wet quantitative pairs from the
canonical Northern Ghana workbook (`data/FieldData/NorthenGhana/
NorthernGhana.xlsx`, Dry/Wet sheets; an earlier revision read a different,
antecedent study's own derived workbook instead, since removed — see
`DECISIONS.md`). That workbook records one dry-season and one wet-season
observation per well and no intra-season sampling-date field, so a real
chronological issue-date sequence cannot be reconstructed. Wells are instead
revealed in a fixed, disclosed, seeded pseudo-random batch order (20
batches); for each batch, predictions are created before that batch's
wet-season observations are assimilated. Dry-season chemistry and static
well features are available for all wells throughout; wet-season labels from
not-yet-revealed batches are inaccessible.

Graph ridge is compared with persistence and an expanding mean-delta baseline.
Uncertainty is calibrated from prior residuals and paired uncertainty uses a
well-block bootstrap. Because coordinates are spatially masked, the data
contain one dry/wet campaign rather than independent years, and the batch
order is arbitrary rather than a real sampling sequence, the permitted claim
is only within-campaign one-step hydrochemistry hold-forward under a
disclosed arbitrary revelation order — not a temporal-sequence validation.

## Interpretation guardrails

- External truth is model-conditioned synthetic truth, not field truth.
- Field data do not validate topology, ages, or reactions.
- The development/test split prevents locked-test tuning but is not
  multi-basin external validation.
- A full M7.2 success strengthens an integration-method paper; it cannot alone
  guarantee Q1 acceptance.
