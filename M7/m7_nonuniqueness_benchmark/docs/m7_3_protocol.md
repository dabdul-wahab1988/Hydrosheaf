# M7.3 conditional-integration benchmark — locked protocol

## Thesis question

M7.3 asks when combining residence-time, connectivity and hydrochemical
evidence reduces interpretive non-uniqueness, when it is redundant, and when a
misspecified evidence stream causes negative transfer.

The benchmark does not define success as universal improvement from adding age
or any other evidence. Reduced uncertainty is counted as beneficial only when
calibration or predictive performance is not materially degraded.

## Independence and analysis split

The external generator is unchanged from M7.2 and imports no HydroSheaf code.
It executes official MODFLOW 6 and MODPATH 7 binaries and separately generates
nonlinear chemistry and tracer observations.

- Development seeds: 5201–5206.
- Locked test seeds: 5301–5312.
- Non-claim-bearing smoke seeds: development 5901–5902 and test 5911–5913.
- Development cases fit all evidence-fusion coefficients.
- Test truth is joined only after truth-blind inference.
- No test result may change a feature, coefficient, condition, threshold or
  stopping rule.

## Experiment 1: evidence integration

Three evidence streams are evaluated:

- `H`: hydraulic log odds;
- `A`: uncertainty-aware age compatibility; and
- `C`: negative constrained-chemistry log objective.

Unweighted logistic models are fitted on development cases for `H`, `A`, `C`,
`HA`, `HC`, `AC` and `HAC`. Unweighted fitting is deliberate because posterior
entropy and Brier/log-loss calibration would be distorted by artificial class
balancing.

The locked test conditions are:

1. `native`: all streams retain their edge association;
2. `age_permuted`: age evidence is permuted within each independent case;
3. `hydraulic_permuted`: hydraulic evidence is permuted within each case; and
4. `joint_misspecified`: both age and hydraulic evidence are independently
   permuted within each case.

Permutation preserves each case's marginal evidence distribution while
destroying its edge-specific meaning. It is a negative control, not a model of
one particular field failure.

Reported metrics are PR-AUC, ROC-AUC, Brier score, log loss, normalized
edgewise Bernoulli entropy, the fraction of edges with posterior probability
between 0.1 and 0.9, expected edge count, calibration error and overconfident
error. Case-block bootstrap intervals use 10,000 replicates.

An entropy reduction is supportable only if the paired mean Brier score and log
loss do not worsen. Confidently wrong predictions are negative transfer, not
resolution.

## Experiment 2: topology assumptions and groundwater age

The same local 3H/39Ar likelihood and lognormal age prior are evaluated under:

- no connectivity assumption;
- a 50% partial true graph;
- the complete externally generated graph; and
- the completely reversed graph.

Two tracer regimes are predeclared:

- `informative`: 3H and 39Ar;
- `tritium_only`: 3H without 39Ar.

Independent local age-posterior particles are reweighted by a soft
downstream-older topology potential with a five-year scale. The topology
potential is an explicit assumption. Its importance-sampling effective sample
size and log mean weight are reported so incompatible graphs cannot silently
produce apparently precise results.

Each full case/regime analysis uses 50,000 common local-posterior particles.

Endpoints are age MAE, bias, 95% coverage, interval width, normalized marginal
age entropy, true-edge ordering violations and importance ESS. Correct versus
no-topology and reversed versus correct contrasts are paired by case.

## Experiment 3: reaction non-uniqueness

Only externally true flow edges are used for reaction-mechanism scoring after
inference. Chemistry is perturbed by 3% lognormal analytical error for 64
bootstrap replicates per case.

Two panels are compared:

- `core`: Ca, Mg, Na, K, HCO3, Cl, SO4 and NO3;
- `enhanced`: core plus F, Fe, PO4 and SiO2.

Both retain PHREEQC direction bounds derived from the available full
speciation/saturation calculation. Outputs are modal family accuracy,
probability assigned to the true family, family-support entropy and effective
number of supported families. Carbonate weathering and precipitation are
reported separately and may not be hidden inside aggregate accuracy.

## Experiment 4: Ghana data-scope audit

The field audit is descriptive, not a truth benchmark. The workbook is checked
for measured chemistry, environmental age tracers, screen intervals,
coordinates, hydraulic observations and independent connectivity labels.

Stable water isotopes are recharge/source tracers, not residence-time tracers.
Elevation minus one static-water-level measurement is a single-occasion head
proxy, not a repeated head field. Masked coordinates and processed graph edges
are not independent site-scale flow truth.

The permitted objective is:

> Apply the framework and its component diagnostics to Ghanaian aquifer
> datasets to determine which integrated interpretations are supportable under
> available data and which remain non-identifiable.

## Decision rules and guardrails

- Evidence complementarity requires improved case-level PR-AUC or log loss with
  a 95% interval excluding zero, or reduced entropy without worsened Brier and
  log loss.
- Redundancy means the added stream changes neither skill nor calibrated
  uncertainty materially.
- Negative transfer means skill/calibration worsens or uncertainty falls while
  overconfident error rises.
- Synthetic truth remains model-conditioned truth.
- Ghana data do not validate age, exact topology or reaction truth.
- M7.3 cannot be described as preregistered; the protocol is Git-locked before
  the declared 5301-series results are generated.
