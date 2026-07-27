# M8 confirmatory protocol

Status: **prospective confirmation after exploratory analysis**. The 2026-07-27
exploratory CSV files and the observations that motivated this protocol have
already been inspected; they are not preregistered evidence and are excluded
from the confirmatory analysis.

## Scope and claim boundary

This benchmark evaluates controlled synthetic parameter recovery through
HydroSheaf's production adapters. It is not field validation, an operational
digital twin, or evidence that one design criterion is universally superior.

## Locked objectives

1. Estimate dispersivity and decay recovery separately across fixed transport
   sampling designs.
2. Test next-measurement choices targeted to dispersivity, decay, D-, A-, and
   E-optimality, and a balanced marginal-variance objective.
3. Verify that the production PHREEQC kinetic law identifies the product
   `k*A`, not `k` and `A` separately, and test an independent surface-area
   measurement as the identifying intervention.

## Transport data-generating model

- Production adapter: `TransportCalibrationAdapter`, analytical 1-D
  advection-dispersion mode.
- Truth: dispersivity = 2.0 m; first-order decay = 0.005 d^-1.
- Distance = 10 m; velocity = 0.1 m d^-1; source concentration = 1.0.
- Error model: independent Gaussian concentration error with
  `sigma_i = max(0.01, 0.03*abs(C_i))` source-concentration units.
- Fisher sensitivities are central finite differences with respect to
  `log10(parameter)` and are whitened by `sigma_i`.
- Primary finite-difference step = 1e-4 log10 units; 1e-3 and 1e-5 are numerical
  stability checks.
- Calibration uses the production `PESTGLM` engine in log-parameter space.
  Adapter-created priors are disabled so random starts do not become random
  regularisation centres.

## Replication and comparisons

- Fresh seed family: 2026072801; bootstrap seed: 2026072802.
- 250 paired replicates per fixed design and per OED strategy.
- Common standard-normal errors and common starts are used across designs or
  strategies within each replicate.
- Primary starts are within +/-0.25 log10 units of truth. A prespecified
  sensitivity analysis uses +/-1.0 log10 units.
- The fixed-design sweep is descriptive: a replicate-wise Spearman coefficient
  across deterministic designs is summarised over stochastic replicates. A
  classical p-value across the 16 hand-selected designs is prohibited.
- Primary recovery endpoint: absolute log10 parameter error, separately for
  dispersivity and decay.
- Secondary endpoints: relative error, linearised 95% coverage, interval width,
  convergence, objective value, marginal local variance, parameter correlation,
  and joint Fisher criteria.
- Uncertainty for median errors and paired median differences is a percentile
  bootstrap with 2,000 resamples.

## Transport decision rule

Full target-dependent trade-off support requires all of the following:

1. The locally predicted dispersivity- and decay-targeted candidate times differ.
2. Dispersivity-targeted sampling has lower realised dispersivity error than
   decay-targeted sampling, with the paired 95% bootstrap interval for the
   median difference excluding zero.
3. Decay-targeted sampling has lower realised decay error than
   dispersivity-targeted sampling, with the paired 95% bootstrap interval for
   the median difference excluding zero.

One successful directional comparison is labelled partial support. Otherwise
the trade-off is not confirmed. D-, A-, and E-optimality are reported as joint
criteria; they are not required to optimise either named parameter.

## Kinetics structural control

- Production adapter: `KineticCalibrationAdapter` using PHREEQC.
- Reaction: calcite; truth `k = 1e-10 mol m^-2 s^-1`, `A = 0.1 m^2 L^-1`.
- Residence times: 10, 30, 100, and 300 d; Ca and HCO3 final concentrations.
- Data-only structural support requires near-collinear log-sensitivity columns,
  numerical rank one, and invariant predictions along a constant-`k*A` ridge.
- A doubled-product off-ridge calculation must produce a non-negligible response,
  ruling out a trivially insensitive experiment.
- An independent `log10(A)` observation with standard deviation 0.10 decades is
  appended as the identifying intervention and must restore numerical rank two.

## Output and stopping rules

Every run receives an immutable run ID and is recorded as RUNNING before model
execution, then PASS or FAIL. A failed run is preserved and never overwritten.
Any code or protocol amendment after locking requires a new lock and a new run
ID. Manuscript claims are revised only after deterministic checks and explicit
scientific review of the generated tables and balanced 2x2 figure.
