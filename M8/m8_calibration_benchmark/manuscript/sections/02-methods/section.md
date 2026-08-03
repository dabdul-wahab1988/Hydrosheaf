## Methods

All benchmarks were executed with the production HydroSheaf package (v0.5.1)
on Windows 11, Python 3.14.6, with numpy 2.4.6, pandas 2.3.3, scipy 1.17.1,
scikit-learn 1.9.0 and flopy 3.10.0; the independent topology benchmark used
MODFLOW 6.7.0 and MODPATH 7.2.001 binaries locked by SHA-256. The three
confirmatory runs were pre-registered: protocol, source hashes and run scripts
were locked before any test-case generation, and each run received an
immutable run ID whose manifest records environment, seed families, file sizes
and artifact checksums.

### 2.1 The calibration engine and information diagnostics

Calibration uses the production `PESTGLM` engine: a Gauss-Levenberg-Marquardt
optimiser solved by trust-region reflective least squares, operating in
log10-parameter space with bound constraints, parallel Jacobian evaluation,
and covariance computed by matrix inverse with an SVD-pseudoinverse fallback
when the inverse fails. Adapter-created priors are disabled in all primary
experiments so that random starting values do not act as random
regularisation centres.

Information diagnostics are computed from the noise-whitened log-parameter
Fisher information matrix. Sensitivities are central finite differences of
the forward model with respect to log10(parameter) (step 1e-4 log10 units,
with 1e-3 and 1e-5 stability checks), divided by the declared Gaussian
observation error standard deviation. From the whitened Jacobian we report
the condition number, minimum eigenvalue, numerical rank, predicted marginal
log-parameter standard deviations, parameter correlation and sensitivity
cosine. Recovery is measured as absolute log10 parameter error between the
estimated and known truth; uncertainty calibration is measured as empirical
coverage of nominal 95% intervals; and optimisation status is recorded per
replicate.

### 2.2 Matched-model transport benchmark

The transport benchmark uses the production `TransportCalibrationAdapter` in
analytical one-dimensional advection-dispersion mode. Truth is dispersivity
2.0 m and first-order decay 0.005 d^-1, over a 10 m domain at 0.1 m d^-1 with
unit source concentration. Observations are concentrations with independent
Gaussian error sigma_i = max(0.01, 0.03 * |C_i|).

Two experiments were run. First, a fixed-design sweep of 16 four-observation
designs spanning early, late, clustered and log-spread sampling. Second, an
optimal-design experiment: starting from a three-point late-time base
[150, 180, 210] d, eight candidate times (50, 70, 90, 105, 120, 150, 180,
240 d) were scored by parameter-marginal-variance (dispersivity-targeted and
decay-targeted), D-, A- and E-optimality, and a balanced marginal-variance
objective in the whitened log-parameter Fisher matrix. Each candidate was
then verified by actually adding it to the design and recalibrating; controls
were a replicate-varying random candidate, the worst-joint candidate, and no
new measurement.

Replication follows a paired design: 250 replicates per fixed design and per
OED strategy, with common standard-normal errors and common starting values
across strategies within each replicate. Primary starts are within +/-0.25
log10 units of truth; a prespecified sensitivity analysis uses +/-1.0 log10
units. Seed families are 2026072801 (data) and 2026072802 (bootstrap).
Median errors and paired median differences carry percentile bootstrap
intervals from 2,000 resamples. The fixed-design sweep is treated as
descriptive: a replicate-wise Spearman coefficient across the deterministic
designs is summarised over replicates, and classical significance testing
across the 16 hand-selected designs is explicitly prohibited.

The prespecified decision rule for the optimal-design experiment: full
target-dependent trade-off support requires (i) locally predicted
dispersivity- and decay-targeted candidate times to differ, (ii)
dispersivity-targeted sampling to have lower realised dispersivity error than
decay-targeted sampling with a paired 95% bootstrap interval excluding zero,
and (iii) the converse for decay. One successful directional comparison is
partial support; otherwise the trade-off is not confirmed.

### 2.3 Kinetic structural control

The kinetic experiment uses the production `KineticCalibrationAdapter`
driving PHREEQC speciation and kinetic reaction calculations
[@parkhurst2013phreeqc]. The reaction is calcite with truth
k = 1e-10 mol m^-2 s^-1 and reactive surface area A = 0.1 m^2 L^-1. Two
data-only designs are considered: a single residence time (two observations,
Ca and HCO3 final concentrations) and four residence times (10, 30, 100,
300 d, eight observations). Sensitivities are whitened by a declared 0.005
mmol L^-1 concentration error.

Structural support for non-identifiability requires near-collinear log
sensitivity columns, numerical rank one, and invariant predictions along a
constant k*A ridge. A doubled-product off-ridge calculation must produce a
non-negligible standardised response, ruling out a trivially insensitive
experiment. The identifying intervention appends an independent log10(A)
observation with standard deviation 0.10 decades and must restore numerical
rank two.

### 2.4 Independent-model robustness

Truth for the independent-model experiment is generated by an implicit
finite-volume/upwind numerical solver of the one-dimensional
advection-dispersion-decay equation, implemented only in the benchmark script
and sharing no code with the calibration forward model. The 240-cell truth is
accepted only if its maximum difference from a 480-cell reference run is no
more than 0.25 declared observation standard deviations (the locked grid gate).
The production analytical adapter then calibrates against this independent
truth.

The development oracle is selected on 80 development replicates (seed family
2026072810) as the candidate with minimum median joint log-parameter error,
then frozen. Verification uses 250 independent locked test replicates (seed
family 2026072820), paired across the strategies no-new-observation,
analytical E-optimal, frozen development oracle, and replicate-wise random,
with 5,000 paired bootstrap resamples. Independent-model robustness is
supported only if the grid gate passes and the E-optimal strategy has
strictly negative 95% paired bootstrap intervals for both dispersivity and
decay absolute log10 error versus no new observation.

### 2.5 Active-learning portability test

The existing `rank_next_measurements` heuristic consumes categorical topology
evidence classes and never touches the Jacobian. Its portability to
transport-time design is tested by supplying candidate times as opaque
candidate IDs with otherwise identical categorical evidence. The heuristic is
actionable only if it returns a non-tied recommendation that identifies a
transport concentration sampling time; no ad hoc translation from
age/isotope/topology recommendations is permitted. Otherwise it is reported
as `NOT_ACTIONABLE_FOR_TRANSPORT_TIME_SELECTION`.

### 2.6 Frontier topology active-learning workflow

A separate, explicit Bayesian value-of-information interface was implemented
for sequential topology-measurement design, following the robust-expected
information-gain formulation of Go and Isaac [@go2022robust] with greedy
joint batch acquisition in the spirit of BatchBALD [@kirsch2019batchbald].
The method contract requires: named joint hypotheses with normalised
posterior weights; actions defined by measurement type, target, feasibility,
relative cost and conditional predictive distributions; scalar expected
information gain by deterministic Gauss-Hermite quadrature; a declared
mean/worst-case robust utility over model-discrepancy scenarios; a combined
utility of normalised expected information gain (weight 0.05) and expected
posterior Brier-risk reduction (weight 0.95) divided by relative cost;
posterior updating with scenario-averaged likelihoods; batch selection by
marginal joint information using a reproducible Sobol estimator so redundant
measurements are not double-counted; and abstention when the maximum robust
information gain is below a locked floor (0.002 nats).

Truth is generated by an independent MODFLOW 6/MODPATH 7 heterogeneous
aquifer generator with nonlinear synthetic geochemistry [@langevin2017modflow6;
@pollock2016modpath]. Twenty-four untouched cases (locked-test seeds
7601-7624) share 256 topology particles, hidden outcomes, a relative-cost
budget of 10 units, a maximum of five sequential actions, and one frozen
prior calibration fitted only on development cases. Action costs are 2
(major-ion chemistry panel), 5 (groundwater-age tracer) and 9 (directed
connectivity tracer) relative units. Predictive likelihoods are
class-conditional ridge regressions fitted on development cases with
coefficient shrinkage of 0.25 and recalculated residual variances. Three
model-discrepancy scenarios (nominal, separation-stress, noise-stress) are
combined with a robust acquisition weight of 0.75 on the worst scenario and
0.25 on the scenario-weighted mean.

Comparators receive identical initial particles, hidden outcomes, budget and
update rule: the proposed robust policy; a mean-information variant without
worst-case weighting; a strong legacy uncertainty-priority policy that
selects the edge of greatest marginal Bernoulli entropy and requests a fixed
chemistry panel; uniformly random affordable action; and an unattainable
realised oracle excluded from superiority gates. Outcomes per case and
strategy are candidate recall, final edge Brier score, PR-AUC and log loss,
joint-hypothesis entropy reduction per cost, cost spent, action count,
abstention and selected measurement types [@brier1950verification;
@davisgoadrich2006prcurves]. All proposed-versus-comparator contrasts use
5,000 paired case-bootstrap resamples with common random numbers.

The frozen claim gate requires: mean candidate recall at least 0.80;
actionability in at least 90% of cases; a strictly negative paired 95%
interval for Brier score versus random and an upper interval versus legacy
below +0.01; a strictly positive paired interval for entropy reduction per
cost versus random and a lower interval versus legacy above -0.01 nats per
cost; PR-AUC lower interval bounds above -0.01 versus both comparators; and
byte-identical regeneration of every registered artifact from the locked
sources. Legacy-superiority, field effectiveness and multimodal usefulness
claims are explicitly outside what any passing result can establish.
