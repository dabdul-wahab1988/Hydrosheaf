# Graph-Constrained TTD Virtual Benchmark Protocol

## Status and scientific scope

This protocol defines a development-stage, controlled-synthetic benchmark for
graph-constrained groundwater transit-time-distribution (TTD) inference.  It
is designed to answer a narrow question:

> Given declared observation, topology, recharge, mixing, and discrepancy
> scenarios, does an inference method recover predeclared TTD functionals or
> abstain more appropriately than matched local and graph controls?

Passing this benchmark supports only controlled-synthetic recovery under its
versioned generators.  It is not field validation, proof of a unique field
flow path, or a universal groundwater-model superiority claim.

The protocol complements, but must not overwrite or re-label, the locked M7
controlled synthetic twin.  Its required new targets are time-indexed tracer
input/output records, local recharge, mixing, and distributed edge kernels.

## Forward-model contract

Each virtual case defines a directed groundwater network with known truth.  A
node response is generated using

\[
g_v(a,t) = \rho_v(t) r_v(a,t) + \sum_{u \to v}\pi_{uv}(t)
             \left[g_u * h_{uv}\right](a,t),
\]

where `r_v` is local recharge/source input, `rho_v` its fraction, `pi_uv` are
non-negative upstream fractions, and `h_uv` is a normalized causal edge TTD.
The simulator must preserve mass and report the truth separately from all
inference-visible observations.

Static cases use stationary `h_uv(a)`.  Dynamic cases introduce declared
nonstationarity in one or more of edge kernels, recharge fractions, mixing
fractions, or input forcing.  Every case includes a branch and merge so that a
serial-only model cannot win by construction.

## Independence, blinding, and provenance

1. The generator must not import the HydroSheaf inverse implementation.
2. An inference function receives only the candidate graph and observation
   view; it must fail if `true_` or `truth_` fields are present.
3. Truth is released only to the scoring stage after inference output is
   persisted.
4. Each run records generator version/source hash, seed, scenario parameters,
   observation mask, code revision, Python environment, input hashes, and
   hashes for every output artifact.
5. At least two structurally distinct generator families are required before a
   headline method claim: the in-repository analytic/particle network family
   and a process-oriented MODFLOW/MODPATH-style family when its executables are
   available.  Neither generator may share inverse code with the evaluated
   method.

## Predeclared estimands

The primary estimands are intentionally coarse and scientifically interpretable:

- young-water fraction `F_y(T)` for predeclared thresholds `T`;
- mean and selected quantile age intervals at observation nodes;
- source/local-recharge contribution intervals;
- edge inclusion/exclusion only where the observation design supports it;
- held-out tracer signal predictions.

Full edge-kernel recovery is a secondary diagnostic.  A wide interval or an
explicit abstention is a valid result when the data do not identify a finer
quantity.

## Observation and stress design

Each case exposes only irregular, possibly sparse tracer observations derived
from seasonal and event-like forcing.  The benchmark varies:

- sampling frequency and missingness;
- tracer noise and censoring;
- correct versus wrong recharge forcing;
- local recharge omitted from the fitted model;
- topology uncertainty, reversed directions, random edges, and edge removal;
- stationary versus transient transport/mixing;
- branches, merges, shortcut paths, and source mixtures.

The observation table must keep tracer-specific missingness.  Missing isotope
or tracer values must never be silently replaced with zero.

## Required comparison arms

Every scored scenario reports, at minimum:

1. a local, no-graph TTD baseline;
2. the declared candidate-graph method;
3. a correct-graph oracle control used only to establish an attainable
   conditional ceiling;
4. reversed, randomised, and edge-removal graph controls;
5. a conventional lumped/mixture TTD comparator where its assumptions apply;
6. a source/mixing-aware graph model.

The correct-graph result does not prove that a field graph is correct.  It is
a conditional virtual reference used to distinguish implementation failure
from topology uncertainty.

## Scoring and decision rules

For every target and scenario, report:

- point error only when a point estimate is declared;
- empirical interval coverage and interval width;
- abstention rate, reason, and whether abstention was appropriate relative to
  declared information content;
- held-out time-series predictive error on dates not used for fitting;
- topology precision/recall and false-edge acceptance/rejection;
- source-fraction and young-water-fraction error;
- sensitivity to discrepancy scenarios.

Metrics are not pooled across incompatible tracer units.  No method may be
called superior from an in-sample residual alone.  A positive claim requires
predeclared performance criteria, confidence intervals across independent
seeds/cases, and no harmful result under the negative controls.

## External WATRES panel

The WATRES release associated with Duchemin et al. is an optional,
no-retuning external temporal-TTD component panel.  It is catchment-level, not
groundwater-network truth.  It can test the temporal-TTD adapter but cannot
validate edge-wise groundwater connectivity or graph recovery.

- Article DOI: `10.1029/2025WR042835`
- Data record: <https://zenodo.org/records/15658651>

Its download is deliberately opt-in and must be checksum-validated.  A missing
or incompatible archive is reported as `ABSTAIN`/unavailable, never synthesized
or treated as a passing external test.

## Field-transfer boundary and future campaign

The Ghana/UE R field records are retained only as a data-readiness and campaign
design panel.  They do not provide independent time-resolved groundwater
input-output TTD truth.  A later prospective field validation would require a
preselected, repeatedly sampled recharge--well--well network with screened
intervals, repeated heads/pumping context, tracer input histories, and
independent connectivity or age evidence.  The virtual benchmark should be
used first to calculate the sampling design needed for that campaign.

## Required run artifacts

Each benchmark run writes a self-contained directory containing:

- a versioned run manifest and generator provenance;
- truth only in a scoring-only artifact;
- an inference-visible observation artifact without truth fields;
- scenario and control definitions;
- predictions, abstentions, and metrics;
- an execution/claim gate describing what passed, failed, or remains deferred.

Any generated report must state the controlled-synthetic evidence ceiling in
its first methods or results paragraph.
