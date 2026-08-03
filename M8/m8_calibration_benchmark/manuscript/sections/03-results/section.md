## Results

All 8,500 calibrations in the matched-model transport experiments (4,000
fixed-design and 4,500 optimal-design calibrations) and all
1,000 calibrations in the independent-model experiments succeeded
(success rate 1.0 in every strategy). All registered artifact tables and the
confirmatory figure are reproduced from the immutable locked runs; the
confirmatory figure is shown in [[FIG:FIG-M8-CONFIRMATORY]]

### 3.1 Sampling placement drives parameter recovery in the fixed design sweep

Across the 16 deterministic four-observation designs, recovery was
parameter-specific and placement-sensitive.

[[TAB:TAB-M8-TRANSPORT]]

The design with the lowest median dispersivity absolute log10 error was `c90_s45`
([45, 75, 105, 135] d) at 0.0172 (95% bootstrap interval 0.0154-0.0185),
while the lowest decay error was achieved by `very_late` ([150, 180, 210,
240] d) at 0.0140 (0.0119-0.0180). The same designs were not jointly optimal:
the tightly clustered early design `c120_s5` had the largest dispersivity
error of the sweep (0.3818) with the worst-conditioned Fisher matrix within
its clustered family (condition number 1025.2), while `very_late`, the design
that was best for decay, was second-worst for dispersivity (0.2404) and had
the overall worst-conditioned Fisher matrix of the sweep (7074.8). Linearised
95% coverage ranged from 0.752 to 0.908 across designs, below the nominal
0.95 in every case. Replicate-wise Spearman correlation between parameter
error and the corrected log-parameter Fisher condition number was positive
for both parameters, with the 95% replicate interval excluding zero for
dispersivity (median 0.587) but not for decay (median 0.490). The fixed
designs are a descriptive panel, not a random population; these correlations
are summarised descriptively as prespecified.

### 3.2 All Fisher criteria select one front observation; the target split is not confirmed

In the optimal-design experiment, every prespecified criterion - local
dispersivity-targeted, decay-targeted, D-, A-, E-optimal and balanced -
selected the same candidate, 50 d, from the three-point late-time base.

[[TAB:TAB-M8-OED]]

The locked decision rule therefore returned
`NOT_SUPPORTED` for the target-specific split: the anticipated divergence
between dispersivity- and decay-targeted next observations did not occur
(registered decision `manuscript/artifacts/m8_claim_decision.json`).

The absence of a split is not a null result. Under local starts, adding the
common 50 d observation reduced median dispersivity absolute log10 error from
0.2154 (no new measurement) to 0.0173, a paired median difference of -0.1952
(95% bootstrap interval -0.2874 to -0.1656); decay error decreased from
0.0182 to 0.0146, a paired difference of -0.0040 (-0.0071 to -0.00003).
Distant starts reproduced essentially the same recovery (dispersivity 0.2118
to 0.0173; decay 0.0171 to 0.0146). The replicate-varying random candidate
was intermediate (dispersivity 0.0883, decay 0.0169 under local starts), and
the worst-joint late candidate was close to no-new-measurement (dispersivity
0.2089). The balanced 2x2 figure summarises predicted marginal uncertainty,
realised recovery and the kinetic structural result
([[FIGREF:FIG-M8-CONFIRMATORY]]).

### 3.3 The kinetic rate law identifies k*A, not k and A separately

In the production PHREEQC kinetic adapter, both the single-time and the
four-residence-time data-only designs produced a local log-sensitivity matrix
of numerical rank one with a zero minimum eigenvalue and infinite condition
number.

[[TAB:TAB-M8-KINETICS]]

Predictions along the constant-k*A ridge
differed by at most 1.33e-6 declared observation standard deviations, while
doubling the product changed predictions by 6.42 standard deviations - so the
rank deficiency is structural, not a dead experiment. Appending an
independent log10(A) observation (standard deviation 0.10 decades) restored
numerical rank two (condition number 45.9, parameter correlation -0.957).
This is controlled evidence about the information required to identify k and
A; it is not field validation of a surface-area measurement programme.

### 3.4 Independent numerical truth: the 50 d choice is parameter-specific, not dual-parameter robust

The independently generated 240-cell numerical truth passed the locked grid
gate (maximum 0.203 declared observation standard deviations from the
480-cell reference, below the 0.25 limit). The analytical E-optimal criterion
and the independent development oracle both selected 50 d, agreeing with the
matched-model optimum.

[[TAB:TAB-M8-INDEPENDENT]]

On the 250 locked test replicates, the 50 d observation improved dispersivity
recovery robustly: median absolute log10 error fell from 0.8262 (no new
measurement) to 0.1674, a paired difference of -0.6690 (95% bootstrap
interval -0.7374 to -0.5757).

[[TAB:TAB-M8-INDEPENDENT-CONTRASTS]]

Decay recovery, however, worsened: error increased from 0.1367 to 0.1541, a paired
difference of +0.0210 (+0.0092 to +0.0276). The prespecified dual-parameter
robustness rule therefore failed: under model-form discrepancy the same
observation that is robustly useful for dispersivity degrades decay.
Linearised coverage collapsed under the independent truth for the 50 d
strategy (0.02 for dispersivity, 0.004 for decay, versus 0.788 and 0.64 for
no new measurement), illustrating that the linearised covariance, not just
the point estimate, is sensitive to model-form discrepancy. The defensible
result is parameter-specific and is stated exactly as such; averaging the two
opposite effects into a joint improvement is not supported.

### 3.5 The existing active-learning heuristic is not actionable for transport-time design

The portability test returned
`NOT_ACTIONABLE_FOR_TRANSPORT_TIME_SELECTION`: the existing
`rank_next_measurements` routine produced eight recommendations with a single
unique priority score and no explicit transport concentration sampling time,
because it consumes categorical topology evidence classes and never evaluates
candidate-time Jacobians (registered record
`provenance/runs/RUN-M8-INDEPENDENT-20260728-01/active_learning_portability.json`).
No transport optimal-design performance claim for the heuristic is made.

### 3.6 A prospectively qualified topology active-learning workflow

Across 24 untouched independent MODFLOW 6/MODPATH 7 cases, the scenario-robust
Bayesian policy was actionable in every case (actionability 1.0) with mean
candidate recall 0.9861.

[[TAB:TAB-M8-FRONTIER-AL-SUMMARY]]

Against random acquisition it improved the median edge Brier score by -0.02483 (95% paired
bootstrap interval -0.02897 to -0.01919, median 0.04810 versus 0.07385) and
joint-hypothesis entropy reduction per cost by +0.03581 nats per relative
cost (0.02206 to 0.04218) - both prespecified superiority gates passed.

[[TAB:TAB-M8-FRONTIER-AL-CONTRASTS]]

Against the strong legacy
uncertainty-chemistry policy it was noninferior, not superior: Brier
difference +0.00148 (interval 0.00000 to 0.00582, below the frozen +0.01
margin) and information efficiency -0.00422 (-0.00955 to 0.00000, above the
frozen -0.01 margin but by only 0.00045 nats per cost - a narrow pass).
PR-AUC differences were exactly zero against every comparator (PR-AUC 1.0 for
all executable strategies). All 120 proposed-policy actions were chemistry
panels: age-tracer and connectivity-tracer actions were available but never
selected under the declared costs and likelihoods, so multimodal usefulness
remains a capability, not an empirical contribution of this run. All seven
registered outputs of the authoritative run regenerated byte-identically from
the locked sources (record
`provenance/runs/RUN-M8-FRONTIER-AL-20260728-01/scientific_workflow_qualification.json`).
Robustness over the mean-information variant did not establish performance
superiority: the two policies chose different action sequences in 15 of 24
cases with no significant outcome difference, and the robust policy remained
below the unattainable realised oracle (+0.00837 Brier, interval 0.00609 to
0.01110), as expected.
