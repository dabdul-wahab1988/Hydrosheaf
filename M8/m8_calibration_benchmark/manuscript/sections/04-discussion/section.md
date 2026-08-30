## Discussion

### 4.1 What a calibration report should contain

The central finding is the dissociation between optimisation success, objective
value and parameter identifiability. In the controlled matched-model
experiments every calibration converged, yet recovery ranged from excellent
(dispersivity error 0.0172) to poor (0.3818) across designs with the same
data-generating model, and linearised 95% intervals covered the truth only
75-91% of the time. A report that contains only convergence flags and residual
statistics therefore cannot distinguish a well-identified calibration from a
poorly identified one. The diagnostics that separated the two in this study
are the noise-whitened log-parameter Fisher condition number, per-parameter
recovery metrics, and coverage reported jointly with interval width. We
recommend that calibration reporting state, alongside the residual, the
condition number of the whitened information matrix, the covariance path used,
and per-parameter uncertainty summaries - and that "converged" and "good fit"
be treated as necessary, not sufficient, evidence.

### 4.2 Measurement placement is parameter- and model-class-specific

The optimal-design result is a common optimum rather than a split: all local
Fisher criteria agreed on the 50 d front observation, and that choice repaired
the late-only design for both transport parameters in the matched model. Under
independently generated numerical truth, the same choice remained the
development oracle and robustly improved dispersivity but degraded decay. The
prespecified dual-parameter robustness claim therefore failed, and the honest
conclusion is model-class-bounded: an observation selected by a linearised
Fisher analysis of one forward model need not be jointly optimal under a
different model form. This is not a failure of the diagnostic; it is the
expected consequence of model-form discrepancy, and it argues for reporting
optimal designs with explicit model-form sensitivity, as in the independent
solver arm of this study.

### 4.3 Structural confounding requires independent parameter information

The kinetic experiment is a clean structural negative result: additional
residence-time sampling cannot separate a rate constant from reactive surface
area because the forward law depends on their product. The experiment is not
insensitive - the off-ridge response is 6.42 standard deviations - so the rank
deficiency is a property of the model parameterisation, not of sampling
density. An independent surface-area observation restored rank two, which
identifies the type of information that would resolve the product ridge. The
generality of this result is limited to the information structure of the rate
law; it demonstrates what sampling can and cannot achieve, and it connects to
the equivalence-class language of the M5/M6 benchmarks and the falsification
framing of M3/M4: identifiability is a first-class property of the model-data
pair, not an incidental attribute of the solver.

### 4.4 Active learning: value bounded by the interface and the likelihoods

The portability test shows that a categorical, Jacobian-free heuristic cannot
express a transport-time design decision; it is topology decision support, and
transport OED claims are removed for it. The separately implemented Bayesian
value-of-information workflow did pass every frozen gate on untouched cases:
superior to random in Brier score and information per cost, noninferior to a
strong legacy policy within tight margins, exactly reproducible. Three
boundaries must be respected. First, the information-efficiency noninferiority
pass is narrow (0.00045 nats per cost above the failure boundary) and the
median differences versus legacy were negative; this is a noninferiority
result, not a superiority result. Second, all selected actions were chemistry
panels, so the usefulness of age and connectivity modalities is undemonstrated
under the declared cost model. Third, the benchmark is controlled-synthetic;
it earns scientific-workflow and bounded-excellence status, not field
effectiveness. These boundaries are recorded in the qualification report
referenced in Section 3.6.

### 4.5 Limitations

This study is a controlled-synthetic evaluation; it is not field validation.
The uncertainty intervals are linearised and should not be interpreted as
calibrated posterior uncertainty; the coverage deficits reported are
descriptive and specific to the Gaussian error model used. The five-adapter
domain envelope (vadose, nitrate source, temporal, isotope mixing and
topology) remains future work, and the canonical-design claims from the
exploratory pilot (C1-C3c) are excluded pending fresh prospective
confirmation. The fixed-design sweep is descriptive by design. Finally, the
frontier active-learning benchmark is conditional on the frozen candidate
graph, likelihood model, discrepancy scenarios and declared relative costs;
different costs or likelihoods could change the modality selection outcome.
