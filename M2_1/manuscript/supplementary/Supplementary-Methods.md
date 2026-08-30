# Supplementary Methods

## S1. Provenance and data contract

The supplementary methods were assembled from the current HydroSheaf source repository and were kept separate from earlier manuscript material. The installed project was identified as HydroSheaf version 0.5.1. Reusable package code, benchmark and validation control scripts, historical analysis records and generated evidence were kept as distinct categories. This separation prevented a benchmark client or an old result from being represented as a production package capability.

The primary controlled-synthetic run was a locked integration run. Its provenance manifest recorded a generating source revision, whether source changes had been committed, generator and executable metadata, and cryptographic hashes for 31 source files. The run was generated before all source changes were committed, so the hashes provide traceability but do not replace a rerun from a committed source tree.

The field-data archive contained a Northern Ghana workbook with wet and dry worksheets and two CSV sources for Lower Anayari and Talensi. The field source material did not contain environmental age-tracer observations, screen intervals, an independently measured directed graph, laboratory covariance metadata or process-truth labels. No additional source document was assumed to be available evidence.

The execution environment was resolved through the project environment. The software-version record included HydroSheaf 0.5.1, the Python dependency lock information, the external MODFLOW 6 executable used by the independent flow generator and the MODPATH 7 executable probe. A successful command was not treated as sufficient by itself: the provenance manifest recorded executable version or hash where available, exit status and output hash. Generated caches, temporary solver folders and bytecode were excluded from the source-hash inventory.

## S2. Evidence representation and candidate units

Each observation was normalised into a node record containing identifiers, measured variables, optional uncertainty and provenance. Candidate directed edges were records connecting two observed nodes. An edge was a hypothesis generated from the declared candidate universe, not a measured flow path. The candidate universe could be generated independently from coordinates and heads, from a specialist comparator, or from a synthetic truth generator; the source of the universe was recorded so that an apparent performance advantage could not be attributed to a hidden candidate-set advantage.

Each evidence channel returned a distribution or a compatible set where possible. Hydraulic evidence was represented as a directional likelihood or support score. Age evidence was represented by candidate age distributions and calibrated intervals. Chemistry evidence was represented by residuals, reaction-family scores and thermodynamic compatibility. Kinetic evidence returned predictive concentrations and an effective-rate or parameter report with identifiability status. Stable isotope fields in the field application were retained as source and mixing diagnostics only. Missing variables were not imputed as positive evidence. Structured missingness, left censoring and permutations were generated only inside the locked synthetic observation contract.

For any response vector (y), the interval coverage estimator was

<!-- EQ:EQ-S1 -->
\[
\widehat{C}_{1-\alpha}=n^{-1}\sum_{i=1}^{n}I\{y_i\in[L_i,U_i]\},
\]

and mean interval width was (n^{-1}\sum_i(U_i-L_i)) over reportable intervals. Abstentions were retained in the denominator for coverage including abstention. Selective risk was calculated only over committed outputs and was accompanied by the commitment rate. Consequently, a method could not obtain a high selective score by silently dropping difficult rows.

For point-error summaries, the locked scorer used the declared target scale and did not mix age units with concentration or probability units. Coverage including abstention used the full eligible locked set; conditional selective risk used only committed predictions and reported the denominator. Censored measurements were represented by their observation model in the synthetic contract and were not silently replaced by a field-style half-detection-limit rule. The one censored fluoride result in Lower Anayari was retained as a field-data limitation rather than used to claim a calibrated synthetic censoring result.

## S3. Specialist age generation, atmospheric histories and calibration

The age candidate generator enumerated declared combinations of mean age, age-distribution shape and mixing parameters. Atmospheric-history candidates represented alternative tracer input histories and excess-air corrections. Multi-tracer likelihoods were formed only for tracers present in the observation record. For a candidate (m) and tracer vector (y), the uncalibrated log score was

<!-- EQ:EQ-S2 -->
\[
\ell(m)=\sum_{j\in\mathcal{O}}\log p\{y_j\mid m,\sigma_j\},
\]

where (mathcal{O}) was the set of observed tracers and (sigma_j) was the declared measurement scale. The age posterior was normalised over the candidate universe and converted to an interval or a set of modes. No sealed synthetic age or pathline field was passed to the candidate selector.

Development observations were divided by case before calibration. Age interval calibration estimated a non-negative scale from development residuals and applied it without re-fitting on locked cases. The locked report recorded coverage, width, relative width, acceptance, mean absolute error (MAE) and selective risk. A result was abstained when no finite likelihood was supported, the required tracer set was absent, the posterior was numerically invalid or the candidate output failed the declared acceptance rule. The USGS TracerLPM literature was used to motivate model-dependent age distributions and multi-tracer comparison, not as field truth for the Ghana data [@jurgens2012tracerlpm; @visser2013multitracer; @mccallum2015limitations].

## S4. Reaction templates, RAPM and calibration

Reaction candidates were generated from a fixed family dictionary with a null option. Candidate features included signed ion changes, scaled chemistry residuals, thermodynamic direction evidence, isotope or tracer support where available, and missing-channel indicators. The candidate generator was independent of the sealed synthetic reaction fields. The sparse inverse component used a bounded objective of the form

<!-- EQ:EQ-S3 -->
\[
\widehat{z}=\arg\min_{z\in\mathcal{Z}}
\left\|Az-r\right\|_2^2+\lambda_1\|z\|_1+\lambda_2\|z\|_2^2,
\]

where (r) was the reactive residual, (A) the applicable reaction dictionary, (mathcal{Z}) the non-negative and thermodynamic constraint set, and (lambda_1,lambda_2) the declared regularisation values. The output was then normalised to reaction families rather than presented as an unrestricted mineral history.

The regularised adjusted plus-minus (RAPM) layer treated evidence channels as on/off contributions to a family score. Coefficients were fitted with case-blocked development data, and a temperature parameter was fitted for probability calibration. Classwise reliability was checked separately because a pooled log loss can conceal a poorly calibrated minority family. The locked report contained multiclass log loss, coverage, classwise expected calibration error, selective risk and false commitment. Bootstrap recurrence was retained as a stability diagnostic, not as evidence that a mineral source was present. External PHREEQC documentation supported the thermodynamic and inverse-geochemical context, but the field application was not labelled as independent PHREEQC validation [@parkhurst2013phreeqc].

The reaction gate thresholds were coverage at least 0.25, maximum classwise expected calibration error no greater than 0.35, false commitment no greater than 0.10 and selective risk no greater than 0.4, with finite output and held-out mechanism requirements. False commitment meant a committed family report incompatible with the sealed outcome under the declared dictionary. The reported coverage of 0.7394 therefore passed the declared minimum, but it was not rounded into a claim of complete classification. The classwise calibration error of 0.2710 remained visible because a pooled log-loss improvement can coexist with uneven reliability among families.

## S5. Discrepancy calibration and discrete model averaging

For each target, component models returned predictive distributions over declared outcomes. Development observations were used to estimate a discrepancy scale sufficient for the observed component residuals. For an interval with estimate (e), lower bound (L), upper bound (U) and fitted scale (s), the calibrated interval was [e-s(e-L), e+s(U-e)]. The scale was the conservative case-blocked order statistic of the per-target factors required to cover development truth at the target coverage. The purpose was not to force all models into agreement; it was to prevent a narrow but misspecified component from dominating the mixture.

Model weights were fitted from case-blocked development log score subject to a simplex and a declared weight floor. The weighted mixture was

<!-- EQ:EQ-S4 -->
\[
q_{\mathrm{mix}}(y)=\sum_{m=1}^{M}w_mq_m(y),
\qquad w_m\geq w_{\min},\quad\sum_mw_m=1.
\]

Weights were fitted by minimising the negative case-equal log score, \( -\sum_c n_c^{-1}\sum_i \log q_{\mathrm{mix},ci}(y_{ci}) \), over the simplex with the declared weight floor. Thus a case with many candidate edges could not dominate the fit merely through its row count. The fitted weights were frozen before locked scoring and were not interpreted as posterior model probabilities.

Pairwise total-variation distance, a measure of distributional disagreement, was used as a diagnostic. At or above 0.25, the target was retained for audit but the compatible outcome set was not promoted to a reportable single interpretation. Convergence was based on the constrained simplex Karush–Kuhn–Tucker (KKT) residual, objective change and weight change. The raw gradient norm was retained because it is a diagnostic of the unconstrained objective and is not, by itself, the correct stationarity measure on a simplex. A model-averaging record therefore included weights, fit hash, convergence status, KKT residual, objective trace, pairwise disagreement and reportability.

## S6. Kinetic forward model and next-measurement policy

The kinetic adapter constructed PHREEQC-compatible kinetic blocks and evaluated a controlled forward response over elapsed time. Predictions were scored against held-out observations using predictive root-mean-square error (RMSE) and calibrated intervals. The kinetic parameter vector included (\log k) and (\log A). When the rate law depended on k×A, the local sensitivity columns were collinear without an independent surface-area observation. Numerical rank, parameter error and prediction error were therefore reported separately. A case was identified only when the declared rank and interval criteria passed; otherwise the parameter report was abstained while the predictive report could remain available.

The next-measurement policy received a finite candidate action list with action IDs, costs and pre-measurement scenario distributions. Expected information gain and declared utility were computed from those distributions. A policy selector was required to be truth-blind and to return an action or abstention before hidden state release. After all selectors returned, the scorer attached the hidden state to the selected action, calculated benefit minus cost, and computed paired contrasts. Uniform random expectation was scored over feasible actions rather than represented by a single unstable random draw. The six locked cases were therefore prospective within the synthetic contract, but not prospective field measurements [@chaloner1995design; @sreekanth2017monitoring].

## S7. Independent generators, held-out cases and adverse controls

The analytic-lattice generator used branching curvilinear paths with terminal merging and an analytic head and tracer forward model. The independent mixing generator used layered branching, a shortcut and terminal merging with a weighted end-member and process-increment chemistry model. The MODFLOW/MODPATH generator used a heterogeneous two-layer flow model with vertical exchange and particle-based branch-to-merge topology, then applied independent nonlinear chemistry and tracer equations. The MODFLOW 6 and MODPATH 7 executables were external to HydroSheaf and their hashes were recorded [@langevin2017modflow6; @pollock2016modpath].

Four development cases were generated for each family and two locked cases for each family. Development cases were used for calibration and model weights. Locked cases were never used to select candidates, set thresholds or choose policies. Observation stress scenarios included structured missingness, left censoring, combined stress, tracer permutation, chemistry permutation and head permutation. Permutation controls preserved marginal values while destroying case-specific association. A generator critic checked finite values, non-negative tracer quantities, edge structure, source independence and stress counts. A truth-blind audit checked that sealed truth fields were not present in observed rows, candidate records or policy inputs.

The performance contract assessed four domains separately. PASS meant that all declared criteria were met, WEAK meant that a finite result fell below at least one performance floor, and ABSTAIN meant that the software withheld a result or claim. Age required calibration, held-out generators, family stratification, target coverage, a width ceiling and selective-risk ceiling. Reaction required calibration, held-out generators, mechanism stratification, log-loss non-inferiority, coverage, expected calibration error, false-commitment and selective-risk limits. Kinetics required held-out regimes, predictive and parameter limits, interval coverage, conditional identification and explicit `k×A` reporting. The integrated gate required discrepancy calibration, converged model averaging, independent complete outcomes, false-commitment control and at least six locked prospective cases. None of the four gates was interpreted as universal superiority.

The locked result table retained the generator family, case identifier, mechanism or stress label, observation subset, baseline identifier and claim status for each row. Separate machine-readable files retained the full row-level records and pooled gate values. This allowed a future review to ask whether a pooled PASS was driven by one easy family or by the no-stress subset. A critic record was blocking if a generator imported HydroSheaf code for truth construction, if a sealed state entered an observed row or candidate list, if a stress transformation changed the declared marginal contract, or if a required field was silently imputed. Non-blocking warnings were retained in the audit but did not disappear from the final report.

## S8. Field-data processing and claim ceiling

Field ion concentrations were harmonised to the declared units and independently screened by charge-balance error. Northern Ghana contained 160 repeated wells and 320 records; the wet and dry labels were retained, but no actual sampling-date sequence was assumed. Lower Anayari had 41 records and one censored fluoride value. Talensi had 63 records and nine missing redox-potential values. Coordinate-neighbour edges from the field chemistry runner were treated as candidate screening constructions. Chemistry R², the coefficient of determination, process-stability or recurrence scores and reaction-family outputs were interpreted as fit or stability diagnostics. They were not interpreted as calibrated field probabilities, measured paths or unique mechanisms.

The field data supported readiness audits, seasonal chemistry prediction, reaction-family plausibility, measurement ablation, alternative-edge sensitivity and non-identifiability reporting. They did not support groundwater-age validation, exact directed connectivity, screen-resolved vertical exchange, a unique reaction history or a field digital twin. Stable isotope values were not promoted to age observations.

The field crosswalk also recorded whether a variable was present, repeated, censored, unit-resolved and suitable for a future claim. The Northern Ghana sheets were treated as repeated seasonal panels rather than a dated time series. Talensi charge-balance quality was sufficiently heterogeneous that its rows were used for screening and exploratory diagnostics as well as quantitative rows. Lower Anayari was retained as a small chemistry example, not as an independent validation sample. These choices preserve the field data's usefulness without manufacturing reference labels that are not present in the files.

## S9. Reproduction and claim decisions

The manuscript assembly process assembled section files, regenerated summary tables from the locked machine-readable gate records, extracted the source-hash inventory and checked the execution gate. Word counts were performed after assembly. A final manuscript-reviewer report was required before the package could be marked complete. Critical and major issues had to be resolved in the text, evidence ledger or claim ledger. The final allowable claim was: HydroSheaf passed bounded controlled-synthetic component and integrated decision diagnostics under the declared generators and observation models, while its field application remained a scope and readiness demonstration. General superiority and field validation therefore remained abstentions.

The exact comparator identifiers, input channels, tuning values and fingerprints are retained in machine-readable baseline and specialist-comparator registries. The minimum regeneration sequence was: activate the project environment; run the locked integration benchmark; verify the execution, generator-quality and external-executable gates; assemble the manuscript and recompute its word counts; run the citation and artifact checks; and then run the manuscript-reviewer audit. A rerun from a committed source tree must replace the provisional provenance before any stronger reproducibility sentence is used. The supplementary package is a method record and audit aid, not a substitute for the machine-readable run records.
