# Supplementary Methods

## S1. Provenance and data contract

The installed project was identified as HydroSheaf 0.5.1. Reusable package code, benchmark control scripts, historical analysis records and generated evidence were treated as distinct categories, and only the first was counted as a package capability. The locked controlled-synthetic run recorded a generating source revision, cryptographic hashes for 31 source files, generator and executable metadata, and a flag indicating that source changes had not been committed at generation time. The hashes therefore establish traceability of the run to a source state, but they do not establish that the current committed tree regenerates the run; that comparison is an outstanding action.

Field sources were three files: a Northern Ghana workbook with Dry and Wet worksheets, and two comma-separated files for Lower Anayari and Talensi. SHA-256 digests of all three are recorded in `field_provenance.json`. None contains environmental age-tracer observations, screen intervals, an independently measured directed graph, laboratory covariance metadata or process-truth labels.

The Northern Ghana seasonal separation was established to be a reconstruction rather than an independent resampling. Static water level, ground elevation, borehole depth and distance to river are identical in 160 of 160 wells between the two worksheets, while pH, temperature, total dissolved solids and electrical conductivity differ in essentially every well. A genuine wet-season resampling cannot return an identical water table at every well. The author confirmed that the chemistry is measured and the seasonal attribute constructed. The Dry worksheet is therefore treated as the primary measured panel and the Wet worksheet is excluded from every inferential result; the diagnostic comparison between the two panels is retained in `field_seasonal_reconstruction_diagnostic.csv` solely to characterise the reconstruction and is never reported as measured seasonality.

## S2. Charge balance and quality tiering

Reported concentrations were converted to milliequivalents per litre by dividing by the equivalent weight, defined as molar mass divided by the absolute charge. The charge-balance error was computed as

$$
\mathrm{CBE}\,(\%) = 100 \times \frac{\sum_i c_i^{+} - \sum_j c_j^{-}}{\sum_i c_i^{+} + \sum_j c_j^{-}},
$$

where $c^{+}$ and $c^{-}$ are cation and anion concentrations in meq/L. Only the species a dataset actually reports enter its own balance, so the Talensi balance omits fluoride, which that dataset does not carry. Samples were tiered as quantitative ($|\mathrm{CBE}| \le 5\%$), screening ($5\% < |\mathrm{CBE}| \le 10\%$) or exploratory ($|\mathrm{CBE}| > 10\%$).

Left-censored values were parsed at the reporting limit and flagged, not substituted by a fraction of the limit; one such value occurs, a fluoride result at Lower Anayari. The Talensi anion excess was tested against two alternative interpretations of the bicarbonate column before being attributed to unmeasured cationic species or an unmatched alkalinity determination: interpreting the column as alkalinity expressed as calcium carbonate gives a median error of −36.0%, and interpreting it as already expressed in meq/L gives −97.5%, both worse than the −29.9% obtained under the stated interpretation.

## S3. Facies and isotope treatment

Dominant-ion facies were assigned from the milliequivalent composition. Cations were grouped as calcium, magnesium and combined sodium plus potassium; anions as bicarbonate, chloride and sulfate. A group was labelled dominant when it exceeded 50% of its own total and mixed otherwise, giving a cation–anion label pair. Samples lacking a complete cation or anion set were labelled unclassified rather than imputed.

Local meteoric water lines were fitted by ordinary least squares of $\delta^2\mathrm{H}$ on $\delta^{18}\mathrm{O}$ for each dataset separately, and are plotted only across the isotopic range that dataset spans. Deuterium excess was computed as $d = \delta^2\mathrm{H} - 8\,\delta^{18}\mathrm{O}$ relative to the global meteoric water line. Isotopes were used only as evidence of recharge source, evaporation and mixing; they were never treated as age tracers.

## S4. Specialist candidate generation and calibration

The age candidate generator enumerated declared combinations of mean age, transit-time distribution shape and mixing parameters, with atmospheric-history candidates representing alternative tracer input histories and excess-air corrections. Multi-tracer likelihoods were formed only over tracers present in the observation record. For a candidate $m$ and observed tracer vector $y$, the uncalibrated log score was

$$
\ell(m) = \sum_{j \in \mathcal{O}} \log p\left(y_j \mid m, \sigma_j\right),
$$

where $\mathcal{O}$ is the set of observed tracers and $\sigma_j$ the declared measurement scale. No sealed synthetic age or pathline field was passed to the candidate selector.

Reaction candidates were generated from a fixed family dictionary including a null option. The sparse inverse component minimised

$$
\hat{z} = \arg\min_{z \in \mathcal{Z}} \left\lVert A z - r \right\rVert_2^2 + \lambda_1 \lVert z \rVert_1 + \lambda_2 \lVert z \rVert_2^2,
$$

where $r$ is the reactive residual, $A$ the applicable reaction dictionary, $\mathcal{Z}$ the non-negative and thermodynamically constrained set, and $\lambda_1, \lambda_2$ the declared regularisation weights. The RAPM layer then treated evidence channels as on-off contributions to a family score, with coefficients fitted on case-blocked development data and a temperature parameter fitted for probability calibration. Classwise reliability was checked separately, because a pooled log loss can conceal a poorly calibrated minority family.

## S5. Discrepancy calibration and model averaging

For an interval with point estimate $e$, lower bound $L$, upper bound $U$ and fitted discrepancy scale $s$, the calibrated interval was

$$
\left[\, e - s\,(e - L),\; e + s\,(U - e) \,\right].
$$

The scale was the conservative case-blocked order statistic of the per-target factors required to cover development truth at the target coverage. Its purpose is not to force agreement between components but to prevent a narrow, misspecified component from dominating the mixture.

Discrete model weights were fitted by minimising the case-equal negative log score

$$
-\sum_{c} n_c^{-1} \sum_{i} \log q_{\mathrm{mix},ci}(y_{ci}),
\qquad
q_{\mathrm{mix}}(y) = \sum_{m=1}^{M} w_m q_m(y),
$$

subject to $w_m \ge w_{\min}$ and $\sum_m w_m = 1$, so that a case contributing many candidate edges cannot dominate the fit through row count alone. Weights were frozen before locked scoring and are not interpreted as posterior model probabilities. Pairwise total-variation distance was used as a disagreement diagnostic; at or above 0.25 the target was retained for audit and not promoted to a single reportable interpretation. Convergence used the constrained simplex Karush–Kuhn–Tucker residual together with objective and weight-change criteria. The unconstrained gradient norm was retained as a separate diagnostic because it is not the correct stationarity measure on a simplex; in the locked run it stood at 0.786 against a nominal raw-gradient tolerance of $10^{-5}$, and the integrated claim is bounded accordingly.

## S6. Kinetic identifiability

The kinetic parameter vector comprised $\log k$ and $\log A$. Where the rate law depends on the product $kA$, the local sensitivity columns are collinear in the absence of an independent surface-area observation, and the numerical rank of the sensitivity matrix falls to one. Parameter error, numerical rank and prediction error were therefore reported separately. A case was declared identified only when the declared rank and interval criteria passed; otherwise the parameter report was abstained while the predictive report remained available. This is the mechanism producing an overall identification rate of 0.167 alongside a conditional rate of 1.0.

## S7. Generators, held-out design and adverse controls

Three independent generator families were used: an analytic lattice with branching curvilinear paths and terminal merging; a layered mixing generator with a shortcut and terminal merging; and a MODFLOW 6 and MODPATH 7 model with heterogeneous two-layer flow, vertical exchange and particle-based branch-to-merge topology, to which independent nonlinear chemistry and tracer equations were applied. None imports HydroSheaf. Four development and two locked cases were generated per family.

Observation stresses comprised structured missingness, left censoring, combined stress, and permutation of the tracer, chemistry and head associations. Permutation preserves marginal distributions while destroying case-specific association, so a method that continues to score well under permutation is drawing on leakage rather than signal. A generator critic checked finite values, non-negative tracer quantities, edge structure, source independence and stress counts; a truth-blind audit checked that sealed fields were absent from observed rows, candidate records and policy inputs. A critic record was blocking if a generator imported package code for truth construction, if a sealed state entered an observed row or candidate list, if a stress transformation altered the declared marginal contract, or if a required field was silently imputed.

## S8. External comparison protocol

Comparison with published lumped-parameter age outputs used leave-one-study-unit-out folds across 20 study units, so that no site contributed to any quantity used in its own prediction. Of 1,272 rows, 675 carry finite paired values on both the estimate and the reference; metrics are computed over those and the denominator is reported. Two variants were computed. The uncalibrated variant compares the framework's own log-space estimate with the published log-space value directly. The calibrated variant applies a ridge correction fitted on the training folds. Because the correction is fitted to reproduce the reference, the calibrated variant measures emulation and not independent agreement, and the two are reported separately throughout.

The directed-topology comparison used a published particle-tracking archive supplying 174 reference edges. Two modes were scored and never pooled: a no-prior mode inferring edges from elevation-as-head and downhill gradient with a nearest-neighbour degree limit, and a prior-assisted mode ingesting the reference structure, which measures ingestion fidelity and contains no independent information. Precision, recall and F1 were recomputed from the confusion counts rather than read from the source record.

## S9. Reporting of the integrated gate

One property of the locked record requires disclosure. In the integrated gate, the declared fields `improves_over_random` and `noninferior_to_strongest_specialist` are null, and the corresponding `recomputed_` fields are true; the same is so for the age gate's pooled `mae` field and the reaction gate's `selective_accuracy`. The gate status is therefore resolved from the recomputed raw-benchmark fields. This paper reports the recomputed values, states the case count, and treats the paired intervals as within-run descriptive summaries. A future release should populate the declared fields directly so that the gate record is self-contained.

## S10. Computation, figures and reproduction

Python owns computation. `derive_field_evidence.py` and `derive_model_evidence.py` read primary data and run records and write tidy read-only CSV exports; `build_manuscript.py` generates the tables and assembles the manuscript. R consumes those exports and owns every figure, recomputing no reported statistic. Figure colours use a colour-vision-deficiency-safe categorical palette verified for lightness band, chroma floor and adjacent-pair separation under deuteranopic and tritanopic simulation, with a worst adjacent-pair separation well above the accepted threshold; identity is always carried by a legend or direct label as well as by colour.

The minimum reproduction sequence is: install the locked environment; run the two derivation scripts; run `make_all_figures.R`; run `build_manuscript.py`; and compare the regenerated exports against the archived copies. Regeneration of the locked programme run itself from a committed source tree remains outstanding and is required before any stronger reproducibility statement is made.
