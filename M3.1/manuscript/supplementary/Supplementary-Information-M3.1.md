# Supplementary Information

**Conditional Value of Graph Priors in Multi-Tracer Groundwater-Age Inference: A Hydrosheaf Benchmark of Public Aquifer Datasets (M3.1)**

Dickson Abdul-Wahab, Dickson Adomako, Gibrilla Abass, Ebenezer Aquisman Asare, Samuel Ganyaglo

## Supplementary Methods

### S1. Residence-time convolution

For a tracer with atmospheric input history $I_j(t)$, decay constant $\lambda_j$, and residence-time distribution $g(\tau \mid \theta)$, the modelled concentration at sample time $t_s$ was:

$$
\widehat{C}_j(t_s, \theta) = \int_0^{\tau_{\max}} I_j(t_s - \tau)\, \exp(-\lambda_j \tau)\, g(\tau \mid \theta)\, \mathrm{d}\tau. \tag{S1}
$$

For stable atmospheric tracers such as SF6 and CFCs, $\lambda_j = 0$. For radioactive parent tracers such as 3H, 14C, 39Ar, and 85Kr, $\lambda_j = \ln(2) / T_{1/2,j}$. For daughter ingrowth in the 3H/3He system, the simplified ingrowth contribution was represented as:

$$
\widehat{C}_{3He_{trit}}(t_s, \theta) = \int_0^{\tau_{\max}} I_{3H}(t_s - \tau)\left[1 - \exp(-\lambda_{3H}\tau)\right] g(\tau \mid \theta)\, \mathrm{d}\tau. \tag{S2}
$$

For carbon-14, the activity prediction was treated as:

$$
\widehat{A}_{14}(\theta) = A_0 \int_0^{\tau_{\max}} \exp(-\lambda_{14}\tau)\, g(\tau \mid \theta)\, \mathrm{d}\tau, \tag{S2b}
$$

where $A_0$ represents the initial or corrected 14C activity. For helium-4, the diagnostic accumulation model was represented as:

$$
\widehat{He}_4(\theta) = He_{4,\mathrm{bg}} + r_{He}\,\bar\tau, \qquad \bar\tau = \int_0^{\tau_{\max}} \tau\, g(\tau \mid \theta)\, \mathrm{d}\tau, \tag{S2c}
$$

where $He_{4,\mathrm{bg}}$ is background helium, $r_{He}$ is the accumulation-rate term, and $\bar\tau$ is mean residence time.

The tritium input function was represented using GNIP/WISER precipitation isotope data from the station nearest each site [@iaeawmo2026gnip]. Most North American WISER stations stop reporting between 1965 and 1999, whereas the benchmark samples were collected between 2004 and 2020. Each station record was therefore continued past its final observation along the regional post-bomb decline curve, rescaled to join the record continuously at the splice year and converging to modern background; all reported results use this continued forcing.

The integration grid over $\tau$ was resolution-graded rather than uniform: 0.25-year steps to 100 years, 1-year steps from 100 to 2,000 years, and 10-year steps beyond 2,000 years, with trapezoidal quadrature weights (`hydrosheaf/nuclear/joint_lpm.py`, `_integration_grid`). An earlier, uniform-step version of this grid, tied to the oldest supported tracer, forced every fit onto 10-year cells regardless of age, which quantised young-water and piston-flow predictions below 10 years; this was corrected prior to the M3 accuracy lock of 2026-07-28 and is unchanged in the present revision (Section S8).

### S2. Transit-time distribution normalisation

All TTDs were normalised so that:

$$
\int_0^{\tau_{\max}} g(\tau \mid \theta)\, \mathrm{d}\tau = 1. \tag{S3}
$$

For binary mixture models, the TTD was expressed as:

$$
g(\tau \mid \theta) = f\, g_1(\tau \mid \theta_1) + (1-f)\, g_2(\tau \mid \theta_2), \qquad 0 \le f \le 1. \tag{S4}
$$

This formulation allowed young and old components to contribute simultaneously to a sampled well, consistent with the Hydrosheaf design emphasis on age-fraction constraints and multi-modal residence-time structure. Supported LPM families were piston-flow (PFM), exponential (EM), dispersion (DM), gamma (GA), exponential-piston (EPM), partial-exponential (PEM), exponential-mixing (EMM), and binary mixtures (BMM-DM-DM, BMM-PEM-PEM), consistent with TracerLPM [@jurgens2012tracerlpm]. Candidate parameter grids used mean ages from 0.01 yr to the model-specific maximum age, dispersion values from 1e-4 to 2.0, gamma shape values from 0.1 to 10, piston fractions from 0.0 to 0.95, capture fractions from 0.05 to 1.0, and binary young-fraction values constrained to 0.001-0.999 after logit transformation. The default benchmark used 90 age-grid steps, deterministic grid search, optional local refinement, and AIC/AICc/BIC ranking.

### S3. Multi-tracer objective and age-fraction penalty

For each sample and scenario, supported tracer observations were fitted jointly:

$$
J_{\mathrm{tr}}(\theta) = \sum_{j=1}^{m} w_j \left[\frac{C_{j,\mathrm{obs}} - \widehat{C}_j(\theta)}{\sigma_j}\right]^2, \tag{S5}
$$

where $C_{j,\mathrm{obs}}$ is the observed tracer value, $\widehat{C}_j(\theta)$ the modelled value, $\sigma_j$ the tracer uncertainty, and $w_j$ the scenario-specific tracer weight. Where age-fraction constraints were activated, the objective was augmented by:

$$
J(\theta) = J_{\mathrm{tr}}(\theta) + \sum_k \left[\frac{F_{k,\mathrm{obs}} - F_{k,\mathrm{pred}}(\theta)}{\sigma_F}\right]^2, \tag{S6}
$$

where $F_k$ denotes a reported or target fraction for age domain $k$ (Anthropocene, Holocene, Pleistocene), evaluated from the TTD using fixed age cut-offs of 70 and 11,700 years.

Reference-free reportability criteria, applied before any reference-age metric was computed, required that: the optimiser converged; a requested reported-model emulation used an exactly supported LPM rather than a fallback; the number of fitted tracer observations exceeded the number of free parameters; standardised tracer-space RMSE did not exceed 2.0; and grid candidates within four objective units of the best candidate spanned no more than 0.5 log10 age units (approximately a factor of three). Reported ages were used only after fitting, to evaluate reportable estimates, and did not participate in the guard.

### S4. Graph-prior construction and penalty

Wells were nodes; hypothesised hydrogeological connections were directed edges. Graph families tested were coordinate-global (one nearest spatial neighbour per node, no uphill movement under an elevation proxy equal to negative depth), study-unit coordinate (the same rule within study-unit groups of at least five samples), depth-constrained (study-unit edges with downstream-minus-upstream depth between 5 and 200 m), hydraulic-proxy-constrained (the depth-constrained subset with edge distance below 25 km), parameter-smoothing (the hydraulic-proxy edge set weighted by similarity of AICc-derived LPM model weights), wrong-direction (every hydraulic-proxy edge reversed), and randomised (within-study-unit node pairs drawn with seed 42, matched in count to the study-unit coordinate graph).

The graph prior was evaluated at four pre-specified strengths, $\alpha \in \{0.00, 0.25, 0.55, 0.85\}$ (none, weak, medium, strong). For a directed edge $u \rightarrow v$ that violated the expected age ordering, the network update over eight iterations was:

$$
a_{v,\mathrm{new}} = (1-\alpha)\, a_v + \alpha\,(a_u + \delta_t), \tag{S7}
$$

where $\delta_t = 0$ in the reported benchmark. The corresponding penalty form is:

$$
J_G = \beta \sum_{(u,v) \in E} \min\left(0,\ a_v - a_u - \frac{L_{uv}}{v_h}\right)^2, \tag{S8}
$$

where $a_u, a_v$ are node ages, $L_{uv}$ is edge length, $v_h$ is inferred hydraulic velocity, and $\beta$ is the graph-prior strength.

### S5. Log-age performance metrics

Age prediction was evaluated on the logarithmic scale, with residuals

$$
e_i = \log_{10}(\widehat a_i) - \log_{10}(a_i), \tag{S9}
$$

median absolute log10 error

$$
\mathrm{MdAE}_{\log} = \mathrm{median}(|e_i|), \tag{S10}
$$

log10 RMSE

$$
\mathrm{RMSE}_{\log} = \sqrt{\frac{1}{n}\sum_{i=1}^n e_i^2}, \tag{S11}
$$

and within-factor agreement

$$
P_f = \frac{1}{n}\sum_{i=1}^n \mathbb{1}\!\left(|e_i| \le \log_{10} f\right), \tag{S12}
$$

with $f = 2$ and $f = 10$.

### S6. Set-valued (identified-TTD) formalism

The exploratory diagnostic reported in the main text (Results 3.5, Figure 7) replaced a single best-fit TTD with the question of which non-negative, unit-mass age-fraction vectors $w$ over a discretised age grid remain consistent with a well's tracer observations. For age bins with response kernel $A$ (Section S1's forward operators, evaluated per unit mass in each bin; `hydrosheaf/nuclear/joint_lpm.py`, `tracer_response_kernel`), a fold was **feasible** at sigma multiplier $k$ if there existed $w \ge 0$, $\sum_i w_i = 1$, satisfying

$$
|A w - C_{\mathrm{obs}}| \le k\,\sigma, \tag{S13}
$$

componentwise, tested as a linear feasibility problem (HiGHS solver). A fold's status was **IDENTIFIED** if every feasible $w$ agreed on a requested age-fractional quantity to within a declared tolerance, **PARTIALLY_IDENTIFIED** if only a bounded range of that quantity was feasible, and **ABSTAIN** where no probability-bearing ensemble was available to report a point estimate; **INFEASIBLE** described a fold with no feasible $w$ at all, i.e. no feasibility problem above admitted a solution.

Two complementary infeasibility diagnostics were run over the same rows and protocol (`configs/identified_ttd_protocol.yaml`, status `development`, claim authority `implementation_only`; `scripts/run_m3_infeasibility_diagnostics.py`). The **envelope** diagnostic tested, for tracer $i$, whether an observation's interval $[\mathrm{obs}_i - k\sigma_i,\ \mathrm{obs}_i + k\sigma_i]$ missed the achievable range $[\min_j A_{ij},\ \max_j A_{ij}]$ over the simplex, which no non-negative unit-mass $w$ could satisfy at any $k$. The **minimal infeasible subset** (MIS) diagnostic tested, for every well with at least three constraints, all singleton and pairwise tracer subsets for feasibility at $k=6.0$, and, for infeasible full panels, searched increasing subset sizes exhaustively for the smallest infeasible combination. Every row and fold excluded from either diagnostic was counted against an explicit, named skip reason and reconciled against the eligible-fold count; no row was silently discarded.

At the baseline sigma multiplier ($k = 1.96$, an independent 95% interval per tracer with no multiplicity adjustment), 975 of 3,501 eligible site-tracer folds (27.85%) were reported infeasible; of these, 974 were confirmed infeasible by the linear-programming solver and one was a solver error rather than a demonstrated infeasibility. A sigma sweep at $k \in \{1.96, 2.5, 3.0, 4.0, 5.0, 6.0\}$ compared observed infeasibility against the multiplicity null $1-(2\Phi(k)-1)^n$ for $n$ independent constraints; the observed rate at $k=6.0$ (14.1% overall; 51.2-65.9% at $n=4$-5) remained far above the null (effectively zero at $k \ge 4$ for $n \ge 4$), excluding both interval multiplicity and uniform interval under-dispersion as explanations. The envelope diagnostic accounted for at most 4.6% of folds with a violation at $k=1.96$ (1.2% at $k=6.0$), excluding individual out-of-range observations as the primary cause. The MIS diagnostic, at $k=6.0$ over 226 infeasible panels with at least three constraints, found 83.2% (188/226) reduced to a minimal infeasible pair, 8.4% (19/226) to a single out-of-envelope tracer, 7.5% (17/226) to a triple, and 0.9% (2/226) to a quadruple; this pairwise concentration argues against a multi-flow-path (no-common-TTD) explanation, which would predominantly produce infeasibility at subset size three or more.

Pairwise infeasibility rates at $k=6.0$ (Figure 7a) were highest among the three chlorofluorocarbons (CFC-11+CFC-12, 32.7%; CFC-11+CFC-113, 29.5%; CFC-113+CFC-12, 19.4%) and between tritium and its decay partner (3H+3H/3He, 17.9%), and lowest between the two most independent tracers (14C+3H, 0.8%).

### S7. Testing candidate explanations for the infeasibility signature

Seven candidate explanations for the pairwise, within-family infeasibility pattern (Section S6) were tested in sequence against the development dataset; none was supported (`docs/m3_identified_ttd_infeasibility_audit_20260731.md`, `docs/m3_cfc_reconciliation_step1_20260731.md`; Supplementary Table S4).

1. **Multiplicity of independent intervals** -- excluded: the observed infeasibility rate at $k=1.96$ (27.8% overall) exceeded the multiplicity null by more than an order of magnitude at every conditioning-set size.
2. **Uniform interval under-dispersion** -- excluded: widening every tracer's interval independently from $k=1.96$ to $k=6.0$ reduced overall infeasibility only from 27.8% to 14.1%, and the rise with conditioning-set size persisted essentially undiminished (5.7% to 56.5% at $k=6.0$, $n=1$ to $n=4$); a uniform under-statement of sigma would have collapsed the curve onto the null at some $k$, which did not occur at any tested value.
3. **Observations outside the achievable envelope** -- excluded: the envelope diagnostic accounted for only 1.2% of folds at $k=6.0$ against 14.1% infeasible overall.
4. **No common TTD across multiple flow paths** -- refuted: 83.2% of infeasible panels reduced to a minimal infeasible pair (Section S6), which a genuine multi-path failure would not predominantly produce.
5. **A refittable per-site recharge-temperature/excess-air correction** -- premise invalid: under the protocol's `gas_correction_mode: usgs_dgm`, corrected chlorofluorocarbon partial pressures are read directly from the USGS release (`Table_5_DissGas_ModOut.txt`); Hydrosheaf's own gas-correction routine is not invoked in this configuration, and none of the 107 sites with at least two measured CFCs carried a finite locally fitted recharge temperature.
6. **A common-mode shared scale error across the CFC family** -- refuted: scanning a single shared multiplicative scale $s \in [0.50, 2.00]$ over the 107 multi-CFC sites at $k=1.96$ reconciled only 19 of 68 already-infeasible sites (27.9%), with the reconciling scale clustering at the grid boundary (median 0.500; only 2 of 19 within +/-20% of unity), consistent with a degenerate old-water solution rather than a plausible correction.
7. **A single culprit tracer (chlorofluorocarbon-11)** -- refuted: among 97 sites measuring all three CFCs, 64 were jointly infeasible; dropping CFC-11 reconciled the remaining pair at 40 of 64 sites (62.5%), more than dropping CFC-113 (48.4%) or CFC-12 (45.3%), consistent with CFC-11's known anoxic degradation, but 37.5% remained infeasible even after removing the single most helpful tracer.

A predeclared eighth step, testing whether CFC-11-implicated infeasibility concentrates under reducing (low-oxygen) conditions, could not be executed: dissolved-oxygen, iron, manganese, and sulphide measurements do not co-occur with multi-CFC tracer measurements at the same wells and sampling times in the evaluated release (`docs/m3_redox_cfc_exclusion_protocol_20260731.md`). A follow-up synthetic control, generated on an independent MODFLOW/MODPATH model with sealed true topology and reused only as a method rather than as M3 evidence, likewise did not support a shared-nuisance-parameter mechanism reproducing the observed CFC-infeasibility signature under controlled conditions. The cause of the 27.85% baseline infeasibility rate is therefore reported as unresolved; no eighth or ninth explanation was available to test with current data, and no correction was applied to any reported scalar result on this basis.

### S8. Reproducibility: regeneration against the current codebase

For this revision, the full result-locking pipeline (`run_m3_manuscript_analysis.py --full --age-steps 90`) was rerun end to end against the current codebase in an isolated copy of the benchmark tree, rather than reusing the previously locked outputs. Every regenerated result file was compared line by line against the prior lock.

The design-matrix scenario summary (Table 2; 13 scenarios, 1,272 wells, 16,536 individual fits) was identical for every headline scenario used in this manuscript (strict reported-configuration parity, N=329; reported-output-constrained fraction sensitivity, N=289; Hydrosheaf model selection, N=309); three non-headline scenarios differed only in the fifteenth significant digit of the reported log10 R2 value, consistent with ordinary floating-point summation-order variation across parallel worker processes rather than a substantive change. All five tracer-withholding cross-validation result files (Figure 6/Supplementary Figure S3) were byte-for-byte identical to the prior lock. In the graph benchmark (Table 5/Figure 5), the hydraulic-proxy-constrained rows, the only family meeting the joint improvement criterion in either version, were byte-for-byte identical; other families showed small numeric drift (single-digit-percent in log10 RMSE and violation counts) attributable to the eight-iteration graph-regularisation update (Equation S7) compounding tiny floating-point ordering differences from the 16-worker parallel design-matrix computation across a connected graph, without changing any family's robust-improvement classification.

A code change identified during this review, the removal of a `reference_age` parameter from the tracer-reliability-weighting routine (`hydrosheaf/nuclear/multi_tracer.py`, `calculate_tracer_reliability_weights`) to close a path by which an evaluation target could otherwise reach tracer weighting, was hypothesised as a candidate source of numerical drift; the line-by-line comparison above shows it did not, in practice, alter any headline scenario on this dataset. This is disclosed as a verified negative finding rather than inferred from the code change alone.

## Supplementary Tables

**Supplementary Table S1.** Reproducibility inventory: source-data tables, variables, filtering rules, and row-count trace used to build the harmonised M3.1 benchmark table.

| Source | Version / DOI | Tables or files used | Variables extracted | Main filters and joins | Row-count trace |
|---|---|---|---|---|---|
| Jurgens et al. (2022) USGS national public-supply aquifers, 2004-2017 | Version 1.1; DOI 10.5066/P9W7T0DN | Table_1_Sites.txt; Table_2_Ages.txt; Table_3_Tracers.txt; Table_5_DissGas_ModOut.txt; Table_7_Carbon14.txt | SampleID, SampleDate, StudyUnit, AqGroup, coordinates, depth, 3H, 3He_trit, SF6, CFCs, 14C, 4He, DGMETA fields, LPM model and age-fraction fields | Inner join by SampleID; duplicate SampleID rows removed; 14C rows deduplicated by corrected and PHREEQC-preferred products; rows without sample year excluded | 1,279 source sites; 1,272 harmonised rows in the principal M3.1 benchmark |
| Faulkner and Jurgens (2019) Western Principal Aquifers, 2004-2018 | DOI 10.5066/P9U9ZSBN | USGS data-release groundwater-age and tracer tables | Tracer concentrations, calculated atmospheric-equivalent gases, age and metadata fields where available | Mapped to Hydrosheaf tracer schema; used for contextual parity/reproducibility checks where compatible | Source release reports 1,353 samples from 1,182 locations; compatible rows retained only when required target fields were finite |
| Gratzer et al. (2025) Mississippi River Valley alluvial aquifer, 2018-2020 | DOI 10.5066/P14DPCXE | Available site/tracer files; reported Table8_LPMs unavailable | Site/tracer context only; no comparable reported-LPM target | Excluded from parity, graph-age, and replication metrics | 0 rows in active performance tables |

**Supplementary Table S2.** Age-stratified reference-workflow agreement for strict reported-configuration emulation on the reportable subset. Age classes derive from model outputs and supports are unequal.

| Age class | N | Median \|log10 error\| | log10 RMSE | Within factor 2 (%) | Within factor 10 (%) | Median proxy coherence |
|---|---:|---:|---:|---:|---:|---:|
| Modern (<=50 yr) | 76 | 0.0890 | 0.3689 | 68.4 | 100.0 | 0.4764 |
| Intermediate (50-1,000 yr) | 85 | 0.0093 | 0.1467 | 90.6 | 100.0 | 0.3466 |
| Old (1,000-30,000 yr) | 158 | 0.0233 | 0.2869 | 94.3 | 97.5 | 0.3060 |
| Very old (>30,000 yr) | 10 | 0.0330 | 0.0626 | 100.0 | 100.0 | 0.3200 |

**Supplementary Table S3.** Diagnostic 14C-4He agreement categories for all harmonised rows. These categories are not independent age validation and do not restore the withdrawn hierarchical prior.

| Old-water status | N | Median 14C age (yr) | Median 4He age (yr) | Median constraint gap (log10) |
|---|---:|---:|---:|---:|
| Agreement | 176 | 5,515 | 6,972 | 0.083 |
| Conflict | 3 | 1,166 | 10,257 | 1.083 |
| Tension | 8 | 2,266 | 1,565 | 0.401 |
| Single tracer | 633 | 9,196 | 42,391 | -- |
| No constraint | 452 | -- | -- | -- |

**Supplementary Table S4.** Status of seven candidate explanations tested for the tracer-infeasibility signature (Section S7). No Supplementary Table S4 was reported in M3; it is added in M3.1 to support the new exploratory finding in Results 3.5.

| Explanation | Verdict | Key evidence |
|---|---|---|
| Multiplicity of independent intervals | Excluded | Observed infeasibility exceeds the multiplicity null by more than an order of magnitude at every conditioning-set size |
| Uniform interval under-dispersion | Excluded | Widening every interval independently from k=1.96 to k=6.0 reduces overall infeasibility only from 27.8% to 14.1%; the rise with conditioning-set size persists undiminished |
| Observations outside the achievable envelope | Excluded | Accounts for only 1.2% of folds at k=6.0 against 14.1% infeasible overall |
| No common TTD across multiple flow paths | Refuted | 83.2% of infeasible panels reduce to a minimal infeasible pair, not to higher-order combinations as a multi-path failure would predict |
| Refittable per-site recharge-temperature / excess-air correction | Premise invalid | Corrected CFC values are read directly from the USGS release under the frozen protocol; Hydrosheaf's own correction routine is not invoked |
| Common-mode shared scale error across the CFC family | Refuted | 72.1% of infeasible multi-CFC sites are unreconcilable at any shared scale; apparent successes cluster at the grid boundary |
| Single culprit tracer (CFC-11) | Refuted | Dropping CFC-11 helps most (62.5% reconciled) but 37.5% remain infeasible after removing the single most helpful tracer |

A predeclared eighth step (redox-stratified analysis) could not be executed because dissolved-oxygen and related redox indicators do not co-occur with multi-CFC tracer measurements in the evaluated release. The cause of the infeasibility rate is reported as unresolved.

## Supplementary Figures

![Supplementary Figure S2. Spatial-distance and vertical-difference distributions for graph edges in the N = 329 strict reportable-node benchmark.](manuscript/artifacts/figures/FIG-S2_edge_geometry.png)

**Supplementary Figure S2.** Spatial-distance and vertical-difference distributions for graph edges in the N = 329 strict reportable-node benchmark. These diagnose hypothesised edge geometry, not verified flow connections.

![Supplementary Figure S3. CFC-11 and CFC-12 target-withheld RMSE after the reportability guard.](manuscript/artifacts/figures/FIG-S3_cfc_withholding.png)

**Supplementary Figure S3.** CFC-11 (N = 4) and CFC-12 (N = 6) target-withheld RMSE after the reportability guard. No reportable CFC node formed an eligible graph edge, so identical bars indicate a non-estimable graph effect, not equivalence.

![Supplementary Figure S4. Controlled-synthetic network-dating ambiguity and recovery diagnostic.](manuscript/artifacts/figures/FIG-S4_network_dating_demo.png)

**Supplementary Figure S4.** Controlled-synthetic network-dating ambiguity and recovery diagnostic. This capability demonstration uses known simulated ages to test implementation behaviour; it is not field validation, does not establish accuracy for the USGS benchmark, and does not justify a general graph-improvement claim.

## References

Jurgens, B., Böhlke, J. K., & Eberts, S. M. (2012). *TracerLPM (Version 1): An Excel workbook for interpreting groundwater age distributions from environmental tracer data* (U.S. Geological Survey Techniques and Methods 4-F3). U.S. Geological Survey. https://doi.org/10.3133/tm4F3

IAEA/WMO. (2026). *Global Network of Isotopes in Precipitation. The GNIP Database.* https://nucleus.iaea.org/wiser
