## Supplementary Information

### S1. Residence-time convolution

For a tracer with atmospheric input history $I_j(t)$, decay constant $\lambda_j$, and residence-time distribution $g(\tau \mid \theta)$, the modelled concentration at sample time $t_s$ was:

<!-- EQ:EQ-S1 -->
$$
\widehat{C}_j(t_s, \theta) = \int_0^{\tau_{\max}} I_j(t_s - \tau)\, \exp(-\lambda_j \tau)\, g(\tau \mid \theta)\, \mathrm{d}\tau. \tag{S1}
$$

For stable atmospheric tracers such as SF~6~ and CFCs, $\lambda_j = 0$. For radioactive parent tracers such as ^3^H, ^14^C, ^39^Ar, and ^85^Kr, $\lambda_j = \ln(2) / T_{1/2,j}$. For daughter ingrowth in the ^3^H/^3^He system, the simplified ingrowth contribution was represented as:

<!-- EQ:EQ-S2 -->
$$
\widehat{C}_{^3He_{trit}}(t_s, \theta) = \int_0^{\tau_{\max}} I_{^3H}(t_s - \tau)\left[1 - \exp(-\lambda_{^3H}\tau)\right] g(\tau \mid \theta)\, \mathrm{d}\tau. \tag{S2}
$$

For carbon-14, the activity prediction was treated as:

<!-- EQ:EQ-S2b -->
$$
\widehat{A}_{14}(\theta) = A_0 \int_0^{\tau_{\max}} \exp(-\lambda_{14}\tau)\, g(\tau \mid \theta)\, \mathrm{d}\tau, \tag{S2b}
$$

where $A_0$ represents the initial or corrected ^14^C activity. For helium-4, the diagnostic accumulation model was represented as:

<!-- EQ:EQ-S2c -->
$$
\widehat{He}_4(\theta) = He_{4,\mathrm{bg}} + r_{He}\,\bar\tau, \qquad \bar\tau = \int_0^{\tau_{\max}} \tau\, g(\tau \mid \theta)\, \mathrm{d}\tau, \tag{S2c}
$$

where $He_{4,\mathrm{bg}}$ is background helium, $r_{He}$ is the accumulation-rate term, and $\bar\tau$ is mean residence time.

The tritium input function was represented using GNIP/WISER precipitation isotope data from the station nearest each site [@iaeawmo2026gnip]. Most North American WISER stations stop reporting between 1965 and 1999, whereas the benchmark samples were collected between 2004 and 2020. Each station record was therefore continued past its final observation along the regional post-bomb decline curve, rescaled to join the record continuously at the splice year and converging to modern background; all reported results use this continued forcing.

The integration grid over $\tau$ was resolution-graded rather than uniform: 0.25-year steps to 100 years, 1-year steps from 100 to 2,000 years, and 10-year steps beyond 2,000 years, with trapezoidal quadrature weights (`hydrosheaf/nuclear/joint_lpm.py`, `_integration_grid`). An earlier, uniform-step version of this grid, tied to the oldest supported tracer, forced every fit onto 10-year cells regardless of age, which quantised young-water and piston-flow predictions below 10 years; this was corrected prior to the M3 accuracy lock of 2026-07-28 and is unchanged in the present revision (Section S8).

### S2. Transit-time distribution normalisation

All TTDs were normalised so that:

<!-- EQ:EQ-S3 -->
$$
\int_0^{\tau_{\max}} g(\tau \mid \theta)\, \mathrm{d}\tau = 1. \tag{S3}
$$

For binary mixture models, the TTD was expressed as:

<!-- EQ:EQ-S4 -->
$$
g(\tau \mid \theta) = f\, g_1(\tau \mid \theta_1) + (1-f)\, g_2(\tau \mid \theta_2), \qquad 0 \le f \le 1. \tag{S4}
$$

This formulation allowed young and old components to contribute simultaneously to a sampled well, consistent with the Hydrosheaf design emphasis on age-fraction constraints and multi-modal residence-time structure. Supported LPM families were piston-flow (PFM), exponential (EM), dispersion (DM), gamma (GA), exponential-piston (EPM), partial-exponential (PEM), exponential-mixing (EMM), and binary mixtures (BMM-DM-DM, BMM-PEM-PEM), consistent with TracerLPM [@jurgens2012tracerlpm]. Candidate parameter grids used mean ages from 0.01 yr to the model-specific maximum age, dispersion values from $1\times10^{-4}$ to 2.0, gamma shape values from 0.1 to 10, piston fractions from 0.0 to 0.95, capture fractions from 0.05 to 1.0, and binary young-fraction values constrained to 0.001-0.999 after logit transformation. The default benchmark used 90 age-grid steps, deterministic grid search, optional local refinement, and AIC/AICc/BIC ranking.

### S3. Multi-tracer objective and age-fraction penalty

For each sample and scenario, supported tracer observations were fitted jointly:

<!-- EQ:EQ-S5 -->
$$
J_{\mathrm{tr}}(\theta) = \sum_{j=1}^{m} w_j \left[\frac{C_{j,\mathrm{obs}} - \widehat{C}_j(\theta)}{\sigma_j}\right]^2, \tag{S5}
$$

where $C_{j,\mathrm{obs}}$ is the observed tracer value, $\widehat{C}_j(\theta)$ the modelled value, $\sigma_j$ the tracer uncertainty, and $w_j$ the scenario-specific tracer weight. Where age-fraction constraints were activated, the objective was augmented by:

<!-- EQ:EQ-S6 -->
$$
J(\theta) = J_{\mathrm{tr}}(\theta) + \sum_k \left[\frac{F_{k,\mathrm{obs}} - F_{k,\mathrm{pred}}(\theta)}{\sigma_F}\right]^2, \tag{S6}
$$

where $F_k$ denotes a reported or target fraction for age domain $k$ (Anthropocene, Holocene, Pleistocene), evaluated from the TTD using fixed age cut-offs of 70 and 11,700 years.

Reference-free reportability criteria, applied before any reference-age metric was computed, required that: the optimiser converged; a requested reported-model emulation used an exactly supported LPM rather than a fallback; the number of fitted tracer observations exceeded the number of free parameters; standardised tracer-space RMSE did not exceed 2.0; and grid candidates within four objective units of the best candidate spanned no more than 0.5 log~10~ age units (approximately a factor of three). Reported ages were used only after fitting, to evaluate reportable estimates, and did not participate in the guard.

### S4. Graph-prior construction and penalty

Wells were nodes; hypothesised hydrogeological connections were directed edges. Graph families tested were coordinate-global (one nearest spatial neighbour per node, no uphill movement under an elevation proxy equal to negative depth), study-unit coordinate (the same rule within study-unit groups of at least five samples), depth-constrained (study-unit edges with downstream-minus-upstream depth between 5 and 200 m), hydraulic-proxy-constrained (the depth-constrained subset with edge distance below 25 km), parameter-smoothing (the hydraulic-proxy edge set weighted by similarity of AICc-derived LPM model weights), wrong-direction (every hydraulic-proxy edge reversed), and randomised (within-study-unit node pairs drawn with seed 42, matched in count to the study-unit coordinate graph).

The graph prior was evaluated at four pre-specified strengths, $\alpha \in \{0.00, 0.25, 0.55, 0.85\}$ (none, weak, medium, strong). For a directed edge $u \rightarrow v$ that violated the expected age ordering, the network update over eight iterations was:

<!-- EQ:EQ-S7 -->
$$
a_{v,\mathrm{new}} = (1-\alpha)\, a_v + \alpha\,(a_u + \delta_t), \tag{S7}
$$

where $\delta_t = 0$ in the reported benchmark. The corresponding penalty form is:

<!-- EQ:EQ-S8 -->
$$
J_G = \beta \sum_{(u,v) \in E} \min\left(0,\ a_v - a_u - \frac{L_{uv}}{v_h}\right)^2, \tag{S8}
$$

where $a_u, a_v$ are node ages, $L_{uv}$ is edge length, $v_h$ is inferred hydraulic velocity, and $\beta$ is the graph-prior strength.

### S5. Log-age performance metrics

Age prediction was evaluated on the logarithmic scale, with residuals

<!-- EQ:EQ-S9 -->
$$
e_i = \log_{10}(\widehat a_i) - \log_{10}(a_i), \tag{S9}
$$

median absolute log~10~ error

<!-- EQ:EQ-S10 -->
$$
\mathrm{MdAE}_{\log} = \mathrm{median}(|e_i|), \tag{S10}
$$

log~10~ RMSE

<!-- EQ:EQ-S11 -->
$$
\mathrm{RMSE}_{\log} = \sqrt{\frac{1}{n}\sum_{i=1}^n e_i^2}, \tag{S11}
$$

and within-factor agreement

<!-- EQ:EQ-S12 -->
$$
P_f = \frac{1}{n}\sum_{i=1}^n \mathbb{1}\!\left(|e_i| \le \log_{10} f\right), \tag{S12}
$$

with $f = 2$ and $f = 10$.

### S6. Set-valued (identified-TTD) formalism

The exploratory diagnostic reported in Results 3.5 and Figure 7 replaced a single best-fit TTD with the question of which non-negative, unit-mass age-fraction vectors $w$ over a discretised age grid remain consistent with a well's tracer observations. For age bins with response kernel $A$ (Section S1's forward operators, evaluated per unit mass in each bin; `hydrosheaf/nuclear/joint_lpm.py`, `tracer_response_kernel`), a fold was **feasible** at sigma multiplier $k$ if there existed $w \ge 0$, $\sum_i w_i = 1$, satisfying

<!-- EQ:EQ-S13 -->
$$
|A w - C_{\mathrm{obs}}| \le k\,\sigma, \tag{S13}
$$

componentwise, tested as a linear feasibility problem (HiGHS solver). A fold's status was **IDENTIFIED** if every feasible $w$ agreed on a requested age-fractional quantity to within a declared tolerance, **PARTIALLY_IDENTIFIED** if only a bounded range of that quantity was feasible, and **ABSTAIN** where no probability-bearing ensemble was available to report a point estimate; **INFEASIBLE** described a fold with no feasible $w$ at all, i.e. no feasibility problem above admitted a solution.

Two complementary infeasibility diagnostics were run over the same rows and protocol (`configs/identified_ttd_protocol.yaml`, status `development`, claim authority `implementation_only`; `scripts/run_m3_infeasibility_diagnostics.py`). The **envelope** diagnostic tested, for tracer $i$, whether an observation's interval $[\mathrm{obs}_i - k\sigma_i,\ \mathrm{obs}_i + k\sigma_i]$ missed the achievable range $[\min_j A_{ij},\ \max_j A_{ij}]$ over the simplex, which no non-negative unit-mass $w$ could satisfy at any $k$. The **minimal infeasible subset** (MIS) diagnostic tested, for every well with at least three constraints, all singleton and pairwise tracer subsets for feasibility at $k=6.0$, and, for infeasible full panels, searched increasing subset sizes exhaustively for the smallest infeasible combination. Every row and fold excluded from either diagnostic was counted against an explicit, named skip reason and reconciled against the eligible-fold count; no row was silently discarded.

At the baseline sigma multiplier ($k = 1.96$, an independent 95% interval per tracer with no multiplicity adjustment), 975 of 3,501 eligible site-tracer folds (27.85%) were reported infeasible; of these, 974 were confirmed infeasible by the linear-programming solver and one was a solver error rather than a demonstrated infeasibility. A sigma sweep at $k \in \{1.96, 2.5, 3.0, 4.0, 5.0, 6.0\}$ compared observed infeasibility against the multiplicity null $1-(2\Phi(k)-1)^n$ for $n$ independent constraints; the observed rate at $k=6.0$ (14.1% overall; 51.2-65.9% at $n=4$-5) remained far above the null (effectively zero at $k \ge 4$ for $n \ge 4$), excluding both interval multiplicity and uniform interval under-dispersion as explanations. The envelope diagnostic accounted for at most 4.6% of folds with a violation at $k=1.96$ (1.2% at $k=6.0$), excluding individual out-of-range observations as the primary cause. The MIS diagnostic, at $k=6.0$ over 226 infeasible panels with at least three constraints, found 83.2% (188/226) reduced to a minimal infeasible pair, 8.4% (19/226) to a single out-of-envelope tracer, 7.5% (17/226) to a triple, and 0.9% (2/226) to a quadruple; this pairwise concentration argues against a multi-flow-path (no-common-TTD) explanation, which would predominantly produce infeasibility at subset size three or more.

Pairwise infeasibility rates at $k=6.0$ (Figure 7a) were highest among the three chlorofluorocarbons (CFC-11+CFC-12, 32.7%; CFC-11+CFC-113, 29.5%; CFC-113+CFC-12, 19.4%) and between tritium and its decay partner (^3^H+^3^H/^3^He, 17.9%), and lowest between the two most independent tracers (^14^C+^3^H, 0.8%).

### S7. Testing candidate explanations for the infeasibility signature

Seven candidate explanations for the pairwise, within-family infeasibility pattern (Section S6) were tested in sequence against the development dataset; none was supported (`docs/m3_identified_ttd_infeasibility_audit_20260731.md`, `docs/m3_cfc_reconciliation_step1_20260731.md`; Supplementary Table S4).

1. **Multiplicity of independent intervals** — excluded: the observed infeasibility rate at $k=1.96$ (27.8% overall) exceeded the multiplicity null by more than an order of magnitude at every conditioning-set size.
2. **Uniform interval under-dispersion** — excluded: widening every tracer's interval independently from $k=1.96$ to $k=6.0$ reduced overall infeasibility only from 27.8% to 14.1%, and the rise with conditioning-set size persisted essentially undiminished (5.7% to 56.5% at $k=6.0$, $n=1$ to $n=4$); a uniform under-statement of sigma would have collapsed the curve onto the null at some $k$, which did not occur at any tested value.
3. **Observations outside the achievable envelope** — excluded: the envelope diagnostic accounted for only 1.2% of folds at $k=6.0$ against 14.1% infeasible overall.
4. **No common TTD across multiple flow paths** — refuted: 83.2% of infeasible panels reduced to a minimal infeasible pair (Section S6), which a genuine multi-path failure would not predominantly produce.
5. **A refittable per-site recharge-temperature/excess-air correction** — premise invalid: under the protocol's `gas_correction_mode: usgs_dgm`, corrected chlorofluorocarbon partial pressures are read directly from the USGS release (`Table_5_DissGas_ModOut.txt`); Hydrosheaf's own gas-correction routine is not invoked in this configuration, and none of the 107 sites with at least two measured CFCs carried a finite locally fitted recharge temperature.
6. **A common-mode shared scale error across the CFC family** — refuted: scanning a single shared multiplicative scale $s \in [0.50, 2.00]$ over the 107 multi-CFC sites at $k=1.96$ reconciled only 19 of 68 already-infeasible sites (27.9%), with the reconciling scale clustering at the grid boundary (median 0.500; only 2 of 19 within $\pm 20\%$ of unity), consistent with a degenerate old-water solution rather than a plausible correction.
7. **A single culprit tracer (chlorofluorocarbon-11)** — refuted: among 97 sites measuring all three CFCs, 64 were jointly infeasible; dropping CFC-11 reconciled the remaining pair at 40 of 64 sites (62.5%), more than dropping CFC-113 (48.4%) or CFC-12 (45.3%), consistent with CFC-11's known anoxic degradation, but 37.5% remained infeasible even after removing the single most helpful tracer.

A predeclared eighth step, testing whether CFC-11-implicated infeasibility concentrates under reducing (low-oxygen) conditions, could not be executed: dissolved-oxygen, iron, manganese, and sulphide measurements do not co-occur with multi-CFC tracer measurements at the same wells and sampling times in the evaluated release (`docs/m3_redox_cfc_exclusion_protocol_20260731.md`). A follow-up synthetic control, generated on an independent MODFLOW/MODPATH model with sealed true topology and reused only as a method rather than as M3 evidence, likewise did not support a shared-nuisance-parameter mechanism reproducing the observed CFC-infeasibility signature under controlled conditions. The cause of the 27.85% baseline infeasibility rate is therefore reported as unresolved; no eighth or ninth explanation was available to test with current data, and no correction was applied to any reported scalar result on this basis.

### S8. Reproducibility: regeneration against the current codebase

For this revision, the full result-locking pipeline (`run_m3_manuscript_analysis.py --full --age-steps 90`) was rerun end to end against the current codebase in an isolated copy of the benchmark tree, rather than reusing the previously locked outputs. Every regenerated result file was compared line by line against the prior lock.

The design-matrix scenario summary (Table 2; 13 scenarios, 1,272 wells, 16,536 individual fits) was identical for every headline scenario used in this manuscript (strict reported-configuration parity, N=329; reported-output-constrained fraction sensitivity, N=289; Hydrosheaf model selection, N=309); three non-headline scenarios differed only in the fifteenth significant digit of the reported log~10~ R² value, consistent with ordinary floating-point summation-order variation across parallel worker processes rather than a substantive change. All five tracer-withholding cross-validation result files (Figure 6/Supplementary Figure S3) were byte-for-byte identical to the prior lock. In the graph benchmark (Table 4/Figure 5), the hydraulic-proxy-constrained rows, the only family meeting the joint improvement criterion in either version, were byte-for-byte identical; other families showed small numeric drift (single-digit-percent in log~10~ RMSE and violation counts) attributable to the eight-iteration graph-regularisation update (Equation S7) compounding tiny floating-point ordering differences from the 16-worker parallel design-matrix computation across a connected graph, without changing any family's robust-improvement classification.

A code change identified during this review, the removal of a `reference_age` parameter from the tracer-reliability-weighting routine (`hydrosheaf/nuclear/multi_tracer.py`, `calculate_tracer_reliability_weights`) to close a path by which an evaluation target could otherwise reach tracer weighting, was hypothesised as a candidate source of numerical drift; the line-by-line comparison above shows it did not, in practice, alter any headline scenario on this dataset. This is disclosed as a verified negative finding rather than inferred from the code change alone.
