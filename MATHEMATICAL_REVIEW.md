# Hydrosheaf — Mathematical & Conceptual Review

**Date:** 2026-05-31
**Scope:** Full `hydrosheaf/` package — inference, sheaf/topology, nuclear/temporal/uncertainty, chemistry/isotopes/causal, vadose/physics/PHREEQC/calibration.
**Method:** Five parallel domain reviews, each reading the actual implementations; the highest-impact claims were independently re-verified against source before inclusion.

> **Headline:** The package is conceptually ambitious and broad, but several of its *flagship* outputs are computed by code that does not implement the math it advertises. The four issues that most undermine scientific credibility are: (1) the transport-model decision ignores the transport fit, (2) the Bayesian topology sampler does not target its stated posterior, (3) the "BCa" confidence intervals are not BCa (and the production path skips bias-correction entirely), and (4) distributed groundwater-age models are biased by truncate-then-renormalize. None of these are cosmetic — they change the numbers users publish.

---

## How to read this report

Each finding has a stable ID (`AREA-n`), a severity, a precise location, what is mathematically wrong, why it matters scientifically, and a concrete fix. Severities:

- **Critical** — produces wrong scientific results that look plausible; fix before any further use in manuscripts.
- **High** — biased/uncalibrated results or invalid diagnostics; fix soon.
- **Medium** — narrower bias, silent failure modes, or misleading labels.
- **Low** — latent traps, dead code, documentation/units hygiene.

Verified findings (re-read against source by the lead reviewer) are marked ✅.

---

## 1. Critical findings

### CORE-1 ✅ — Transport-model selection ignores the transport fit entirely
**`inference/edge_fit.py:323-325, 361, 366-373`**

The per-edge `objective` (line 323) is built **only** from the chemistry LASSO fit:
`chem_fit.residual_norm + λ·l1 + penalties`. It never adds the *transport* misfit. Yet this objective decides the winning transport model (`best_result`, line 361) and the reported `transport_probabilities` (lines 366-373). Because evaporation fits a free scalar γ over `[-0.99, 999]` while mixing fits a constrained `f∈[0,1]`, evaporation can absorb common-mode variance "for free" and leave a smaller chemistry residual — so it wins on a quantity that never charged it for a bad transport fit. (Line 325 compounds the confusion: it stores `chem_fit.residual_norm` under the key `transport_residual_norm`.)

**Impact:** The package's central per-edge claim — *was this evaporation or mixing?* — is decided by a biased score. Evap is systematically favored; the reported probabilities are not posterior model probabilities.

**Fix:** Use one objective that includes the transport residual and charges transport degrees of freedom (AIC/BIC term for evap's extra free parameter): `objective = chem_residual + transport_residual + λ·l1 + penalties + complexity(k)`. Use the *same* objective for selection and for the probability weights.

---

### CORE-2 ✅ — "Transport probabilities" are an un-scaled softmax → unit-dependent, uncalibrated
**`inference/edge_fit.py:366-373`**

`weights = exp(-(score - min_score))` with no noise/variance scaling. `score` is a weighted RSS in concentration² units. In mg/L the differences are huge (one model → prob ≈ 1 always); in mmol/L they're tiny (≈ uniform). The "probability" is an artifact of the unit system and the absolute weight magnitudes.

**Fix:** Scale by an estimated noise variance, `exp(-(Δscore)/(2σ²))`, or use proper AIC weights `w_i ∝ exp(-ΔAIC_i/2)` (the `ReactionFit.aic` field already exists).

---

### CORE-3 ✅ — Evaporation model permits dilution (γ<1) and applies one conservative-ion factor to all ions
**`inference/edge_fit.py:227-246`**

Evap is `Δ = δγ·x_u` fit with `lb=[-0.99]`, so `γ = 1+δγ` can be 0.01 (i.e. **dilution**, not evaporation). The single δγ is fit against the conservative-weight vector (≈ Cl only) and then force-subtracted from **every** ion's target, treating reactive ions as if they evaporate proportionally. A correct, clamped `fit_evaporation` exists in `models/mixing.py:26-42` but is **not** used on this path.

**Fix:** Enforce `γ ≥ 1` (`lb=[0.0]`), fit γ on conservative tracers only, and reuse/route through `models/mixing.fit_evaporation`. Encode which ion is the conservative tracer rather than relying on a hard-coded weight vector.

---

### TOPO-1 ✅ — Metropolis–Hastings omits the Hastings proposal-asymmetry correction → chain does not target the stated posterior
**`inference/topology_posterior.py:171-199` (propose), `283-285` (accept)**

`_propose_flip` is state-dependent and asymmetric: P(remove e) = `0.5/|G|`, P(add e) = `0.5/(|U|−|G|)`, and the split collapses to 1.0 at empty/full boundaries. The acceptance test uses only `proposed_logp − current_logp`, assuming a symmetric proposal. The Hastings ratio `q(G|G')/q(G'|G)` is missing.

**Impact:** Detailed balance is violated; the stationary distribution is **not** `p(G|data)`. Every headline output — edge-inclusion probabilities, log-odds, MAP topology, entropy, edge-count CI — is biased.

**Fix:** Add `log_ratio += log|current| − log(|U|−|current|+1)` (with boundary care), **or** switch to a symmetric single-bit-flip proposal (pick a uniform edge from the universe, toggle it) which needs no correction. *(Flagged independently by both the inference and topology reviewers.)*

---

### TOPO-2 ✅ — Posterior marginals divided by the wrong sample count (noop/rejected handling)
**`inference/topology_posterior.py:273-304, 324`**

Post-burn-in iterations that hit `label == "noop"` execute `continue` (lines 276-277) **before** recording, so fewer than `n_samples` states are tallied — but `edge_counts`, `n_edges_trace`, and `acceptance_count` are all normalized by the fixed `n_samples` (lines 301-303, 324). Inclusion probabilities are biased low whenever the chain visits the empty/full graph. Separately, a rejected M-H proposal must still re-count the *current* state as a sample; `continue` drops it.

**Fix:** Maintain an explicit `recorded` counter incremented exactly when a post-burn-in state is tallied, and divide all marginals/acceptance by it. Do **not** `continue` on noop/reject — record the unchanged current state.

---

### NUC-1 ✅ — "BCa" intervals are not BCa; the production path skips bias-correction entirely
**`uncertainty/bootstrap.py:277-334` (esp. 315-318); `:150-168`**

`compute_bca_ci` hard-codes acceleration `a = 0.0` ("requires jackknife, not implemented"), reducing BCa to bias-corrected percentile. The jackknife `a = Σ(θ̄−θ_i)³ / [6(Σ(θ̄−θ_i)²)^{3/2}]` is never computed — exactly the correction that matters for the skewed, bounded LASSO/WLS estimators here. Worse, `bootstrap_edge_fit` (lines 150-168) doesn't even call `compute_bca_ci`; it uses raw `np.percentile(...,[2.5,97.5])`, so production CIs have **no** bias correction at all.

**Fix:** Implement jackknife acceleration and route `bootstrap_edge_fit`'s intervals through the real BCa function. Handle degenerate `z0` (NUC-M1) with a continuity correction and a flag.

---

### NUC-2 — Distributed-age TTDs are truncate-then-renormalized → ages biased young
**`nuclear/joint_lpm.py:287-358`; `nuclear/lpm.py:79-113`**

Exponential/EPM/PEM/dispersion transit-time PDFs are evaluated on a finite grid (`max_tau = max(100, 5·T)`) and **renormalized after** decay/survival are applied. Renormalizing the truncated, already-decayed mass over-weights short transit times, so inferred concentrations come out high and apparent ages young. 5 mean lives captures only ~99.3% of an exponential — and the missing tail still carries decay weight that renormalization misallocates.

**Fix:** Normalize the **undecayed** TTD once over the full support (integrate analytically or to ≥10–15 mean lives), *then* apply decay. Use the published Kreft–Zuber dispersion TTD (see NUC-M5).

---

### NUC-3 — Reactive-TTD convolution ignores the time-varying input (collapses to a constant)
**`temporal/convolution.py:6-43`**

`convolve_reactive_ttd` computes `Σ_i w_i · c_in · e^{−kτ_i}` using a **single scalar** `c_in` for all lags — the lagged input history is never sampled. A TTD convolution requires `Σ_i w_i · C_in(t−τ_i) · e^{−kτ_i}`. The result is just "constant × attenuation," defeating the purpose. `lpm.py:114-124` already does the correct lagged sampling.

**Fix:** Pass the input time series and sample `C_in(t−τ_i)` per lag, mirroring `lpm.py`.

---

### CHEM-1 — "Optimal transport" is not OT: incommensurable cost scales, non-conservative cross-species flow, non-metric output
**`models/optimal_transport.py:29-70, 101-319`**

Off-diagonal transport costs are hand-assigned integers (0,1,2,5,10) by charge-sign heuristics; creation/destruction costs are a *separate* arbitrary scale (0.5/2/5) summed directly with them (line 279). The LP lets any upstream ion convert to any downstream ion (e.g. Ca→Cl) at finite cost — no element/charge/mass conservation. Marginals are raw concentrations, normalized by `max(total_up,total_down)` (neither balanced nor standard unbalanced-OT). The output is a heuristic penalty, not a distance, and is asymmetric. Additionally `c_create == c_destroy` (lines 178-179) lets the solver invent or annihilate species at equal cost.

**Fix:** Either restrict flows to chemically meaningful (stoichiometric, element/charge-conserving) transformations, or implement a true (unbalanced) OT with a single regularizer and a physically grounded ground metric on a conserved quantity. Don't sum two unrelated cost scales; make creation ≠ destruction.

---

### CHEM-2 — Nitrate mixing variance uses `Σfᵢσᵢ²` (should be `Σfᵢ²σᵢ²`); analytical and MCMC models disagree
**`models/nitrate_isotopes_mcmc.py:230-241, 430-443`; mirrored in `models/nitrate_isotopes.py:115-147`**

For a deterministic mixture mean `Σfᵢμᵢ` with uncertain source means, the predictive variance is `Σfᵢ²σᵢ²` (squared weights), not `Σfᵢσᵢ²`. Since `Σfᵢσᵢ² ≥ Σfᵢ²σᵢ²`, the likelihood is too flat → source fractions under-identified, manure-vs-fertilizer discrimination weakened, posterior-predictive QC desensitized. The analytical path uses yet another covariance, so the two "models" answer different questions.

**Fix:** Use `mix_var = Σfᵢ²σᵢ² + σ_resid²` consistently in both the MCMC and analytical likelihoods.

---

### PHYS-1 ✅ — PHREEQC dissolution is gated on *both* endpoints undersaturated (excludes dissolution-to-saturation)
**`phreeqc/constraints.py:80-87`**

Dissolution (`lb=0, ub=∞`) is permitted only when `si_u < −τ AND si_v < −τ`. The common, physically expected case — undersaturated upstream dissolving *to* saturation downstream (`si_u<0, si_v≈0`) — is wrongly forbidden, forcing `free`/`precipitation_only`. Precipitation is gated only on `si_v > τ` (downstream), an inconsistent asymmetry.

**Fix:** Gate dissolution on the **upstream** state (`si_u < −τ → dissolution_only`); gate precipitation on downstream supersaturation (`si_v > τ`); leave `free` when `si_u<−τ` and `si_v>τ`.

---

### PHYS-2 ✅ — Darcy residence time mixes three inconsistent gradient definitions; can divide head drop by vertical separation
**`inference/network_fit.py:82-97`; `graph3d/build_3d.py:313-329`**

`horizontal_gradient = Δh/d_xy` and `vertical_gradient = Δh/d_z` use the *same* head drop over *different* denominators — they are not orthogonal components, just three scalars. Selecting `vertical_gradient` for a near-horizontal edge (tiny `d_z`) inflates `i` by orders of magnitude, collapsing `τ` toward zero and corrupting every downstream kinetic/TTD prior.

**Fix:** Use one definition, `i = |Δh| / distance_3d_m` (head drop over actual flow-path length). Remove the `d_z`-only fallback.

---

## 2. High-severity findings

### Inference / chemistry objective
- **CORE-4 ✅ `edge_fit.py:319-321`** — `chemistry_r2` compares the post-transport anomaly RSS against the *raw downstream* SST (different response variables, mixed weighting, unweighted `mean_v`); can be <0 or >1. Compute R² on one consistent, identically-weighted response.
- **CORE-5 `edge_fit.py:220-225`** — `compositional_weighting` sets WLS weights to `1/x_v²` (the **response**), making weights endogenous and biasing toward small-concentration ions. The "relative error" comment doesn't match the implementation. Weight by `1/σ²` from a measurement-error model, or iterate weights on the fitted value.
- **CORE-6 `models/gibbs.py:79-96,132-150`** — Gibbs class "probabilities" are products of independent sigmoids that don't sum to 1; `1−p_evap` then penalizes the objective with an uncalibrated amount. Normalize over the three classes (softmax) before using as a penalty.
- **CORE-7 `models/reactions.py:127-137`, `ec_tds.py:13-19`** — Reaction stoichiometry assumes mmol/L (and meq for 1:2 exchange `Ca:1,Na:−2`) but EC/TDS/Gibbs thresholds assume mg/L-scale sums; no unit guard and no electroneutrality constraint in the fit (only a post-hoc flag). Assert/convert units in `Config`, use equivalents for charge-balanced reactions, add an electroneutrality residual.

### Sheaf / topology
- **SHEAF-1 ✅ `sheaf/cohomology.py:182-184`** — `h1_dim = n_edge_dofs − rank_D` equals the true first Betti number only when α=1; with affine α≠1 maps it is data-dependent slack, not a topological invariant. Compute the cycle-space dimension `dim·(E−V+C)` structurally and report affine inconsistency separately.
- **SHEAF-2 ✅ `sheaf/cohomology.py:92-103,186-199`** — "Obstruction energy" is **not orientation-invariant** when α≠1 (flipping an edge scales the residual by 1/α²), and it folds the *prior confidence* weight into the least-squares **precision**, so low-confidence wrong edges are rewarded with lower inconsistency. Normalize rows / use a direction-invariant (holonomy or log-space) form; separate measurement precision from prior confidence.
- **SHEAF-3 `sheaf/cohomology.py:209-296`** — Leverage and cycle-obstruction recompute least-squares over *different* node universes (removing a degree-1 edge drops a node) and solve cycles in isolation (under-determined → understates closure failure). Fix the node universe across subproblems; compute cycle closure as affine holonomy around the loop.
- **TOPO-3 ✅ `topology_posterior.py:149-159` + `graph/build.py:256-269`** — The flow-direction probability `Φ(Δh/σ)` (a P(downhill)) is reused as the Bernoulli **prior of edge existence**, and edges are treated as independent despite sharing nodes/heads. Distinguish direction confidence from inclusion prior; use a structured (per-source multinomial) prior.
- **TOPO-4 ✅ `topology_posterior.py:309-313`** — "Posterior entropy over graphs" is actually `Σ` marginal Bernoulli entropies, which upper-bounds (overstates) the true joint entropy under correlation. Rename, or estimate joint entropy from sampled bitsets.

### Nuclear / temporal / uncertainty
- **NUC-4 `temporal/residence_time.py:280,300-317`** — Cross-correlation lag is a one-sided, correlation-weighted **centroid** (not the peak) over lags ≥0; asymmetric correlation curves bias τ high, and negative lags can't be detected. Use parabolic-interpolated `argmax`, scan negative lags, derive uncertainty from peak curvature.
- **NUC-5 `nuclear/joint_lpm.py:409`** — Pre-1950 tritium input set to 0 TU (should be ~3–5 TU natural background); three different defaults (0.0 / 0.5 / 5.0) exist across the codebase. Standardize one documented pre-bomb background.
- **NUC-6 `nuclear/network_aging.py:177-178`** — The "network flow constraint" is `Potential(−20.0·Σ relu(...)²)` with a fixed unitless coefficient unrelated to the age scale — an uncalibrated tuning knob masquerading as physics. Replace with a hierarchical prior `age_v ~ Normal(age_u + L/v, σ_process)`.
- **NUC-7 `uncertainty/bayesian.py:417-458`** — Hand-rolled ESS truncates the autocorrelation sum at the first `ρ_k<0.05`, overestimating ESS and giving false convergence confidence. Use `arviz.ess` (Geyer initial-positive-sequence).
- **NUC-8 `nuclear/joint_lpm.py:411-414,471-480`** — 3H/³He "daughter" ingrowth is convolved separately then fed to a ratio-age formula valid only under piston flow; for distributed flow the ratio of integrals ≠ TTD-weighted age. Fit 3H and ³He as separate observables, or restrict the ratio-age to PFM.

### Chemistry / isotopes / causal
- **CHEM-3 `models/latent.py:92-128`** — PCA "latent endmembers" are built from a single right-singular vector (not a real sample projection or vertex), overshoot by an ad-hoc ×1.1, then a closed composition is rescaled by an unrelated sample's raw ion-sum; PCA sign/scale ambiguity unresolved. Result: endmembers can be chemically impossible and need not bracket the data. Use actual extreme samples, or reconstruct `clr = mean + U[k,i]·S[i]·Vt[i]` with consistent TDS; enforce nonnegativity/closure.
- **CHEM-4 `models/nitrate_isotopes.py:337-345`** — Nitrification expected δ¹⁸O hard-codes a 2:1 water:O₂ ratio and then **up-weights whichever source mean sits near the nitrification line**, injecting source-specific bias unrelated to the sample. Make the ratio/O₂-δ¹⁸O configurable with uncertainty; apply the constraint to the sample likelihood, not per-source mean alignment.
- **CHEM-5 `null_models/{chemistry,endmembers,lithology}.py`** — None of the "null models" builds a null distribution; they return weighted sums of magic constants (e.g. `+0.3`, `+0.2`, fixed 0.5/0.7) with hard thresholds and no p-value or multiple-testing control. Build genuine permutation nulls with the `(1+#≥obs)/(1+N)` p-value and BH/FDR across edges.
- **CHEM-6 `causal/discovery.py:125-230`** — The "partial correlation" computes correlations **across ion species** (composition shape), not temporal precedence; the returned p-value is the **mean of per-species p-values** (not a valid combined test); negative lag correlations are discarded. Use residualized partial correlation / Granger or transfer entropy with a real significance test; combine p-values with Fisher's method.
- **CHEM-7 `models/nitrate_isotopes.py:115-147,181-192`** — Posterior-predictive QC uses the mixture-distribution variance (`Σpᵢσᵢ² + Var(μ)`) instead of the mixture-mean variance (`Σpᵢ²σᵢ²`), so the tail probability is too large and `qc_posterior_predictive_mismatch` rarely fires (bad fits pass). Same root cause as CHEM-2.

### Vadose / physics / calibration
- **PHYS-3 `vadose/richards1d.py:436-443`** — Column storage `Σθ·dz` uses full `dz` for the two boundary half-cells, an O(dz) mass-balance error reported as solver error (flags good runs). Use trapezoidal weights (½·dz at ends) consistent with the flux discretization.
- **PHYS-4 `vadose/ttd.py:154-213`** — `summary` moments (→ `tt_mean_days` → edge residence time) are computed from the advective matrix only, while the returned kernel includes the preferential fast path; the two disagree whenever preferential flow > 0. Compute moments from the final mixed PDF.
- **PHYS-5 `vadose/ttd.py:35-83` + `network_fit.py:97`** — Advective piston transit time is labeled "residence time"/"mean age" interchangeably and uses **total** porosity (not effective/mobile) with no tortuosity, overestimating τ and relaxing the kinetic feasibility filter. Use effective porosity + tortuosity; distinguish transit time from mean age.
- **PHYS-6 `physics/kinetic_limit.py:23-77`** — `β_max = k·A·τ` uses the far-from-equilibrium rate for the whole residence time, ignoring the affinity term `(1−10^SI)` that → 0 near saturation, so the feasibility filter is far too permissive near equilibrium. Integrate the affinity-dependent rate or cap at the equilibrium-attainable mass.
- **CAL-1 `calibration/glm.py:96-116,279-298`** — Covariance double-counts the weighting: residuals are pre-weighted by `1/σ` (so `JᵀJ` ≈ inverse covariance), then multiplied **again** by `σ²=φ/dof`; and `dof = n_obs − n_par` ignores the Tikhonov rows present in the Jacobian. Pick one convention; compute covariance from the measurement block (or add the prior precision explicitly).

---

## 3. Medium-severity findings (condensed)

| ID | Location | Issue | Fix |
|---|---|---|---|
| CORE-M1 | `models/reactions.py:201-223` | AIC/BIC use **weighted** RSS inside `log(RSS/n)` and count k by a `1e-6` magnitude threshold (not LASSO DOF); n≈11 makes AICc unstable | Use unweighted/σ-normalized RSS; document DOF approximation |
| CORE-M2 | `models/reactions.py:276-282` | Non-negativity applied **before** user `lb/ub`; L1 stays centered at 0 even when feasible set excludes 0 | Derive effective `[lb,ub]` once; apply L1 relative to boundary |
| CORE-M3 | `inference/network_fit.py:101-110` | τ-uncertainty drops the `CV_L` term, but `τ ∝ L²/Δh`, so the largest term is omitted | Propagate `(2·CV_L)² + CV_Δh²` |
| TOPO-M1 | `topology_posterior.py:316-320` | `n_edges_ci95` uses floor-truncated index for both bounds (interval too narrow) | `np.percentile(...,[2.5,97.5])` |
| TOPO-M2 | `topology_posterior.py:387-443` | Sheaf-fail penalty `1e6` is unitful vs `β`; additive costs scale with `|G|`, an implicit sparsity penalty competing with `edge_penalty` | Relative/explicit rejection; normalize size-scaling |
| TOPO-M3 | `sheaf/topology_refine.py:276-303` | `pi_evap` multiplied by ad-hoc `×0.5` gates then clamped to `[0,0.8]` — not a probability calculus | Combine adjustments in logit space |
| NUC-M1 | `uncertainty/bootstrap.py:315` | BCa `z0` collapses to 0 when all replicates one-sided (degenerate bootstrap hidden) | Continuity correction `(n_less+0.5)/(n_boot+1)` + flag |
| NUC-M2 | `uncertainty/bootstrap.py:80-90` | Bootstrap fixes model selection to the original fit → CIs exclude transport/endmember selection uncertainty (overconfident) | Re-select per replicate or document conditionality |
| NUC-M3 | `nuclear/old_groundwater.py:224,308-317` | "Hierarchical" A₀ is a fixed 75/25 blend, not precision-weighted pooling | Weight by precisions or sample A₀ hierarchically |
| NUC-M4 | `temporal/temporal_edge_fit.py:209-265` | Seasonal model = 1 harmonic; `residual_std` uses `ddof=0` with 4 params; aliasing on irregular sampling | ≥2 harmonics; `ddof=n−p`; check condition number |
| NUC-M5 | `nuclear/lpm.py:28-40`, `joint_lpm.py:291-298` | Dispersion PDF not the standard Kreft–Zuber Green's function; relies on post-hoc renorm | Use published DM TTD, normalize undecayed PDF analytically |
| CHEM-M1 | `nitrate_isotopes_mcmc.py:232,434-458` | `mix_cov` not guaranteed PD; invalid MvNormal silently falls back to **uniform 50/50** posing as a result | Clip cov to ±0.999·√(varₙvarₒ); surface fallback as non-converged |
| CHEM-M2 | `nitrate_isotopes.py:311-322` | Denitrification slope measured from source mean (conflates mixing displacement with fractionation); contested 0.65; hard `dx_n≤0` cutoff | Anchor at pre-denit signature; configurable slope; soft directional weight |
| CHEM-M3 | `isotopes.py:160-188` | `craig_gordon_enrichment` is a no-op (`return delta_L`) | Implement or raise `NotImplementedError` |
| CHEM-M4 | `nitrate_source_v2.py:169-235` | Naive-Bayes LLR sums correlated ratios (Cl shared across NO3/Cl, PO4/Cl) → overconfident | Use ilr/balance coordinates or joint covariance |
| PHYS-M1 | `reactive_transport/rate_laws.py:35-67` | Carbonate/sulfate TST allows precipitation at the **dissolution** rate constant (no `SI>0` guard like silicates) | Separate precipitation rate / inhibition factor |
| PHYS-M2 | `reactive_transport/metrics.py:50-96` | RMSE/NSE/R² computed across the ion vector → dominated by the largest ion | Normalize per ion before aggregating |
| PHYS-M3 | `reactive_transport/metrics.py:297-324` | Damköhler divides by equilibrium **concentration** not distance-to-equilibrium | Use driving-force scale `|C₀−C_eq|` or `Da=k·τ` |
| PHYS-M4 | `vadose/nitrate.py:179-233` | Convolution assumes lag index == τ-grid index (breaks if TTD grid not zero-origin/uniform) | Assert uniform-from-zero or map `round(τⱼ/dt)` |
| CAL-M1 | `calibration/glm.py:170-180` | Forward-difference Jacobian can probe outside bounds; small `|x|` → noisy derivatives | Central differences, reflect probe inward at bounds |

---

## 4. Low-severity / hygiene

- **CORE-L1 ✅ `edge_fit.py:317`** — `residence_time_days or residence_time_days` is a copy-paste no-op (intended a fallback). Verify intended fallback value.
- **NUC-L1 `nuclear/lpm.py:8-17`** — `piston_flow_model` returns all zeros (not a delta); dead stub that could be mistaken for a PDF.
- **NUC-L2 `nuclear/nuclides.py:18,28`** — `decay_constant_per_day` stored but unused; latent per-day/per-year unit trap. Remove or document.
- **CHEM-L1 `isotopes.py:7-8,100-116`** — `d_excess` hard-codes slope 8 while the module fits a local MWL slope ≠ 8; reference lines inconsistent across modules.
- **CHEM-L2 `nitrate_isotopes.py:431-468`** — MVN density classifier lets endmembers with artificially small reported σ dominate apportionment; enforce a σ floor.
- **TOPO-L1 `sheaf/isotope_metrics.py:102-108`** — Evap "probability" weights sum to 1.1 and double-apply the depth gate; relabel as heuristic score.
- **PHYS-L1 `physics/kinetic_limit.py:46`** — Single constant reactive surface area for all minerals/edges (the dominant kinetic uncertainty).
- **PHYS-L2 `physics/kinetic_limit.py:45`** — Feasibility filter assumes isothermal 25 °C; the (correct) Arrhenius correction is never wired in.
- **CAL-L1 `calibration/benchmark_bootstrap.py:60-78`** — Precision/F1 computed over labeled edges only; label clearly as conditional.

### Confirmed correct (no action)
- CoDA `coda_sbp.py`: valid SBP tree, correct `√(rs/(r+s))` ilr scaling, correct clr geometric mean and softmax inverse. *(The problem is that other modules abandon this correct machinery — see CHEM-4, CHEM-M4.)*
- Arrhenius **sign** `rate_laws.py:267`, R=8.314, Kelvin conversion — correct (reviewer self-corrected an initial false flag).
- Gamma-kernel TTD parameterization `ttd.py:99-117`; free-drainage Darcy flux sign `richards1d.py:423-433`; PHREEQC SI sign convention in `metrics.check_thermodynamic_consistency`.

---

## 5. Cross-cutting root causes

Most findings trace to four recurring patterns. Fixing the *patterns* prevents regressions better than patching sites individually.

1. **"Scores" dressed as "probabilities."** CORE-2/6, TOPO-3/4, CHEM-5/6, TOPO-M3, NUC-L… repeatedly take an unnormalized heuristic, sigmoid/exp it, and treat it as a calibrated probability. **Remedy:** a single `probabilities.py` with audited helpers (softmax-with-temperature, logit-combine, AIC-weights, permutation-p-value) and a rule that no field named `*_probability`/`*_prob` may be set except through them.

2. **Mixture-variance algebra (`Σfσ²` vs `Σf²σ²`).** CHEM-2/7, CHEM-M1. **Remedy:** one `mixture_moments(f, mu, sigma)` function used by every mixing model; unit-test against a Monte-Carlo reference.

3. **Truncate-then-renormalize TTDs.** NUC-2/M5, PHYS-4. **Remedy:** a shared `ttd.normalize(undecayed_pdf)` enforced *before* decay; integrate to ≥10 mean lives or analytically.

4. **Unit/convention drift.** CORE-3/7, NUC-5, PHYS-2/5, NUC-L2. **Remedy:** enforce mmol/L (+ meq for charge-balanced reactions) at the `Config` boundary with an assertion, and a single canonical gradient/porosity definition.

A fifth, organizational issue: several flagship modules (OT, null-models, causal, latent endmembers) are **conceptually under-specified** — they ship a heuristic where a defined statistical/physical model is implied. These need a written model spec (generative story, estimand, validation) *before* re-implementation, not just a code patch.

---

## 6. Proposed remediation plan

### Phase 0 — Stop-the-bleeding (days)
Make wrong-but-plausible outputs either correct or clearly gated, so nothing new is published on broken math.
- CORE-1, CORE-2, CORE-3 — fix the transport objective, probability scaling, and γ bound. *(Smallest change with the largest credibility impact.)*
- PHYS-1 — correct PHREEQC dissolution gating.
- PHYS-2 — single Darcy gradient definition.
- CHEM-M1 / CHEM-M3 — make silent uniform/no-op fallbacks **raise or flag** instead of returning a fake result.
- Add `WARNING`-level "experimental / not calibrated" guards on OT, null-models, and causal outputs until Phase 2.

### Phase 1 — Correct the Bayesian & uncertainty core (1–2 weeks)
These produce the headline CIs and posteriors.
- TOPO-1, TOPO-2 — Hastings correction (or symmetric proposal) + correct sample normalization. Add a Geweke/`arviz` convergence check and a unit test that recovers a known posterior on a 3-node toy graph.
- NUC-1, NUC-M1, NUC-M2 — real jackknife BCa, route production CIs through it, handle degeneracy.
- NUC-7 — replace hand-rolled ESS with `arviz.ess`.
- NUC-6 — hierarchical network-age prior replacing the `−20·relu` potential.

### Phase 2 — Correct the physics/age forward models (2–3 weeks)
- NUC-2, NUC-M5, PHYS-4 — shared `ttd.normalize` (undecayed, full-support) used everywhere; Kreft–Zuber DM.
- NUC-3 — lagged-input reactive convolution.
- NUC-4 — peak-based cross-correlation lag.
- PHYS-5, PHYS-6 — effective porosity/tortuosity; affinity-dependent kinetic limit.
- CHEM-2, CHEM-7 — shared `mixture_moments`; reconcile analytical/MCMC nitrate likelihoods.

### Phase 3 — Re-specify the under-defined modules (3–4 weeks, design-first)
For each, write a one-page model spec (estimand, generative model, identifiability, validation dataset) **before** coding:
- CHEM-1 — OT: conserved-quantity ground metric + true unbalanced-OT, or replace with explicit stoichiometric reaction-distance.
- CHEM-5 — null models: permutation framework + FDR.
- CHEM-6 — causal: residualized partial correlation / transfer entropy with significance.
- CHEM-3 — latent endmembers: real vertices or proper PCA reconstruction with nonnegativity/closure.
- SHEAF-1/2/3 — orientation-invariant obstruction (holonomy), structural Betti number, fixed-universe leverage.

### Phase 4 — Hardening (ongoing)
- The four shared utilities from §5 (`probabilities`, `mixture_moments`, `ttd.normalize`, unit-assertion at `Config`).
- Property-based tests: probabilities ∈ [0,1] and sum to 1; obstruction energy invariant to edge re-orientation; mixture variance ≤ Monte-Carlo; MCMC recovers a known posterior; TTD integrates to 1 before decay.
- A `docs/MODEL_SPECS.md` linking each estimator to its math reference (the README already claims "all math examples are computationally verified" — extend that contract to the estimators here).

### Conceptual / organizational improvements (beyond bug-fixes)
- **Separate "calibrated" from "heuristic" outputs** in the result schema (e.g. a `provenance`/`calibrated: bool` field), so manuscripts can't accidentally cite a heuristic score as a probability.
- **Single source of truth for units and conventions** at the API boundary; reject ambiguous inputs early.
- **A small benchmark with known ground truth** (synthetic network with known topology, ages, and reactions) run in CI to catch regressions in the estimators, complementing the existing M1–M6 reference benchmarks.

---

## 7. Top-10 fix priority (impact × tractability)

1. CORE-1 — transport objective ignores transport fit
2. TOPO-1 — M-H missing Hastings correction
3. NUC-1 — BCa not implemented; production CIs uncorrected
4. PHYS-1 — PHREEQC dissolution gating
5. CORE-2 — uncalibrated transport "probabilities"
6. NUC-2 — truncate-then-renormalize ages
7. CHEM-2 — mixture variance `Σfσ²` vs `Σf²σ²`
8. PHYS-2 — inconsistent Darcy gradient
9. TOPO-2 — posterior normalization / noop handling
10. CHEM-M1/CHEM-M3 — silent fake-result fallbacks → make them fail loudly

---

*Reviewer note:* CORE-1/2/3, PHYS-1, PHYS-2, TOPO-1/2, SHEAF-1/2, TOPO-3/4, NUC-1, and the CORE-L1 no-op were re-read against source and confirmed (✅). The remaining findings are reported as read by the domain reviewers from the actual code with cited line numbers; each should be confirmed at the cited location before implementing its fix.
