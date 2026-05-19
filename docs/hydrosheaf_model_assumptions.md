# Hydrosheaf Model Assumptions Catalogue

## 1. Lumped Parameter Models (9 types)

| Model | PDF Form | Parameters | Key Assumptions |
|-------|----------|------------|-----------------|
| **PFM** | δ(t − τ) | `mean_age_years` | Zero dispersion; all parcels travel at identical speed |
| **EM** | (1/τ)·e^(−t/τ) | `mean_age_years` | Well-mixed reservoir; exponential residence times |
| **DM** | 1D advection-dispersion flux solution | `mean_age_years`, `dispersion` | 1D transport; Pd ∈ [1e−4, 2.0]; grid: 0.01–1.0 |
| **GA** | Gamma(shape, scale=τ/shape) | `mean_age_years`, `shape` | Shape ∈ [0.1, 10.0]; grid: 0.5–5.0 |
| **EPM** | Piston delay + exponential tail | `mean_age_years`, `piston_fraction` | Piston fraction ∈ [0, 0.95]; remainder is EM |
| **PEM** | Truncated exponential (short-circuited) | `mean_age_years`, `capture_fraction` | Capture fraction ∈ [0.05, 1.0]; cutoff = −τ·ln(1−η) |
| **EMM** | Two-component EM mixture | `mean_age_years`, `young_fraction`, `old_age_ratio` | Old/young ratio ∈ [2, 10]; only 2 components |
| **BMM-DM-DM** | Binary mix of two DMs | `mean_age_1_years`, `mean_age_2_years`, `binary_fraction`, `dispersion`, `dispersion_2` | Age₂ > Age₁ always; fraction ∈ [0.001, 0.999] |
| **BMM-PEM-PEM** | Binary mix of two PEMs | Same pattern + `capture_fraction`, `capture_fraction_2` | Same binary constraints |

### Common LPM assumptions

- Integration: dτ = 0.25 yr (age ≤ 250), 1.0 yr (≤ 2000), 10 yr (otherwise)
- Upper bound: max(100, 6·τ_mean, max_age)
- PDFs always re-normalized to unit area
- PFM pre-bomb fallback: 0.5 TU when recharge date precedes history

---

## 2. Tracer Nuclide Assumptions

| Nuclide | Half-life | Units | Decay Model |
|---------|-----------|-------|-------------|
| ³H | 12.32 yr | TU | exp(−λt), year = 365.25 days |
| ¹⁴C | 5730 yr (Libby) | pmc | exp(−λt), mean life = 8033 yr |
| ³⁹Ar | 269 yr | pmc | exp(−λt), initial = 100 pmc |
| ⁸⁵Kr | 10.76 yr | dpm/cc | exp(−λt) |

### Tritium input function

- Ottawa-style Northern Continental (peak 3000 TU in 1963–64)
- Tropical scaled by 0.2
- Cosmogenic baseline: 3.0 TU
- Pre-1952 background: 5.0 TU

### Atmospheric gas histories

Compact Northern Hemisphere defaults:

| Tracer | Range | Peak |
|--------|-------|------|
| SF6 | 0 → 12 pptv | monotonic increase |
| CFC11 | ~260 pptv | ~1990 |
| CFC12 | ~540 pptv | ~2000 |
| CFC113 | ~85 pptv | ~1995 |
| ⁸⁵Kr | 0 → 110 pptv | monotonic increase |

Assumed replaceable with station-specific data.

### 3H/³He age

- Closed-system ingrowth: `age = ln(1 + He³/³H) / λ`
- No degassing, no dispersion

### 4He age

- Formula: `(He − 4.6e−8) / 2.0e−11`
- Default background: 4.6e−8 ccSTP/g
- Default accumulation rate: 2.0e−11 ccSTP/g/year
- "Highly uncertain unless site-calibrated"

---

## 3. Observation Uncertainty Assumptions

| Tracer | Relative σ | Floor |
|--------|-----------|-------|
| ³H | 15% | 0.1 TU |
| ³H/³He | 20% | 0.1 TU |
| ¹⁴C | 5% | 1.0 pmc |
| ³⁹Ar | 10% | 5.0 pmc |
| ⁸⁵Kr | 10% | 0.05 |
| SF6/CFCs | 10% | 0.05 pptv |
| ⁴He | 50% | 1.0e−9 cc/g |

### Likelihood types

- **gaussian** (default): squared standardized residual
- **upper_censored**: supersaturated gas → no penalty when pred ≤ obs
- **lower_censored**: no penalty when pred ≥ obs
- **contaminated_mixture**: Cauchy-like robust loss = log1p(z²)

---

## 4. ¹⁴C Correction Strategies

| Mode | Assumption |
|------|-----------|
| `raw` | A₀ = 100 pmc, no correction |
| `fixed_a0` | A₀ = 100 pmc fixed |
| `selected` | Use pre-computed corrected ¹⁴C from USGS table |
| `ensemble` | Weighted average: PHREEQC(8), NETPATH(7), F&G(5), Tamers(3), scaling(1), other(0.5) |
| `hierarchical` | Study-unit/aquifer-level A₀ prior from grouped data |

Valid pmc range: (0, 130].

---

## 5. ⁴He Uncertainty Modes

| Mode | Assumption |
|------|-----------|
| `calibrated` | Uses provided σ |
| `calibrated_uncertainty` | Inflates σ = max(σ_meas, excess·fraction, 1e−12); fraction = 0.25 (4+ fields), 0.35 (2+), 0.50 (default) |
| `hierarchical` | Study-unit/aquifer prior for background and accumulation rate |
| `disabled` | ⁴He excluded entirely |

---

## 6. Dissolved Gas Correction

### Excess-air models

- **UA** (Unfractionated Air): single excess-air parameter
- **CE** (Closed-system Equilibration): excess air + fractionation

### Solubility model

- Exponential temperature correction from 20°C reference
- Pressure and salinity corrections applied multiplicatively
- Default solubilities at 20°C (ccSTP/g): He=4.6e−8, Ne=1.9e−7, Ar=3.2e−4, Kr=7.3e−8, Xe=1.05e−8

### Fitting

- Grid search: T ∈ [0, 35]°C, excess air ∈ [0, 0.02] ccSTP/g, fractionation ∈ {0, 0.25, 0.5, 0.75, 1.0}
- Observation σ: 3% relative
- Objective: least squares
- Model ranking: AIC

### Monte Carlo

- 100 samples, seed=42
- Each realization samples from N(obs, σ)
- Summary: p05, p50, p95

---

## 7. Model Selection & BMA

- **Criterion:** AICc (corrected for small n)
- AICc = ∞ when n ≤ k+1
- **Gated BMA:**
  - Only converged fits within ΔAICc ≤ 4.0
  - Fallback to top model if log₁₀(age span) > 0.7
  - Fallback to top model if fewer than 2 models pass gates
- **Refinement:** Nelder-Mead (300 iter) then L-BFGS-B fallback
- **Effective n_params:** PFM/EM=1, DM/GA/EPM/PEM=2, EMM=3, BMM=5

---

## 8. Geochemical Reaction Models

### Sparse fitting

- Elastic net coordinate descent (L1 + L2 penalties)
- Non-negativity enforced for dissolution reactions
- Convergence: 300 iter, tol=1e−6
- Model selection: AICc, BIC

### Geologic bias scales

| Setting | Favored (×0.5) | Penalized (×10–50) |
|---------|---------------|-------------------|
| Crystalline | albite, anorthite, feldspar, biotite, chlorite, pyroxene | calcite, dolomite, gypsum, anhydrite, halite |
| Sedimentary | calcite, dolomite, magnesite, aragonite | albite, anorthite, feldspar, biotite, chlorite, pyroxene |

### Indicator ion gating

Mandatory ions must be present for each reaction:

| Reaction | Required Ions |
|----------|--------------|
| pyrite_oxidation_aerobic | Fe, SO₄ |
| pyrite_oxidation_denit | Fe, SO₄ |
| biotite | K, Mg, Fe |
| chlorite | Mg, Fe |
| fluorite | F, Ca |
| sylvite | K, Cl |

### Concentration gates

| Condition | Threshold |
|-----------|----------|
| Gypsum SO₄ minimum | 0.208 mmol/L |
| Halite Cl minimum | 0.564 mmol/L |
| Fluorite F minimum | 0.026 mmol/L |
| Denitrification NO₃ pruning | < 0.16 mmol/L |
| High NO₃ source forcing | > 0.8 mmol/L |

### Ion exchange

- 1:2 divalent:monovalent stoichiometry
- CaNa_exch: Ca=+1, Na=−2 (freshening)
- NaCa_exch: Ca=−1, Na=+2 (salinization)

### Redox classification

- Oxic: NO₃ > 0.05 mmol/L
- Reducing: Fe > 0.01 mmol/L
- In reducing conditions: pyrite_oxidation_aerobic forced to [0, 0]

---

## 9. Graph / Network Topology

### Head inference

- Bayesian linear-Gaussian hierarchical model
- 3-tier observations: direct measurement, depth-to-water, topographic proxy
- Topographic prior: hᵢ + μ_dtw ~ N(elevᵢ, σ_topo²)
- MCMC alternative via PyMC/NUTS (falls back to closed-form)

### Edge confidence

- Prior p_uv from edge attribute, default = 1.0
- Clamped to [1e−6, 1.0]

### Topology refinement

- Edge penalty = −log(p_uv)·head_weight + isotope cost + Cl cost + age cost
- pi_evap refined by depth, d-excess, and thermodynamic heuristics
- pi_evap capped to [0, 0.8]
- Softmax selection with β = 2.0
- Max 3 refinement iterations

### Network aging

- Hierarchical Bayesian (PyMC/NUTS)
- log_age prior ~ N(log(20), 2.0)
- Network constraint: age_v ≥ age_u + distance/velocity
- ¹⁴C A₀: hierarchical with global mean ~ N(85, 10)
- 4 chains, target_accept=0.95

---

## 10. Residence Time Estimation

Five methods in priority order:

| Method | Approach | Assumptions |
|--------|----------|------------|
| Cross-correlation | Multi-tracer (Cl, δ¹⁸O, δ²H) consensus | Minimum 3 points; agreement tolerance 0.4 spread |
| Bayesian lag | Jeffreys prior on log(τ) | Truncated normal prior; 5-point minimum |
| TTD convolution | NNLS non-negative weights | Toeplitz lag matrix; smoothness penalty |
| Darcy gradient | τ = L·φ/(K·i) | 50% uncertainty |
| Tracer decay | PFM age difference | Negative differences clamped to 0 |

### Quality gates

- Isotope evaporation: LMWL RMSE > 6 → weight ×0.35; d-excess < 5 → weight ×0.35
- Cl non-conservative: stepiness ratio > 2.0 → weight ×0.4
- GMWL default: δ²H = 8.0·δ¹⁸O + 10.0

---

## 11. Age Fraction / Geological Epoch Boundaries

| Epoch | Transit Time Range | Rationale |
|-------|-------------------|-----------|
| Anthropocene | τ ≤ 70 yr | Post-1950 (~present minus 70) |
| Holocene | 70 < τ ≤ 11,700 yr | Standard Holocene boundary |
| Pleistocene | τ > 11,700 yr | Pre-Holocene |

- Fractions normalized to sum to 1.0
- Default σ_fraction = 0.10 for fitting constraints

### ¹⁴C/⁴He constraint agreement

| Gap (log₁₀) | Status |
|-------------|--------|
| ≤ 0.30 | Agreement |
| 0.30–0.70 | Tension |
| > 0.70 | Conflict |

---

## 12. Reactive Transport Kinetics

### Rate laws (PHREEQC-style)

| Reaction | Rate Expression | Condition |
|----------|----------------|-----------|
| calcite, dolomite, gypsum, halite, fluorite | k·A·(1 − 10^SI) | Dissolution + precipitation |
| albite, anorthite | k·A·(1 − 10^SI) | Dissolution only (SI < 0) |
| pyrite_oxidation_aerobic | k·A·[O₂]^0.5 | O₂ > 1e−6 |
| denitrification | −k·[NO₃] | NO₃ > 1e−6 |

### Default kinetic parameters

| Reaction | k (mol/m²/s) | Surface Area (m²/L) | Ea (J/mol) |
|----------|-------------|--------------------|------------|
| calcite | 1e−6 | 0.1 | 41,840 |
| dolomite | 1e−8 | 0.05 | 52,000 |
| gypsum | 1e−4 | 0.1 | 15,000 |
| halite | 1e−3 | 0.1 | 10,000 |
| fluorite | 1e−7 | 0.05 | 45,000 |
| albite | 1e−12 | 0.1 | 70,000 |
| anorthite | 1e−11 | 0.1 | 65,000 |
| pyrite_ox | 1e−9 | 0.05 | 56,000 |
| denitrification | 1e−7 (1/s) | 1.0 | 69,000 |

### Temperature correction (Arrhenius)

- k(T) = k_ref · exp(−Ea/R · (1/T − 1/T_ref))
- R = 8.314 J/(mol·K), T_ref = 298.15 K

---

## 13. Key Hardcoded Constants

### Nuclear / Decay

- Tritium half-life: 12.32 years
- ¹⁴C half-life: 5730.0 years (Libby)
- ³⁹Ar half-life: 269.0 years
- ⁸⁵Kr half-life: 10.76 years
- ¹⁴C mean life: 8033 years
- Year length: 365.25 days

### Atmospheric / Dissolved Gas

- Pre-bomb tritium: 5 TU
- Cosmogenic tritium base: 3.0 TU
- Cosmogenic ¹⁴C initial: 100 pmc
- ¹⁴C valid range: (0, 130] pmc
- ⁴He background: 4.6e−8 ccSTP/g
- ⁴He accumulation rate: 2.0e−11 ccSTP/g/year
- Molar volume STP: 22,414 cc/mol
- Gas constant R: 0.082057366 L·atm/(mol·K)

### Age Boundaries

- "Modern" water: age < 60 years
- Anthropocene: τ ≤ 70 years
- Holocene: 70 < τ ≤ 11,700 years
- Pleistocene: τ > 11,700 years
- Old groundwater: ¹⁴C age > 1,000 years
- Young gas tracer max age: 85 years

### Geochemical

- Denitrification NO₃ pruning: < 0.16 mmol/L
- High nitrate forcing: > 0.8 mmol/L
- Gypsum SO₄ threshold: < 0.208 mmol/L
- Halite Cl threshold: < 0.564 mmol/L
- Fluorite F threshold: < 0.026 mmol/L

### Transport

- Default porosity: 0.25
- Default dispersivity: 1.0 m
- Default denitrification rate: 0.001 day⁻¹
- Default hydraulic gradient: 0.01

### Numerical

- Age grid search: 90 steps
- Refinement: Nelder-Mead (300 iter) → L-BFGS-B
- PDF integration steps: 0.25, 1.0, or 10.0 years
- Picard tolerance: 1e−6
- Reaction convergence tolerance: 1e−6
- Least-squares tolerance: 1e−8
