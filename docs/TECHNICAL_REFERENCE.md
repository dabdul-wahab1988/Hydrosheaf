# Hydrosheaf Technical Reference

## 1. Mathematical Framework

Hydrosheaf solves the inverse geochemical modeling problem by finding the optimal combination of **transport** (mixing/evaporation) and **reaction** processes that evolve waters from an upstream source ($U$) to a downstream target ($V$).

The core optimization problem minimizes the following objective function for each edge $(u, v)$:

$$
\mathcal{L} = \| \mathbf{x}_v^{obs} - \mathbf{x}_v^{pred} \|_W^2 + \lambda_{L1} \sum_{k} |z_k| + \mathcal{P}_{SI} + \mathcal{P}_{Iso} + \mathcal{P}_{Gibbs}
$$

Where:

* $\mathbf{x}$ is the vector of ion concentrations (mmol/L).
* $z_k$ is the molar transfer (extent) of reaction $k$.
* $\|\cdot\|_W^2$ is the weighted squared Euclidean distance.

---

## 2. Transport Models

The model evaluates multiple transport hypotheses and selects the best fit based on the AIC/BIC-like objective score.

### 2.1 Evaporation Model (`evap`)

Assumes water is concentrated by a factor $\gamma$ due to evapotranspiration.
$$
\mathbf{x}_{transport} = \gamma \mathbf{x}_u \quad \text{where } \gamma = \frac{\mathbf{x}_u^T W \mathbf{x}_v}{\mathbf{x}_u^T W \mathbf{x}_u}
$$

* Constraint: $\gamma \ge 1.0$ (Net water loss only).

### 2.2 Mixing Model (`mix`)

Assumes simple binary mixing with a known end-member $\mathbf{x}_{end}$ (e.g., seawater).
$$
\mathbf{x}_{transport} = (1-f)\mathbf{x}_u + f\mathbf{x}_{end}
$$

* $f$ is the mixing fraction ($0 \le f \le 1$).

---

## 3. Reaction Optimization (LASSO)

After accounting for transport, the residual mass $\mathbf{r} = \mathbf{x}_v^{obs} - \mathbf{x}_{transport}$ must be explained by mineral reactions.

We solve for reaction extents $\mathbf{z}$ using **Coordinate Descent** with soft-thresholding:

$$
z_k^{(t+1)} = \frac{\mathcal{S}(\rho_k, \lambda/2)}{(\mathbf{R}_k^T W \mathbf{R}_k)}
$$

where $\mathcal{S}$ is the soft-thresholding operator:
$$
\mathcal{S}(x, \tau) = \text{sign}(x) \max(|x| - \tau, 0)
$$

This induces **sparsity**, forcing the model to explain the chemistry with the fewest possible minerals (Occam's Razor).

---

## 4. Constraints & Penalties

### 4.1 Thermodynamic Constraints (PHREEQC)

The model integrates PHREEQC to calculate Saturation Indices (SI) for the downstream water.

* **Dissolution Constraint**: If $SI_{mineral} > 0$ (oversaturated), dissolution is forbidden ($z_{min} \le 0$).
* **Precipitation Constraint**: If $SI_{mineral} < 0$ (undersaturated), precipitation is forbidden ($z_{min} \ge 0$).
This hard constraint corresponds to setting $LB=0$ or $UB=0$ in the optimization.

### 4.2 Isotope-Chemistry Consistency

To prevent unphysical decoupling of water isotopes and salinity, we impose a consistency penalty:

$$
\mathcal{P}_{consistency} = \lambda_{iso} \cdot | \gamma_{iso} - \gamma_{Cl} |
$$

* $\gamma_{iso}$: Concentration factor predicted by Rayleigh distillation of $\delta^{18}O$.
* $\gamma_{Cl}$: Concentration factor observed from Chloride ratios.

### 4.3 Redox Constraints

The model automatically classifies the redox environment:

* **Oxic**: $NO_3 > 0.05$ mmol/L.
* **Reducing**: $Fe > 0.01$ mmol/L.

In **Reducing** zones, aerobic oxidation reactions (e.g., $Pyrite + O_2$) are strictly disabled ($z = 0$).

---

## 5. Isotope Fractionation

Isotope evolution follows Rayleigh distillation kinetics or simple mixing.

### 5.1 Enrichment Slope

For each flow path, we calculate the local enrichment slope in $\delta^2H$-$\delta^{18}O$ space:
$$
S_{local} = \frac{\delta^2H_v - \delta^2H_u}{\delta^{18}O_v - \delta^{18}O_u}
$$

* $S \approx 8.0$: Meteoric water (Rain).
* $S \approx 4.0 - 6.0$: Evaporation signal.
* $S \to \infty$ (Vertical): Chemical alteration without water loss (Dissolution).

### 5.2 LMWL Fitting

The Local Meteoric Water Line is auto-calibrated using Orthogonal Distance Regression (ODR) or simple least squares:
$$
\delta^2H = a \cdot \delta^{18}O + b
$$

---

## 6. Nitrate Source Discrimination

This module implements a probabilistic classifier to distinguish between **Manure** and **Inorganic Fertilizer** sources of nitrate.

### 6.1 Methodology (Bayesian Accumulator)

The system calculates a Manure Probability ($P_m$) by accumulating "evidence" ($\phi_k$) in log-odds space:

$$
\text{Logit} = \ln\left(\frac{P_{prior}}{1-P_{prior}}\right) + \sum_{k} w_k \cdot \phi_k
$$

$$
P(\text{Manure}) = \frac{1}{1 + e^{-\text{Logit}}}
$$

### 6.2 Evidence Terms ($\phi$)

Features are z-scored using robust statistics (Median / MAD) to handle outliers:

1. **NO3/Cl Ratio ($\phi_1$)**: High ratios indicate Fertilizer; Low ratios indicate Manure/Sewage.
2. **NO3/K Ratio ($\phi_2$)**: Manure is rich in Potassium. High NO3/K strongly suggests Fertilizer.
3. **Denitrification ($\phi_5$)**: Derived from the LASSO reaction model. Strong denitrification supports Manure (organic carbon donor).
4. **Alkalinity Coupling ($\phi_6$)**: Correlation between NO3 loss and HCO3 gain indicates organic mineralization (Manure).

### 6.3 Contextual Gating

to prevent false positives in evaporative environments:

* **Evaporation Gate**: If Deuterium Excess ($d < 10$) or Transport Model = `evap`, the weight of ratio-based evidence is reduced by 50%.
* **CoDA Context**: A 7-ion Compositional Data Analysis (ilr coordinates) is used to validate that samples belong to the fresh/brackish conceptual model.

---

## 7. Advanced Extensions (v0.2.0)

### 7.1 Reactive Transport

Integrates kinetic rate laws to validate the "equilibrium" assumption of inverse models.

* **Reference**: `docs/REACTIVE_TRANSPORT_IMPLEMENTATION.md`
* **Verification**: `tests/test_accuracy_reactive.py` confirms Arrhenius temperature dependence and Damköhler number logic.

### 7.2 3D Flow Networks

Extends graph inference to layered aquifer systems.

* **Reference**: `docs/3D_FLOW_NETWORKS_IMPLEMENTATION.md`
* **Features**: Vertical anisotropy ($\alpha_v$), aquitard connectivity, and Bayesian topographic priors.
* **Verification**: `tests/test_accuracy_3d.py` and `tests/test_constraints.py` (probabilistic fallback).

### 7.3 Temporal Dynamics

Resolves time-variant signals and estimates residence times.

* **Reference**: `docs/TEMPORAL_IMPLEMENTATION.md`
* **Features**: Time-series interpolation, cross-correlation residence time (Center-of-Mass), and seasonal decomposition.
* **Verification**: `tests/test_accuracy_temporal.py` confirms retrieval of 10-day lag in synthetic signals.

### 7.4 Uncertainty Quantification

Provides statistical bounds on model outputs.

* **Reference**: `docs/UNCERTAINTY_IMPLEMENTATION.md`
* **Methods**: Bayesian MCMC (NUTS sampler) for rigorous posteriors; Bias-Corrected Bootstrap (BCa) for confidence intervals.
* **Verification**: `tests/test_accuracy_uncertainty.py` confirms coverage probabilities for Gaussian benchmarks.

---

## 8. Dual Isotope Nitrate Source Apportionment (v0.3.0)

Version 0.3.0 introduces a rigorous Bayesian mixing model for analyzing $\delta^{15}\text{N}_{\text{NO}_3}$ and $\delta^{18}\text{O}_{\text{NO}_3}$ data. this module works in tandem with the hydrochemical evidence accumulator (Section 6).

### 8.1 Process Logic

The system follows a strict priority hierarchy:

1. **Tier 1 (Isotopes)**: If dual isotope data is present, a 2D Bayesian mixing model calculates source probabilities directly.
2. **Tier 2 (Hydrochemistry)**: If isotope data is missing, the system falls back to the robust evidence accumulation method (Section 6), using NO3/Cl ratios, Potassium, and Denitrification signals.

### 8.2 Endmembers

Default endmembers are derived from *Kendall (1998)* and *Xue et al. (2009)*:

* **Manure**: High $\delta^{15}\text{N}$ (+8 to +25 permil).
* **Fertilizer**: Low $\delta^{15}\text{N}$ (-3 to +3 permil).
* **Soil N**: Intermediate (+3 to +8 permil).
* **Precipitation**: High $\delta^{18}\text{O}$ (+25 to +70 permil).

### 8.3 Verification

* **Test File**: `tests/test_nitrate_isotopes.py`
* **Validation**: Confirms that a sample with Manure isotopic signature is correctly identified even if its NO3/Cl ratio suggests Fertilizer, proving the priority of the isotopic signal.
