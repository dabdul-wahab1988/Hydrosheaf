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
