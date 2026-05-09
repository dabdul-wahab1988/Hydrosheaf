# Hydrosheaf Functional Specification Document (FSD)

**Version:** 3.0.0
**Status:** Active Development
**Target Audience:** Scientific Supervisors, Hydrogeologists, and Research Engineers.

---

## 1. Executive Summary

**Hydrosheaf** is a comprehensive Python framework for **Inverse Hydrogeochemical Modeling**, developed to solve the "Sheaf" problem: reconstructing the optimal network of flow paths and biogeochemical processes that govern groundwater evolution.

Unlike traditional forward models (PHREEQC/MODFLOW) that require exhaustive parameterization, Hydrosheaf is **inference-driven**. It works backwards from sparse water quality data to identify:

1. **Topology**: Which wells are hydraulically connected?
2. **processes**: What reactions (Mixing, Evaporation, Mineral Dissolution) occur along these connections?
3. **Sources**: Are nitrate contaminants originating from Manure or Fertilizer?
4. **Timescales**: How long does water reside in the aquifer?

The system integrates **Sparse Optimization (LASSO)**, **Bayesian MCMC Inference**, and **Physical Process Models (Richards Eq, TTDs)** into a unified pipeline.

---

## 2. Core Module 1: Inverse Modeling Engine (The "Sheaf" Solver)

The mathematical core (`hydrosheaf.inference`) explains the chemical evolution between an upgradient node $\mathbf{u}$ and downgradient node $\mathbf{v}$ by solving a constrained optimization problem.

### 2.1 The Inverse Equation

We model the observed concentration vector $\mathbf{x}_v$ as:
$$ \mathbf{x}_v \approx T(\mathbf{x}_u \mid \theta) + S \cdot \xi + \epsilon $$

Where:

* $T(\cdot)$: Physical transport operator (Mixing or Evaporation).
* $S$: Stoichiometric matrix of candidate phases (minerals, gases).
* $\xi$: Vector of reaction mass transfers (mol/L).
* $\epsilon$: Residual error.

### 2.2 Transport Operators

The framework evaluates competing transport hypotheses:

1. **Evaporation (Rayleigh-Type)**:
    $$ \mathbf{x}_{pred} = \gamma \cdot \mathbf{x}_u $$
    * $\gamma \ge 1$: Concentration factor ($1 / \text{Leaching Fraction}$).
2. **Conservative Mixing**:
    $$ \mathbf{x}_{pred} = (1 - f) \mathbf{x}_u + f \cdot \mathbf{x}_{end} $$
    * $f \in [0, 1]$: Mixing fraction of a known endmember (e.g., Seawater).

### 2.3 Sparse Reaction Solver (LASSO)

To ensure geochemical plausibility, we enforce **Sparsity** (Occam's Razor) using **L1-Regularization**. We solve for $\xi$ via:

$$ \min_{\xi} \frac{1}{2} \| W^{1/2} (\mathbf{x}_{res} - S \xi) \|_2^2 + \lambda \|\xi\|_1 $$

**Algorithm**: We implement **Coordinate Descent** with Soft-Thresholding:
$$ \xi_j^{(k+1)} = \frac{\text{SoftThreshold}(\rho_j, \lambda)}{ (S^T W S)_{jj} } $$
$$ \text{SoftThreshold}(z, \tau) = \text{sgn}(z) \max(|z| - \tau, 0) $$

### 2.4 Physical Constraints & Penalties

The objective function $J$ includes penalty terms used to guide the solver towards physically valid solutions:

* **Gibbs Penalty**: Penalizes "Evaporation" models if the sample plots in the "Rock Dominance" or "Precipitation Dominance" fields of a Gibbs Diagram.
* **Isotope Consistency**: Checks if $\delta^{18}O$ and $\delta^2H$ enrichment aligns with the computed chloride concentration factor $\gamma$.
* **EC/TDS Penalty**: $\eta (Obs_{EC} - Est_{EC}(\mathbf{x}_{pred}))^2$, ensuring the major ion reconstruction matches bulk salinity.

---

## 3. Core Module 2: Nitrate Source Discrimination

Located in `hydrosheaf.nitrate_source_v2`, this module uses a hierarchical Bayesian approach to distinguish Manure vs. Fertilizer sources.

### 3.1 Dual Isotope Mixing (MCMC)

When $\delta^{15}N$ and $\delta^{18}O_{NO3}$ data are available, we use **PyMC** to sample the posterior source probabilities:

* **Model**: $\mathbf{Obs} \sim \mathcal{N}(\sum f_i \mu_i, \sigma^2)$
* **Priors**: $f \sim \text{Dirichlet}(\alpha)$, $\sigma \sim \text{HalfNormal}$.
* **Output**: Full posterior distribution $P(\text{Manure} \mid \text{Isotopes})$.

### 3.2 Hydrochemical Evidence Fusion

When isotopes are missing, we use **Bayesian Evidence Fusion** of hydrochemical proxies.

* **Method**: Weighted Sigmoid Fusion of Robust Z-Scores.
* **Robust Z-Score**: $Z = \frac{x - \text{Median}(x)}{1.4826 \cdot \text{MAD}(x)}$ (Insensitive to outliers).
* **Evidence Nodes**:
    1. **$NO_3/Cl$ Ratio**: Low $\to$ Manure.
    2. **$NO_3/K$ Ratio**: High $\to$ Fertilizer (Potash).
    3. **High $PO_4$**: Low mobility, high in manure.
    4. **Redox State**: $Fe > 0$ implies reducing conditions (often manure-linked).

### 3.3 Compositional Data Analysis (CoDA)

To avoid spurious correlations in closed datasets (mg/L sums to TDS), we transform data into **Isometric Log-Ratio (ilr)** coordinates using a specific **Sequential Binary Partition (SBP)** tree (`coda_sbp.py`):

1. Cations ($Ca, Mg, Na, K$) vs Anions ($HCO_3, Cl, SO_4$).
2. Alkaline Earths ($Ca, Mg$) vs Alkalis ($Na, K$).
3. $Ca$ vs $Mg$.
4. $Na$ vs $K$.
5. $HCO_3$ vs Salinity Ions ($Cl, SO_4$).
6. $Cl$ vs $SO_4$.

---

## 4. Core Module 3: Temporal Dynamics (Residence Time)

The `hydrosheaf.temporal` module estimates groundwater travel times ($\tau$).

### 4.1 Cross-Correlation Consensus

Estimates lag between time-series $u(t)$ and $v(t)$ for multiple tracers ($Cl$, $\delta^{18}O$):

1. **Normalize**: Standardize signals ($Z$-score).
2. **Cross-Correlate**: Calculate $r(\tau) = (u * v)(\tau)$.
3. **Consensus**: Weighted average of $\tau_{peak}$ from all valid tracers.
4. **Gating**:
    * **LMWL Gate**: Rejects isotope lags if samples deviate significantly from the Local Meteoric Water Line (indicating local evaporation, not transport lag).
    * **Chloride Step Gate**: Rejects lags if $Cl$ shows non-conservative step-changes.

### 4.2 TTD Convolution (NNLS)

Deconvolves the Transit Time Distribution (TTD) from input/output signals:
$$ v(t) = \int_0^{\infty} u(t - \tau) h(\tau) d\tau $$

* **Solver**: Non-Negative Least Squares (NNLS) to solve for $h(\tau)$.
* **Constraints**: Smoothness penalty (Tikhonov on 1st derivative) and exponential decay prior.

---

## 5. Core Module 4: Vadose Zone Physics

The `hydrosheaf.vadose` module simulates vertical recharge using the **1D Richards Equation**:
$$ C(\psi) \frac{\partial \psi}{\partial t} = \frac{\partial}{\partial z} \left[ K(\psi) \left( \frac{\partial \psi}{\partial z} + 1 \right) \right] - S(\psi) $$

* **Solver**: Finite Volume Method (Cell-Centered).
* **Non-Linearity**: Solved via implicit **Picard Iteration**.
* **Hydraulic Functions**: Van Genuchten-Mualem ($K(\psi), \theta(\psi)$).
* **Root Uptake**: Feddes reduction function for transpiration sink $S(\psi)$.

---

## 6. Core Module 5: Calibration (PESTGLM)

The `hydrosheaf.calibration` module implements a custom **Gauss-Levenberg-Marquardt** solver, mimicking PEST++ functionality in pure Python.

* **Objective**: Minimize weighted sum of squared residuals (Phi).
    $$ \Phi = \sum w_i (Obs_i - Sim_i)^2 + \Phi_{reg} $$
* **Jacobian**: Calculated via Parallel Finite Differences (`ThreadPoolExecutor`).
* **Regularization**: Tikhonov regularization ($w_{reg} (p_{curr} - p_{prior}) = 0$) to stabilize ill-posed inversion.
* **Uncertainty**: Approximates posterior covariance as $\Sigma = \sigma^2 (J^T W J)^{-1}$.

---

## 7. Software Architecture

### 7.1 Directory Map

* `hydrosheaf/inference/`: Core LASSO and Graph solvers (`edge_fit.py`, `network_fit.py`).
* `hydrosheaf/models/`: Geochemical logic (`reactions.py`, `mixing.py`, `nitrate_isotopes.py`, `gibbs.py`).
* `hydrosheaf/temporal/`: Time-series analysis (`residence_time.py`, `time_series.py`).
* `hydrosheaf/vadose/`: Physics-based recharge (`richards1d.py`, `soil.py`).
* `hydrosheaf/calibration/`: Optimization engines (`glm.py`).
* `hydrosheaf/graph/`: Topology construction and Head Inference (`head_inference.py`).

### 7.2 Key Dependencies

* **NumPy/SciPy**: Linear algebra, NNLS, Sparse matrices.
* **Pandas**: Data handling.
* **PyMC**: Bayesian MCMC sampling (optional).
* **NetworkX**: Graph theory operations.

---
