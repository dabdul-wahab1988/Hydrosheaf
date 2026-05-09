<!-- Markdown companion for Milestone 6 (M6) manuscript planning. -->

# A) Quantifying Information-Theoretic Parsimony and Cascading Uncertainty in the Hydrosheaf Groundwater Discovery Framework

# B) Research Questions & Guardrails

- **Primary Question:** How sensitive is Hydrosheaf’s geochemical discovery to the choice of mathematical penalties and the inherent uncertainty in residence-time physics?
- **Sub-question 1:** Does the L1-regularization path exhibit a stable "Discovery Plateau" where mineral identification is invariant to the exact penalty weight ($\lambda$)?
- **Sub-question 2:** Can Information-Theoretic metrics (AIC/BIC/AICc) successfully resolve system states (e.g., Open vs. Closed Carbonate evolution) without manual intervention?
- **Sub-question 3:** How does residence-time uncertainty (±σ years) propagate through the "Triplet" into the stability of the final reaction assemblage?
- **Required wording guardrail:** M6 must be presented as a "Mathematical Stress Test" of the framework. It should not claim to find new field results, but rather to prove the statistical reliability of the results found in M2 and M5.

# C) Detailed Outline with word counts

## Abstract — [350] words

- Rationale for robustness testing in inverse modeling, summary of sensitivity methods (L1-Paths, AICc, Monte Carlo), key findings on Phase Stability Indices, and implications for automated discovery.

## 1. Introduction — [1,200] words

### 1.1 The "Ill-Posed" Problem of Geochemical Inversion — [300] words
- High non-uniqueness and the critique of manual mineral selection.

### 1.2 Moving from Calibration to Regularization — [300] words
- The need for formal statistical constraints (L1, Information Theory).

### 1.3 The Problem of Cascading Uncertainty — [300] words
- How errors in age (M3) and topology (M4) "pollute" the reaction discovery (M5).

### 1.4 Aim and Objectives of the Robustness Milestone — [300] words

## 2. Materials and Methods — [4,500] words

### 2.1 The Regularization Path Algorithm — [500] words
- Mathematical definition of the L1-gradient and the automated $\lambda$-search.

### 2.2 Information-Theoretic Selection (AIC/BIC/AICc) — [600] words
- Application of AIC and Corrected AIC (AICc) to sparse hydrogeochemical mass-balance.
- Formal penalty for model complexity and handling small sample sizes (few ions).

### 2.3 Global Sensitivity Analysis (GSA) Framework — [750] words
- Design of the Monte Carlo perturbation suite for chemical and physical inputs.
- Definition of the **Phase Stability Index (PSI)**: the probability that a mineral is included in the sparse solution under Gaussian noise.

### 2.4 Quantifying the Uncertainty Cascade — [850] words
- Step-by-step propagation: Tracer Error $\rightarrow$ Age Variance $\rightarrow$ Kinetic Bound Shift $\rightarrow$ Mineral Selection change.

### 2.5 Theoretical Limits of Reduced-Order discovery — [600] words
- Defining the "Detection Limit" of the framework (when is data too sparse to trust?).

### 2.6 Field Application: Stress-Testing the Ghana (M2) Results — [600] words
- Sub-sampling and perturbation tests on the `manu` and `talensi` datasets.

### 2.7 Software implementation and reproducibility — [600] words
- Integration of `sensitivity.py` and `kinetic_limit.py` into the Hydrosheaf API.

## 3. Results and Discussion — [2,500] words

### 3.1 Regularization Path and Discovery Plateaus — [600] words
- Figure: Mineral extents vs. $\lambda$. Demonstration of the "Optimal Lambda" identification via AICc.

### 3.2 AICc-Based System Resolution — [500] words
- Evidence of AICc correctly distinguishing Open from Closed system carbonate evolution in deep vs. shallow samples.

### 3.3 The Uncertainty Cascade Matrix — [700] words
- Heatmap showing the elasticity of minerals to input parameters.
- Identifying "Robust" minerals (insensitive to age) vs. "Kinetic" minerals (sensitive to age).

### 3.4 Benchmarking the Phase Stability Index (PSI) — [700] words
- Quantifying the reliability of the Ghana discovery results under field-noise conditions.
- Demonstration of high PSI (>60%) for dominant processes even under 5% measurement error.

## 4. Conclusion — [500] words
- Summary of the framework's mathematical defense.
- Final statement on the reliability of the "Hydrosheaf Triplet" as a scientific standard.

# D) Proposed Tables and Figures

- **Table 1: Information Criterion Summary**: Comparison of AIC/BIC/AICc across model versions.
- **Table 2: Phase Stability Matrix**: Percent probability of mineral recovery (PSI) under ±5% chemical noise.
- **Figure 1: The Regularization Path**: Visualization of mineral selection "switching" as $\lambda$ changes.
- **Figure 2: The Uncertainty Cascade**: Flow diagram showing error propagation from tracer to reaction.
- **Figure 3: Global Sensitivity Heatmap**: Elasticity of reaction extents to input concentration and age variance.
