# HYDROSHEAF TECHNICAL MANUAL STRUCTURE
Based on the structure of HYDRUS_Technical_Manual_2D3D_V5.pdf

## FRONT MATTER
- **DISCLAIMER**
- **ABSTRACT**
- **TABLE OF CONTENTS**
- **LIST OF FIGURES**
- **LIST OF TABLES**
- **LIST OF MAIN VARIABLES** (Derived from existing "Notation" section)

## 1. INTRODUCTION
- **1.1. Overview** (Problem Motivation, Solution Strategy)
- **1.2. System Architecture** (Sheaf-Theoretic Perspective)
- **1.3. Input-Output Contract** (Brief high-level overview)

## 2. GOVERNING MATHEMATICAL FRAMEWORK
- **2.1. Graph-Theoretic Foundation**
    - 2.1.1. Network Representation
    - 2.1.2. Probabilistic Edge Inference
    - 2.1.3. Sheaf-Theoretic Topology Refinement
- **2.2. Transport Processes**
    - 2.2.1. Evaporation Model
    - 2.2.2. Mixing Model
    - 2.2.3. Transport Model Selection
- **2.3. Geochemical Reactions**
    - 2.3.1. Mass Balance and Stoichiometry
    - 2.3.2. Sparse Reaction Formulation (LASSO)
    - 2.3.3. Thermodynamic Constraints (Saturation Indices)
- **2.4. Isotope Hydrogeology**
    - 2.4.1. Local Meteoric Water Line (LMWL) & d-excess
    - 2.4.2. Isotope-Based Penalties
    - 2.4.3. Nitrate Source Discrimination (Dual Isotopes & Ratios)
- **2.5. Vadose Zone Physics** (from Extensions)
    - 2.5.1. Richards Equation and ET
    - 2.5.2. Vadose-to-Groundwater Travel Time
    - 2.5.3. Nitrate Breakthrough Convolution

## 3. NUMERICAL SOLUTION STRATEGY
- **3.1. Network Optimization Algorithm**
    - 3.1.1. Nested Model Selection Strategy
    - 3.1.2. Global Algorithm Structure
- **3.2. Solution of Reaction Equations**
    - 3.2.1. Coordinate Descent for Bounded LASSO
    - 3.2.2. Convergence Properties
    - 3.2.3. Stability-Based Hyperparameter Tuning
- **3.3. Graph Inference Algorithms**
    - 3.3.1. Bayesian Head Estimation
    - 3.3.2. Soft Selection and MAP Inference
- **3.4. Temporal Deconvolution**
    - 3.4.1. Cross-Correlation Residence Time
    - 3.4.2. TTD Convolution Solver

## 4. INVERSE MODELING AND PARAMETER ESTIMATION
- **4.1. Definition of the Objective Function**
    - 4.1.1. Weighted Least Squares with Regularization
    - 4.1.2. Penalty Functions (EC, TDS, Isotope, Gibbs)
- **4.2. Optimization Algorithms**
    - 4.2.1. Native GLM Engine
    - 4.2.2. External PEST++ Integration
- **4.3. Uncertainty Quantification**
    - 4.3.1. Bayesian MCMC (NUTS)
    - 4.3.2. Bootstrap Confidence Intervals
    - 4.3.3. Statistics of the Inverse Solution

## 5. PROBLEM DEFINITION AND PROGRAM IMPLEMENTATION
- **5.1. Implementation Map** (Math to Code mapping)
- **5.2. Programmatic API**
- **5.3. Computational Complexity**
- **5.4. Extension Capabilities** (3D, Reactive Transport, etc.)

## 6. EXAMPLE PROBLEMS AND VALIDATION
- **6.1. Validation Methodology** (Charge Balance, Cross-Validation)
- **6.2. Representative Case Studies**
    - 6.2.1. Salinization in Irrigated Aquifers
    - 6.2.2. Nitrate Contamination and Denitrification
    - 6.2.3. Aquifer Mixing Analysis
- **6.3. Benchmarks**

## 7. INPUT DATA REFERENCE
- **7.1. Command Line Interface (CLI)**
- **7.2. Data Schemas**
    - 7.2.1. Samples Input
    - 7.2.2. Temporal Data
    - 7.2.3. Physics Priors

## 8. OUTPUT DATA REFERENCE
- **8.1. Output Files Description**
- **8.2. Result Tables Structure** (JSON/CSV columns)

## 9. REFERENCES
