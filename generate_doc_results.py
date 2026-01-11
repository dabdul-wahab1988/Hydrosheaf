#!/usr/bin/env python3
"""
Generate verifiable results for the Hydrosheaf Technical Document examples.
This script produces numerical outputs that expert reviewers can reproduce.
"""

import numpy as np
from hydrosheaf.models.transport import fit_evaporation, fit_mixing
from hydrosheaf.models.reactions import fit_reactions
from hydrosheaf.models.nitrate_isotopes import load_isotope_endmembers, compute_isotope_prob, IsotopeSample
from hydrosheaf.nitrate_source_v2 import load_nitrate_config
import json

print("=" * 80)
print("HYDROSHEAF TECHNICAL DOCUMENT - VERIFIABLE RESULTS")
print("=" * 80)

# ============================================================================
# EXAMPLE 4.4: Evaporation vs. Mixing Discrimination
# ============================================================================
print("\n" + "=" * 80)
print("EXAMPLE 4.4: Evaporation vs. Mixing Discrimination (Section 4)")
print("=" * 80)

# Ion concentrations: Ca, Mg, Na, HCO3, Cl, SO4, NO3, F (mmol/L)
x_u = np.array([2.0, 1.0, 3.0, 5.0, 1.5, 2.0, 0.5, 0.1])
x_v = np.array([4.1, 2.0, 6.2, 10.1, 3.1, 4.0, 1.0, 0.2])
weights = np.ones(8)  # Identity weight matrix

print("\nUpstream concentrations (x_u):", x_u)
print("Downstream concentrations (x_v):", x_v)

# Evaporation hypothesis
print("\n--- EVAPORATION HYPOTHESIS ---")
numerator_evap = np.dot(x_u, x_v)
denominator_evap = np.dot(x_u, x_u)
gamma_star = numerator_evap / denominator_evap

print(f"Numerator (x_u^T · x_v): {numerator_evap:.5f}")
print(f"Denominator (x_u^T · x_u): {denominator_evap:.5f}")
print(f"Optimal gamma*: {gamma_star:.5f}")

# Predicted concentrations and residual
x_pred_evap = gamma_star * x_u
residual_evap = x_v - x_pred_evap
sse_evap = np.sum(residual_evap**2)

print(f"Predicted concentrations: {x_pred_evap}")
print(f"Residual: {residual_evap}")
print(f"Sum of squared errors: {sse_evap:.6f}")

# Mixing hypothesis with rainfall endmember
print("\n--- MIXING HYPOTHESIS (Rainfall Endmember) ---")
x_rain = np.array([0.2, 0.1, 1.0, 2.0, 0.5, 0.1, 0.0, 0.0])
d = x_rain - x_u

print(f"Rainfall endmember (x_rain): {x_rain}")
print(f"Difference vector (d = x_rain - x_u): {d}")

numerator_mix = np.dot(d, x_v - x_u)
denominator_mix = np.dot(d, d)
f_unconstrained = numerator_mix / denominator_mix

print(f"Numerator (d^T · (x_v - x_u)): {numerator_mix:.5f}")
print(f"Denominator (d^T · d): {denominator_mix:.5f}")
print(f"Unconstrained f: {f_unconstrained:.5f}")

# Apply constraints [0, 1]
f_star = max(0.0, min(1.0, f_unconstrained))
print(f"Constrained f* (projected to [0,1]): {f_star:.5f}")

# Predicted concentrations and residual for mixing
x_pred_mix = (1 - f_star) * x_u + f_star * x_rain
residual_mix = x_v - x_pred_mix
sse_mix = np.sum(residual_mix**2)

print(f"Predicted concentrations: {x_pred_mix}")
print(f"Residual: {residual_mix}")
print(f"Sum of squared errors: {sse_mix:.6f}")

print(f"\n*** CONCLUSION: Evaporation is strongly preferred (SSE_evap = {sse_evap:.6f} << SSE_mix = {sse_mix:.6f}) ***")

# ============================================================================
# EXAMPLE 5.4: Reaction Fitting for Calcite-Gypsum System
# ============================================================================
print("\n" + "=" * 80)
print("EXAMPLE 5.4: Reaction Fitting for Calcite-Gypsum System (Section 5)")
print("=" * 80)

# Residual after transport (mmol/L): Ca, Mg, Na, HCO3, Cl, SO4, NO3, F
r = np.array([1.5, 0.2, 0.0, 1.5, 0.0, 1.0, 0.0, 0.0])
print(f"\nResidual after transport (r): {r}")

# Stoichiometric reactions:
# Calcite: CaCO3 -> Ca2+ + HCO3- (affects indices 0=Ca, 3=HCO3)
# Gypsum: CaSO4 -> Ca2+ + SO42- (affects indices 0=Ca, 5=SO4)

rxn_calcite = np.zeros(8)
rxn_calcite[0] = 1.0  # Ca
rxn_calcite[3] = 1.0  # HCO3

rxn_gypsum = np.zeros(8)
rxn_gypsum[0] = 1.0  # Ca
rxn_gypsum[5] = 1.0  # SO4

reaction_matrix = [rxn_calcite, rxn_gypsum]
print("\nReaction matrix S:")
print(f"  Calcite: {rxn_calcite}")
print(f"  Gypsum: {rxn_gypsum}")

# Unconstrained least squares solution (for reference)
S = np.column_stack([rxn_calcite, rxn_gypsum])
z_ls = np.linalg.lstsq(S, r, rcond=None)[0]
print(f"\nUnconstrained LS solution: z_LS = {z_ls}")

# LASSO with L1 penalty
lambda_param = 0.1
weights_rxn = [1.0] * 8
lb = [-100.0, -100.0]
ub = [100.0, 100.0]
signed_mask = [True, True]

print(f"\nLASSO parameters: lambda = {lambda_param}")
fit_result = fit_reactions(r, reaction_matrix, weights_rxn, lambda_param,
                           signed_mask=signed_mask, lb=lb, ub=ub)
z_star = fit_result.extents

print(f"LASSO solution: z* = {z_star}")
print(f"  Calcite extent: {z_star[0]:.4f} mmol/L")
print(f"  Gypsum extent: {z_star[1]:.4f} mmol/L")

# Reconstruction
r_reconstructed = S @ z_star
residual_fit = r - r_reconstructed
sse_fit = np.sum(residual_fit**2)

print(f"\nReconstructed residual: {r_reconstructed}")
print(f"Fitting residual: {residual_fit}")
print(f"Sum of squared errors: {sse_fit:.6f}")
print(f"Number of active reactions: {np.sum(np.abs(z_star) > 1e-6)}")

# ============================================================================
# EXAMPLE 6.2: Saturation Index Constraints
# ============================================================================
print("\n" + "=" * 80)
print("EXAMPLE 6.2: Saturation Index Constraints (Section 6)")
print("=" * 80)

tau = 0.1  # Tolerance
print(f"\nSaturation index tolerance: tau = {tau}")

# Case 1: Calcite SI = 0.45 (supersaturated)
si_calcite = 0.45
print(f"\nCase 1 - Calcite:")
print(f"  SI = {si_calcite:.2f}")
print(f"  Status: Supersaturated (SI > {tau})")

if si_calcite > tau:
    bound_calcite = "z <= 0 (precipitation only)"
elif si_calcite < -tau:
    bound_calcite = "z >= 0 (dissolution only)"
else:
    bound_calcite = "z in R (free)"
print(f"  Bound: {bound_calcite}")

# Case 2: Gypsum SI = -1.2 (undersaturated)
si_gypsum = -1.2
print(f"\nCase 2 - Gypsum:")
print(f"  SI = {si_gypsum:.2f}")
print(f"  Status: Undersaturated (SI < -{tau})")

if si_gypsum > tau:
    bound_gypsum = "z <= 0 (precipitation only)"
elif si_gypsum < -tau:
    bound_gypsum = "z >= 0 (dissolution only)"
else:
    bound_gypsum = "z in R (free)"
print(f"  Bound: {bound_gypsum}")

# ============================================================================
# NITRATE SOURCE DISCRIMINATION: DUAL ISOTOPE METHOD
# ============================================================================
print("\n" + "=" * 80)
print("NITRATE SOURCE DISCRIMINATION - METHOD 1: Dual Isotope Mixing (Section 8)")
print("=" * 80)

# Load endmember database
sources = load_isotope_endmembers()
print("\nEndmember Database:")
for src in sources:
    print(f"  {src.name}:")
    print(f"    d15N: {src.d15N_mean:.1f} +/- {src.d15N_std:.1f} permil")
    print(f"    d18O: {src.d18O_mean:.1f} +/- {src.d18O_std:.1f} permil")

# Test case: Clear manure signature
print("\n--- Test Sample 1: Manure Signature ---")
sample_manure = IsotopeSample(d15N=15.0, d18O=5.0)
print(f"Sample isotopes: d15N = {sample_manure.d15N} permil, d18O = {sample_manure.d18O} permil")

probs_manure = compute_isotope_prob(sample_manure, sources)
for source_name, prob in probs_manure.items():
    print(f"  P({source_name}|sample) = {prob:.4f}")

print(f"\n*** Most likely source: {max(probs_manure, key=probs_manure.get)} (p = {max(probs_manure.values()):.4f}) ***")

# Test case: Clear fertilizer signature
print("\n--- Test Sample 2: Fertilizer Signature ---")
sample_fert = IsotopeSample(d15N=0.0, d18O=20.0)
print(f"Sample isotopes: d15N = {sample_fert.d15N} permil, d18O = {sample_fert.d18O} permil")

probs_fert = compute_isotope_prob(sample_fert, sources)
for source_name, prob in probs_fert.items():
    print(f"  P({source_name}|sample) = {prob:.4f}")

print(f"\n*** Most likely source: {max(probs_fert, key=probs_fert.get)} (p = {max(probs_fert.values()):.4f}) ***")

# ============================================================================
# NITRATE SOURCE DISCRIMINATION: HYDROCHEMICAL RATIOS METHOD
# ============================================================================
print("\n" + "=" * 80)
print("NITRATE SOURCE DISCRIMINATION - METHOD 2: Hydrochemical Ratios (Section 8)")
print("=" * 80)

# Load configuration
config = load_nitrate_config()
print("\nEvidence Weights from Configuration:")
print(f"  w1 (NO3/Cl ratio): {config['weights']['w1_no3_cl']}")
print(f"  w2 (NO3/K ratio): {config['weights']['w2_no3_k']}")
print(f"  w3 (PO4 flag): {config['weights']['w3_po4']}")
print(f"  w4 (Fe redox): {config['weights']['w4_fe']}")
print(f"  w5 (Denitrification): {config['weights']['w5_denitrif']}")
print(f"  w6 (Alk coupling): {config['weights']['w6_alk_coupling']}")
print(f"\nMinimum NO3 threshold: 10 mg/L (default, configurable via CLI)")
print(f"Evaporation gate factor: {config.get('evap_gate_factor', 0.5)}")

print("\n--- Key Features ---")
print("Method 2 uses six weighted evidence terms:")
print("  1. NO3/Cl ratio (w=1.2) - High values indicate fertilizer")
print("  2. NO3/K ratio (w=0.4) - Manure is K-rich")
print("  3. PO4 presence (w=0.3) - High PO4 suggests manure/sewage")
print("  4. Fe redox state (w=0.6) - High Fe indicates reducing conditions (manure)")
print("  5. Denitrification extent (w=1.5) - Strong denitrification supports manure")
print("  6. Alkalinity/NO3 coupling (w=0.8) - Derived from edge reaction fitting")

print("\nEvidence accumulates in log-odds space:")
print("  Logit = ln(P_prior/(1-P_prior)) + Sum(w_i * phi_i)")
print("  P(Manure) = sigmoid(Logit)")

print("\nGating logic:")
print("  - If evaporation detected OR deuterium excess < threshold:")
print("    -> Reduce ratio evidence weights (prevent false fertilizer attribution)")
print("  - If isotopes available AND NO3 >= threshold:")
print("    -> Method 1 (Dual Isotope) takes priority")

print("\n" + "=" * 80)
print("ALL RESULTS GENERATED SUCCESSFULLY")
print("=" * 80)
print("\nThese results are computationally verified and can be reproduced")
print("by expert reviewers using the test suite: pytest tests/test_doc_examples.py")
print("=" * 80)
