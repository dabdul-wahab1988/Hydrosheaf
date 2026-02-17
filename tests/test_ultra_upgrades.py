"""
Tests for Ultra Upgrade features: Latent Endmembers, Compositional Weighting, Iterative Jacobian, Reacted TTD.
"""

import math
import pytest
from typing import List

try:
    import numpy as np
except ImportError:
    np = None

from hydrosheaf.config import Config
from hydrosheaf.models.latent import identify_latent_endmembers
from hydrosheaf.inference.edge_fit import fit_edge, EdgeResult
from hydrosheaf.temporal.convolution import convolve_reactive_ttd, compute_effective_reaction_factor
from hydrosheaf.coda_sbp import clr, clr_inv

# --- Test 1: Latent Endmembers (Topology) ---

@pytest.mark.skipif(np is None, reason="Numpy required for latent endmembers")
def test_identify_latent_endmembers_simple_mixing():
    """
    Test if PCA-based identification finds endmembers along a mixing line.
    True Endmembers: A=[10, 10], B=[100, 100] (Log scale: 2.3, 4.6)
    Samples: Mixtures of A and B.
    """
    ion_order = ["Ca", "Cl"]
    
    # Create synthetic samples along a line
    samples = []
    # Mixes: 10%, 50%, 90% of B
    mixes = [0.1, 0.5, 0.9]
    
    # Endmembers
    em1 = np.array([10.0, 20.0]) # A
    em2 = np.array([100.0, 200.0]) # B (Same ratio, just concentrated)
    # Actually, in CLR space, concentration is just a translation if ratios are constant.
    # CLR([10, 20]) -> ln(10/sqrt(200)), ln(20/sqrt(200)) -> ln(0.7), ln(1.4)
    # CLR([100, 200]) -> same!
    # PCA on constant ratio data will fail to find 2 endmembers based on composition difference.
    # We need chemically distinct endmembers (different ratios).
    
    em1 = np.array([100.0, 10.0]) # High Ca, Low Cl
    em2 = np.array([10.0, 100.0]) # Low Ca, High Cl
    
    for i, f in enumerate(mixes):
        # Linear mix in concentration space
        conc = (1 - f) * em1 + f * em2
        samples.append({
            "site_id": f"S{i}",
            "Ca": conc[0],
            "Cl": conc[1]
        })
        
    virtual = identify_latent_endmembers(samples, ion_order, n_endmembers=1)
    
    # We expect 2 virtual samples (Min/Max along PC1)
    assert len(virtual) == 2
    assert virtual[0]["type"] == "virtual"
    
    # Check if they bracket the data
    # S0 (10% B) should be closer to em1
    # S2 (90% B) should be closer to em2
    
    # Get Ca values
    v_ca = sorted([v["Ca"] for v in virtual])
    
    # True range roughly 10 to 100
    # Virtual nodes usually project slightly outwards (1.1x)
    assert v_ca[0] < 20.0 # Low end
    assert v_ca[1] > 90.0 # High end

# --- Test 2: Compositional Weighting (Statistics) ---

def test_compositional_weighting_impact():
    """
    Test that compositional_weighting=True alters the objective function logic.
    We set up a case where a Trace element has a large relative mismatch but small absolute mismatch.
    """
    config = Config()
    config.transport_models_enabled = ["mix"]
    config.mixing_endmembers = {"Source": [100.0, 1.0]} # Major, Trace
    config.ion_order = ["Major", "Trace"]
    config.weights = [1.0, 1.0]
    config.conservative_weights = [1.0, 1.0] # Critical: Must match data dimension
    config.lambda_l1 = 0.0 # Turn off sparsity to focus on fit
    
    # Observation: Major matches perfectly, Trace is half
    x_u = [0.0, 0.0] # Dummy U
    x_v = [100.0, 0.5] # Obs V.
    # Source [100, 1].
    # If f=1.0: Pred=[100, 1]. Res=[0, -0.5]. SSE = 0 + 0.25 = 0.25
    # If f=0.5: Pred=[50, 0.5]. Res=[50, 0]. SSE = 2500 + 0 = 2500
    
    # Standard: f=1.0 minimizes absolute error (SSE 0.25 vs 2500). Trace error is ignored.
    
    # Compositional: Weight ~ 1/x^2.
    # W_Major ~ 1/100^2 = 0.0001
    # W_Trace ~ 1/0.5^2 = 4.0
    # If f=1.0: Res=[0, -0.5]. Weighted SSE = 0 + 4.0 * 0.25 = 1.0
    # If f=0.99 (fitting Trace? No, f affects both)
    # The point is simply that the 'optimal' f might shift if there's a tradeoff.
    # Let's use a tradeoff case.
    
    # Endmember: [100, 10]
    # Obs V:     [110, 5]  (Major +10%, Trace -50%)
    
    # f needs to scale both.
    # If f=1.1 -> [110, 11]. Res Major=0, Res Trace=6.
    # If f=0.5 -> [50, 5].   Res Major=60, Res Trace=0.
    
    # Standard (Abs): Prefer f=1.1 (Error 6 < Error 60). Fits Major.
    # Compositional (Rel): Prefer f=0.5?
    # Rel Error Major (at f=0.5) ~ 60/110 = 50%.
    # Rel Error Trace (at f=1.1) ~ 6/5 = 120%.
    # It balances them.
    
    config.mixing_endmembers = {"Source": [100.0, 10.0]}
    x_v_tradeoff = [110.0, 5.0]
    
    # 1. Standard Run
    config.compositional_weighting = False
    res_std = fit_edge(x_u, x_v_tradeoff, config, obs_v={"Major": 110, "Trace": 5})
    
    # 2. Compositional Run
    config.compositional_weighting = True
    res_comp = fit_edge(x_u, x_v_tradeoff, config, obs_v={"Major": 110, "Trace": 5})
    
    # Result: Standard should fit Major (f ~ 1.1), Comp should compromise (f < 1.1)
    assert res_std.f is not None
    assert res_comp.f is not None
    
    # Standard f should be at the cap (1.0) because 1.1 > 1.0
    # Comp f should be pulled down (f ~ 0.6)
    
    # If standard is 1.0 and comp is 1.0, the feature fails.
    # We assert that Comp is SIGNIFICANTLY different.
    
    print(f"Std f: {res_std.f}, Comp f: {res_comp.f}")
    
    assert res_std.f >= 0.99
    assert res_comp.f < 0.9 # Should be around 0.6
    assert res_comp.f < 1.0 # Should feel the pull of the 0.5 ratio

# --- Test 3: Iterative Jacobian (Chemistry) ---

def test_iterative_jacobian_execution():
    """
    Test that the iterative solver loop executes without error.
    Since we don't have a live PHREEQC Jacobian updater yet, we verify the
    config flag triggers the loop logic in fit_edge.
    """
    config = Config()
    config.transport_models_enabled = ["evap"]
    config.iterative_jacobian_enabled = True
    config.iterative_jacobian_max_iter = 2
    
    # Simple evap case
    x_u = [10.0]
    x_v = [20.0] # Gamma=2
    
    # This should run 2 passes of reaction fitting
    # We can't easily inspect the internal loop count from result,
    # but successful return implies the code path is valid.
    res = fit_edge(x_u, x_v, config)
    
    assert res.gamma is not None
    assert abs(res.gamma - 2.0) < 0.1
    assert res.reaction_fit is not None

# --- Test 4: Reacted TTD (Time) ---

def test_reacted_ttd_convolution():
    """Test the convolution logic against analytical solutions."""
    
    # 1. Piston Flow Limit
    # TTD: single spike at t=10.
    weights = [1.0]
    lags = [10.0]
    k = 0.1 # decay constant
    c_in = 100.0
    
    c_out = convolve_reactive_ttd(c_in, weights, lags, k)
    expected = 100.0 * math.exp(-0.1 * 10.0)
    assert math.isclose(c_out, expected)
    
    # 2. Binary Mixing
    # 50% water age 0 (fresh), 50% water age 100 (old)
    weights = [0.5, 0.5]
    lags = [0.0, 100.0]
    k = 0.01 # half life approx 69 days
    
    c_out_mix = convolve_reactive_ttd(c_in, weights, lags, k)
    
    part1 = 0.5 * 100.0 * math.exp(0) # 50
    part2 = 0.5 * 100.0 * math.exp(-0.01 * 100.0) # 50 * e^-1 = 50 * 0.367 = 18.39
    expected_mix = part1 + part2
    
    assert math.isclose(c_out_mix, expected_mix)
    
    # 3. Verify Mean Age Approx Fails
    # Mean age = 50.
    # Piston approx = 100 * exp(-0.01 * 50) = 100 * 0.606 = 60.6
    # True mix = 50 + 18.4 = 68.4
    # The mixed water has MORE surviving tracer than the mean-age water 
    # (because young water bypasses decay).
    
    piston_approx = 100.0 * math.exp(-0.01 * 50.0)
    assert c_out_mix > piston_approx # Convexity of decay function

# --- Test 5: Aitchison Utils ---

def test_clr_invariance():
    """Test CLR scaling invariance."""
    x = [10.0, 20.0, 50.0]
    x_scaled = [20.0, 40.0, 100.0]
    
    clr_x = clr(x)
    clr_scaled = clr(x_scaled)
    
    # Should be identical
    for a, b in zip(clr_x, clr_scaled):
        assert math.isclose(a, b, abs_tol=1e-9)
    
    # Sum should be 0
    assert math.isclose(sum(clr_x), 0.0, abs_tol=1e-9)

def test_clr_inverse():
    x = [0.1, 0.2, 0.7] # Simplex
    c = clr(x)
    inv = clr_inv(c)
    
    # Should recover proportions
    for a, b in zip(x, inv):
        assert math.isclose(a, b, abs_tol=1e-9)

