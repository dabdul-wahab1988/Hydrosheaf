import numpy as np
from hydrosheaf.vadose.ttd import mixture_ttd_from_series

def test_dual_domain_ttd_effect():
    # Setup: Constant advective travel time of 100 days
    tau_series = [100.0] * 10
    weights = [1.0] * 10
    
    # 1. Base case (Single domain)
    grid1, pdf1, _ = mixture_ttd_from_series(
        tau_series, weights,
        grid_dt_days=1.0,
        max_lag_days=200.0,
        cv=0.1,
        preferential_flow_fraction=0.0
    )
    
    # 2. Dual domain case (20% preferential flow, 10x faster)
    # Fast path should be around 100/10 = 10 days
    grid2, pdf2, _ = mixture_ttd_from_series(
        tau_series, weights,
        grid_dt_days=1.0,
        max_lag_days=200.0,
        cv=0.1,
        preferential_flow_fraction=0.2,
        preferential_velocity_factor=10.0
    )
    
    # Check early arrival
    # In single domain (mean 100, cv 0.1 -> std 10), prob at t=10 should be near zero.
    idx_10 = 10 
    prob1_early = sum(pdf1[:20]) # Sum prob up to 20 days
    prob2_early = sum(pdf2[:20])
    
    print(f"Prob early (single): {prob1_early}")
    print(f"Prob early (dual): {prob2_early}")
    
    assert prob1_early < 1e-3, "Single domain should have negligible early arrival"
    assert prob2_early > 0.1, "Dual domain should have significant early arrival (approx 20%)"
    
    # Check that total mass is preserved (should sum to ~1/dt if PDF, or sum*dt = 1)
    # mixture_ttd_from_series returns PDF values.
    # sum(pdf) * dt should be 1.0
    integral1 = sum(pdf1) * 1.0
    integral2 = sum(pdf2) * 1.0
    
    assert abs(integral1 - 1.0) < 1e-2
    assert abs(integral2 - 1.0) < 1e-2
