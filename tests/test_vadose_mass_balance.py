
import pytest
import numpy as np
from datetime import datetime, timedelta
from hydrosheaf.vadose.richards1d import run_richards_column
from hydrosheaf.vadose.contracts import VadoseProfile, VadoseForcingSample, VadoseLayer, VadoseRunConfig

def test_richards_mass_balance_closure():
    """
    PhD-Level Verification: Ensure mass balance closes in 1D Richards solver.
    storage_change = (Inflow - Outflow - Sink) * dt
    """
    # 1. Setup a simple profile
    layer = VadoseLayer(
        thickness_m=2.0,
        texture="silt_loam",
        ks_m_day=0.1,
        theta_s=0.45,
        theta_r=0.06,
        alpha_1_m=1.9,
        n=1.3,
        pore_connectivity=0.5
    )
    profile = VadoseProfile(
        profile_id="test",
        layers=[layer],
        depth_m=2.0,
        root_depth_m=0.5
    )
    
    # 2. Forcing: Constant rain for 5 days
    t0 = datetime(2023, 1, 1)
    forcing = [
        VadoseForcingSample(timestamp=t0 + timedelta(days=i), 
                            precipitation_mm_day=10.0, 
                            et0_mm_day=2.0)
        for i in range(6)
    ]
    
    # 3. Run simulation
    config = VadoseRunConfig(dz_m=0.1, max_picard_iter=50)
    sim = run_richards_column(profile, forcing, config=config, initial_psi_m=-1.0)
    
    # 4. Verify Mass Balance across all steps
    # mass_balance_error_m is computed internally as:
    # (storage_np1 - storage_n) - (inflow - outflow - sink) * dt
    errors = [d.mass_balance_error_m for d in sim.diagnostics]
    
    # Tolerance for mass balance: should be very small (numerical precision)
    # Given Picard iteration tolerance, we expect < 1e-4 m
    max_error = np.max(np.abs(errors))
    
    print(f"Max Mass Balance Error: {max_error:.2e} m")
    assert max_error < 1e-4, f"Mass balance error too high: {max_error:.2e}"

def test_richards_bottom_flux_dirichlet():
    """
    Verify that bottom flux is correctly computed with constant head BC.
    """
    layer = VadoseLayer(thickness_m=1.0, texture="sand", ks_m_day=1.0)
    profile = VadoseProfile(profile_id="bc_test", layers=[layer], depth_m=1.0)
    
    t0 = datetime(2023, 1, 1)
    forcing = [
        VadoseForcingSample(
            timestamp=t0 + timedelta(days=i),
            precipitation_mm_day=0.0,
            et0_mm_day=0.0
        ) for i in range(3)
    ]
    
    # Prescribe water table at 1.5m (bottom is at 1.0m)
    # Expected suction at bottom: psi ≈ -0.5m
    wt_data = [(t0 + timedelta(days=i), 1.5) for i in range(3)]
    
    config = VadoseRunConfig(lower_boundary="constant_head_from_wt", dz_m=0.1)
    sim = run_richards_column(profile, forcing, config=config, water_table_depth_m=wt_data)
    
    # Bottom flux should NOT be exactly K (which is 1.0) because of head gradient
    # With unit gradient (gravity only) and saturated K=1, q=1.
    # Here we have suction pulling UP (psi=-0.5 at bot, likely higher above).
    # Recharge should be less than K if table is deep.
    for q in sim.recharge_m_day:
        assert q != 1.0, "Bottom flux failed to account for Dirichlet BC gradient"
        assert q < 1.0, "Expected reduced recharge due to capillary suction from deep WT"

if __name__ == "__main__":
    pytest.main([__file__])
