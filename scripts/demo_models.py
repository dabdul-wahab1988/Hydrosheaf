
def coupled_model(params):
    # 1. Physics: Calculate time
    v = params.get("velocity", 0.1)
    dist = 100.0
    time_days = dist / v
    
    # 2. Chemistry: Calculate Dissolution over that time
    # We also calibrate the rate constant 'log_k_calcite'
    # PEST passes log_k, we need 10^val
    log_k = params.get("log_k_calcite", -6.0)
    k = 10**log_k
    A = 1.0 # Fixed surface area
    
    # Simple dissolution model: C = C0 + rate * time
    # True k = 1e-5. True v = 0.5 (Time=200).
    # C_final = 0 + 1e-5 * 1.0 * 200 * scaling
    
    scaling = 1e3 # Make units visible
    c_final = k * A * time_days * scaling
    
    return {
        "obs_travel_time": time_days,
        "obs_calcium": c_final
    }
