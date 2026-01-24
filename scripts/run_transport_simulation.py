import numpy as np
import matplotlib.pyplot as plt
import json
import os
from scipy.special import erfc


def ogata_banks_decay(C0, x, t, v, D, k):
    """
    Analytical solution for 1D advection-dispersion with first-order decay.
    Ref: Ogata & Banks (1961), Bear (1979)

    C0: Initial concentration at source
    x: Distance from source
    t: Time
    v: Pore water velocity (q/porosity)
    D: Dispersion coefficient (dispersivity * v + diffusion)
    k: First order decay constant (denitrification)
    """
    if t <= 0:
        return 0.0

    # Check for zero dispersion to avoid division by zero
    if D <= 1e-9:
        if x < v * t:
            return C0 * np.exp(-k * t)
        else:
            return 0.0

    # Transformation for decay
    # For steady state or transport with decay, standard solution involves
    # exponential terms.
    # Simplified approximation for transport + reaction:
    # C(x,t) = (C0/2) * exp(v*x/2D) * [ exp(-U*x/2D) * erfc((x - U*t)/(2*sqrt(Dt))) + ... ]
    # where U = sqrt(v^2 + 4*k*D)

    U = np.sqrt(v**2 + 4 * k * D)

    term1 = np.exp((v - U) * x / (2 * D)) * erfc((x - U * t) / (2 * np.sqrt(D * t)))
    term2 = np.exp((v + U) * x / (2 * D)) * erfc((x + U * t) / (2 * np.sqrt(D * t)))

    return (C0 / 2) * (term1 + term2)


def main():
    # Paths
    base_dir = "analysis_results_complete/objective4_transport"
    params_file = os.path.join(base_dir, "flopy_parameters.json")
    output_plot = os.path.join(base_dir, "transport_simulation_analytical.png")

    # Load parameters
    with open(params_file, "r") as f:
        params = json.load(f)

    print("Loading parameters from:", params_file)
    print(json.dumps(params, indent=2))

    # Extract physical parameters
    hk = params.get("hk", 1.0)
    porosity = params.get("porosity", 0.25)
    alpha_L = params.get("dispersivity", 10.0)
    k_decay = params.get("denitrification_k", 0.001)
    C0 = params.get("source_conc", 50.0)

    # Hydraulic gradient (assumed based on typical network analysis or set to generic)
    # Since we computed a lag time of ~75 days for some edges, let's estimate velocity
    # v = distance / time. If delr=2m and we have ~50 cells, length=100m.
    # We'll assume a hydraulic gradient i that makes physical sense or use the derived v from lag time if available
    # For this simulation, we'll calculate v based on Darcy's law with an assumed gradient i
    i = 0.01  # Assumed gradient
    v_darcy = hk * i
    v_pore = v_darcy / porosity

    # Dispersion coefficient D = alpha_L * v_pore + D_molecular (ignore molecular)
    D = alpha_L * v_pore

    print(f"\nSimulation Parameters:")
    print(f"  Pore Velocity (v): {v_pore:.4f} m/d")
    print(f"  Dispersion Coeff (D): {D:.4f} m^2/d")
    print(f"  Decay Rate (k): {k_decay:.4f} /d")
    print(f"  Half-life: {np.log(2)/k_decay:.1f} days")

    # Setup domains
    x_max = 100.0
    x_vals = np.linspace(0, x_max, 100)

    times = [50, 100, 365, 730]  # Days

    # Plotting
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

    # Subplot 1: Spatial Profiles at different times
    colors = plt.cm.viridis(np.linspace(0, 1, len(times)))
    for t, color in zip(times, colors):
        c_vals = [ogata_banks_decay(C0, x, t, v_pore, D, k_decay) for x in x_vals]
        ax1.plot(x_vals, c_vals, lw=2, label=f"t = {t} days", color=color)

    # Also plot conservative tracer (no decay) for the last time step for comparison
    c_vals_no_decay = [
        ogata_banks_decay(C0, x, times[-1], v_pore, D, 0) for x in x_vals
    ]
    ax1.plot(
        x_vals,
        c_vals_no_decay,
        lw=2,
        linestyle="--",
        color="gray",
        alpha=0.6,
        label=f"Conservative (t={times[-1]})",
    )

    ax1.set_title("Spatial Concentration Profiles", fontsize=14)
    ax1.set_xlabel("Distance (m)", fontsize=12)
    ax1.set_ylabel("Nitrate Concentration (mg/L)", fontsize=12)
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax1.set_ylim(0, C0 * 1.05)

    # Subplot 2: Breakthrough Curves at specific observation points
    obs_points = [20, 50, 80]  # meters
    t_vals = np.linspace(0, 1000, 200)
    colors_obs = ["#e74c3c", "#f39c12", "#2ecc71"]

    for x_obs, color in zip(obs_points, colors_obs):
        c_vals = [ogata_banks_decay(C0, x_obs, t, v_pore, D, k_decay) for t in t_vals]
        ax2.plot(t_vals, c_vals, lw=2, label=f"x = {x_obs} m", color=color)

    ax2.set_title("Breakthrough Curves", fontsize=14)
    ax2.set_xlabel("Time (days)", fontsize=12)
    ax2.set_ylabel("Nitrate Concentration (mg/L)", fontsize=12)
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_plot, dpi=300)
    print(f"\nSimulation plot saved to: {output_plot}")


if __name__ == "__main__":
    main()
