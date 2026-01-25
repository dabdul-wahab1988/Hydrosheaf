"""Scientific visualization module for advanced hydrogeologic concepts.

This module provides publication-quality figures for:
1. Transit Time Distributions (TTDs) - The 'memory' of the aquifer.
2. Contaminant Breakthrough Curves - Convolution of input history.
3. Bayesian Posterior Ridges - Uncertainty in parameter estimates.
4. Reactive Transport Validation - Modeled vs Observed comparisons.
"""

from typing import List, Optional
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats

from ..inference.edge_fit import EdgeResult
from .utils import PlotConfig, save_with_metadata


def _check_mpl():
    pass


def plot_ttd_kernel(
    edge_results: List[EdgeResult], 
    max_tau: float = 365.0, 
    path: Optional[str] = None,
    config: Optional[PlotConfig] = None
) -> None:
    """Plot the Transit Time Distribution (TTD) kernel g(tau).

    Visualizes the probability density of travel times for edges with temporal data.
    """
    if config is None:
        config = PlotConfig(figsize=(10, 6))
    config.apply()

    # Filter for edges with temporal TTD data
    valid_edges = [
        r
        for r in edge_results
        if r.temporal_residence_time_details
        and "ttd_g" in r.temporal_residence_time_details
    ]

    if not valid_edges:
        print("Warning: No edges with TTD data found for plotting.")
        return

    fig, ax = plt.subplots(figsize=config.figsize)

    # Plot each edge's kernel
    try:
         colors = plt.cm.viridis(np.linspace(0, 1, len(valid_edges)))
    except Exception:
         colors = plt.cm.get_cmap("viridis")(np.linspace(0, 1, len(valid_edges)))

    for i, res in enumerate(valid_edges):
        details = res.temporal_residence_time_details
        tau_grid = np.array(details.get("ttd_tau_days", []))
        g = np.array(details.get("ttd_g", []))

        if len(tau_grid) == 0 or len(g) == 0:
            continue

        label = f"{res.edge_id} (Mean={res.temporal_residence_time_days:.1f}d)"
        ax.plot(tau_grid, g, lw=2, label=label, color=colors[i], alpha=0.8)
        ax.fill_between(tau_grid, g, alpha=0.1, color=colors[i])

    ax.set_xlabel("Transit Time $\\tau$ (days)", fontsize=12*config.font_scale)
    ax.set_ylabel("Probability Density $g(\\tau)$", fontsize=12*config.font_scale)
    ax.set_title("Aquifer Memory: Transit Time Distributions", fontsize=14*config.font_scale)
    ax.set_xlim(0, max_tau)
    ax.grid(True, alpha=0.3)

    if len(valid_edges) < 15:
        ax.legend(frameon=True, fontsize=9*config.font_scale)

    plt.tight_layout()

    if path:
        save_with_metadata(fig, path, config, extra_metadata={"n_edges": len(valid_edges)})
        plt.close(fig)
    else:
        plt.show()


def plot_breakthrough(
    edge_results: List[EdgeResult], 
    tracer: str = "Cl", 
    path: Optional[str] = None,
    config: Optional[PlotConfig] = None
) -> None:
    """Plot input vs. output breakthrough curves (convolution visualization).

    Shows how the input signal (Upstream) is smoothed and delayed to match the
    output signal (Downstream).
    """
    if config is None:
        config = PlotConfig() # Default figsize is handled per subplot logic below if needed, or overridden
    config.apply()

    # Filter for edges with temporal TTD data
    valid_edges = [
        r
        for r in edge_results
        if r.temporal_residence_time_details
        and "convolved_series" in r.temporal_residence_time_details
    ]

    if not valid_edges:
        print("Warning: No edges with temporal convolution data found.")
        return

    # Limit to first few edges to avoid clutter, or create subplots
    n_plots = min(len(valid_edges), 4)
    
    # Calculate figure size based on plots
    figsize = (config.figsize[0], 3 * n_plots)
    
    fig, axes = plt.subplots(n_plots, 1, figsize=figsize, sharex=False)
    if n_plots == 1:
        axes = [axes]

    for i in range(n_plots):
        res = valid_edges[i]
        ax = axes[i]
        details = res.temporal_residence_time_details

        # Extract series
        timestamps = details.get("timestamps_str", [])  # Assuming stored as ISO strings
        u_series = details.get("u_series_norm", [])  # Normalized usually
        v_series = details.get("v_series_norm", [])
        pred_series = details.get("convolved_series", [])

        if not timestamps:
            continue

        # Parse dates
        try:
            dates = [np.datetime64(ts) for ts in timestamps]
        except (TypeError, ValueError):
            dates = np.arange(len(u_series))

        ax.plot(dates, u_series, "k--", alpha=0.4, label="Input (Upstream)", lw=1)
        ax.plot(dates, v_series, "bo", alpha=0.5, ms=4, label="Observed (Downstream)")
        ax.plot(dates, pred_series, "r-", lw=2, label="Model (TTD Convolution)")

        ax.set_title(f"Edge: {res.edge_id} ($R^2$={details.get('ttd_r2', 0):.2f})", fontsize=10*config.font_scale)
        ax.set_ylabel(f"Normalized {tracer}")
        ax.grid(True, alpha=0.3)
        if i == 0:
            ax.legend(loc="best", fontsize=9*config.font_scale)

    plt.tight_layout()

    if path:
        save_with_metadata(fig, path, config, extra_metadata={"edges": [e.edge_id for e in valid_edges[:n_plots]]})
        plt.close(fig)
    else:
        plt.show()


def plot_posterior_ridges(
    edge_result: EdgeResult, 
    param_type: str = "extents", 
    path: Optional[str] = None,
    config: Optional[PlotConfig] = None
) -> None:
    """Plot Bayesian posterior distributions (Ridge Plot).

    Args:
        edge_result: The edge result containing uncertainty samples.
        param_type: 'extents' (reaction rates) or 'transport' (gamma/f).
    """
    if config is None:
        config = PlotConfig()
    config.apply()

    if not edge_result.uncertainty or not getattr(edge_result.uncertainty, "extents_samples", None):
        print(f"Warning: No Bayesian samples found for edge {edge_result.edge_id}")
        return

    samples = np.array(
        edge_result.uncertainty.extents_samples
    )  # (n_samples, n_reactions)
    labels = edge_result.reaction_labels

    if param_type == "transport":
        # Plot gamma or f
        samples = np.array(edge_result.uncertainty.gamma_samples).reshape(-1, 1)
        labels = ["Gamma (Evap Factor)"]

    n_params = samples.shape[1]

    # Filter only active params (non-zero mean or high variance)
    active_indices = []
    for i in range(n_params):
        if np.abs(np.mean(samples[:, i])) > 1e-6 or np.std(samples[:, i]) > 1e-6:
            active_indices.append(i)

    if not active_indices:
        print("No active parameters to plot.")
        return

    n_active = len(active_indices)
    
    # Dynamic height
    figsize = (8, max(4, n_active * 1.2))
    
    fig, axes = plt.subplots(
        n_active, 1, figsize=figsize, sharex=False
    )
    if n_active == 1:
        axes = [axes]

    try:
        colors = plt.cm.plasma(np.linspace(0, 0.8, n_active))
    except Exception:
        colors = plt.cm.get_cmap("plasma")(np.linspace(0, 0.8, n_active))


    for idx, param_idx in enumerate(active_indices):
        ax = axes[idx]
        data = samples[:, param_idx]
        lbl = labels[param_idx]

        # Kernel Density Estimate
        density = stats.gaussian_kde(data)
        xs = np.linspace(min(data), max(data), 200)
        ys = density(xs)

        ax.plot(xs, ys, color=colors[idx], lw=2)
        ax.fill_between(xs, ys, color=colors[idx], alpha=0.3)

        # Stats annotations
        mean_val = np.mean(data)
        ci_low = np.percentile(data, 2.5)
        ci_high = np.percentile(data, 97.5)

        ax.axvline(mean_val, color="k", ls="--", alpha=0.5, lw=1)
        ax.text(
            0.02, 0.8, f"{lbl}", transform=ax.transAxes, fontweight="bold", fontsize=10*config.font_scale
        )
        ax.text(
            0.98,
            0.8,
            f"Mean: {mean_val:.4f}\n95% CI: [{ci_low:.4f}, {ci_high:.4f}]",
            transform=ax.transAxes,
            ha="right",
            fontsize=9*config.font_scale,
            va="top",
        )

        ax.set_yticks([])
        ax.spines["left"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)

    ax.set_xlabel("Parameter Value (mmol/L or factor)")
    fig.suptitle(
        f"Bayesian Posterior Distributions: {edge_result.edge_id}", fontsize=14*config.font_scale
    )
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.5)

    if path:
        save_with_metadata(fig, path, config, extra_metadata={"edge_id": edge_result.edge_id})
        plt.close(fig)
    else:
        plt.show()


def plot_reactive_transport_validation(
    edge_results: List[EdgeResult], 
    path: Optional[str] = None,
    config: Optional[PlotConfig] = None
) -> None:
    """Plot modeled vs observed concentrations from forward RT validation."""
    if config is None:
        config = PlotConfig(figsize=(12, 5))
    config.apply()
    
    # Filter validated edges
    valid_edges = [
        r
        for r in edge_results
        if r.rt_validation and r.rt_validation.get("rmse") is not None
    ]

    if not valid_edges:
        print("Warning: No reactive transport validation data found.")
        return

    rmses = [r.rt_validation["rmse"] for r in valid_edges]
    nses = [r.rt_validation.get("nse", -999) for r in valid_edges]
    ids = [r.edge_id for r in valid_edges]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=config.figsize)

    # RMSE Bar Plot
    ax1.bar(ids, rmses, color="indianred", alpha=0.7)
    ax1.set_title("Model Error (RMSE)")
    ax1.set_ylabel("RMSE (mmol/L)")
    ax1.tick_params(axis="x", rotation=45)
    ax1.axhline(0, color="k", lw=0.5)

    # NSE Scatter Plot
    ax2.scatter(ids, nses, color="teal", s=100, edgecolor="k")
    ax2.set_title("Model Efficiency (NSE)")
    ax2.set_ylabel("Nash-Sutcliffe Efficiency")
    ax2.axhline(0, color="r", ls="--", label="Mean Predictor")
    ax2.axhline(0.5, color="g", ls="--", label="Good Fit")
    ax2.set_ylim(-1, 1)
    ax2.tick_params(axis="x", rotation=45)
    ax2.legend()

    plt.tight_layout()

    if path:
        save_with_metadata(fig, path, config, extra_metadata={"n_edges": len(valid_edges)})
        plt.close(fig)
    else:
        plt.show()
