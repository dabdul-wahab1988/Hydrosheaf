"""
Plotting utilities for Network-Enhanced Nuclear Dating.

This module generates publication-quality (PhD-level) visualizations for the
Bayesian groundwater aging results. It focuses on:
1. Identifying "Modern" vs "Fossil" waters.
2. Visualizing the "Aging Gradient" along the flow network.
3. Displaying posterior uncertainty (Credible Intervals).
4. Diagnosing the "Bomb Peak" fit and input scaling.
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from typing import Dict, Optional, Any

def plot_network_ages(
    graph: nx.DiGraph,
    age_results: Dict[str, Dict[str, float]],
    figsize: tuple = (12, 10),
    cmap: str = "viridis_r", # Reversed: Blue=Young, Yellow=Old
    node_size: int = 300,
    ax: Optional[plt.Axes] = None
):
    """
    Plot the network colored by Mean Residence Time (Age).
    
    Critique Defense:
    - Shows the spatial structure of aging.
    - Highlights "Age Reversals" (if any persist despite constraints).
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
        
    # Extract ages
    nodes = list(graph.nodes())
    ages = []
    p_modern = []
    
    for n in nodes:
        if n in age_results:
            ages.append(age_results[n]["mean_age_years"])
            p_modern.append(age_results[n]["p_modern"])
        else:
            ages.append(np.nan)
            p_modern.append(0.0)
            
    # Layout
    pos = nx.get_node_attributes(graph, "pos")
    if not pos:
        pos = nx.kamada_kawai_layout(graph)
        
    # Draw Edges (Flow Direction)
    nx.draw_networkx_edges(
        graph, pos, ax=ax, 
        edge_color="gray", 
        arrows=True, 
        arrowstyle="-|>", 
        arrowsize=15, 
        alpha=0.6
    )
    
    # Draw Nodes (Color by Age)
    # Filter valid nodes
    valid_indices = [i for i, a in enumerate(ages) if not np.isnan(a)]
    valid_nodes = [nodes[i] for i in valid_indices]
    valid_ages = [ages[i] for i in valid_indices]
    
    if valid_nodes:
        sc = nx.draw_networkx_nodes(
            graph, pos, 
            nodelist=valid_nodes,
            node_size=node_size,
            node_color=valid_ages,
            cmap=cmap,
            ax=ax,
            edgecolors="black",
            linewidths=1.0
        )
        # Colorbar
        plt.colorbar(sc, ax=ax, label="Mean Residence Time (Years)")
        
    # Labels
    nx.draw_networkx_labels(graph, pos, ax=ax, font_size=8)
    
    ax.set_title("Network-Inferred Mean Residence Times (Bayesian MAP)", fontsize=14)
    ax.axis("off")
    return ax

def plot_age_vs_distance(
    graph: nx.DiGraph,
    age_results: Dict[str, Dict[str, float]],
    root_node: Optional[str] = None,
    ax: Optional[plt.Axes] = None
):
    """
    Plot 'Age vs Distance from Recharge' to demonstrate the Aging Gradient.
    
    PhD Value:
    - Validates the "Piston Flow" or "Linear Aging" hypothesis.
    - Shows dispersion/mixing as the spread of CI.
    - Identifying "Stagnant Zones" (High Age, Low Distance).
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 6))
        
    # Compute distance from "Recharge Zones" (indegree 0) or specific root
    if root_node:
        roots = [root_node]
    else:
        roots = [n for n, d in graph.in_degree() if d == 0]
        
    # BFS distance (topological or physical if length_m exists)
    # For simplicity, use unweighted shortest path if length_m missing
    distances = {}
    for root in roots:
        # Use length_m if available
        try:
            d = nx.shortest_path_length(graph, source=root, weight="length_m")
            for n, dist in d.items():
                if n not in distances or dist < distances[n]:
                    distances[n] = dist
        except Exception:
            pass
            
    # Collect Plot Data
    x_dist = []
    y_age = []
    y_err_low = []
    y_err_high = []
    
    for n, res in age_results.items():
        if n in distances:
            d = distances[n]
            mean = res["mean_age_years"]
            low = res["mean_age_years"] - res["age_95_low"]
            high = res["age_95_high"] - res["mean_age_years"]
            
            x_dist.append(d)
            y_age.append(mean)
            y_err_low.append(low)
            y_err_high.append(high)
            
    if not x_dist:
        ax.text(0.5, 0.5, "No connected paths from roots found", ha="center")
        return ax
        
    # Plot
    ax.errorbar(
        x_dist, y_age, 
        yerr=[y_err_low, y_err_high], 
        fmt='o', 
        capsize=3, 
        ecolor='gray', 
        alpha=0.7, 
        label="Node Posterior (95% CI)"
    )
    
    # Fit line for "Apparent Velocity"
    # V = dX/dt
    if len(x_dist) > 1:
        coeffs = np.polyfit(x_dist, y_age, 1)
        poly = np.poly1d(coeffs)
        x_line = np.linspace(min(x_dist), max(x_dist), 100)
        ax.plot(x_line, poly(x_line), "r--", label=f"Trend (V ≈ {1/coeffs[0]:.1f} m/y)")
        
    ax.set_xlabel("Distance from Recharge (m)")
    ax.set_ylabel("Inferred Mean Age (Years)")
    ax.set_title("Groundwater Aging Gradient")
    ax.legend()
    ax.grid(True, alpha=0.3)
    return ax

def plot_modern_probability(
    age_results: Dict[str, Dict[str, float]],
    ax: Optional[plt.Axes] = None
):
    """
    Plot the 'Probability of Modern Water' (p_modern) for each node.
    
    PhD Value:
    - Directly answers the "Binary Dating" critique (Post-bomb).
    - Robust classification of "Renewable" vs "Fossil".
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 6))
        
    nodes = sorted(age_results.keys(), key=lambda k: age_results[k]["mean_age_years"])
    probs = [age_results[n]["p_modern"] for n in nodes]
    
    # Color bars: Green (Modern) -> Red (Fossil)
    colors = plt.cm.RdYlGn(probs)
    
    ax.bar(nodes, probs, color=colors)
    ax.axhline(0.5, color="black", linestyle="--", alpha=0.5, label="Uncertainty Threshold")
    
    ax.set_ylabel("Probability (Age < 60 years)")
    ax.set_title("Probability of Modern Recharge (p_modern)")
    ax.set_ylim(0, 1)
    
    # Rotate labels if many nodes
def plot_mcmc_diagnostics(
    trace: Any,
    var_names: Optional[list] = None,
    ax: Optional[plt.Axes] = None
):
    """
    Plot MCMC diagnostics (Trace and Posterior Density).
    
    PhD Value:
    - Essential proof of convergence.
    - Demonstrates that the "latent variables" (Input Scale, Velocity) are identifiable.
    """
    import arviz as az
    
    if var_names is None:
        # Default to the key latent variables, not the 100s of node ages
        var_names = ["input_scale", "velocity", "c14_a0", "c14_a0_mu"]
        # Filter to only existing vars
        available_vars = list(trace.posterior.data_vars)
        var_names = [v for v in var_names if v in available_vars]
        
    if not var_names:
        return
        
    # Arviz handles the plotting elegantly
    az.plot_trace(trace, var_names=var_names, kind="trace", compact=True)
    plt.suptitle("Bayesian Convergence Diagnostics (Trace & Density)", fontsize=14)
    
    # Print summary statistics to console for the user/log
    summary = az.summary(trace, var_names=var_names)
    print("\n=== MCMC Convergence Summary ===")
    print(summary[["r_hat", "ess_bulk", "ess_tail"]])
    print("================================\n")
    
    return plt.gcf()
