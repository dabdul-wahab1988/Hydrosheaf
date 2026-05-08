import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

BENCHMARK_ROOT = Path(__file__).resolve().parent.parent

def plot_e3_diagnostics():
    print("Generating Figure S3: PHREEQC Forward Diagnostics...")
    
    forward_path = BENCHMARK_ROOT / "results" / "phreeqc_forward_validation.csv"
    if not forward_path.exists():
        print(f"Error: {forward_path} not found.")
        return
        
    forward = pd.read_csv(forward_path)
    
    # Filter out NaNs for NSE plot
    forward_clean = forward.dropna(subset=['nse'])
    
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    
    # A. RMSE Histogram
    axes[0].hist(forward["rmse_mmolL"], bins=25, color="#2563eb", alpha=0.75, edgecolor='white')
    axes[0].axvline(forward["rmse_mmolL"].median(), color="red", linestyle="--", label=f"Median: {forward['rmse_mmolL'].median():.2f}")
    axes[0].set_title("A. Forward residual RMSE", fontsize=12, weight='bold')
    axes[0].set_xlabel("RMSE (mmol/L)")
    axes[0].set_ylabel("Frequency")
    axes[0].legend(frameon=False)
    
    # B. Model Efficiency (NSE)
    axes[1].scatter(forward_clean["rmse_mmolL"], forward_clean["nse"], s=15, alpha=0.3, color="#0f766e")
    axes[1].axhline(0.5, color="red", linestyle="--", lw=1.2, label="Acceptable fit (NSE=0.5)")
    axes[1].set_title("B. Forward model efficiency", fontsize=12, weight='bold')
    axes[1].set_xlabel("RMSE (mmol/L)")
    axes[1].set_ylabel("NSE (Nash-Sutcliffe)")
    axes[1].set_ylim(-0.5, 1.05)
    axes[1].legend(loc='lower right', frameon=False)
    
    # C. Mineral Saturation Indices
    si_cols = [c for c in forward.columns if c.startswith('si_')]
    if si_cols:
        si_data = forward[si_cols].dropna()
        # Rename for plot labels
        si_data.columns = [c.replace('si_', '').capitalize() for c in si_data.columns]
        
        si_data.boxplot(ax=axes[2], patch_artist=True, 
                        boxprops=dict(facecolor="#f1f5f9", color="#334155"),
                        medianprops=dict(color="red"))
        axes[2].axhline(0, color="black", linestyle="-", lw=0.8)
        axes[2].axhline(0.2, color="gray", linestyle=":", lw=0.8)
        axes[2].axhline(-0.2, color="gray", linestyle=":", lw=0.8)
        axes[2].set_title("C. Mineral Saturation Indices", fontsize=12, weight='bold')
        axes[2].set_ylabel("SI")
        axes[2].tick_params(axis='x', rotation=45)
    
    fig.tight_layout()
    output_path = BENCHMARK_ROOT / "figures" / "figure_s3_phreeqc_forward_diagnostics_live.png"
    fig.savefig(output_path, dpi=300)
    print(f"Figure saved to {output_path}")

if __name__ == "__main__":
    plot_e3_diagnostics()
