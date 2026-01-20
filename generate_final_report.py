import json
import pandas as pd
import os
import datetime


def generate_html_report():
    base_dir = "analysis_results_complete"
    summary_file = os.path.join(base_dir, "master_summary.json")
    output_file = os.path.join(base_dir, "final_project_report.html")

    # Load summary
    with open(summary_file, "r") as f:
        summary = json.load(f)

    # HTML Header
    html = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>Hydrosheaf Groundwater Analysis Report</title>
        <style>
            body {{ font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; line-height: 1.6; color: #333; max_width: 1200px; margin: 0 auto; padding: 20px; }}
            h1 {{ color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 10px; }}
            h2 {{ color: #2980b9; margin-top: 30px; border-bottom: 1px solid #eee; padding-bottom: 5px; }}
            h3 {{ color: #16a085; }}
            .container {{ display: flex; flex-wrap: wrap; gap: 20px; }}
            .box {{ flex: 1; min-width: 300px; background: #f9f9f9; padding: 15px; border-radius: 5px; box-shadow: 0 2px 5px rgba(0,0,0,0.05); }}
            .stat-value {{ font-size: 1.2em; font-weight: bold; color: #2c3e50; }}
            img {{ max_width: 100%; height: auto; border: 1px solid #ddd; border-radius: 5px; margin: 10px 0; }}
            table {{ border-collapse: collapse; width: 100%; margin: 15px 0; }}
            th, td {{ text-align: left; padding: 8px; border-bottom: 1px solid #ddd; }}
            th {{ background-color: #f2f2f2; }}
            .alert {{ background-color: #ffebee; color: #c62828; padding: 10px; border-radius: 5px; }}
            .success {{ background-color: #e8f5e9; color: #2e7d32; padding: 10px; border-radius: 5px; }}
            footer {{ margin-top: 50px; font-size: 0.8em; color: #777; text-align: center; border-top: 1px solid #eee; padding-top: 20px; }}
        </style>
    </head>
    <body>
    
    <h1>Hydrosheaf Groundwater Analysis: Final Report</h1>
    <p><strong>Date:</strong> {datetime.datetime.now().strftime("%Y-%m-%d")}</p>
    <p><strong>Project:</strong> Synthetic Groundwater Dataset Analysis</p>
    
    <div class="success">
    <h3>Executive Summary</h3>
    <p>This comprehensive analysis addresses the four key research objectives regarding nitrate transport and fate in the study area. 
    Key findings indicate a strong link between agricultural inputs and groundwater nitrate, with a mean transport lag time of {summary['objective3_recharge']['estimated_lag_days']} days.
    Isotopic evidence confirms synthetic fertilizer as the dominant source (~{summary['objective2_sources']['fertilizer_contribution_pct']:.1f}%), 
    while transport modeling highlights specific vulnerability zones.</p>
    </div>

    <h2>1. Vadose Zone Nitrate Transport</h2>
    <div class="container">
        <div class="box">
            <h3>Key Statistics</h3>
            <ul>
                <li><strong>30cm Depth Mean NO3:</strong> <span class="stat-value">{summary['objective1_vadose']['depth_30cm_mean_no3']:.1f} mg/kg</span></li>
                <li><strong>60cm Depth Mean NO3:</strong> <span class="stat-value">{summary['objective1_vadose']['depth_60cm_mean_no3']:.1f} mg/kg</span></li>
                <li><strong>Wet Season Mean:</strong> <span class="stat-value">{summary['objective1_vadose']['wet_season_mean_no3']:.1f} mg/kg</span></li>
                <li><strong>Dry Season Mean:</strong> <span class="stat-value">{summary['objective1_vadose']['dry_season_mean_no3']:.1f} mg/kg</span></li>
            </ul>
            <p>Significant attenuation is observed between 30cm and 60cm depths. Wet season concentrations are markedly higher, suggesting flush events drive transport.</p>
        </div>
        <div class="box">
            <h3>Visual Analysis</h3>
            <img src="objective1_vadose/vadose_no3_dynamics.png" alt="Vadose Zone Dynamics">
        </div>
    </div>

    <h2>2. Nitrate Source Identification</h2>
    <div class="container">
        <div class="box">
            <h3>Source Apportionment</h3>
            <table>
                <tr><th>Source</th><th>Contribution (%)</th></tr>
                <tr><td>Synthetic Fertilizer</td><td>{summary['objective2_sources']['fertilizer_contribution_pct']:.1f}%</td></tr>
                <tr><td>Soil Organic N</td><td>{summary['objective2_sources']['soil_n_contribution_pct']:.1f}%</td></tr>
                <tr><td>Manure</td><td>{summary['objective2_sources']['manure_contribution_pct']:.1f}%</td></tr>
                <tr><td>Atmospheric</td><td>{summary['objective2_sources']['atmospheric_contribution_pct']:.1f}%</td></tr>
            </table>
            <p>Denitrification was detected in <strong>{summary['objective2_sources']['denitrification_fraction']*100:.1f}%</strong> of samples.</p>
        </div>
        <div class="box">
            <h3>Isotopic Evidence</h3>
            <img src="objective2_sources/nitrate_source_analysis.png" alt="Nitrate Isotopes">
        </div>
    </div>

    <h2>3. Recharge Pathways & Lag Times</h2>
    <div class="container">
        <div class="box">
            <h3>Recharge Characteristics</h3>
            <ul>
                <li><strong>Mean d-excess:</strong> {summary['objective3_recharge']['mean_d_excess']:.2f} ‰</li>
                <li><strong>Est. Transport Lag:</strong> <span class="stat-value">{summary['objective3_recharge']['estimated_lag_days']:.0f} days</span></li>
            </ul>
            <p>The ~2.5 month lag time suggests a relatively fast transport pathway from the vadose zone to the shallow groundwater monitoring points.</p>
        </div>
        <div class="box">
            <h3>Recharge Analysis</h3>
            <img src="objective3_recharge/recharge_analysis.png" alt="Recharge Analysis">
        </div>
    </div>

    <h2>4. Transport Model Framework</h2>
    <div class="container">
        <div class="box">
            <h3>Model Parameters</h3>
            <ul>
                <li><strong>Dispersivity:</strong> {summary['objective4_transport']['flopy_params']['dispersivity']:.2f} m</li>
                <li><strong>Denitrification Rate:</strong> {summary['objective4_transport']['mean_denitrification_rate']:.3f} mmol/L/d</li>
                <li><strong>High Vulnerability Stations:</strong> {summary['objective4_transport']['high_vulnerability_stations']}</li>
            </ul>
            <p>Parameters derived from inverse modeling were used to simulate 1D transport (below).</p>
        </div>
        <div class="box">
            <h3>Vulnerability Map</h3>
            <img src="objective4_transport/transport_model_framework.png" alt="Vulnerability Map">
        </div>
    </div>
    
    <h3>Simulation Results</h3>
    <div class="container">
        <div class="box" style="flex: 2;">
            <p><strong>Scenario:</strong> 1D Advection-Dispersion with First-Order Decay (Ogata-Banks solution) using derived parameters.</p>
            <img src="objective4_transport/transport_simulation_analytical.png" alt="Transport Simulation">
        </div>
    </div>

    <footer>
        <p>Generated by Hydrosheaf Analysis Pipeline | {datetime.datetime.now().strftime("%Y")}</p>
    </footer>
    </body>
    </html>
    """

    with open(output_file, "w") as f:
        f.write(html)

    print(f"Report generated: {output_file}")


if __name__ == "__main__":
    generate_html_report()
