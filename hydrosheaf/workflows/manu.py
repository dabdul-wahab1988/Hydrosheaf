import pandas as pd
import sys
from pathlib import Path
from typing import List, Dict

# Attempt absolute import if package installed, else relative
try:
    from hydrosheaf import Config
    from hydrosheaf.api import fit_network_pipeline
    from hydrosheaf.data.units import mgL_to_mmolL, mmolL_to_mgL
    from hydrosheaf.graph.types import Edge
    from hydrosheaf.outputs.export import export_edge_results_json, export_edge_results_csv
    from hydrosheaf.log import setup_logging
    from hydrosheaf.inference.network_fit import infer_edges
except ImportError:
    # Fallback for running as script from root
    sys.path.append(str(Path(__file__).parent.parent.parent))
    from hydrosheaf import Config, fit_network_pipeline
    from hydrosheaf.data.units import mgL_to_mmolL, mmolL_to_mgL
    from hydrosheaf.graph.types import Edge
    from hydrosheaf.outputs.export import export_edge_results_json, export_edge_results_csv
    from hydrosheaf.log import setup_logging
    from hydrosheaf.inference.network_fit import infer_edges

class ManuWorkflow:
    """
    Advanced Workflow for Ghana Groundwater Analysis (Manu et al. Dataset).
    Features:
    - High-Resolution Topology (Probabilistic Geography-based)
    - Expanded Mineralogy (Silicates + Fluoride)
    - Water Quality Compliance Checks (WHO)
    - Detailed Geochemical Reporting
    """
    
    def __init__(self):
        self.config = Config(
            ion_order=['Ca', 'Mg', 'Na', 'K', 'HCO3', 'Cl', 'SO4', 'NO3', 'F', 'Fe', 'PO4'],
            weights=[1.0]*11,
            # Local Meteoric Water Line (Ghana - Gibrilla et al. 2022)
            lmwl_a=7.22, 
            lmwl_b=8.66,
            isotope_enabled=True,
            latent_endmembers_enabled=True,
            phreeqc_enabled=False, 
            transport_models_enabled=['evap', 'mix'],
            active_minerals=[
                "calcite", "dolomite", "gypsum", "albite", "halite", 
                "fluorite", "k_feldspar", "biotite", "CO2(g)", 
                "NO3src", "CaNa_exch", "pyrite_oxidation_aerobic"
            ],
            exchange_enabled=True,
            uncertainty_method="bayesian",
            bayesian_n_samples=500,
            # Edge Inference Settings
            edge_radius_km=15.0,
            edge_max_neighbors=2,
            edge_p_min=0.3,  # Allow some uncertainty in flow direction
            edge_head_key="head_meas",
            edge_elevation_key="elevation"
        )

    def load_data(self, file_path: str) -> List[Dict]:
        path = Path(file_path)
        if not path.exists():
            print(f"Error: {file_path} not found.")
            return []
            
        try:
            df = pd.read_excel(path)
        except Exception as e:
            print(f"Error reading Excel: {e}")
            return []
            
        samples = []
        print(f"Columns found: {list(df.columns)}")
        
        for idx, row in df.iterrows():
            # Robust ID extraction
            raw_station = row.get('Station')
            raw_id = row.get('Sample ID')
            
            # Combine Station and ID for uniqueness
            station_str = str(raw_station) if pd.notna(raw_station) else "UNK"
            id_str = str(raw_id) if pd.notna(raw_id) else f"S{idx}"
            site_id = f"{station_str}_{id_str}"

            s = {
                'site_id': site_id,
                'sample_id': id_str,
                # Map to standard Hydrosheaf coordinate keys (lon/lat)
                'lon': row.get('X coordinate', 0.0),
                'lat': row.get('Y coordinate', 0.0),
                'elevation': row.get('Elevation', 0.0),
                'head_meas': row.get('Elevation', 0.0), # Use Elevation as Head proxy for flow direction
                'pH': row.get('pH', 7.0),
                'temp_C': row.get('Temp', 25.0),
                '18O': row.get('d18O'), 
                '2H': row.get('d2H'),
                'type': 'sample'
            }
            
            # Map ions
            for ion in self.config.ion_order:
                val = row.get(ion)
                if pd.isna(val): 
                    s[ion] = None
                else:
                    try:
                        s[ion] = mgL_to_mmolL(float(val), ion)
                    except:
                        s[ion] = None
            
            samples.append(s)
        
        return samples

    def build_topology(self, samples: List[Dict]) -> List[Edge]:
        """
        Constructs flow network based on Probabilistic Head Inference.
        """
        print(f"Building geography-based topology for {len(samples)} samples...")
        
        # Try probabilistic first
        edges = infer_edges(
            samples,
            method="probabilistic",
            config=self.config,
            head_key="head_meas",
            elevation_key="elevation",
            allow_uphill=True
        )
        
        # Fallback to simple distance-based if probabilistic fails (e.g. flat area)
        if not edges:
            print("  Note: Probabilistic inference found no edges. Using distance-based fallback.")
            edges = infer_edges(
                samples,
                method="simple",
                max_neighbors=2,
                head_key="head_meas",
                elevation_key="elevation",
                allow_uphill=True
            )
        
        return edges

    def run(self, file_path: str, output_dir: str = "manu_results"):
        out_path = Path(output_dir)
        out_path.mkdir(exist_ok=True, parents=True)
        
        # Initialize Logging to File
        setup_logging(verbose=False, log_file=str(out_path / "analysis.log"))
        
        print("="*60)
        print("   GHANA GROUNDWATER ADVANCED ANALYSIS")
        print("="*60)
        
        samples = self.load_data(file_path)
        if not samples: return

        print(f"Loaded {len(samples)} samples.")
        
        edges = self.build_topology(samples)
        print(f"Identified {len(edges)} likely flow paths based on topography.")
        
        if not edges:
            print("No flow paths found.")
            return

        print("\nRunning Geochemical Inverse Modeling...")
        try:
            results, extras = fit_network_pipeline(samples, edges, self.config)
            
            # Export Detailed Results (Crucial for diagnostics/analysis)
            print(f"Saving detailed results to {output_dir}/...")
            export_edge_results_json(results, str(out_path / "results.json"))
            export_edge_results_csv(results, str(out_path / "results.csv"))
            
            self.generate_report(results, samples)
            
            # Generate Plots
            self.plot_results(results, samples, extras, output_dir)
            
        except Exception as e:
            print(f"Analysis Failed: {e}")
            import traceback
            traceback.print_exc()

    def plot_results(self, results, samples, extras, output_dir: str):
        """
        Generates comprehensive visualization suite.
        """
        from hydrosheaf.outputs.plots import plot_gibbs, plot_ilr, plot_edge_anomalies
        from hydrosheaf.outputs.plots_3d import plot_network_3d
        from hydrosheaf.outputs.utils import PlotConfig
        
        out_path = Path(output_dir)
        out_path.mkdir(exist_ok=True, parents=True)
        
        print(f"\nGenerating plots in {output_dir}/...")
        
        config = PlotConfig(font_scale=1.2)
        
        # 1. Gibbs Plot (Mechanism Discovery)
        plot_gibbs(samples, path=str(out_path / "gibbs_plot.png"), config=config)
        
        # 2. ILR Plot (Facies Classification)
        plot_ilr(samples, path=str(out_path / "facies_ilr.png"), config=config)
        
        # 3. Model Accuracy (Anomalies)
        plot_edge_anomalies(results, path=str(out_path / "model_fit_anomalies.png"), config=config)
        
        # 4. 3D Network Visualization
        try:
            plot_network_3d(
                samples, 
                edges=extras.get('edges', []), 
                output_path=str(out_path / "network_3d.png"),
                z_exaggeration=15.0,
                show=False
            )
        except Exception as e:
            print(f"  * 3D Plotting skipped: {e}")

        # 5. Bayesian Uncertainty (Top 3 active edges)
        from hydrosheaf.outputs.science_plots import plot_posterior_ridges
        # Sort results by Nitrate source loading if present
        try:
            interesting_results = sorted(results, key=lambda x: x.objective_score, reverse=True)[:3]
            for i, res in enumerate(interesting_results):
                if res.uncertainty and hasattr(res.uncertainty, "extents_samples"):
                    plot_posterior_ridges(
                        res, 
                        path=str(out_path / f"uncertainty_edge_{res.edge_id.replace('->','_')}.png"),
                        config=config
                    )
        except Exception as e:
            print(f"  * Uncertainty plotting skipped: {e}")

        print(f"Successfully generated diagnostic plots in {output_dir}/")

    def generate_report(self, results, samples):
        print("\n" + "="*60)
        print("   INSIGHTS REPORT")
        print("="*60)
        
        # 1. Water Quality Check (WHO Standards)
        print("\n[1] WATER QUALITY ALERTS (WHO Guidelines)")
        alerts = 0
        for s in samples:
            issues = []
            # NO3
            no3_mg = mmolL_to_mgL(s.get('NO3', 0), 'NO3')
            if no3_mg > 50: issues.append(f"NO3 {no3_mg:.1f} mg/L (>50)")
            # F
            f_mg = mmolL_to_mgL(s.get('F', 0), 'F')
            if f_mg > 1.5: issues.append(f"F {f_mg:.1f} mg/L (>1.5)")
            # Cl (Taste)
            cl_mg = mmolL_to_mgL(s.get('Cl', 0), 'Cl')
            if cl_mg > 250: issues.append(f"Cl {cl_mg:.1f} mg/L (>250)")
            
            if issues:
                print(f"  * {s['site_id']}: {', '.join(issues)}")
                alerts += 1
        
        if alerts == 0: print("  All samples within standard limits for NO3, F, Cl.")

        # 2. Geochemical Processes
        print("\n[2] DOMINANT GEOCHEMICAL PROCESSES")
        process_map = {}
        total_paths = len(results)
        
        for res in results:
            if res.transport_model == 'evap':
                process_map['Evaporation'] = process_map.get('Evaporation', 0) + 1
            
            if res.z_labels:
                for lbl, ext in zip(res.z_labels, res.z_extents):
                    if abs(ext) > 0.05:
                        process_map[lbl] = process_map.get(lbl, 0) + 1

        if process_map:
            sorted_p = sorted(process_map.items(), key=lambda x: x[1], reverse=True)
            for p, count in sorted_p:
                pct = (count / total_paths) * 100
                print(f"  * {p}: {pct:.1f}% of flow paths")
        else:
            print("  No dominant reactions detected.")

        # 3. Source Identification
        virtual = [s for s in samples if s.get('type') == 'virtual']
        if virtual:
            print("\n[3] INFERRED SOURCE CANDIDATES (Potential Pollution)")
            for vn in virtual:
                no3_mg = mmolL_to_mgL(vn.get('NO3', 0), 'NO3')
                cl_mg = mmolL_to_mgL(vn.get('Cl', 0), 'Cl')
                so4_mg = mmolL_to_mgL(vn.get('SO4', 0), 'SO4')
                
                print(f"  * {vn['site_id']} (End-member):")
                print(f"      Nitrate: {no3_mg:.1f} mg/L")
                print(f"      Chloride: {cl_mg:.1f} mg/L")
                if so4_mg > 10: print(f"      Sulfate: {so4_mg:.1f} mg/L")
        else:
            print("\n[3] No hidden contamination sources detected.")

        # 4. Summary Recommendation
        print("\n[4] SUMMARY & RECOMMENDATIONS")
        if process_map.get('Evaporation', 0) > total_paths * 0.3:
            print("  - High evaporation signature: Suggests shallow circulation or arid conditions.")
        if process_map.get('fluorite', 0) > total_paths * 0.1:
            print("  - Fluorite dissolution detected: Monitor Fluoride levels closely.")
        if process_map.get('NO3src', 0) > 0:
            print("  - Nitrate inputs detected: Possible agricultural or septic contamination.")
        if process_map.get('CaNa_exch', 0) > total_paths * 0.2:
            print("  - Cation Exchange active: Indicates freshening or hardening (check flow direction).")

def run_manu(file_path: str = "manu.xlsx"):
    wf = ManuWorkflow()
    wf.run(file_path)

if __name__ == "__main__":
    if len(sys.argv) > 1:
        run_manu(sys.argv[1])
    else:
        run_manu()
