"""
Automated Intelligent Workflows.
Encapsulates the "Universal Pipeline" for end-users.
"""

import pandas as pd
import yaml
import os
import math
from pathlib import Path
from typing import Optional, List, Dict, Any

from ..config import Config
from ..api import fit_network_pipeline
from ..data.units import mgL_to_mmolL, mmolL_to_mgL
from ..graph.types import Edge
from ..log import get_logger

logger = get_logger(__name__)

class HydrosheafAnalyzer:
    def __init__(self, config_path: Optional[str] = None):
        self.config = Config()
        
        # 1. Set Universal Defaults (The "Broad Spectrum" Logic)
        # We enable a representative mineral from every major geological class
        # so the model can handle any terrain (Sedimentary, Crystalline, Coastal).
        self.config.active_minerals = [
            # Carbonates (Limestone/Dolomite terrains)
            "calcite", "dolomite", 
            # Evaporites (Saline/Gypsiferous terrains)
            "gypsum", "halite", "fluorite",
            # Silicates (Granite/Gneiss/Volcanic terrains)
            "albite",       # Na-Feldspar
            "anorthite",    # Ca-Feldspar
            "k_feldspar",   # K-Feldspar
            "biotite",      # Mica (Source of F, Fe, Mg)
            # Redox & Pollution
            "pyrite_oxidation_aerobic", # Acid mine drainage / natural oxidation
            "NO3src",       # Anthropogenic Nitrate input
            # Critical Processes
            "CO2(g)",       # Open system degassing/recharge
            "CaNa_exch"     # Cation Exchange (Clay interactions)
        ]
        
        self.config.transport_models_enabled = ['evap', 'mix']
        self.config.latent_endmembers_enabled = True
        self.config.exchange_enabled = True
        self.config.gibbs_enabled = True

        # 2. Override with User Config if provided
        if config_path and os.path.exists(config_path):
            self.load_config(config_path)
        
    def load_config(self, path: str):
        with open(path, 'r') as f:
            cfg_dict = yaml.safe_load(f)
        valid_keys = Config.__init__.__code__.co_varnames
        filtered = {k: v for k, v in cfg_dict.items() if k in valid_keys}
        for k, v in filtered.items():
            setattr(self.config, k, v)

    def load_data(self, file_path: str) -> List[Dict[str, Any]]:
        path = Path(file_path)
        if not path.exists():
            raise FileNotFoundError(f"{file_path} not found.")
            
        if path.suffix in ['.xls', '.xlsx']:
            try:
                df = pd.read_excel(path)
            except ImportError:
                raise ImportError("Please install openpyxl to read Excel files: pip install openpyxl")
        else:
            df = pd.read_csv(path)
            
        samples = []
        # Robust Ion Mapping
        ion_map = {
            'Ca': ['Ca', 'Ca_mgL', 'Calcium'],
            'Mg': ['Mg', 'Mg_mgL', 'Magnesium'],
            'Na': ['Na', 'Na_mgL', 'Sodium'],
            'K': ['K', 'K_mgL', 'Potassium'],
            'HCO3': ['HCO3', 'HCO3_mgL', 'Bicarbonate', 'Alkalinity'],
            'Cl': ['Cl', 'Cl_mgL', 'Chloride'],
            'SO4': ['SO4', 'SO4_mgL', 'Sulfate'],
            'NO3': ['NO3', 'NO3_mgL', 'Nitrate'],
            'F': ['F', 'F_mgL', 'Fluoride'],
            'Fe': ['Fe', 'Fe_mgL', 'Iron'],
            'PO4': ['PO4', 'PO4_mgL', 'Phosphate']
        }
        
        col_map = {}
        for ion, candidates in ion_map.items():
            for c in candidates:
                if c in df.columns:
                    col_map[ion] = c
                    break
                    
        # Metadata map
        meta_map = {}
        for key in ['Station', 'Sample ID', 'site_id', 'id']:
            if key in df.columns:
                meta_map['site_id'] = key
                break
        
        for key in ['X coordinate', 'X', 'Lon', 'Longitude']:
            if key in df.columns: meta_map['x'] = key
        for key in ['Y coordinate', 'Y', 'Lat', 'Latitude']:
            if key in df.columns: meta_map['y'] = key
        for key in ['Elevation', 'Z', 'Elev', 'Depth']:
            if key in df.columns: meta_map['z'] = key
            
        # Isotopes
        iso_map = {}
        for key in ['d18O', '18O', 'd18O_permil']:
            if key in df.columns: iso_map['18O'] = key
        for key in ['d2H', '2H', 'd2H_permil', 'D']:
            if key in df.columns: iso_map['2H'] = key
            
        for _, row in df.iterrows():
            sid = str(row.get(meta_map.get('site_id', df.columns[0])))
            s = {'site_id': sid, 'type': 'sample'}
            
            if 'x' in meta_map: s['x'] = row.get(meta_map['x'], 0.0)
            if 'y' in meta_map: s['y'] = row.get(meta_map['y'], 0.0)
            if 'z' in meta_map: s['z'] = row.get(meta_map['z'], 0.0)
            
            for ion in self.config.ion_order:
                src_col = col_map.get(ion)
                val = 0.0
                if src_col:
                    raw = row.get(src_col, 0.0)
                    if not pd.isna(raw):
                        try:
                            val = mgL_to_mmolL(float(raw), ion)
                        except:
                            val = 0.0
                s[ion] = val
            
            if '18O' in iso_map: s['18O'] = row.get(iso_map['18O'])
            if '2H' in iso_map: s['2H'] = row.get(iso_map['2H'])
            
            samples.append(s)
            
        return samples

    def build_topology(self, samples: List[Dict]):
        has_xyz = all('x' in s for s in samples[:5])
        edges = []
        if has_xyz and len(samples) > 2:
            print("[AUTO] Building 3D Hydraulic Topology...")
            for i, u in enumerate(samples):
                neighbors = []
                for j, v in enumerate(samples):
                    if i == j: continue
                    dx = u.get('x',0)-v.get('x',0)
                    dy = u.get('y',0)-v.get('y',0)
                    dz = u.get('z',0)-v.get('z',0)
                    dist = math.sqrt(dx*dx + dy*dy)
                    if dz > -5.0:
                        neighbors.append((v['site_id'], dist))
                neighbors.sort(key=lambda x: x[1])
                for vid, dist in neighbors[:2]:
                    edges.append(Edge(u=u['site_id'], v=vid, edge_id=f"{u['site_id']}->{vid}"))
        else:
            print("[AUTO] No spatial data. Building sequential chain for demo.")
            for i, u in enumerate(samples):
                # Fallback: connect to next sample
                if i < len(samples) - 1:
                    edges.append(Edge(u=u['site_id'], v=samples[i+1]['site_id'], edge_id="Chain"))
        return edges

    def run(self, file_path: str):
        print(f"Analyzing {file_path}...")
        samples = self.load_data(file_path)
        print(f"Loaded {len(samples)} samples.")
        
        edges = self.build_topology(samples)
        print(f"Identified {len(edges)} flow paths.")
        
        try:
            results, extras = fit_network_pipeline(samples, edges, self.config)
            self.print_report(results, samples)
            return results
        except Exception as e:
            print(f"Analysis Error: {e}")
            import traceback
            traceback.print_exc()
            return None

    def print_report(self, results, samples):
        print("\n" + "="*60)
        print("   HYDROSHEAF INTELLIGENT REPORT")
        print("="*60)
        
        virtual = [s for s in samples if s.get('type') == 'virtual']
        if virtual:
            print(f"\n[SOURCE DISCOVERY] Found {len(virtual)} inferred source candidates:")
            for vn in virtual:
                no3_mg = mmolL_to_mgL(vn.get('NO3', 0), 'NO3')
                cl_mg = mmolL_to_mgL(vn.get('Cl', 0), 'Cl')
                print(f"  * {vn['site_id']}: NO3={no3_mg:.1f} mg/L, Cl={cl_mg:.1f} mg/L")
        
        if not results:
            print("No results generated.")
            return
            
        n_evap = sum(1 for r in results if r.transport_model == 'evap')
        print(f"\n[SYSTEM] Dominant Regime: {'Evaporation' if n_evap > len(results)/2 else 'Mixing'}")
        
        print("\n[HOTSPOTS] Active Contamination Loading:")
        hotspots = []
        for res in results:
            if "NO3src" in res.z_labels:
                val = res.z_extents[res.z_labels.index("NO3src")] * 62.0
                if val > 10:
                    hotspots.append((res.v, val))
        
        if hotspots:
            hotspots.sort(key=lambda x: x[1], reverse=True)
            # Deduplicate
            seen = set()
            for s, v in hotspots:
                if s not in seen:
                    print(f"  * {s}: +{v:.1f} mg/L Nitrate Added")
                    seen.add(s)
                if len(seen) >= 5: break
        else:
            print("  None detected.")

def analyze_dataset(data_path: str, config_path: Optional[str] = None):
    """
    Main entry point for one-line analysis.
    """
    analyzer = HydrosheafAnalyzer(config_path)
    return analyzer.run(data_path)
