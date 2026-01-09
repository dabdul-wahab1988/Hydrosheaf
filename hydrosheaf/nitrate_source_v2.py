import math
import yaml
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from .coda_sbp import ilr_from_sbp, robust_zscore
from .config import Config
from .data.units import mgL_to_mmolL
from .models import nitrate_isotopes

# Default config path
DEFAULT_CONFIG_PATH = Path(__file__).parent / "config" / "nitrate_source_v2.yaml"


@dataclass
class NitrateSourceResult:
    """Result of nitrate source discrimination."""
    p_manure: Optional[float]
    p_fertilizer: Optional[float]
    logit_score: Optional[float]
    top_evidence: List[str]
    gating_flags: List[str]
    ilr_valid: bool
    reason_code: Optional[str] = None


@dataclass
class NitrateStats:
    """Robust statistics for z-scoring."""
    ln_no3_cl_median: float = 0.0
    ln_no3_cl_mad: float = 1.0
    ln_no3_k_median: float = 0.0
    ln_no3_k_mad: float = 1.0
    ln_po4_cl_median: float = 0.0
    ln_po4_cl_mad: float = 1.0
    ln_fe_median: float = 0.0
    ln_fe_mad: float = 1.0
    d_excess_p25: float = 10.0
    po4_p90: float = 0.1
    # Denitrification extent (reaction) stats
    denit_median: float = 0.0
    denit_mad: float = 1.0
    # Alkalinity ratio stats
    alk_ratio_median: float = 0.0
    alk_ratio_mad: float = 1.0


def load_nitrate_config(path: Path = DEFAULT_CONFIG_PATH) -> dict:
    if not path.exists():
        return {}
    with open(path, "r") as f:
        return yaml.safe_load(f)


def _safe_log_ratio(num: float, den: float, eps: float = 1e-12) -> float:
    return math.log(max(num, eps) / max(den, eps))


def _sigmoid(x: float) -> float:
    return 1.0 / (1.0 + math.exp(-x))


def fit_robust_stats(
    samples_df: pd.DataFrame, 
    edge_results: Optional[List[Any]] = None
) -> NitrateStats:
    """Compute global robust statistics (Median/MAD) from data."""
    
    stats = NitrateStats()
    eps = 1e-12
    
    # Helper to get series
    def get_series(col):
        return samples_df[col] if col in samples_df.columns else pd.Series(dtype=float)

    # 1. NO3/Cl
    no3 = get_series("NO3")
    cl = get_series("Cl")
    if not no3.empty and not cl.empty:
        # Filter out NaN or negative
        valid = (no3 > 0) & (cl > 0)
        if valid.any():
            ratio = np.log((no3[valid] + eps) / (cl[valid] + eps))
            stats.ln_no3_cl_median = ratio.median()
            stats.ln_no3_cl_mad = (ratio - stats.ln_no3_cl_median).abs().median()

    # 2. NO3/K
    k_ion = get_series("K")
    if not no3.empty and not k_ion.empty:
        valid = (no3 > 0) & (k_ion > 0)
        if valid.any():
            ratio = np.log((no3[valid] + eps) / (k_ion[valid] + eps))
            stats.ln_no3_k_median = ratio.median()
            stats.ln_no3_k_mad = (ratio - stats.ln_no3_k_median).abs().median()

    # 3. PO4/Cl and P90
    po4 = get_series("PO4")
    if not po4.empty:
        stats.po4_p90 = po4.quantile(0.9)
        if not cl.empty:
            valid = (po4 > 0) & (cl > 0)
            if valid.any():
                ratio = np.log((po4[valid] + eps) / (cl[valid] + eps))
                stats.ln_po4_cl_median = ratio.median()
                stats.ln_po4_cl_mad = (ratio - stats.ln_po4_cl_median).abs().median()

    # 4. Fe
    fe = get_series("Fe")
    if not fe.empty:
        valid = fe > 0
        if valid.any():
            val = np.log(fe[valid] + eps)
            stats.ln_fe_median = val.median()
            stats.ln_fe_mad = (val - stats.ln_fe_median).abs().median()

    # 5. d-excess
    # d_excess = d2H - 8 * d18O. Check if computed in df
    if "d_excess" in samples_df.columns:
        stats.d_excess_p25 = samples_df["d_excess"].quantile(0.25)

    # 6. Edge stats (Denitrification & Alkalinity)
    if edge_results:
        denit_vals = []
        alk_ratios = []
        
        for e in edge_results:
            # Denitrification extent
            if hasattr(e, "z_labels") and hasattr(e, "z_extents"):
                try:
                    idx = e.z_labels.index("denit")
                    z_val = e.z_extents[idx]
                    denit_vals.append(z_val)
                except ValueError:
                    pass

        if denit_vals:
            arr = np.array(denit_vals)
            stats.denit_median = np.median(arr)
            stats.denit_mad = np.median(np.abs(arr - stats.denit_median))

    return stats


def compute_evidence(
    sample: Dict[str, float],
    stats: NitrateStats,
    weights: Dict[str, float],
    gate_factor: float = 1.0,
    eps: float = 1e-12,
) -> Tuple[float, List[str]]:
    """Compute node-level evidence logits."""
    logit_add = 0.0
    evidence = []

    # w1: NO3/Cl (High -> Fertilizer -> Negative logit)
    if "NO3" in sample and "Cl" in sample:
        val = _safe_log_ratio(sample["NO3"], sample["Cl"], eps)
        z = robust_zscore(val, stats.ln_no3_cl_median, stats.ln_no3_cl_mad)
        # Sign: Negative. High ratio -> High val -> High Z -> Negative contribution
        term = -1.0 * z * weights.get("w1_no3_cl", 1.2) * gate_factor
        logit_add += term
        if term < -0.5: evidence.append("NO3/Cl_high_fert")
        if term > 0.5: evidence.append("NO3/Cl_low_manure")

    # w2: NO3/K (High -> Fertilizer -> Negative logit)
    if "NO3" in sample and "K" in sample:
        val = _safe_log_ratio(sample["NO3"], sample["K"], eps)
        z = robust_zscore(val, stats.ln_no3_k_median, stats.ln_no3_k_mad)
        term = -1.0 * z * weights.get("w2_no3_k", 0.4) * gate_factor
        logit_add += term

    # w3: PO4 High (High -> Fertilizer -> Negative logit)
    if "PO4" in sample and "Cl" in sample:
        po4_high = sample["PO4"] > stats.po4_p90
        flag_val = 1.0 if po4_high else 0.0
        val = _safe_log_ratio(sample["PO4"], sample["Cl"], eps)
        z = robust_zscore(val, stats.ln_po4_cl_median, stats.ln_po4_cl_mad)
        
        direction = 2.0 * flag_val - 1.0
        term = -1.0 * direction * z * weights.get("w3_po4", 0.3) * gate_factor
        logit_add += term
        if po4_high and term < -0.2: evidence.append("PO4_high_fert")

    # w4: Fe (High -> Reducing -> Manure -> Positive logit)
    if "Fe" in sample:
        val = math.log(max(sample["Fe"], eps))
        z = robust_zscore(val, stats.ln_fe_median, stats.ln_fe_mad)
        term = 1.0 * z * weights.get("w4_fe", 0.6)
        logit_add += term
        if term > 0.5: evidence.append("Fe_high_manure")

    return logit_add, evidence


def compute_edge_evidence(
    edge_result: Any,
    u_vals: Dict[str, float],
    v_vals: Dict[str, float],
    stats: NitrateStats,
    weights: Dict[str, float],
    eps: float = 1e-12,
) -> Tuple[float, List[str]]:
    """Compute edge-level evidence (Denitrification, Alkalinity)."""
    logit_add = 0.0
    evidence = []

    # w5: Denitrification Extent
    denit_extent = 0.0
    if hasattr(edge_result, "z_labels") and hasattr(edge_result, "z_extents"):
        try:
            idx = edge_result.z_labels.index("denit")
            denit_extent = edge_result.z_extents[idx]
        except ValueError:
            pass
            
    # Logic: if denitrification active (removal > 0)
    if denit_extent > 1e-9: # Positive z = removal
        mad = stats.denit_mad if stats.denit_mad > 1e-9 else 1.0
        z = (denit_extent - stats.denit_median) / (1.4826 * mad)
        term = 1.0 * z * weights.get("w5_denitrif", 1.5)
        # Only apply if it votes manure (positive)
        if term > 0:
            logit_add += term
            if term > 0.5: evidence.append("denitrif_strong")
            
    # w6: Alkalinity Coupling
    dn = v_vals.get("NO3", 0.0) - u_vals.get("NO3", 0.0)
    dh = v_vals.get("HCO3", 0.0) - u_vals.get("HCO3", 0.0)
    
    # Check if nitrate removed (dn < 0)
    if dn < -eps:
        r_alk = dh / (-dn)
        r_alk = min(r_alk, 2.0)
        z = r_alk # Simplified if no stats
        term = 1.0 * z * weights.get("w6_alk_coupling", 0.8)
        logit_add += term
        if term > 0.5: evidence.append("alk_coupling_manure")
        
    return logit_add, evidence


def infer_node_posteriors(
    nodes_df: pd.DataFrame,
    edge_results: List[Any],
    config_overrides: Optional[dict] = None,
) -> Dict[str, NitrateSourceResult]:
    """Main inference function."""
    
    # 0. Load Config
    file_conf = load_nitrate_config()
    if config_overrides:
        file_conf.update(config_overrides)
        
    weights = file_conf.get("weights", {})
    prior_p = file_conf.get("prior_prob", 0.5)
    prior_logit = math.log(prior_p / (1.0 - prior_p))
    evap_gate = file_conf.get("evap_gate_factor", 0.5)
    
    # Isotope Config
    iso_enabled = file_conf.get("nitrate_isotope_mixing_enabled", True)
    n15_col = file_conf.get("nitrate_isotope_n15_col", "d15N")
    o18_col = file_conf.get("nitrate_isotope_o18_col", "d18O_NO3")
    
    iso_sources = []
    if iso_enabled:
        iso_sources = nitrate_isotopes.load_isotope_endmembers()

    # Threshold for Background
    min_mg_L = float(file_conf.get("nitrate_source_min_mg_L", 10.0))
    # Convert to internal units (likely mol/L if using mgL_to_mmolL)
    min_conc = mgL_to_mmolL(min_mg_L, "NO3")

    # 1. Compute Stats
    stats = fit_robust_stats(nodes_df, edge_results)
    
    # Apply Config Overrides
    if file_conf.get("nitrate_source_d_excess_p25") is not None:
        stats.d_excess_p25 = float(file_conf["nitrate_source_d_excess_p25"])
    if file_conf.get("nitrate_source_po4_p90") is not None:
        stats.po4_p90 = float(file_conf["nitrate_source_po4_p90"])
    
    # 2. Node Processing
    results = {}
    
    # Pre-index edges by target v
    incoming_edges = {}
    for e in edge_results:
        if hasattr(e, "v"):
            incoming_edges.setdefault(e.v, []).append(e)
            
    for idx, row in nodes_df.iterrows():
        sample = row.to_dict()
        node_id = str(sample.get("site_id", idx))
        
        # Check CoDA validity
        ilr, ilr_valid = ilr_from_sbp(sample)
        
        # Check Nitrate Threshold
        no3_val = sample.get("NO3", 0.0)
        if no3_val < min_conc:
            # Skip inference for background
            results[node_id] = NitrateSourceResult(
                p_manure=None,
                p_fertilizer=None,
                logit_score=None,
                top_evidence=[],
                gating_flags=["below_detection_threshold"],
                ilr_valid=ilr_valid,
                reason_code=f"Low Nitrate (Background < {min_mg_L} mg/L)"
            )
            continue
            
        # --- Dual Isotope Logic (v0.3.0) ---
        used_isotope_model = False
        if iso_enabled and iso_sources:
            try:
                d15_val = sample.get(n15_col)
                d18_val = sample.get(o18_col)
                
                # Check valid float
                if (d15_val is not None and not math.isnan(d15_val) and 
                    d18_val is not None and not math.isnan(d18_val)):
                    
                    iso_s = nitrate_isotopes.IsotopeSample(float(d15_val), float(d18_val))
                    probs = nitrate_isotopes.compute_isotope_prob(iso_s, iso_sources)
                    
                    p_man = probs.get("Manure", 0.0)
                    # p_fert = probs.get("Fertilizer", 0.0) # Could be multiple non-manure sources
                    
                    results[node_id] = NitrateSourceResult(
                        p_manure=p_man,
                        p_fertilizer=1.0 - p_man,
                        logit_score=None, # Not applicable for mixing model
                        top_evidence=[f"d15N={d15_val:.1f}", f"d18O={d18_val:.1f}"],
                        gating_flags=["dual_isotope_priority"],
                        ilr_valid=ilr_valid,
                        reason_code="Dual Isotope Mixing"
                    )
                    used_isotope_model = True
            except Exception:
                # Fallback on error
                pass
                
        if used_isotope_model:
            continue
        
        # Check Evap Gate
        # Flag if d_excess < P25
        d_excess = sample.get("d_excess")
        is_evap = False
        gates = []
        
        if d_excess is not None and d_excess < stats.d_excess_p25:
            is_evap = True
            gates.append("low_d_excess")
            
        # Also check incoming edges for evap model
        edges = incoming_edges.get(node_id, [])
        for e in edges:
            if hasattr(e, "transport_model") and e.transport_model == "evap":
                is_evap = True
                if "transport_evap" not in gates:
                    gates.append("transport_evap")
                    
        gate_factor = evap_gate if is_evap else 1.0
        
        # Base Logit
        logit = prior_logit
        all_evidence = []
        
        # Node Evidence
        node_logit, node_ev = compute_evidence(sample, stats, weights, gate_factor)
        logit += node_logit
        all_evidence.extend(node_ev)
        
        # Edge Evidence (incoming)
        if edges:
            edge_logits = []
            for e in edges:
                u_id = getattr(e, "u", None)
                if u_id in nodes_df.index: # Assuming index is site_id
                     u_row = nodes_df.loc[u_id].to_dict()
                     elogit, eev = compute_edge_evidence(e, u_row, sample, stats, weights)
                     edge_logits.append(elogit)
                     all_evidence.extend(eev)
                elif hasattr(e, "u") and e.u in nodes_df["site_id"].values:
                     # fallback search
                     u_rows = nodes_df[nodes_df["site_id"] == e.u]
                     if not u_rows.empty:
                        u_row = u_rows.iloc[0].to_dict()
                        elogit, eev = compute_edge_evidence(e, u_row, sample, stats, weights)
                        edge_logits.append(elogit)
                        all_evidence.extend(eev)
            
            if edge_logits:
                logit += sum(edge_logits) / len(edge_logits)

        # Final Sigmoid
        p_manure = _sigmoid(logit)
        
        results[node_id] = NitrateSourceResult(
            p_manure=p_manure,
            p_fertilizer=1.0 - p_manure,
            logit_score=logit,
            top_evidence=list(set(all_evidence)),
            gating_flags=gates,
            ilr_valid=ilr_valid,
            reason_code="Hydrochemical Ratios (No Isotopes)"
        )
        
    return results
