
import math
import copy
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from .coda_sbp import ilr_from_sbp, robust_zscore
from .config import Config
from .data.units import mgL_to_mmolL
from .models import nitrate_isotopes
import yaml

if TYPE_CHECKING:
    from .models.nitrate_isotopes_mcmc import MCMCMixingResult

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
    # MCMC-specific fields
    mcmc_result: Optional["MCMCMixingResult"] = None
    source_fractions: Optional[Dict[str, float]] = None
    ci_lower: Optional[Dict[str, float]] = None
    ci_upper: Optional[Dict[str, float]] = None
    diagnostics: Optional[Dict[str, float]] = None


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

@lru_cache(maxsize=8)
def _load_nitrate_config_cached(path_str: str) -> Dict[str, Any]:
    path = Path(path_str)
    if not path.exists():
        return {}
    with open(path, "r", encoding="utf-8") as handle:
        loaded = yaml.safe_load(handle)
    return loaded if isinstance(loaded, dict) else {}


def load_nitrate_config(path: Path = DEFAULT_CONFIG_PATH) -> dict:
    return copy.deepcopy(_load_nitrate_config_cached(str(path)))




def _safe_log_ratio(num: float, den: float, eps: float = 1e-12) -> float:
    return math.log(max(num, eps) / max(den, eps))


def _sigmoid(x: float) -> float:
    return 1.0 / (1.0 + math.exp(-x))


def fit_robust_stats(
    samples_df: pd.DataFrame, edge_results: Optional[List[Any]] = None
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
    """
    Compute node-level evidence using Naive Bayes assumption.
    
    Returns:
        log_likelihood_ratio: ln( P(Data|Manure) / P(Data|Fertilizer) )
        evidence_list: List of triggered evidence strings
    """
    llr = 0.0
    evidence = []

    # 1. NO3/Cl (Manure has Low Ratio, Fertilizer has High Ratio)
    # Using Z-scores to estimate likelihood strength
    if "NO3" in sample and "Cl" in sample:
        val = _safe_log_ratio(sample["NO3"], sample["Cl"], eps)
        z = robust_zscore(val, stats.ln_no3_cl_median, stats.ln_no3_cl_mad)
        # Low ratio (Negative Z) -> Manure -> Positive LLR
        # High ratio (Positive Z) -> Fertilizer -> Negative LLR
        # weight w1 corresponds to the strength of this evidence
        w1 = weights.get("w1_no3_cl", 1.2) * gate_factor
        term = -1.0 * z * w1
        llr += term
        if term > 0.5:
            evidence.append("NO3/Cl_low_manure")
        elif term < -0.5:
            evidence.append("NO3/Cl_high_fert")

    # 2. NO3/K (Fertilizer (Potash) has High K, so Low Ratio? No, Fertilizer is N+K usually.
    # Actually, manure has very high K. So Low NO3/K -> Manure.
    # Existing logic was: High NO3/K -> Fertilizer. 
    if "NO3" in sample and "K" in sample:
        val = _safe_log_ratio(sample["NO3"], sample["K"], eps)
        z = robust_zscore(val, stats.ln_no3_k_median, stats.ln_no3_k_mad)
        w2 = weights.get("w2_no3_k", 0.4) * gate_factor
        term = -1.0 * z * w2
        llr += term

    # 3. PO4 (Manure has High PO4)
    if "PO4" in sample:
        val = _safe_log_ratio(sample["PO4"], sample.get("Cl", 1.0), eps) # Normalized to Cl
        z = robust_zscore(val, stats.ln_po4_cl_median, stats.ln_po4_cl_mad)
        # High PO4 -> Positive Z -> Manure
        w3 = weights.get("w3_po4", 0.3) * gate_factor
        term = 1.0 * z * w3
        llr += term
        if term > 0.3:
            evidence.append("PO4_high_manure")

    # 4. Fe (Reducing Conditions)
    # Reducing conditions often associated with manure (high DOC), but not strictly source.
    # Treat as weak evidence.
    if "Fe" in sample:
        val = math.log(max(sample["Fe"], eps))
        z = robust_zscore(val, stats.ln_fe_median, stats.ln_fe_mad)
        w4 = weights.get("w4_fe", 0.6)
        term = 1.0 * z * w4
        llr += term
        if term > 0.5:
            evidence.append("Fe_high_manure")

    return llr, evidence


def compute_edge_evidence(
    edge_result: Any,
    u_vals: Dict[str, float],
    v_vals: Dict[str, float],
    stats: NitrateStats,
    weights: Dict[str, float],
    eps: float = 1e-12,
) -> Tuple[float, List[str]]:
    """
    Compute edge-level evidence.
    
    REMOVED: Circular 'denitrification extent' logic.
    Instead, we look for 'alkalinity coupling' which suggests organic carbon oxidation (manure).
    """
    llr = 0.0
    evidence = []

    # w6: Alkalinity Coupling
    # Reaction: 5CH2O + 4NO3- -> 2N2 + 4HCO3- + H2CO3 + H2O
    # Stoichiometry: Consumes NO3, Produces HCO3. Ratio ~ 1:1 to 1:2.
    
    dn = v_vals.get("NO3", 0.0) - u_vals.get("NO3", 0.0)
    dh = v_vals.get("HCO3", 0.0) - u_vals.get("HCO3", 0.0)

    # If NO3 decreases and HCO3 increases
    if dn < -eps and dh > eps:
        ratio = dh / (-dn)
        # Ideal heterotrophic denitrification ratio is near 1.0 - 1.25
        # Sulfide-driven is different.
        # If ratio is close to expected, it supports organic C source (Manure)
        if 0.5 < ratio < 2.5:
            w6 = weights.get("w6_alk_coupling", 0.8)
            llr += w6
            evidence.append("alk_coupling_manure")

    # w5: Denitrification Boost (Physical Evidence of Nitrate Removal)
    # Strong denitrification extent in reactive model supports manure/organic C presence.
    if hasattr(edge_result, "z_labels") and hasattr(edge_result, "z_extents"):
        try:
            idx = edge_result.z_labels.index("denit")
            z_val = edge_result.z_extents[idx]
            z_score = robust_zscore(z_val, stats.denit_median, stats.denit_mad)
            
            if z_score > 1.0:
                w5 = weights.get("w5_denitrif", 1.0)
                llr += w5
                evidence.append("denitrif_strong")
        except (ValueError, IndexError):
            pass

    return llr, evidence


def infer_node_posteriors(
    nodes_df: pd.DataFrame,
    edge_results: List[Any],
    config_overrides: Optional[dict] = None,
    config: Optional[Config] = None,
) -> Dict[str, NitrateSourceResult]:
    """Main inference function.

    Parameters
    ----------
    nodes_df : pd.DataFrame
        DataFrame with sample data indexed by site_id
    edge_results : List[Any]
        List of edge fitting results
    config_overrides : dict, optional
        Override config values from YAML
    config : Config, optional
        Hydrosheaf Config object for MCMC settings

    Returns
    -------
    Dict[str, NitrateSourceResult]
        Results keyed by node/site ID
    """

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
    water_o18_col = file_conf.get("nitrate_isotope_water_o18_col", "d18O")
    process_constraints_enabled = file_conf.get(
        "nitrate_isotope_process_constraints_enabled", True
    )
    process_constraints_conf = file_conf.get("isotope_process_constraints", {})
    if not isinstance(process_constraints_conf, dict):
        process_constraints_conf = {}
    isotope_qc_enabled = bool(file_conf.get("nitrate_isotope_qc_enabled", True))
    isotope_qc_conf = file_conf.get("isotope_qc", {})
    if not isinstance(isotope_qc_conf, dict):
        isotope_qc_conf = {}
    process_mcmc_alpha_scale = float(
        process_constraints_conf.get("mcmc_alpha_scale", 8.0)
    )

    # MCMC Config
    mcmc_enabled = False
    mcmc_n_samples = 2000
    mcmc_n_chains = 4
    mcmc_target_accept = 0.9
    mcmc_warmup = 1000
    mcmc_hierarchical_enabled = bool(
        file_conf.get("isotope_mcmc_hierarchical_enabled", False)
    )

    if config is not None:
        mcmc_enabled = config.isotope_mcmc_enabled
        mcmc_n_samples = config.isotope_mcmc_n_samples

        mcmc_n_chains = config.isotope_mcmc_n_chains
        mcmc_target_accept = config.isotope_mcmc_target_accept
        mcmc_warmup = config.isotope_mcmc_warmup
        mcmc_hierarchical_enabled = config.isotope_mcmc_hierarchical_enabled
        water_o18_col = config.nitrate_isotope_water_o18_col
        process_constraints_enabled = (
            config.nitrate_isotope_process_constraints_enabled
        )
        isotope_qc_enabled = config.nitrate_isotope_qc_enabled

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
    incoming_edges: Dict[str, List[Any]] = {}
    for e in edge_results:
        if hasattr(e, "v"):
            incoming_edges.setdefault(str(e.v), []).append(e)

    node_rows: List[Tuple[str, Dict[str, Any]]] = []
    node_lookup: Dict[Any, Dict[str, Any]] = {}
    samples_records = nodes_df.to_dict(orient="records")
    node_index_values = list(nodes_df.index)
    for idx, sample in zip(node_index_values, samples_records):
        site_value = sample.get("site_id", idx)
        node_id = str(site_value)
        node_rows.append((node_id, sample))
        node_lookup[idx] = sample
        node_lookup[node_id] = sample
        if site_value is not None:
            node_lookup[site_value] = sample

    def _compute_process_priors(
        node_id: str,
        sample: Dict[str, Any],
        iso_sample: nitrate_isotopes.IsotopeSample,
    ) -> Tuple[Optional[Dict[str, float]], List[str]]:
        if not process_constraints_enabled:
            return None, []
        edges_local = incoming_edges.get(node_id, [])
        water_o18 = sample.get(water_o18_col)
        water_o18_float = (
            float(water_o18)
            if water_o18 is not None and not pd.isna(water_o18)
            else None
        )
        denit_vals: List[float] = []
        for edge in edges_local:
            if hasattr(edge, "z_labels") and hasattr(edge, "z_extents"):
                try:
                    denit_idx = edge.z_labels.index("denit")
                    denit_val = float(edge.z_extents[denit_idx])
                    if denit_val > 0.0:
                        denit_vals.append(denit_val)
                except ValueError:
                    continue
        denit_extent = (
            float(sum(denit_vals) / len(denit_vals)) if denit_vals else None
        )
        process_probs, process_flags, _ = nitrate_isotopes.compute_process_prior_probs(
            sample=iso_sample,
            sources=iso_sources,
            water_d18O=water_o18_float,
            denitrification_extent=denit_extent,
            process_config=process_constraints_conf,
        )
        return process_probs, process_flags

    hierarchical_results: Dict[str, Any] = {}
    hierarchical_flags: Dict[str, List[str]] = {}
    if mcmc_enabled and mcmc_hierarchical_enabled and iso_enabled and iso_sources:
        try:
            from .models.nitrate_isotopes_mcmc import run_mcmc_mixing_hierarchical

            hier_samples: List[nitrate_isotopes.IsotopeSample] = []
            hier_node_ids: List[str] = []
            pooled_process_priors: List[Dict[str, float]] = []
            for node_id, sample in node_rows:
                no3_val = sample.get("NO3", 0.0)
                if no3_val < min_conc:
                    continue
                d15_val = sample.get(n15_col)
                d18_val = sample.get(o18_col)
                if d15_val is None or d18_val is None:
                    continue
                try:
                    d15_float = float(d15_val)
                    d18_float = float(d18_val)
                except (TypeError, ValueError):
                    continue
                if math.isnan(d15_float) or math.isnan(d18_float):
                    continue
                iso_sample = nitrate_isotopes.IsotopeSample(d15_float, d18_float)
                process_prior_probs, process_flags = _compute_process_priors(
                    node_id, sample, iso_sample
                )
                hier_samples.append(iso_sample)
                hier_node_ids.append(node_id)
                if process_prior_probs:
                    pooled_process_priors.append(process_prior_probs)
                hierarchical_flags[node_id] = list(process_flags)
            if len(hier_samples) >= 2:
                prior_alpha = None
                if pooled_process_priors:
                    avg_prior: Dict[str, float] = {}
                    for source in iso_sources:
                        values = [
                            p.get(source.name, 0.0) for p in pooled_process_priors
                        ]
                        avg_prior[source.name] = (
                            float(sum(values) / len(values)) if values else 0.0
                        )
                    prior_alpha = [
                        max(
                            avg_prior.get(source.name, 0.0) * process_mcmc_alpha_scale,
                            1e-3,
                        )
                        for source in iso_sources
                    ]
                hier_result = run_mcmc_mixing_hierarchical(
                    samples=hier_samples,
                    sources=iso_sources,
                    sample_ids=hier_node_ids,
                    n_samples=mcmc_n_samples,
                    n_chains=mcmc_n_chains,
                    target_accept=mcmc_target_accept,
                    warmup=mcmc_warmup,
                    prior_alpha=prior_alpha,
                )
                hierarchical_results = hier_result.sample_results
                fallback_used = any(
                    "fallback" in str(msg).lower() for msg in hier_result.warnings
                )
                for node_id in hier_result.sample_results:
                    flags = hierarchical_flags.setdefault(node_id, [])
                    if "hierarchical_mcmc" not in flags:
                        flags.append("hierarchical_mcmc")
                    if (
                        fallback_used
                        and "hierarchical_mcmc_fallback" not in flags
                    ):
                        flags.append("hierarchical_mcmc_fallback")
        except Exception:
            hierarchical_results = {}
            hierarchical_flags = {}

    for node_id, sample in node_rows:

        edges = incoming_edges.get(node_id, [])

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
                reason_code=f"Low Nitrate (Background < {min_mg_L} mg/L)",
            )
            continue

        # --- Dual Isotope Logic (v0.3.0) ---
        used_isotope_model = False
        if iso_enabled and iso_sources:
            try:
                d15_val = sample.get(n15_col)
                d18_val = sample.get(o18_col)

                d15_float = float(d15_val) if d15_val is not None else None
                d18_float = float(d18_val) if d18_val is not None else None

                if (
                    d15_float is not None
                    and d18_float is not None
                    and not math.isnan(d15_float)
                    and not math.isnan(d18_float)
                ):

                    iso_s = nitrate_isotopes.IsotopeSample(d15_float, d18_float)
                    process_flags: List[str] = list(hierarchical_flags.get(node_id, []))
                    process_prior_probs: Optional[Dict[str, float]] = None

                    if process_constraints_enabled and node_id not in hierarchical_results:
                        process_prior_probs, computed_flags = _compute_process_priors(
                            node_id, sample, iso_s
                        )
                        for flag in computed_flags:
                            if flag not in process_flags:
                                process_flags.append(flag)

                    # Use MCMC if enabled, otherwise use analytical method
                    if mcmc_enabled:
                        if node_id in hierarchical_results:
                            mcmc_result = hierarchical_results[node_id]
                        else:
                            from .models.nitrate_isotopes_mcmc import run_mcmc_mixing

                            prior_alpha = None
                            if process_prior_probs:
                                prior_alpha = [
                                    max(
                                        process_prior_probs.get(source.name, 0.0)
                                        * process_mcmc_alpha_scale,
                                        1e-3,
                                    )
                                    for source in iso_sources
                                ]

                            # Run MCMC Bayesian mixing
                            mcmc_result = run_mcmc_mixing(
                                sample=iso_s,
                                sources=iso_sources,
                                n_samples=mcmc_n_samples,
                                n_chains=mcmc_n_chains,
                                target_accept=mcmc_target_accept,
                                warmup=mcmc_warmup,
                                prior_alpha=prior_alpha,
                            )

                        p_man = mcmc_result.source_fractions.get("Manure", 0.0)
                        evidence_list = [f"d15N={d15_float:.1f}", f"d18O={d18_float:.1f}"]
                        if mcmc_result.warnings:
                            evidence_list.extend(mcmc_result.warnings[:2])
                        qc_diagnostics: Optional[Dict[str, float]] = None
                        qc_flags: List[str] = []
                        if isotope_qc_enabled:
                            qc_diagnostics, qc_flags = (
                                nitrate_isotopes.compute_isotope_qc_diagnostics(
                                    sample=iso_s,
                                    sources=iso_sources,
                                    source_probs=mcmc_result.source_fractions,
                                    qc_config=isotope_qc_conf,
                                    prior_probs=process_prior_probs,
                                    posterior_samples=mcmc_result.posterior_samples,
                                    ci_lower=mcmc_result.ci_lower,
                                    ci_upper=mcmc_result.ci_upper,
                                )
                            )
                        isotope_flags = list(
                            dict.fromkeys(
                                ["mcmc_isotope_mixing"] + process_flags + qc_flags
                            )
                        )

                        results[node_id] = NitrateSourceResult(
                            p_manure=p_man,
                            p_fertilizer=mcmc_result.source_fractions.get(
                                "Fertilizer", 1.0 - p_man
                            ),
                            logit_score=None,
                            top_evidence=evidence_list,
                            gating_flags=isotope_flags,
                            ilr_valid=ilr_valid,
                            reason_code="MCMC Bayesian Isotope Mixing",
                            mcmc_result=mcmc_result,
                            source_fractions=mcmc_result.source_fractions,
                            ci_lower=mcmc_result.ci_lower,
                            ci_upper=mcmc_result.ci_upper,
                            diagnostics=qc_diagnostics,
                        )
                        used_isotope_model = True
                    else:
                        # Analytical mixing model
                        probs = nitrate_isotopes.compute_isotope_prob(
                            iso_s, iso_sources, prior_probs=process_prior_probs
                        )

                        p_man = probs.get("Manure", 0.0)
                        qc_diagnostics = None
                        qc_flags: List[str] = []
                        if isotope_qc_enabled:
                            qc_diagnostics, qc_flags = (
                                nitrate_isotopes.compute_isotope_qc_diagnostics(
                                    sample=iso_s,
                                    sources=iso_sources,
                                    source_probs=probs,
                                    qc_config=isotope_qc_conf,
                                    prior_probs=process_prior_probs,
                                )
                            )
                        isotope_flags = list(
                            dict.fromkeys(
                                ["dual_isotope_priority"] + process_flags + qc_flags
                            )
                        )

                        results[node_id] = NitrateSourceResult(
                            p_manure=p_man,
                            p_fertilizer=1.0 - p_man,
                            logit_score=None,
                            top_evidence=[f"d15N={d15_float:.1f}", f"d18O={d18_float:.1f}"],
                            gating_flags=isotope_flags,
                            ilr_valid=ilr_valid,
                            reason_code="Dual Isotope Mixing",
                            source_fractions=probs,
                            diagnostics=qc_diagnostics,
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
        for e in edges:
            if hasattr(e, "transport_model") and e.transport_model == "evap":
                is_evap = True
                if "transport_evap" not in gates:
                    gates.append("transport_evap")

        gate_factor = evap_gate if is_evap else 1.0

        # Base Logit
        logit = prior_logit
        all_evidence_set = set()

        # Node Evidence
        node_logit, node_ev = compute_evidence(sample, stats, weights, gate_factor)
        logit += node_logit
        all_evidence_set.update(node_ev)

        # Edge Evidence (incoming)
        if edges:
            edge_logits = []
            for e in edges:
                u_id = getattr(e, "u", None)
                u_row = node_lookup.get(u_id)
                if u_row is None and u_id is not None:
                    u_row = node_lookup.get(str(u_id))
                if u_row is None:
                    continue
                elogit, eev = compute_edge_evidence(
                    e, u_row, sample, stats, weights
                )
                edge_logits.append(elogit)
                all_evidence_set.update(eev)

            if edge_logits:
                logit += sum(edge_logits) / len(edge_logits)

        # Final Sigmoid
        p_manure = _sigmoid(logit)

        results[node_id] = NitrateSourceResult(
            p_manure=p_manure,
            p_fertilizer=1.0 - p_manure,
            logit_score=logit,
            top_evidence=list(all_evidence_set),
            gating_flags=gates,
            ilr_valid=ilr_valid,
            reason_code="Hydrochemical Ratios (No Isotopes)",
        )

    return results
