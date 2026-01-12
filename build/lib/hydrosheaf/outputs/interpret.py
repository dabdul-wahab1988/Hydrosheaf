"""Geochemical interpretation and report generation for Hydrosheaf results."""

import collections
import json
import math
from typing import Any, Dict, List


def print_interpretation_report(results: List[Dict[str, Any]]):
    """Print a comprehensive geochemical interpretation report based on model results."""
    if not results:
        print("No results found.")
        return

    n_edges = len(results)
    print(f"ANALYSIS OF {n_edges} GROUNDWATER FLOW PATHS")
    print("=" * 60)

    # 1. Transport Mechanisms
    transport_counts = collections.Counter(r.get("transport_model") for r in results)
    print("TRANSPORT PHYSICS:")
    for model, count in transport_counts.items():
        if model:
            print(
                f"  - {model.capitalize()}-Dominant: {count} edges ({100*count/n_edges:.1f}%)"
            )

    # 2. Gibbs Classification
    gibbs_counts: Dict[str, int] = collections.Counter()
    gibbs_penalties = []
    for r in results:
        gm_data = r.get("gibbs_metrics")
        if gm_data:
            if isinstance(gm_data, str):
                try:
                    gm = json.loads(gm_data)
                except json.JSONDecodeError:
                    gm = {}
            else:
                gm = gm_data
            gibbs_counts[gm.get("gibbs_dominance", "unknown")] += 1
        gibbs_penalties.append(r.get("gibbs_penalty", 0) or 0)

    if gibbs_counts:
        print("\nGIBBS DIAGRAM CLASSIFICATION:")
        for cls, count in gibbs_counts.items():
            print(f"  - {cls.capitalize()}: {count} edges ({100*count/n_edges:.1f}%)")
        print(
            f"  - Mean Gibbs Penalty: {sum(gibbs_penalties)/n_edges:.4f} (Lower = better match with isotopes)"
        )

    # 2.5 Isotope Consistency
    iso_slopes = []
    iso_penalties = []
    for r in results:
        m_data = r.get("isotope_metrics")
        if m_data:
            if isinstance(m_data, str):
                try:
                    metrics = json.loads(m_data)
                except json.JSONDecodeError:
                    metrics = {}
            else:
                metrics = m_data
            slope = metrics.get("enrichment_slope")
            if slope is not None and not (isinstance(slope, float) and math.isnan(slope)):
                iso_slopes.append(slope)
        iso_penalties.append(r.get("isotope_consistency_penalty", 0) or 0)

    if iso_slopes:
        print("\nISOTOPE CONSISTENCY:")
        print(
            f"  - Mean Enrichment Slope: {sum(iso_slopes)/len(iso_slopes):.2f} (Typical Evap: 4-6)"
        )
        print(f"  - Mean Consistency Penalty: {sum(iso_penalties)/n_edges:.4f}")
        active_penalties = sum(1 for p in iso_penalties if p > 0.1)
        if active_penalties > 0:
            print(
                f"  - [WARNING] {active_penalties} flow paths show significant Isotope-Chloride mismatch."
            )

    # 3. Geochemical Reactions
    # Identify reaction keys (starting with 'reaction_')
    reaction_keys = [
        k
        for k in results[0].keys()
        if k.startswith("reaction_")
        and k not in ("reaction_iterations", "reaction_converged")
    ]

    print("\nGEOCHEMICAL PROCESSES (Mass Transfer > 1e-4 mmol/L):")
    print(
        f"{'Mineral/Reaction':<25} | {'Activity (%)':<12} | {'Direction':<20} | {'Mean Flux (mmol/L)':<15}"
    )
    print("-" * 80)

    for key in sorted(reaction_keys):
        name = key.replace("reaction_", "")
        values = [r.get(key, 0) or 0 for r in results]

        # Filter significant reactions
        active_values = [v for v in values if v is not None and abs(v) > 1e-4]
        if not active_values:
            continue

        pct_active = 100 * len(active_values) / n_edges
        mean_flux = sum(v for v in values if v is not None) / n_edges

        # Determine dominant direction
        n_pos = sum(1 for v in active_values if v > 0)
        n_neg = sum(1 for v in active_values if v < 0)
        if n_pos > n_neg * 2:
            direction = "Dissolution/Release"
        elif n_neg > n_pos * 2:
            direction = "Precip/Uptake"
        else:
            direction = "Reversible/Mixed"

        print(f"{name:<25} | {pct_active:<12.1f} | {direction:<20} | {mean_flux:<15.4f}")

    print("-" * 80)
    print("\nINTERPRETATION HIGHLIGHTS:")

    # Highlight specific stories
    # Pyrite check
    pyr_keys = [k for k in reaction_keys if "pyrite" in k]
    pyr_active = sum(
        1 for r in results if any(abs(r.get(k, 0) or 0) > 1e-4 for k in pyr_keys)
    )
    if pyr_active > 0:
        print(f"* Oxidative Weathering: Pyrite oxidation detected in {pyr_active} flow paths.")
        # Check if Siderite acts as sink
        sid_flux = sum(r.get("reaction_siderite", 0) or 0 for r in results)
        if sid_flux < -1e-4:
            print(
                "  - Iron precipitation (as Siderite) suggests pH buffering or reducing conditions."
            )
        goe_flux = sum(r.get("reaction_goethite", 0) or 0 for r in results)
        if goe_flux < -1e-4:
            print(
                f"  - Iron precipitation (as Goethite/Fe(OH)3) detected. Mean Flux: {goe_flux/len(results):.4f}"
            )

    # Apatite vs Fluorite check
    apatite_active = sum(1 for r in results if abs(r.get("reaction_apatite", 0) or 0) > 1e-4)
    fluorite_active = sum(
        1 for r in results if abs(r.get("reaction_fluorite", 0) or 0) > 1e-4
    )

    if apatite_active > 0:
        print(
            f"* Apatite Weathering: Detected in {apatite_active} paths. Linking Fluoride to Phosphate source."
        )
    if fluorite_active > 0:
        print(f"* Fluorite Dissolution: Detected in {fluorite_active} paths.")

    # Silicate check
    sil_keys = [
        k
        for k in reaction_keys
        if any(x in k for x in ["albite", "anorthite", "feldspar", "NaSil", "biotite", "chlorite"])
    ]
    sil_flux = sum(sum(r.get(k, 0) or 0 for r in results) for k in sil_keys)
    if sil_flux > 0:
        print("* Silicate Weathering: Net cation release from silicates (Albite/Feldspars/Biotite).")
        if any(abs(r.get("reaction_biotite", 0) or 0) > 1e-4 for r in results):
            print("  - Biotite weathering identified as a potential source of K, Mg, Fe, and Fluoride.")

    # Redox check
    redox_constrained = 0
    for r in results:
        ca_data = r.get("constraints_active")
        if ca_data:
            if isinstance(ca_data, str):
                try: ca = json.loads(ca_data)
                except: ca = {}
            else: ca = ca_data
            if ca.get("redox") == "active":
                redox_constrained += 1
    
    if redox_constrained > 0:
        print(f"* Redox Evolution: {redox_constrained} paths were identified as 'Reducing'. Aerobic oxidation was automatically constrained.")

    # Exchange check
    exch_keys = ["reaction_CaNa_exch", "reaction_MgNa_exch"]
    exch_active = sum(
        1 for r in results if any(abs(r.get(k, 0) or 0) > 1e-4 for k in exch_keys)
    )
    if exch_active > 0:
        print(
            f"* Cation Exchange: Ca-Na/Mg-Na exchange active in {exch_active} paths, modifying final water type."
        )

    # Consistency highlight
    if any(p > 0.5 for p in iso_penalties):
        print(
            "* Salinization Mechanism: Isotope-Chloride mismatch detected. High salt input is likely DISCONNECTED from water loss."
        )

    # 4. Detailed Example
    print("\nDETAILED EXAMPLE (High Reaction Path):")
    # Find edge with highest total reaction mass transfer
    best_edge = max(
        results, key=lambda r: sum(abs(r.get(k, 0) or 0) for k in reaction_keys)
    )
    print(f"path: {best_edge.get('edge_id', 'Unknown')}")
    print(
        f"Model: {best_edge.get('transport_model')} (gamma={best_edge.get('gamma', 1.0):.2f})"
    )
    print("Reactions:")
    for key in reaction_keys:
        val = best_edge.get(key, 0) or 0
        if abs(val) > 1e-4:
            print(f"  - {key.replace('reaction_', '')}: {val:.4f} mmol/L")
