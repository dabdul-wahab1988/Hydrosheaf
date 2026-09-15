"""Bounded Objective 5 (O5) execution for Northern Ghana New (UER).

Implements:
1. Multi-physics cellular sheaf assembly, ungated vs. gated Sheaf Laplacian regularisation.
2. The Tripartite Boundary diagnostics on real UER field data:
   - Candidate-network screening and explicit edge-level missingness
   - Conditional minimax frontier and diminishing model value
   - Inconsistency Trap & Conflict Isolation: Diagnostic sheaf re-solve with explicit evidence limits
3. Conditional Minimax Measurement Design from UER tritium cohort summaries and 38 individual boreholes.
4. Generates 600 DPI publication figures, audit tables, and comprehensive report with exact metrics.
"""

from __future__ import annotations

import json
import hashlib
import math
import os
from pathlib import Path
import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[4] / "July_2026" / "NeutroProject" / "Groundwater" / "Hydrosheaf"
if not ROOT.exists():
    # Fallback to local repo resolution
    for p in Path(__file__).resolve().parents:
        if (p / "hydrosheaf").exists():
            ROOT = p
            break
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from hydrosheaf.config import Config
from hydrosheaf.graph.types import Edge
from hydrosheaf.nuclear.joint_lpm import tracer_response_kernel
from hydrosheaf.nuclear.ttd_certified_design import (
    CertifiedCandidateTracer,
    evaluate_worst_case_ambiguity,
    solve_budgeted_minimax_design,
    solve_certified_measurement_design,
)
from hydrosheaf.nuclear.ttd_identified import AgeFunctional, TracerConstraint
from hydrosheaf.sheaf.directed_section import (
    DirectedEdgeMap,
    build_edge_maps,
    compute_edge_section_residuals,
    solve_directed_section,
)

# Output paths.  An override keeps corrected/audit reruns beside the original
# package so the previously generated artifacts remain recoverable.
_OUT_OVERRIDE = os.environ.get("HYDROSHEAF_O5_OUT_DIR")
OUT_DIR = Path(_OUT_OVERRIDE).expanduser() if _OUT_OVERRIDE else ROOT / "outputs" / "uer_objective5"
FIG_DIR = OUT_DIR / "figures"
TAB_DIR = OUT_DIR / "tables"
DATA_CSV = ROOT / "data" / "FieldData" / "derived" / "uer_field_integration_dataset.csv"
EDGES_CSV = ROOT / "outputs" / "uer_field_integration" / "uer_network_edges.csv"
COST_MANIFEST = ROOT / "provenance" / "uer_objective5_cost_manifest_2026-09-15.json"


def load_cost_manifest() -> dict:
    """Load and validate the dated, source-linked Objective 5 cost basis."""
    if not COST_MANIFEST.exists():
        raise FileNotFoundError(f"Objective 5 cost manifest is missing: {COST_MANIFEST}")
    with COST_MANIFEST.open("r", encoding="utf-8") as handle:
        manifest = json.load(handle)

    if manifest.get("target_currency") != "USD":
        raise ValueError("Objective 5 cost manifest must define USD as its target currency")
    candidates = manifest.get("candidates")
    if not isinstance(candidates, list) or not candidates:
        raise ValueError("Objective 5 cost manifest contains no candidate costs")

    option_ids = []
    for row in candidates:
        option_id = row.get("option_id")
        cost = row.get("usd_per_sample")
        if not option_id or option_id in option_ids:
            raise ValueError(f"Objective 5 cost manifest has an invalid or duplicate option_id: {option_id!r}")
        if not isinstance(cost, (int, float)) or not math.isfinite(float(cost)) or float(cost) <= 0.0:
            raise ValueError(f"Objective 5 cost manifest has an invalid USD cost for {option_id!r}: {cost!r}")
        for required in ("lab", "source_url", "price_date", "quote_status"):
            if not row.get(required):
                raise ValueError(f"Objective 5 cost manifest is missing {required!r} for {option_id!r}")
        option_ids.append(option_id)
    return manifest


def setup_matplotlib():
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman", "DejaVu Serif"],
        "mathtext.fontset": "cm",
        "font.size": 9.5,
        "axes.titlesize": 11,
        "axes.labelsize": 10,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.fontsize": 8.5,
        "figure.dpi": 600,
        "savefig.dpi": 600,
    })


def run_sheaf_and_tripartite_boundary(df: pd.DataFrame, df_edges: pd.DataFrame):
    print("\n--- Pillar 1: Sheaf Assembly, Multi-Physics Screening & Gating Verification ---")
    
    ion_cols = ["ca_mmol_L", "mg_mmol_L", "na_mmol_L", "k_mmol_L", "hco3_mmol_L", "cl_mmol_L", "so4_mmol_L", "no3_mmol_L"]
    samples_dict = {}
    for _, row in df.iterrows():
        sid = str(row["node_id"])
        cbe = float(row["cbe_percent"]) if pd.notna(row["cbe_percent"]) else 0.0
        samples_dict[sid] = {
            "site_id": sid,
            "sample_type": str(row["sample_type"]),
            "elevation": float(row["elevation_m"]),
            "chem_vector": [float(row[col]) if pd.notna(row[col]) else 0.0 for col in ion_cols],
            "cbe_percent": cbe,
            "cbe_pass": abs(cbe) <= 10.0,
            "d18O": float(row["d18O_permil"]) if pd.notna(row["d18O_permil"]) else None,
            "d2H": float(row["d2H_permil"]) if pd.notna(row["d2H_permil"]) else None,
            "tritium": float(row["tritium_TU"]) if pd.notna(row["tritium_TU"]) else None,
            "d15N": float(row["d15N_NO3_permil_air"]) if pd.notna(row["d15N_NO3_permil_air"]) else None,
            "d18O_NO3": float(row["d18O_NO3_permil_VSMOW"]) if pd.notna(row["d18O_NO3_permil_VSMOW"]) else None,
            "nitrate_class": str(row["nitrate_source_candidate"]),
            "community": str(row["community"]),
        }

    edge_objs = []
    for _, erow in df_edges.iterrows():
        edge_objs.append(
            Edge(
                edge_id=str(erow["edge_id"]),
                u=str(erow["u"]),
                v=str(erow["v"]),
                attrs={
                    "distance_km": float(erow["distance_km"]),
                    "delta_h": float(erow["delta_h_m"]),
                    "p_uv": float(erow["p_uv"]) if pd.notna(erow["p_uv"]) else 0.5,
                }
            )
        )

    cfg = Config()
    cfg.geologic_bias = "crystalline"
    cfg.active_minerals = ["calcite", "dolomite", "albite", "halite", "pyrite_oxidation_aerobic"]
    cfg.transport_models_enabled = ["evap"]
    cfg.weights = [1.0] * len(ion_cols)

    sample_chem_vectors = {sid: s["chem_vector"] for sid, s in samples_dict.items()}
    edge_maps = build_edge_maps(edge_objs, sample_chem_vectors, cfg)
    print(f"Built {len(edge_maps)} Sheaf DirectedEdgeMaps across {len(sample_chem_vectors)} nodes.")

    # -------------------------------------------------------------
    # STAGE 1: UNGATED / FORCED INVERSION (ALL 1,600 CANDIDATE EDGES)
    # -------------------------------------------------------------
    node_ids = list(samples_dict.keys())
    solved_states_ungated = solve_directed_section(
        node_ids,
        edge_maps.values(),
        sample_chem_vectors,
        obs_weight=1.0,
        diag_eps=0.01,
        non_negative=True,
    )
    residuals_ungated = compute_edge_section_residuals(edge_maps, solved_states_ungated, cfg.weights)
    dirichlet_energy_ungated = sum(r ** 2 for r in residuals_ungated.values())
    print(f"[Ungated Solve] Solved across {len(edge_maps)} edges. Total Dirichlet Energy: {dirichlet_energy_ungated:.2f}")

    # -------------------------------------------------------------
    # MULTI-CHANNEL CONFLICT SCREENING & MISSINGNESS TRACKING
    # -------------------------------------------------------------
    audit_rows = []
    # Keep missingness separate from a passed measurement.  Chemistry and
    # CBE are available for every edge; isotope channels are tested only when
    # both endpoints contain the relevant observations.
    synergy_counts = {"heads_only": 0, "heads_plus_chem": 0, "heads_chem_isotopes": 0, "fully_coherent": 0}
    channel_denominators = {
        "cbe_tested": 0, "cbe_pass": 0, "cbe_fail": 0,
        "tritium_dual_tested": 0, "tritium_pass": 0, "tritium_fail": 0, "tritium_untested": 0,
        "nitrate_dual_tested": 0, "nitrate_pass": 0, "nitrate_fail": 0, "nitrate_untested": 0,
        "chem_tested": 0, "chem_pass": 0, "chem_fail": 0,
    }

    for e in edge_objs:
        u_data = samples_dict[e.u]
        v_data = samples_dict[e.v]
        res_e_ungated = residuals_ungated.get(e.edge_id, 0.0)
        dist_km = e.attrs.get("distance_km", 1.0)
        delta_h = e.attrs.get("delta_h", 0.0)
        sample_pair_type = f"{u_data['sample_type']}->{v_data['sample_type']}"

        # 1. Chemical QC & Stoichiometry
        channel_denominators["chem_tested"] += 1
        chem_conflict = res_e_ungated >= 1.5
        if chem_conflict:
            channel_denominators["chem_fail"] += 1
        else:
            channel_denominators["chem_pass"] += 1

        channel_denominators["cbe_tested"] += 1
        cbe_ok = u_data["cbe_pass"] and v_data["cbe_pass"]
        if cbe_ok:
            channel_denominators["cbe_pass"] += 1
        else:
            channel_denominators["cbe_fail"] += 1

        # 2. Radioactive Tritium Inversion
        # Tested ONLY if both endpoints possess measured tritium
        tritium_tested = (u_data["tritium"] is not None) and (v_data["tritium"] is not None)
        tritium_inversion = False
        tritium_status = "UNTESTED_MISSING_DATA"
        if tritium_tested:
            channel_denominators["tritium_dual_tested"] += 1
            if v_data["tritium"] > u_data["tritium"] + 0.5:
                tritium_inversion = True
                tritium_status = "FAIL_TRITIUM_INVERSION"
                channel_denominators["tritium_fail"] += 1
            else:
                tritium_status = "TESTED_COHERENT"
                channel_denominators["tritium_pass"] += 1
        else:
            channel_denominators["tritium_untested"] += 1

        # 3. Nitrate Contamination Spike
        # The source-candidate label is derived and remains populated for
        # missing isotope values.  It is therefore not a valid missingness
        # test; require both d15N and d18O at both endpoints.
        nitrate_tested = all(
            u_data[key] is not None and v_data[key] is not None
            for key in ("d15N", "d18O_NO3")
        )
        nitrate_conflict = False
        nitrate_status = "UNTESTED_MISSING_DATA"
        if nitrate_tested:
            channel_denominators["nitrate_dual_tested"] += 1
            if "Manure" in v_data["nitrate_class"] and "Low Nitrate" in u_data["nitrate_class"]:
                nitrate_conflict = True
                nitrate_status = "FAIL_MANURE_CONTAMINATION_SPIKE"
                channel_denominators["nitrate_fail"] += 1
            else:
                nitrate_status = "TESTED_COHERENT"
                channel_denominators["nitrate_pass"] += 1
        else:
            channel_denominators["nitrate_untested"] += 1

        # Overall conflict flag. CBE failure is a quality-control conflict,
        # separate from a physical mechanism diagnosis, but it is included in
        # the admitted quantitative section.
        cbe_conflict = not cbe_ok
        is_conflict = chem_conflict or cbe_conflict or tritium_inversion or nitrate_conflict

        # Progressive filter count
        synergy_counts["heads_only"] += 1
        if not chem_conflict:
            synergy_counts["heads_plus_chem"] += 1
            if not cbe_conflict and not tritium_inversion:
                synergy_counts["heads_chem_isotopes"] += 1
                if not nitrate_conflict:
                    synergy_counts["fully_coherent"] += 1

        conflict_types = []
        if chem_conflict:
            conflict_types.append("Stoichiometric Breakdown")
        if cbe_conflict:
            conflict_types.append("Charge Balance QC Failure")
        if tritium_inversion:
            conflict_types.append("Tritium Vertical Inversion Candidate")
        if nitrate_conflict:
            conflict_types.append("Anthropogenic Nitrate Spike")

        audit_rows.append({
            "edge_id": e.edge_id,
            "u": e.u,
            "v": e.v,
            "sample_pair_type": sample_pair_type,
            "u_community": u_data["community"],
            "v_community": v_data["community"],
            "distance_km": dist_km,
            "delta_h_m": delta_h,
            "res_ungated": res_e_ungated,
            "cbe_ok": cbe_ok,
            "cbe_conflict": cbe_conflict,
            "u_tritium": u_data["tritium"],
            "v_tritium": v_data["tritium"],
            "tritium_status": tritium_status,
            "nitrate_status": nitrate_status,
            "chem_conflict": chem_conflict,
            "tritium_conflict": tritium_inversion,
            "nitrate_conflict": nitrate_conflict,
            "is_conflict": is_conflict,
            "conflict_diagnosis": " + ".join(conflict_types) if conflict_types else "None (Coherent)",
            "gating_action": "REJECT_EDGE_ISOLATE_CONFLICT" if is_conflict else "ADMIT_COHERENT_SECTION"
        })

    df_audit = pd.DataFrame(audit_rows)

    # -------------------------------------------------------------
    # STAGE 2: GATED INVERSION (RE-SOLVING ON THE ADMITTED SUBGRAPH)
    # -------------------------------------------------------------
    coherent_edge_ids = set(df_audit[~df_audit["is_conflict"]]["edge_id"].values)
    edge_maps_gated = [em for em in edge_maps.values() if em.edge.edge_id in coherent_edge_ids]
    
    solved_states_gated = solve_directed_section(
        node_ids,
        edge_maps_gated,
        sample_chem_vectors,
        obs_weight=1.0,
        diag_eps=0.01,
        non_negative=True,
    )
    residuals_gated = compute_edge_section_residuals(
        {em.edge.edge_id: em for em in edge_maps_gated},
        solved_states_gated,
        cfg.weights
    )
    dirichlet_energy_gated = sum(r ** 2 for r in residuals_gated.values())

    # Map gated residuals back to audit dataframe
    df_audit["res_gated"] = df_audit["edge_id"].map(residuals_gated)
    
    # Calculate energy & residual statistics
    coherent_mask = ~df_audit["is_conflict"]
    mean_res_ungated_coherent = df_audit.loc[coherent_mask, "res_ungated"].mean()
    mean_res_gated_coherent = df_audit.loc[coherent_mask, "res_gated"].mean()
    mean_res_conflicts = df_audit.loc[df_audit["is_conflict"], "res_ungated"].mean()
    coherent_energy_ungated = float((df_audit.loc[coherent_mask, "res_ungated"] ** 2).sum())
    coherent_energy_gated = float((df_audit.loc[coherent_mask, "res_gated"] ** 2).sum())
    coherent_rms_ungated = float(np.sqrt(np.mean(df_audit.loc[coherent_mask, "res_ungated"] ** 2)))
    coherent_rms_gated = float(np.sqrt(np.mean(df_audit.loc[coherent_mask, "res_gated"] ** 2)))

    print(f"[Gated Re-Solve] Coherent Edges: {len(edge_maps_gated)} (Dirichlet Energy: {dirichlet_energy_gated:.2f}).")
    print(f"  Mean Coherent Residual (Ungated): {mean_res_ungated_coherent:.4f}")
    print(f"  Mean Coherent Residual (Gated):   {mean_res_gated_coherent:.4f}")
    print(f"  Same-edge coherent energy:         {coherent_energy_ungated:.2f} -> {coherent_energy_gated:.2f}")
    print(f"  Mean Conflict Residual (Ungated): {mean_res_conflicts:.4f}")

    # Export audit files
    df_audit.to_csv(TAB_DIR / "uer_objective5_tripartite_boundary_audit.csv", index=False)
    
    df_conflicts = df_audit[df_audit["is_conflict"]].copy()
    df_conflicts.to_csv(TAB_DIR / "uer_objective5_sheaf_conflict_localisation.csv", index=False)

    # Save Gating Comparison CSV
    df_gating_comp = pd.DataFrame([{
        "solve_mode": "Ungated (Forced All)",
        "edge_count": len(edge_objs),
        "dirichlet_energy": dirichlet_energy_ungated,
        "mean_residual_all": df_audit["res_ungated"].mean(),
        "mean_residual_coherent_edges": mean_res_ungated_coherent,
        "rms_residual_coherent_edges": coherent_rms_ungated,
        "coherent_edge_energy": coherent_energy_ungated,
        "edge_set_comparison": "same coherent edges evaluated under all-edge solve",
    }, {
        "solve_mode": "Sheaf Gated (Isolated)",
        "edge_count": len(edge_maps_gated),
        "dirichlet_energy": dirichlet_energy_gated,
        "mean_residual_all": df_audit.loc[coherent_mask, "res_gated"].mean(),
        "mean_residual_coherent_edges": mean_res_gated_coherent,
        "rms_residual_coherent_edges": coherent_rms_gated,
        "coherent_edge_energy": coherent_energy_gated,
        "edge_set_comparison": "same coherent edges evaluated under gated solve",
    }])
    df_gating_comp.to_csv(TAB_DIR / "uer_objective5_sheaf_gating_comparison.csv", index=False)

    # Markdown Summary
    md_lines = [
        "# Table O5.3: Sheaf Conflict Isolation, Gating Verification, and Evidence Coverage\n",
        f"- Total candidate topographic flow edges evaluated: **{len(df_audit)}** (Tier-C DEM proxy)",
        f"- Edges flagged and isolated by chemistry/QC/isotope screens: **{len(df_conflicts)} ({len(df_conflicts)/len(df_audit):.1%})**",
        f"- Candidate edges admitted to the regularised sheaf section: **{len(edge_maps_gated)} ({len(edge_maps_gated)/len(df_audit):.1%})**\n",
        "### 1. Evidence Channel Coverage & Tested Denominators",
        "| Evidence Channel | Dual-Endpoint Tested | Passed | Flagged / Conflict | Untested (Missing Data) |",
        "| :--- | :---: | :---: | :---: | :---: |",
        f"| **Stoichiometric Chemistry** | {channel_denominators['chem_tested']} | {channel_denominators['chem_pass']} | {channel_denominators['chem_fail']} | 0 |",
        f"| **Charge Balance QC** ($|\\text{{CBE}}| \\le 10\\%$) | {channel_denominators['cbe_tested']} | {channel_denominators['cbe_pass']} | {channel_denominators['cbe_tested'] - channel_denominators['cbe_pass']} | 0 |",
        f"| **Radioactive Tritium** ($^3\\text{{H}}$) | {channel_denominators['tritium_dual_tested']} | {channel_denominators['tritium_pass']} | {channel_denominators['tritium_fail']} | {channel_denominators['tritium_untested']} |",
        f"| **Dual Nitrate Isotopes** ($\\delta^{{15}}\\text{{N}}/\\delta^{{18}}\\text{{O}}$) | {channel_denominators['nitrate_dual_tested']} | {channel_denominators['nitrate_pass']} | {channel_denominators['nitrate_fail']} | {channel_denominators['nitrate_untested']} |\n",
        "*Note: Only 17 candidate edges possess measured tritium at both endpoints. Only edges with both nitrate isotope endpoints are tested on that channel; missing values remain 'Untested'. Candidate edges are not verified groundwater flowpaths.*\n",
        "### 2. Inconsistency Diagnostic: Gated vs. Ungated Sheaf Re-solve",
        "| Inversion Solve Mode | Network Edges | Total Sheaf Dirichlet Energy | Mean Residual on Same Admitted Edges | Same-edge Energy | Interpretation |",
        "| :--- | :---: | :---: | :---: | :---: | :--- |",
        f"| **Ungated (Forced All Edges)** | {len(edge_objs)} | {dirichlet_energy_ungated:.2f} | {mean_res_ungated_coherent:.4f} | {coherent_energy_ungated:.2f} | Reference fit; includes flagged edges |",
        f"| **Sheaf Gated (Admitted Re-solve)** | {len(edge_maps_gated)} | {dirichlet_energy_gated:.2f} | {mean_res_gated_coherent:.4f} | {coherent_energy_gated:.2f} | Diagnostic only; total energy is not an apples-to-apples comparison |",
        f"\n*The total energy changes with the number of included edges. No independent flow or age truth is available here, so this re-solve does not prove error-propagation protection.*",
    ]
    with open(TAB_DIR / "uer_objective5_sheaf_conflict_localisation.md", "w", encoding="utf-8") as f:
        f.write("\n".join(md_lines))

    return df_audit, synergy_counts, channel_denominators, df_gating_comp


def run_minimax_design_on_empirical_uer(df: pd.DataFrame):
    print("\n--- Pillar 2 & 3: Conditional Minimax Design from UER Tritium Cohort Summaries ---")
    
    # 1. Load and verify empirical tritium samples
    df_3h = df[df["tritium_TU"].notna()].copy()
    if len(df_3h) == 0:
        raise ValueError("Critical failure: No empirical tritium records found in dataset!")
    
    n_3h = len(df_3h)
    med_3h = float(df_3h["tritium_TU"].median())
    mean_3h = float(df_3h["tritium_TU"].mean())
    min_3h = float(df_3h["tritium_TU"].min())
    max_3h = float(df_3h["tritium_TU"].max())
    std_3h = float(df_3h["tritium_TU"].std())

    print(f"Empirical UER Tritium Dataset: n={n_3h} samples, median={med_3h:.3f} TU, mean={mean_3h:.3f} TU, range={min_3h:.2f}..{max_3h:.2f} TU, std={std_3h:.3f} TU.")

    # 2. Define three cohort summaries directly from the UER distribution.
    # These are summary priors, not 42 simultaneous observation constraints.
    cohort_low = df_3h[df_3h["tritium_TU"] <= 1.5]
    cohort_med = df_3h[(df_3h["tritium_TU"] > 1.5) & (df_3h["tritium_TU"] <= 3.0)]
    cohort_high = df_3h[df_3h["tritium_TU"] > 3.0]

    prior_cohorts = [
        ("Sub-modern Basement Cohort", float(cohort_low["tritium_TU"].mean()), float(cohort_low["tritium_TU"].std()) if len(cohort_low) > 1 else 0.20, len(cohort_low)),
        ("Regional Median Screening Cohort", med_3h, float(cohort_med["tritium_TU"].std()) if len(cohort_med) > 1 else 0.35, len(cohort_med)),
        ("Active Modern Recharge Cohort", float(cohort_high["tritium_TU"].mean()), float(cohort_high["tritium_TU"].std()) if len(cohort_high) > 1 else 0.45, len(cohort_high)),
    ]

    age_grid = np.linspace(0.5, 80.0, 40)
    sample_year = 2016.0

    target_mtt = AgeFunctional(
        name="mean_transit_time",
        coefficients=age_grid / 80.0,
        maximum_reportable_width=0.20,
        units="fraction_of_80yr"
    )

    resp_3h = tracer_response_kernel("3H", age_grid, sample_year=sample_year)

    # The frontier cost is loaded from a dated provenance manifest. The
    # manifest distinguishes public list prices from vendor quotations and
    # keeps missing field/logistics components visible.
    cost_manifest = load_cost_manifest()
    cost_rows = {row["option_id"]: row for row in cost_manifest["candidates"]}
    cost_by_option = {option_id: float(row["usd_per_sample"]) for option_id, row in cost_rows.items()}
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    TAB_DIR.mkdir(parents=True, exist_ok=True)
    snapshot_path = OUT_DIR / "cost_manifest_snapshot.json"
    with snapshot_path.open("w", encoding="utf-8") as handle:
        json.dump(cost_manifest, handle, indent=2, ensure_ascii=False)
    snapshot_hash = hashlib.sha256(snapshot_path.read_bytes()).hexdigest()
    (OUT_DIR / "cost_manifest_snapshot.sha256").write_text(
        f"{snapshot_hash}  {snapshot_path.name}\n", encoding="utf-8"
    )
    pd.DataFrame(cost_manifest["candidates"]).to_csv(
        TAB_DIR / "uer_objective5_cost_basis.csv", index=False
    )
    print(
        "Cost basis: "
        f"{cost_manifest['manifest_id']} ({cost_manifest['quote_status']}; "
        f"{cost_manifest.get('evidence_classification', 'UNCLASSIFIED')})"
    )

    def cost_metadata(option_id: str) -> dict:
        row = cost_rows[option_id]
        return {
            "cost_usd_per_sample": float(row["usd_per_sample"]),
            "cost_basis": row["price_basis"],
            "cost_source": row["lab"],
            "cost_price_date": row["price_date"],
            "cost_quote_status": row["quote_status"],
            "cost_evidence_classification": cost_manifest.get(
                "evidence_classification", "UNCLASSIFIED"
            ),
            "cost_manifest_id": cost_manifest["manifest_id"],
            "cost_scope": cost_manifest["frontier_unit"],
        }

    candidates = [
        CertifiedCandidateTracer(
            option_id="SF6",
            tracer="SF6",
            sample_year=sample_year,
            error_bound=0.04,
            cost=cost_by_option["SF6"],
            response=tracer_response_kernel("SF6", age_grid, sample_year=sample_year),
            metadata={
                "unit": "pptv",
                "model": "exponential_lpm",
                "target": "Modern water (<35 yr)",
                "assumption_status": "hypothetical_candidate",
                "history": "default_approximate_northern_hemisphere",
                **cost_metadata("SF6"),
            }
        ),
        CertifiedCandidateTracer(
            option_id="CFC12",
            tracer="CFC12",
            sample_year=sample_year,
            error_bound=0.05,
            cost=cost_by_option["CFC12"],
            response=tracer_response_kernel("CFC12", age_grid, sample_year=sample_year),
            metadata={
                "unit": "pptv",
                "model": "exponential_lpm",
                "target": "Industrial 1950-2000 plateau",
                "assumption_status": "hypothetical_candidate",
                "history": "default_approximate_northern_hemisphere",
                **cost_metadata("CFC12"),
            }
        ),
        CertifiedCandidateTracer(
            option_id="3H_resample",
            tracer="3H_high_prec",
            sample_year=sample_year,
            error_bound=0.10,
            cost=cost_by_option["3H_resample"],
            response=resp_3h,
            metadata={
                "unit": "TU",
                "model": "decay_lpm",
                "target": "Confirmatory Tritium",
                **cost_metadata("3H_resample"),
            }
        ),
        CertifiedCandidateTracer(
            option_id="14C",
            tracer="14C",
            sample_year=sample_year,
            error_bound=1.0,
            cost=cost_by_option["14C"],
            response=tracer_response_kernel("14C", age_grid, sample_year=sample_year),
            metadata={
                "unit": "pMC",
                "model": "decay_lpm",
                "target": "Sub-modern / fossil fraction",
                **cost_metadata("14C"),
            }
        ),
        CertifiedCandidateTracer(
            option_id="3H_3He",
            tracer="3H/3He",
            sample_year=sample_year,
            error_bound=0.02,
            cost=cost_by_option["3H_3He"],
            response=tracer_response_kernel("3H/3He", age_grid, sample_year=sample_year),
            metadata={
                "unit": "TU_equivalent",
                "model": "linear_lpm",
                "target": "High-precision age dating",
                "assumption_status": "hypothetical_candidate",
                "history": "default_tritium_input",
                **cost_metadata("3H_3He"),
            }
        ),
    ]

    # Include every subset-cost breakpoint so the discrete frontier cannot
    # skip a portfolio transition, plus common caps retained for comparison.
    subset_costs = {0.0}
    for mask in range(1, 1 << len(candidates)):
        subset_costs.add(round(
            sum(
                candidates[i].cost
                for i in range(len(candidates))
                if mask & (1 << i)
            ),
            2,
        ))
    requested_caps = {
        float(value) for value in cost_manifest.get("frontier_budget_caps_usd", [])
    }
    budgets = sorted(subset_costs | requested_caps)
    all_frontier_records = []

    for label, tu_val, tu_sig, n_samples in prior_cohorts:
        prior_constraint = TracerConstraint(
            tracer="3H",
            observed=tu_val,
            sigma=tu_sig,
            response=resp_3h,
            units="TU"
        )
        b0_amb = None
        for b in budgets:
            cert = solve_budgeted_minimax_design(
                age_grid,
                (prior_constraint,),
                candidates,
                target_mtt,
                budget=b,
            )
            mtt_amb_yr = cert.achieved_ambiguity * 80.0
            if b == 0.0:
                b0_amb = mtt_amb_yr
            red_pct = (1.0 - mtt_amb_yr / b0_amb) * 100.0 if b0_amb and b0_amb > 0 else 0.0

            all_frontier_records.append({
                "prior_regime": label,
                "cohort_size": n_samples,
                "baseline_tu": tu_val,
                "baseline_sigma": tu_sig,
                "budget": b,
                "selected_options": "+".join(cert.selected_option_ids) if cert.selected_option_ids else "None (Prior 3H only)",
                "total_cost": cert.total_cost,
                "mtt_ambiguity_yr": mtt_amb_yr,
                "mtt_ambiguity_norm": cert.achieved_ambiguity,
                "status": cert.status,
                "ambiguity_reduction_pct": red_pct,
                "cost_manifest_id": cost_manifest["manifest_id"],
                "cost_basis_status": cost_manifest["quote_status"],
            })

    df_frontiers = pd.DataFrame(all_frontier_records)
    df_frontiers.to_csv(TAB_DIR / "uer_objective5_minimax_pareto_frontier.csv", index=False)

    # 3. Individual Borehole Minimax Designs (Row-level provenance)
    borehole_records = []
    df_bh_3h = df_3h[df_3h["sample_type"] == "BH"].copy()
    for _, row in df_bh_3h.iterrows():
        tu_val = float(row["tritium_TU"])
        tu_sig = max(0.15, tu_val * 0.10)
        p_con = TracerConstraint(tracer="3H", observed=tu_val, sigma=tu_sig, response=resp_3h, units="TU")
        cert_0 = solve_budgeted_minimax_design(age_grid, (p_con,), candidates, target_mtt, budget=0.0)
        cert_800 = solve_budgeted_minimax_design(age_grid, (p_con,), candidates, target_mtt, budget=800.0)
        
        borehole_records.append({
            "node_id": row["node_id"],
            "site_id": row["site_id"],
            "community": row["community"],
            "measured_3H_TU": tu_val,
            "baseline_ambiguity_yr": cert_0.achieved_ambiguity * 80.0,
            "opt_800_cost": cert_800.total_cost,
            "opt_800_portfolio": "+".join(cert_800.selected_option_ids),
            "opt_800_ambiguity_yr": cert_800.achieved_ambiguity * 80.0,
            "opt_800_status": cert_800.status,
            "cost_manifest_id": cost_manifest["manifest_id"],
            "cost_basis_status": cost_manifest["quote_status"],
        })
    df_bh_designs = pd.DataFrame(borehole_records)
    df_bh_designs.to_csv(TAB_DIR / "uer_objective5_well_specific_minimax_designs.csv", index=False)
    print(f"Generated sample-specific minimax designs for {len(df_bh_designs)} UER boreholes.")

    # 4. Mathematical Impossibility Witness (14C on modern water < 25 yr)
    target_modern = AgeFunctional(
        name="modern_water_fraction_25yr",
        coefficients=(age_grid <= 25.0).astype(float),
        maximum_reportable_width=0.05,
        units="fraction"
    )
    median_sigma = float(cohort_med["tritium_TU"].std()) if len(cohort_med) > 1 else 0.35
    p_median = TracerConstraint(tracer="3H", observed=med_3h, sigma=median_sigma, response=resp_3h, units="TU")
    c14_cand = [c for c in candidates if c.option_id == "14C"]

    cert_imposs = solve_certified_measurement_design(
        age_grid,
        (p_median,),
        c14_cand,
        target_modern,
        target_tolerance=0.05,
    )
    print(f"Computed Impossibility Witness: Status={cert_imposs.status}, Achieved Ambiguity={cert_imposs.achieved_ambiguity:.4f} (Target <= {cert_imposs.target_tolerance}).")

    # 5. Save Markdown Pareto Frontier Table
    md_front_lines = [
        "# Table O5.1: Conditional Minimax Measurement-Design Frontier from UER Tritium Cohorts\n",
        f"- Target Decision Functional: **Mean Transit Time (MTT)** over finite grid $[0.5, 80.0\\text{{ yr}}]$",
        f"- Certification Target Tolerance: $\\delta = 0.20 \\times 80\\text{{ yr}} = **16.0\\text{{ yr}}**$",
        f"- Equivalence Margin for Redundancy Plateau: $\\Delta W < 0.05\\text{{ yr}}$",
        f"- Three summary priors derived from **{n_3h} empirical UER samples** (median ${med_3h:.3f}\\text{{ TU}}$); each frontier row uses one cohort-level 3H constraint",
        f"- Candidate tracer responses use the declared default histories at an assumed sample year of **{sample_year:.0f}**.",
        f"- Cost manifest: **{cost_manifest['manifest_id']}**. Evidence class: **{cost_manifest.get('evidence_classification', 'UNCLASSIFIED')}**. Rates below are provisional public-list analytical prices, not vendor quotations; field, shipping, taxes and unpriced QA/QC remain excluded.\n",
        "| Option | Laboratory | Posted rate | USD rate used | Price date | Quote status | Evidence class |",
        "| :--- | :--- | ---: | ---: | :--- | :--- | :--- |",
    ]
    for row in cost_manifest["candidates"]:
        md_front_lines.append(
            f"| {row['option_id']} | {row['lab']} | {row['original_currency']} {row['original_price']:,.2f} | ${row['usd_per_sample']:,.2f} | {row['price_date']} | `{row['quote_status']}` | `{cost_manifest.get('evidence_classification', 'UNCLASSIFIED')}` |"
        )
    md_front_lines.extend([
        "",
        "| Prior Cohort | Allocation Budget (USD) | Actual Cost (USD) | Selected Tracer Portfolio | Worst-Case MTT Ambiguity | Reduction (%) | Mathematical Status |",
        "| :--- | :---: | :---: | :--- | :---: | :---: | :--- |",
    ])
    for _, r in df_frontiers.iterrows():
        md_front_lines.append(
            f"| {r['prior_regime']} | ${r['budget']:,.0f} | ${r['total_cost']:,.0f} | {r['selected_options']} | **{r['mtt_ambiguity_yr']:.3f} yr** | {r['ambiguity_reduction_pct']:.1f}% | `{r['status']}` |"
        )
    with open(TAB_DIR / "uer_objective5_minimax_pareto_frontier.md", "w", encoding="utf-8") as f:
        f.write("\n".join(md_front_lines))

    return df_frontiers, cert_imposs, age_grid, df_bh_designs


def generate_publication_figures(df_audit: pd.DataFrame, df_frontiers: pd.DataFrame, cert_imposs, age_grid, synergy_counts):
    print("\n--- Pillar 4: Generating Publication-Grade 600 DPI Figures ---")
    setup_matplotlib()

    # ------------------------------------------------------------------
    # FIGURE O5.1: Tripartite Boundary (Canonical 2x2 with 3 plots)
    # ------------------------------------------------------------------
    fig, axes = plt.subplots(2, 2, figsize=(12.0, 9.4), dpi=600, gridspec_kw={'hspace': 0.34, 'wspace': 0.28})
    ax1 = axes[0, 0]
    ax2 = axes[0, 1]
    ax3 = axes[1, 0]
    axes[1, 1].axis('off')

    # Subplot A: Synergistic Compression of Admissible Flow Network
    stages = ["Heads Only\n(Tier-C DEM)", "+ Major\nChemistry", "+ CBE +\nTritium", "Admitted\nCandidate Set"]
    counts = [
        synergy_counts["heads_only"],
        synergy_counts["heads_plus_chem"],
        synergy_counts["heads_chem_isotopes"],
        synergy_counts["fully_coherent"]
    ]
    colors = ["#94a3b8", "#38bdf8", "#3b82f6", "#1e3a8a"]
    bars = ax1.bar(stages, counts, color=colors, edgecolor="black", linewidth=0.8, width=0.55)
    for bar, count in zip(bars, counts):
        ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 25, f"{count:,}\n({count/counts[0]:.1%})",
                 ha="center", va="bottom", fontsize=8.2, fontweight="bold")
    ax1.set_ylabel("Candidate Flow Edges in Network", fontsize=10.0, fontweight="bold")
    ax1.set_title("(a) Synergistic Network Compression", fontsize=11.0, fontweight="bold", loc="left", pad=8)
    ax1.set_ylim(0, 1950)
    ax1.grid(True, ls="--", alpha=0.35, axis="y")

    # Subplot B: the costed frontier. Use the actual budget grid generated
    # from the manifest's subset breakpoints instead of a stale hard-coded
    # price grid.
    med_curve = df_frontiers[df_frontiers["prior_regime"] == "Regional Median Screening Cohort"]
    budget_pts = med_curve["budget"].tolist()
    mtt_pts = med_curve["mtt_ambiguity_yr"].tolist()
    opt_row = med_curve[med_curve["budget"] == 800.0].iloc[0]
    opt_val = float(opt_row["mtt_ambiguity_yr"])
    opt_portfolio = str(opt_row["selected_options"])
    opt_cost = float(opt_row["total_cost"])
    max_budget = max(budget_pts)
    
    ax2.plot(budget_pts, mtt_pts, "o-", color="#1e40af", lw=2.2, ms=6, label=r"Conditional UER cohort frontier $W(S^*)$")
    ax2.axvspan(800, max_budget, color="#f1f5f9", alpha=0.8, label="Higher-budget region\n(model conditional)")
    ax2.annotate(
        f"$800 budget:\n{opt_portfolio} (${opt_cost:,.2f}, {opt_val:.1f} yr)",
        xy=(800, opt_val), xytext=(850, 30),
        fontsize=8.2, fontweight="bold",
        arrowprops=dict(arrowstyle="->", color="#b91c1c", lw=1.2),
        bbox=dict(boxstyle="round,pad=0.25", fc="#fef2f2", ec="#b91c1c", alpha=0.9)
    )
    ax2.set_xlabel("Tracer Investment Budget (USD)", fontsize=10.0, fontweight="bold")
    ax2.set_ylabel("Worst-Case MTT Ambiguity (Years)", fontsize=10.0, fontweight="bold")
    ax2.set_title("(b) The Redundancy Plateau", fontsize=11.0, fontweight="bold", loc="left", pad=8)
    ax2.set_xlim(-0.05 * max_budget, max_budget * 1.02)
    ax2.set_ylim(0, 75)
    ax2.grid(True, ls="--", alpha=0.35)
    ax2.legend(loc="upper right", fontsize=8.0)

    # Subplot C: Sheaf Residual Conflict Localisation & Gating
    residuals_ungated = df_audit["res_ungated"].values
    conf_mask = df_audit["is_conflict"].values
    
    ax3.hist(residuals_ungated[~conf_mask], bins=30, color="#10b981", alpha=0.7, label=f"Coherent Edges (n={(~conf_mask).sum()})", edgecolor="black", lw=0.4)
    ax3.hist(residuals_ungated[conf_mask], bins=30, color="#ef4444", alpha=0.8, label=f"Gated Conflicts (n={conf_mask.sum()})", edgecolor="black", lw=0.4)
    ax3.axvline(1.5, color="#7f1d1d", ls="--", lw=1.8, label=r"Sheaf Gate Threshold ($r_e=1.5$)")
    ax3.set_xlabel(r"Sheaf Coboundary Residual $r_e = \|\rho_{u \to e}(x_u) - \rho_{v \to e}(x_v)\|_2$", fontsize=9.2, fontweight="bold")
    ax3.set_ylabel("Edge Count", fontsize=10.0, fontweight="bold")
    ax3.set_title("(c) Inconsistency Trap Localisation", fontsize=11.0, fontweight="bold", loc="left", pad=8)
    ax3.grid(True, ls="--", alpha=0.35)
    ax3.legend(loc="upper right", fontsize=8.0)

    fig1_path = FIG_DIR / "figure_o5_1_tripartite_boundary.png"
    plt.savefig(fig1_path, dpi=600, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"Saved Figure O5.1 Canonical 2x2 to {fig1_path.name}")

    # Spanned companion layout
    fig_s = plt.figure(figsize=(12.0, 9.4), dpi=600)
    gs = fig_s.add_gridspec(2, 2, hspace=0.34, wspace=0.28)
    s_ax1 = fig_s.add_subplot(gs[0, 0])
    s_ax2 = fig_s.add_subplot(gs[0, 1])
    s_ax3 = fig_s.add_subplot(gs[1, :])

    bars = s_ax1.bar(stages, counts, color=colors, edgecolor="black", linewidth=0.8, width=0.55)
    for bar, count in zip(bars, counts):
        s_ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 25, f"{count:,}\n({count/counts[0]:.1%})",
                   ha="center", va="bottom", fontsize=8.2, fontweight="bold")
    s_ax1.set_ylabel("Candidate Flow Edges in Network", fontsize=10.0, fontweight="bold")
    s_ax1.set_title("(a) Synergistic Network Compression", fontsize=11.0, fontweight="bold", loc="left", pad=8)
    s_ax1.set_ylim(0, 1950)
    s_ax1.grid(True, ls="--", alpha=0.35, axis="y")

    s_ax2.plot(budget_pts, mtt_pts, "o-", color="#1e40af", lw=2.2, ms=6, label=r"Conditional UER cohort frontier $W(S^*)$")
    s_ax2.axvspan(800, max_budget, color="#f1f5f9", alpha=0.8, label="Higher-budget region\n(model conditional)")
    s_ax2.annotate(
        f"$800 budget:\n{opt_portfolio} (${opt_cost:,.2f}, {opt_val:.1f} yr)",
        xy=(800, opt_val), xytext=(850, 30),
        fontsize=8.2, fontweight="bold",
        arrowprops=dict(arrowstyle="->", color="#b91c1c", lw=1.2),
        bbox=dict(boxstyle="round,pad=0.25", fc="#fef2f2", ec="#b91c1c", alpha=0.9)
    )
    s_ax2.set_xlabel("Tracer Investment Budget (USD)", fontsize=10.0, fontweight="bold")
    s_ax2.set_ylabel("Worst-Case MTT Ambiguity (Years)", fontsize=10.0, fontweight="bold")
    s_ax2.set_title("(b) The Redundancy Plateau", fontsize=11.0, fontweight="bold", loc="left", pad=8)
    s_ax2.set_xlim(-0.05 * max_budget, max_budget * 1.02)
    s_ax2.set_ylim(0, 75)
    s_ax2.grid(True, ls="--", alpha=0.35)
    s_ax2.legend(loc="upper right", fontsize=8.0)

    s_ax3.hist(residuals_ungated[~conf_mask], bins=45, color="#10b981", alpha=0.7, label=f"Coherent Edges (n={(~conf_mask).sum()})", edgecolor="black", lw=0.4)
    s_ax3.hist(residuals_ungated[conf_mask], bins=45, color="#ef4444", alpha=0.8, label=f"Gated Conflicts (n={conf_mask.sum()})", edgecolor="black", lw=0.4)
    s_ax3.axvline(1.5, color="#7f1d1d", ls="--", lw=1.8, label=r"Sheaf Gate Threshold ($r_e=1.5$)")
    s_ax3.set_xlabel(r"Sheaf Coboundary Residual $r_e = \|\rho_{u \to e}(x_u) - \rho_{v \to e}(x_v)\|_2$", fontsize=9.8, fontweight="bold")
    s_ax3.set_ylabel("Edge Count", fontsize=10.0, fontweight="bold")
    s_ax3.set_title("(c) Inconsistency Trap Localisation", fontsize=11.0, fontweight="bold", loc="left", pad=8)
    s_ax3.grid(True, ls="--", alpha=0.35)
    s_ax3.legend(loc="upper right", fontsize=8.5)

    fig1_span_path = FIG_DIR / "figure_o5_1_tripartite_boundary_spanned.png"
    plt.savefig(fig1_span_path, dpi=600, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"Saved Figure O5.1 Spanned Companion to {fig1_span_path.name}")

    # ------------------------------------------------------------------
    # FIGURE O5.2: Empirical Minimax Pareto Frontier and Impossibility Witness
    # ------------------------------------------------------------------
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.4), constrained_layout=True)

    styles = [
        ("Sub-modern Basement Cohort", "#64748b", "s--"),
        ("Regional Median Screening Cohort", "#2563eb", "o-"),
        ("Active Modern Recharge Cohort", "#059669", "^-."),
    ]
    for label, col, fmt in styles:
        sub = df_frontiers[df_frontiers["prior_regime"] == label]
        ax1.plot(sub["budget"], sub["mtt_ambiguity_yr"], fmt, color=col, lw=2.0, ms=5, label=f"{label} (n={sub['cohort_size'].iloc[0]})")

    # True certification threshold band: 0.20 * 80 = 16 yr
    ax1.axhspan(0, 16.0, color="#dcfce7", alpha=0.4, label="Certified Quantitative Target (<= 16.0 yr)")
    max_frontier_budget = float(df_frontiers["budget"].max())
    ax1.set_xlabel("Provisional analytical budget per borehole (USD; logistics excluded)", fontsize=10.5, fontweight="bold")
    ax1.set_ylabel("Worst-Case MTT Ambiguity Width (Years)", fontsize=10.5, fontweight="bold")
    ax1.set_title("(a) Conditional UER Minimax Frontiers", fontsize=11.5, fontweight="bold", loc="left", pad=8)
    ax1.set_xlim(-0.05 * max_frontier_budget, max_frontier_budget * 1.02)
    ax1.set_ylim(0, 75)
    ax1.grid(True, ls="--", alpha=0.4)
    ax1.legend(loc="upper right", fontsize=8.0, framealpha=0.95)

    if cert_imposs.status == "IMPOSSIBILITY_WITNESS" and cert_imposs.lower_witness is not None:
        ax2.plot(age_grid, cert_imposs.lower_witness, color="#047857", lw=2.2, label=r"Lower witness $x^*$")
        ax2.plot(age_grid, cert_imposs.upper_witness, color="#be185d", lw=2.2, ls="--", label=r"Upper witness $(x')^*$")
        ax2.fill_between(age_grid, cert_imposs.lower_witness, cert_imposs.upper_witness, color="#94a3b8", alpha=0.3, label="Unresolvable Ambiguity Gap (100%)")
        ax2.axvline(25.0, color="#0f172a", ls=":", lw=1.5, label="Modern Water Cutoff (25 yr)")
        ax2.set_xlabel("Groundwater Transit Time (years)", fontsize=10.5, fontweight="bold")
        ax2.set_ylabel(r"Probability Mass Density $x_k$", fontsize=10.5, fontweight="bold")
        ax2.set_title(r"(b) Minimax LP Impossibility Witness ($^{14}\mathrm{C}$ on Modern Water)", fontsize=11.5, fontweight="bold", loc="left", pad=8)
        ax2.set_xlim(0, 80)
        ax2.set_ylim(-0.03, 1.08)
        ax2.grid(True, ls="--", alpha=0.4)
        ax2.legend(loc="upper right", fontsize=8.0, framealpha=0.95)

    fig2_path = FIG_DIR / "figure_o5_2_minimax_frontier.png"
    plt.savefig(fig2_path, dpi=600)
    plt.close()
    print(f"Saved Figure O5.2 to {fig2_path.name}")


def generate_full_report(df_audit: pd.DataFrame, df_frontiers: pd.DataFrame, df_gating_comp: pd.DataFrame, channel_denominators: dict):
    print("\n--- Generating Full Objective 5 Technical Synthesis Report ---")
    rep_path = OUT_DIR / "UER_OBJECTIVE_5_FULL_REPORT.md"
    cost_manifest = load_cost_manifest()
    try:
        out_rel = OUT_DIR.relative_to(ROOT).as_posix()
    except ValueError:
        out_rel = OUT_DIR.as_posix()
    
    n_total = len(df_audit)
    n_conflicts = int(df_audit["is_conflict"].sum())
    n_coherent = n_total - n_conflicts
    
    med_curve = df_frontiers[df_frontiers["prior_regime"] == "Regional Median Screening Cohort"]
    b0_row = med_curve[med_curve["budget"] == 0.0].iloc[0]
    opt_row = med_curve[med_curve["budget"] == 800.0].iloc[0]
    high_row = med_curve[med_curve["budget"] == 2500.0].iloc[0]

    delta_plat = opt_row["mtt_ambiguity_yr"] - high_row["mtt_ambiguity_yr"]
    opt_portfolio = str(opt_row["selected_options"])
    high_portfolio = str(high_row["selected_options"])
    cost_statuses = sorted({row["quote_status"] for row in cost_manifest["candidates"]})
    cost_status_text = ", ".join(f"`{status}`" for status in cost_statuses)
    campaign_cost_model = cost_manifest.get("campaign_cost_model", {})
    campaign_cost_status = campaign_cost_model.get("status", "missing")
    campaign_cost_component_count = len(campaign_cost_model.get("components", []))
    cost_table_lines = [
        "| Option | Laboratory | Method/matrix | Detection limit / uncertainty | Posted rate and date | Shipping/tax/field status |",
        "| :--- | :--- | :--- | :--- | :--- | :--- |",
    ]
    for row in cost_manifest["candidates"]:
        cost_table_lines.append(
            "| {option} | {lab} | {method} / {matrix} | {dl}; {unc} | {currency} {price:,.2f} "
            "= USD {usd:,.2f} ({date}) | {shipping}; {taxes} Field costs excluded. |".format(
                option=row["option_id"],
                lab=row["lab"],
                method=row["method"],
                matrix=row["matrix"],
                dl=row["detection_limit"],
                unc=row["analytical_uncertainty"],
                currency=row["original_currency"],
                price=row["original_price"],
                usd=row["usd_per_sample"],
                date=row["price_date"],
                shipping=row["shipping"],
                taxes=row["taxes"],
            )
        )
    cost_table = "\n".join(cost_table_lines)
    unpriced_lines = "\n".join(f"- {item}" for item in cost_manifest["unpriced_components"])

    content = f"""# Conditional Objective 5 (O5) UER Application: Evidence Integration and Minimax Design
## Upper East Region (UER) Data Audit and Model-Based Screening

---

### Executive Summary

This report applies the Objective 5 (Hypothesis 3) diagnostics to the harmonized Upper East Region (UER) dataset (237 sampling locations, 1,600 candidate directed edges). The field package supports data screening and conditional model calculations; it does not contain independent flow-path, age, or reaction truth.

Objective 5 investigates the conditions under which integrating multi-source hydrogeological data **improves**, **adds no value to**, or **weakens** groundwater inference. 

The current run establishes:
1. **Candidate-network screening**: Chemistry, CBE and available isotope checks are applied to 1,600 Tier-C DEM candidate edges. The admitted count is an algorithmic screening result, not a validated flow-path count.
2. **Missingness semantics**: Only 17 edges have measured tritium at both endpoints (13 coherent, 4 inversion candidates); the remaining 1,583 are **Untested** on that channel. Only edges with both nitrate isotope endpoints are tested on the nitrate channel.
3. **Gating diagnostic**: The gated re-solve changes both the edge set and the fitted section. Total energy falls from {df_gating_comp.loc[0, 'dirichlet_energy']:.2f} to {df_gating_comp.loc[1, 'dirichlet_energy']:.2f}, but this is not an apples-to-apples protection proof. On the same admitted edges, the mean residual changes from {df_gating_comp.loc[0, 'mean_residual_coherent_edges']:.4f} to {df_gating_comp.loc[1, 'mean_residual_coherent_edges']:.4f}.
4. **Conditional minimax frontier**: Three cohort-level priors are derived from **42 empirical UER tritium samples** (median ${med_curve['baseline_tu'].iloc[0]:.3f}\\text{{ TU}}$). Candidate responses and error bounds remain model assumptions at an assumed sample year of 2016. Analytical costs are loaded from the dated manifest **{cost_manifest['manifest_id']}** as **{cost_manifest.get('evidence_classification', 'UNCLASSIFIED')}** public list prices; they are not vendor quotations and exclude field/logistics costs. The frontier is not a measured UER field optimum.

---

### Core Pillars Mapped to Objective 5 Hypotheses

#### 1. Conditional Network Screening (Objective 5, Condition 1: Potentially Improves Inference)
- **Topographic proxy**: The 1,600 directed edges are generated from DEM elevation gradients. They are candidate adjacency hypotheses, not measured hydraulic flowpaths.
- **Progressive screening**:
  - The chemical residual threshold flags {channel_denominators['chem_fail']} candidate edges.
  - CBE QC flags {channel_denominators['cbe_fail']} edge pairs involving a sample with $|\\text{{CBE}}| > 10\\%$.
  - The dual-endpoint tritium screen flags {channel_denominators['tritium_fail']} inversion candidates among {channel_denominators['tritium_dual_tested']} tested edges.
  - The dual nitrate-isotope screen tests {channel_denominators['nitrate_dual_tested']} edges and flags {channel_denominators['nitrate_fail']} under the current rule; {channel_denominators['nitrate_untested']} remain untested.
  - The admitted candidate set contains **{n_coherent} edges** ({n_coherent/n_total:.1%}).
- **Evidence boundary**: This shows constraint-based screening. It does not establish that the admitted edges are physically correct or that inference improved without independent flow truth.

#### 2. Model-Defined Redundancy Plateau (Objective 5, Condition 2: Adds No Value)
- **Baseline Ambiguity**: Across the empirical UER Regional Median Screening Cohort ($^3\\text{{H}} = {b0_row['baseline_tu']:.2f} \\pm {b0_row['baseline_sigma']:.2f}\\text{{ TU}}$), initial unconstrained worst-case MTT ambiguity width is **{b0_row['mtt_ambiguity_yr']:.3f} years**.
- **At the $800 budget cap**: The costed solver selects **{opt_portfolio}** at an analytical cost of **${opt_row['total_cost']:.2f}** and gives **{opt_row['mtt_ambiguity_yr']:.3f} years** of model-defined worst-case ambiguity (a **{opt_row['ambiguity_reduction_pct']:.1f}% reduction**).
- **Higher-budget comparison**: At the $2,500 cap the solver selects **{high_portfolio}** at an analytical cost of **${high_row['total_cost']:.2f}** and gives **{high_row['mtt_ambiguity_yr']:.3f} years**. The incremental model reduction from the $800 row is **{delta_plat:.5f} years** (the sign and magnitude are model outputs, not field performance).
- **Cost evidence boundary**: The manifest contains {len(cost_manifest['candidates'])} posted laboratory rates with statuses {cost_status_text}. They are suitable for a provisional sensitivity run only. A dated vendor quotation is still required before a field budget or procurement claim.

#### 3. Inconsistency Screening (Objective 5, Condition 3: Weakens Inference)
- The UER run does not execute a naive-versus-gated predictive comparison against independent truth, so it cannot establish parameter distortion or improved regional inference.
- **Re-solve diagnostic**:
  - **Ungated solve (all {len(df_audit)} edges)**: total Dirichlet energy = {df_gating_comp.loc[0, 'dirichlet_energy']:.2f}.
  - **Admitted-edge re-solve ({n_coherent} edges)**: total Dirichlet energy = {df_gating_comp.loc[1, 'dirichlet_energy']:.2f}.
  - On the same admitted edges, mean residual changes from {df_gating_comp.loc[0, 'mean_residual_coherent_edges']:.4f} to {df_gating_comp.loc[1, 'mean_residual_coherent_edges']:.4f}; the total-energy decrease alone is not evidence of protection.
- **Screened mechanisms**: chemical residuals, CBE QC failures and tritium inversion candidates are diagnostic flags. Nitrate source labels are candidate classifications and are not proof of anthropogenic contamination without measured dual-isotope support.

#### 4. Finite-Grid LP Impossibility Witnesses
- Solved as a paired-primal linear program with SciPy HiGHS, the calculation generates **Impossibility Witnesses** within the declared model:
  - Measuring $^{{14}}\\text{{C}}$ alone cannot resolve modern drinking-water fractions ($<25\\text{{ yr}}$) below $\\delta = 5\\%$.
  - Generates extremal witness distributions ($x^*, (x')^*$) spanning an unresolvable ambiguity gap of **100 percentage points (1.0)**.

---

### Evidence Coverage & Denominators

| Evidence Channel | Tested Edges | Passed | Conflicted | Untested (Missing Data) | Data Source & Quality Tier |
| :--- | :---: | :---: | :---: | :---: | :--- |
| **Topographic Gradient** | 1,600 | 1,600 | 0 | 0 | 30m SRTM/Copernicus DEM (Tier-C candidate proxy) |
| **Major Ion Chemistry** | 1,600 | {channel_denominators['chem_pass']} | {channel_denominators['chem_fail']} | 0 | 237 complete 8-ion vectors |
| **Charge Balance Error** | 1,600 | {channel_denominators['cbe_pass']} | {channel_denominators['cbe_fail']} | 0 | 6 samples exceed 10% CBE; used as a QC gate |
| **Radioactive Tritium** | 17 | {channel_denominators['tritium_pass']} | {channel_denominators['tritium_fail']} | {channel_denominators['tritium_untested']} | 42 measured samples (38 boreholes) |
| **Nitrate Isotopes** | {channel_denominators['nitrate_dual_tested']} | {channel_denominators['nitrate_pass']} | {channel_denominators['nitrate_fail']} | {channel_denominators['nitrate_untested']} | 74 samples with both isotope values |

---

### Cost Basis and Quotation Gap

The frontier reads the following dated public-list analytical rates from `{COST_MANIFEST.relative_to(ROOT).as_posix()}`. The evidence class is **{cost_manifest.get('evidence_classification', 'UNCLASSIFIED')}**. The rates are converted to USD using the nearest prior Bank of Canada observation (2026-09-14, 1 CAD = 0.7190 USD). They are **not vendor quotations**.

{cost_table}

QA/QC, duplicate/blank policy, sample rejection rules and turnaround must be confirmed in a written quotation. The following costs remain outside the frontier because no dated, site-specific amount was available:

{unpriced_lines}

The Ghana quote target is the National Isotope Hydrology Laboratory at GAEC/NNRI. Its public service page confirms isotope-analysis capability but publishes no rates; it must be contacted for a dated quotation before the frontier is used for procurement or a field budget.

The full campaign ledger is currently **{campaign_cost_status}** and defines {campaign_cost_component_count} cost components. Once the quotation fields are populated, the solver's optional subset-cost function can charge shared visits, batch freight and customs once rather than repeating them for every tracer. The present run deliberately uses the analytical-only default because the required amounts are not yet quoted.

---

### Summary of Objective 5 Artifacts

| Category | Artifact | Path |
| :--- | :--- | :--- |
| **Figure** | Figure O5.1: The Tripartite Boundary | `{out_rel}/figures/figure_o5_1_tripartite_boundary.png` (600 DPI) |
| **Figure** | Figure O5.1: Spanned Companion | `{out_rel}/figures/figure_o5_1_tripartite_boundary_spanned.png` (600 DPI) |
| **Figure** | Figure O5.2: Minimax Pareto Frontier | `{out_rel}/figures/figure_o5_2_minimax_frontier.png` (600 DPI) |
| **Table** | Table O5.1: Minimax Pareto Frontier | `{out_rel}/tables/uer_objective5_minimax_pareto_frontier.md` |
| **Table** | Table O5.2: Tripartite Boundary Audit | `{out_rel}/tables/uer_objective5_tripartite_boundary_audit.csv` |
| **Table** | Table O5.3: Conflict Localisation Table | `{out_rel}/tables/uer_objective5_sheaf_conflict_localisation.md` |
| **Table** | Table O5.4: Gating Comparison | `{out_rel}/tables/uer_objective5_sheaf_gating_comparison.csv` |
| **Table** | Table O5.5: Well-Specific Designs | `{out_rel}/tables/uer_objective5_well_specific_minimax_designs.csv` |
| **Table** | Cost basis with source dates and quote status | `{out_rel}/tables/uer_objective5_cost_basis.csv` |
| **Manifest** | Cost provenance snapshot | `{out_rel}/cost_manifest_snapshot.json` |
| **Hash** | SHA-256 for cost provenance snapshot | `{out_rel}/cost_manifest_snapshot.sha256` |
| **Template** | Dated laboratory quotation request | `provenance/uer_objective5_quote_request_2026-09-15.md` |

---

### Scientific Conclusions
The clean UER dataset and this bounded computational audit support the following:
1. Multi-source evidence integration should be treated as conditional on compatibility, measurement coverage and the declared model.
2. The reported redundancy plateau is a finite-grid result under assumed tracer histories/error bounds and provisional public-list analytical costs; it is not a field-validated optimum or a procurement budget.
3. The screens identify candidate inconsistencies and produce an admitted candidate subgraph. Gating benefit for groundwater inference remains unverified without independent flow/age truth or held-out prediction.
"""
    with open(rep_path, "w", encoding="utf-8") as f:
        f.write(content)
    print(f"Wrote full technical report to {rep_path.name}")


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    TAB_DIR.mkdir(parents=True, exist_ok=True)

    print(f"Loading field integration dataset from {DATA_CSV.name}...")
    df = pd.read_csv(DATA_CSV)
    print(f"Loading network edges from {EDGES_CSV.name}...")
    df_edges = pd.read_csv(EDGES_CSV)

    # 1. Sheaf Assembly, Gating & Tripartite Boundary
    df_audit, synergy_counts, channel_denominators, df_gating_comp = run_sheaf_and_tripartite_boundary(df, df_edges)

    # 2. Empirical Minimax Design on UER Wells & Cohorts
    df_frontiers, cert_imposs, age_grid, df_bh_designs = run_minimax_design_on_empirical_uer(df)

    # 3. Publication Figures (600 DPI)
    generate_publication_figures(df_audit, df_frontiers, cert_imposs, age_grid, synergy_counts)

    # 4. Full Technical Report
    generate_full_report(df_audit, df_frontiers, df_gating_comp, channel_denominators)

    print("\n=======================================================")
    print("OBJECTIVE 5 (O5) BOUNDED AUDIT EXECUTION COMPLETED SUCCESSFULLY!")
    print("=======================================================")


if __name__ == "__main__":
    main()
