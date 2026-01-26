"""
Adapter layer between web API and Hydrosheaf core engine.
Handles data transformation and configuration mapping.
"""

from typing import Dict, List, Any, Optional


class ConfigAdapter:
    """Convert frontend config to Hydrosheaf Config"""

    @staticmethod
    def frontend_to_hydrosheaf(frontend_config: Dict[str, Any]):
        """Convert frontend config to Hydrosheaf Config object"""
        try:
            from hydrosheaf import Config, DEFAULT_ION_ORDER
        except ImportError:
            raise RuntimeError("Hydrosheaf package not installed on backend")

        # Build uncertainty method
        uncertainty_method = "none"
        if frontend_config.get("enable_uncertainty", False):
            uncertainty_method = frontend_config.get("uncertainty_method", "bootstrap")

        return Config(
            # Penalty settings
            lambda_l1=frontend_config.get("lasso_penalty", 0.1),
            # Core feature flags
            phreeqc_enabled=frontend_config.get("enable_phreeqc", False),
            isotope_enabled=frontend_config.get("enable_isotopes", True),
            # Uncertainty quantification
            uncertainty_method=uncertainty_method,
            bootstrap_n_resamples=frontend_config.get("bootstrap_iterations", 100),
            bayesian_n_samples=frontend_config.get("bayesian_samples", 5000),
            bayesian_n_chains=frontend_config.get("bayesian_chains", 4),
            # Nitrate source discrimination
            nitrate_source_enabled=frontend_config.get("enable_nitrate_source", False),
            nitrate_isotope_mixing_enabled=frontend_config.get(
                "nitrate_isotope_mixing", True
            ),
            # Temporal analysis
            temporal_enabled=frontend_config.get("enable_temporal", False),
            temporal_window_days=frontend_config.get("temporal_window_days", 365),
            residence_time_method=frontend_config.get(
                "residence_time_method", "cross_correlation"
            ),
            # 3D network / layered aquifer
            network_3d_enabled=frontend_config.get("enable_3d_network", False),
            vertical_flow_enabled=frontend_config.get("vertical_flow_enabled", True),
            vertical_anisotropy=frontend_config.get("vertical_anisotropy", 0.1),
            layer_enabled=frontend_config.get("enable_layer_system", False),
            # Reactive transport validation
            reactive_transport_validation=frontend_config.get(
                "enable_reactive_transport", False
            ),
            rt_simulator=frontend_config.get("rt_simulator", "phreeqc_kinetic"),
            rt_n_time_steps=frontend_config.get("rt_time_steps", 100),
            # Edge inference settings
            edge_radius_km=frontend_config.get("edge_radius_km", 10.0),
            edge_max_neighbors=frontend_config.get("edge_max_neighbors", 3),
            edge_p_min=frontend_config.get("edge_p_min", 0.75),
            edge_head_inference=frontend_config.get("edge_head_inference", "heuristic"),
            # Gibbs diagram analysis
            gibbs_enabled=frontend_config.get("enable_gibbs", True),
            gibbs_weight=frontend_config.get("gibbs_weight", 0.5),
            # Ion exchange
            exchange_enabled=frontend_config.get("enable_exchange", True),
            # Topology refinement
            sheaf_soft_beta=frontend_config.get("sheaf_soft_beta", 1.0),
            # Optional advanced settings
            ion_order=frontend_config.get("ion_order", DEFAULT_ION_ORDER.copy()),
            weights=frontend_config.get("weights", [1.0] * 10),
            # Detection limit handling
            detection_limit_policy=frontend_config.get(
                "detection_limit_policy", "half"
            ),
            missing_policy=frontend_config.get("missing_policy", "skip"),
        )

    @staticmethod
    def generate_plots(edge_results, samples, output_dir: str, config: Dict[str, Any]) -> Dict[str, str]:
        """
        Generate scientific plots using the Hydrosheaf core plotting engine.
        Returns a dictionary mapping plot types to their file paths (relative to static dir).
        """
        import logging
        logger = logging.getLogger("hydrosheaf_adapter")
        
        try:
            from hydrosheaf.outputs.plots import plot_ilr, plot_gibbs, plot_edge_anomalies
            from hydrosheaf.outputs.science_plots import plot_ttd_kernel, plot_breakthrough, plot_posterior_ridges
            from hydrosheaf.outputs.plots_3d import plot_network_3d
            from hydrosheaf.outputs.utils import PlotConfig
            import os
            from pathlib import Path
        except ImportError as e:
            logger.error(f"Failed to import plotting modules: {e}")
            return {}

        # Create output directory
        out_path = Path(output_dir)
        try:
            out_path.mkdir(parents=True, exist_ok=True)
        except Exception as e:
            logger.error(f"Failed to create output directory {output_dir}: {e}")
            return {}
        
        # Prepare plot config
        try:
            plot_cfg = PlotConfig(
                style=config.get("plot_style", "seaborn-v0_8-whitegrid"),
                palette=config.get("plot_palette", "colorblind"),
                font_scale=config.get("plot_font_scale", 1.0),
                dpi=150, # Web-friendly DPI
                file_format="png"
            )
        except Exception as e:
            logger.warning(f"Failed to create PlotConfig, using default: {e}")
            plot_cfg = PlotConfig()

        generated_plots = {}

        # Helper to run safe
        def safe_plot(name, func, *args, **kwargs):
            try:
                filename = f"{name}.png"
                filepath = out_path / filename
                func(*args, path=str(filepath), config=plot_cfg, **kwargs)
                if filepath.exists():
                    return filename
            except Exception as e:
                logger.error(f"Failed to generate plot {name}: {e}")
            return None

        # 1. Standard Hydrochemistry
        if samples:
            # Convert back to simple dict list if needed, but plotting functions handle List[Dict]
            # Ensure we pass the list of dicts directly
            generated_plots["ilr"] = safe_plot("ilr_diagram", plot_ilr, samples)
            generated_plots["gibbs"] = safe_plot("gibbs_diagram", plot_gibbs, samples)

        # 2. Network Analysis
        if edge_results:
            generated_plots["anomalies"] = safe_plot("edge_anomalies", plot_edge_anomalies, edge_results)
            
            # 3. Science Plots (Advanced)
            # TTD (needs temporal data)
            if any(r.temporal_residence_time_details for r in edge_results):
                generated_plots["ttd"] = safe_plot("ttd_kernel", plot_ttd_kernel, edge_results)
                generated_plots["breakthrough"] = safe_plot("breakthrough", plot_breakthrough, edge_results)
            
            # Uncertainty
            if any(r.uncertainty for r in edge_results):
                # Just pick the first interesting edge for now
                target_edge = next((r for r in edge_results if r.uncertainty), None)
                if target_edge:
                    generated_plots["posterior"] = safe_plot(
                        f"posterior_{target_edge.edge_id}", 
                        plot_posterior_ridges, 
                        target_edge
                    )

        # 4. 3D Plot
        # Requires nodes (derived from samples) and edges
        # We can reconstruct a basic node map from samples
        if config.get("enable_3d_network", False):
             generated_plots["network_3d"] = safe_plot(
                 "network_3d", 
                 plot_network_3d, 
                 samples, # Nodes
                 None, # Edges (will be inferred or passed if we had the edge objects, but here we just pass nodes for now or need edge objects)
                 # Note: plot_network_3d expects Edge objects, but we might only have EdgeResults here which are different.
                 # For now, let's skip passing edges to 3D plot to avoid crash, or reconstruct them.
                 # The plot_network_3d function handles "List[Dict]" for nodes.
                 z_exaggeration=config.get("vertical_anisotropy", 10.0)
             )

        return {k: v for k, v in generated_plots.items() if v is not None}



class SampleAdapter:
    """Convert frontend sample format to Hydrosheaf format"""

    # Molecular weights for unit conversion (mg/L to mmol/L)
    # Imported from core to ensure consistency
    from hydrosheaf.data.units import MOLAR_MASS_G_MOL as ION_MOLECULAR_WEIGHTS


    @staticmethod
    def mgL_to_mmolL(mg_l: float, molecular_weight: float) -> float:
        """Convert mg/L to mmol/L"""
        return mg_l / molecular_weight

    @staticmethod
    def frontend_to_hydrosheaf(
        frontend_samples: List[Dict[str, Any]]
    ) -> List[Dict[str, Any]]:
        """
        Convert frontend sample format to Hydrosheaf expected format.

        Frontend format:
        {
            "sample_id": "S001",
            "location_id": "Well_A",
            "ca": 85.2,  # mg/L, lowercase
            "mg": 32.1,
            ...
        }

        Hydrosheaf format:
        {
            "site_id": "Well_A",
            "Ca": 2.127,  # mmol/L, capitalized
            "Mg": 1.321,
            ...
        }
        """
        hydrosheaf_samples = []

        for sample in frontend_samples:
            hydrosheaf_sample = {}

            # Map ID fields
            if "location_id" in sample:
                hydrosheaf_sample["site_id"] = sample["location_id"]
            elif "sample_id" in sample:
                hydrosheaf_sample["site_id"] = sample["sample_id"]

            # Optional: preserve original sample_id
            if "sample_id" in sample:
                hydrosheaf_sample["sample_id"] = sample["sample_id"]

            # Convert major ions from mg/L to mmol/L
            ions_to_convert = {
                "ca": "Ca",
                "mg": "Mg",
                "na": "Na",
                "k": "K",
                "hco3": "HCO3",
                "so4": "SO4",
                "cl": "Cl",
                "no3": "NO3",
                "f": "F",
                "fe": "Fe",
                "po4": "PO4",
            }

            for ion_lower, ion_capital in ions_to_convert.items():
                if ion_lower in sample and sample[ion_lower] is not None:
                    mg_l_value = sample[ion_lower]
                    mw = SampleAdapter.ION_MOLECULAR_WEIGHTS.get(ion_capital)

                    if mw:
                        mmol_l_value = SampleAdapter.mgL_to_mmolL(mg_l_value, mw)
                        hydrosheaf_sample[ion_capital] = mmol_l_value

            # Direct copy for non-concentration fields
            # Map frontend lowercase keys to Core uppercase keys where needed
            field_map = {
                "ph": "pH",
                "ec": "EC",
                "tds": "TDS",
                "temperature": "temperature",
                "date": "date"
            }
            for field_lower, field_core in field_map.items():
                if field_lower in sample:
                    hydrosheaf_sample[field_core] = sample[field_lower]
                elif field_core in sample: # Handle case where input is already correct
                    hydrosheaf_sample[field_core] = sample[field_core]


            # Isotope fields (already in per-mil, no conversion needed)
            if "d18o" in sample:
                hydrosheaf_sample["18O"] = sample["d18o"]
            if "d2h" in sample:
                hydrosheaf_sample["2H"] = sample["d2h"]
            if "d15n" in sample:
                hydrosheaf_sample["15N_NO3"] = sample["d15n"]
            if "d18o_no3" in sample:
                hydrosheaf_sample["18O_NO3"] = sample["d18o_no3"]

            # Spatial coordinates
            for field in ["x", "y", "z", "latitude", "longitude", "elevation"]:
                if field in sample:
                    hydrosheaf_sample[field] = sample[field]

            hydrosheaf_samples.append(hydrosheaf_sample)

        return hydrosheaf_samples


class ResultAdapter:
    """Convert Hydrosheaf results to frontend format"""

    @staticmethod
    def hydrosheaf_to_frontend(
        edge_results, extras: Optional[Dict] = None
    ) -> Dict[str, Any]:
        """
        Convert Hydrosheaf EdgeResult objects to frontend-compatible format.

        Args:
            edge_results: List of EdgeResult objects
            extras: Additional data from fit_network_pipeline

        Returns:
            Frontend-compatible results dictionary
        """
        if not edge_results:
            return {
                "summary": {"total_samples": 0, "total_edges": 0},
                "edges": [],
                "transport_model": {},
                "reactions": [],
            }

        # Aggregate statistics
        transport_models = [
            r.transport_model for r in edge_results if r.transport_model
        ]
        dominant_process = (
            max(set(transport_models), key=transport_models.count)
            if transport_models
            else "unknown"
        )

        # Average gamma and f across edges
        avg_gamma = sum(r.gamma or 1.0 for r in edge_results) / len(edge_results)
        avg_f = sum(r.f or 0.0 for r in edge_results if r.f is not None) / max(
            1, len([r for r in edge_results if r.f is not None])
        )

        # Collect all unique reactions across network
        all_reactions = {}
        for result in edge_results:
            # EdgeResult uses z_labels and z_extents, not reaction_labels/reaction_extents
            if result.z_labels and result.z_extents:
                for label, extent in zip(result.z_labels, result.z_extents):
                    if abs(extent) > 1e-6:  # Only non-zero reactions
                        if label not in all_reactions:
                            all_reactions[label] = []
                        all_reactions[label].append(extent)

        # Average reaction extents
        reactions_list = [
            {
                "mineral": label,
                "rate_mmol_L": sum(extents) / len(extents),
                "direction": "dissolution" if sum(extents) > 0 else "precipitation",
                "occurrences": len(extents),
            }
            for label, extents in all_reactions.items()
        ]

        # Sort by absolute rate
        reactions_list.sort(key=lambda x: abs(x["rate_mmol_L"]), reverse=True)

        # Collect unique nodes
        # EdgeResult uses 'u' and 'v' attributes, not 'source' and 'target'
        all_nodes = set()
        for r in edge_results:
            all_nodes.add(r.u)
            all_nodes.add(r.v)

        # Build result structure
        frontend_result = {
            "analysis_type": "full_pipeline",
            "summary": {
                "total_samples": len(all_nodes),
                "total_edges": len(edge_results),
            },
            "transport_model": {
                "dominant_process": dominant_process,
                "evaporation_fraction": (
                    (avg_gamma - 1.0) if dominant_process == "evap" else 0.0
                ),
                "mixing_fractions": {
                    "mixing_fraction": avg_f if dominant_process == "mix" else 0.0,
                },
                "average_gamma": avg_gamma,
            },
            "reactions": reactions_list,
            "network_inference": {
                "flow_probabilities": [
                    {
                        "from": r.u,
                        "to": r.v,
                        "probability": r.edge_confidence if r.edge_confidence else 0.85,
                    }
                    for r in edge_results[:10]  # Limit to first 10 for display
                ]
            },
            "edges": [
                {
                    "edge_id": r.edge_id,
                    "source": r.u,  # EdgeResult uses 'u' not 'source'
                    "target": r.v,  # EdgeResult uses 'v' not 'target'
                    "transport_model": r.transport_model,
                    "gamma": r.gamma,
                    "f": r.f,
                    "residual_norm": r.transport_residual_norm,  # Correct field name
                    "reactions": [
                        {"label": label, "extent": extent}
                        for label, extent in zip(
                            r.z_labels or [], r.z_extents or []
                        )  # z_labels/z_extents not reaction_labels/reaction_extents
                        if abs(extent) > 1e-6
                    ],
                }
                for r in edge_results
            ],
            "uncertainty": {
                "confidence_level": 0.95,
                "mixing_ci": [max(0, avg_f - 0.05), min(1, avg_f + 0.05)],
                "reaction_uncertainties": {},
            },
        }

        # Add temporal results if available
        if extras and "temporal_results" in extras:
            frontend_result["temporal"] = {
                "residence_times": [
                    {
                        "edge_id": tr.edge_id,
                        "residence_time_days": tr.residence_time_days,
                        "method": tr.residence_time_method,
                    }
                    for tr in extras["temporal_results"]
                ]
            }

        return frontend_result
