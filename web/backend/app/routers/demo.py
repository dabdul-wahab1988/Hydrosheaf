from fastapi import APIRouter, HTTPException
from typing import Dict, Any, List, Optional
from pathlib import Path
import pandas as pd
import numpy as np
import json

# Hydrosheaf imports
from hydrosheaf.config import Config
from hydrosheaf.data.units import mgL_to_mmolL
from hydrosheaf.vadose import (
    VadoseProfile,
    VadoseLayer,
)
from hydrosheaf.models.nitrate_isotopes import (
    IsotopeSample,
    load_isotope_endmembers,
    compute_isotope_prob,
)
from hydrosheaf.models.nitrate_isotopes_mcmc import (
    check_pymc_available,
    run_mcmc_mixing,
)
from hydrosheaf.isotopes import (
    compute_d_excess,
    evaporation_index,
)
from hydrosheaf.transport import (
    check_flopy_available,
    couple_vadose_saturated,
)

router = APIRouter()


# Helper to load synthetic data
def get_synthetic_data_path():
    # Relative to web/backend/app/routers/demo.py -> ... -> ... -> root -> hydrosheaf_synthetic_csv
    base_path = Path(__file__).resolve().parents[4]
    return base_path / "hydrosheaf_synthetic_csv"


@router.get("/objectives/1/vadose")
async def get_vadose_analysis():
    """Objective 1: Vadose Zone Nitrate Transport Dynamics"""
    try:
        data_dir = get_synthetic_data_path()
        chem_df = pd.read_csv(data_dir / "water_chem_full.csv")

        # Filter lysimeter data
        lysimeter_data = chem_df[chem_df["station_type"] == "lysimeter"].copy()

        # Seasonal stats
        lysimeter_data["date"] = pd.to_datetime(lysimeter_data["collection_date"])
        lysimeter_data["month"] = lysimeter_data["date"].dt.month
        lysimeter_data["season"] = lysimeter_data["month"].apply(
            lambda m: "Dry" if m in [11, 12, 1, 2, 3] else "Wet"
        )

        seasonal_stats = (
            lysimeter_data.groupby("season")["NO3_mg_L"]
            .agg(["mean", "std", "count"])
            .reset_index()
        )

        # Prepare time series data for frontend
        time_series = []
        for _, row in lysimeter_data.sort_values("date").iterrows():
            time_series.append(
                {
                    "date": row["collection_date"],
                    "station": row["station_code"],
                    "no3": row["NO3_mg_L"],
                    "event": row["event_code"],
                }
            )

        return {
            "title": "Vadose Zone Transport",
            "seasonal_stats": seasonal_stats.to_dict(orient="records"),
            "time_series": time_series,
            "profile_config": {
                "profile_id": "L1_profile",
                "depth_m": 2.0,
                "layers": 3,
                "root_depth_m": 0.6,
            },
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/objectives/2/sources")
async def get_source_discrimination():
    """Objective 2: Nitrate Source Discrimination"""
    try:
        data_dir = get_synthetic_data_path()
        chem_df = pd.read_csv(data_dir / "water_chem_full.csv")

        # Load endmembers (simulated loading since we might not have the file handy or want to hardcode for stability)
        sources = load_isotope_endmembers()
        # If load fails or returns empty, we can mock it, but load_isotope_endmembers usually has defaults

        source_data = []
        for s in sources:
            source_data.append(
                {
                    "name": s.name,
                    "d15N_mean": s.d15N_mean,
                    "d15N_std": s.d15N_std,
                    "d18O_mean": s.d18O_mean,
                    "d18O_std": s.d18O_std,
                }
            )

        samples = []
        for station in chem_df["station_code"].unique():
            station_data = chem_df[chem_df["station_code"] == station].iloc[0]
            d15n = station_data.get("d15N_NO3_permil")
            d18o = station_data.get("d18O_NO3_permil")

            if pd.notna(d15n) and pd.notna(d18o):
                # Analytical probabilities
                sample = IsotopeSample(d15N=d15n, d18O=d18o)
                probs = compute_isotope_prob(sample, sources)

                samples.append(
                    {
                        "station": station,
                        "d15N": d15n,
                        "d18O": d18o,
                        "probabilities": probs,
                    }
                )

        return {
            "sources": source_data,
            "samples": samples,
            "mcmc_available": check_pymc_available(),
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/objectives/3/recharge")
async def get_recharge_tracing():
    """Objective 3: Groundwater Recharge Pathways"""
    try:
        data_dir = get_synthetic_data_path()
        chem_df = pd.read_csv(data_dir / "water_chem_full.csv")

        results = []
        for _, row in chem_df.iterrows():
            d18o = row.get("d18O_H2O_permil")
            d2h = row.get("d2H_H2O_permil")

            if pd.notna(d18o) and pd.notna(d2h):
                d_excess = compute_d_excess(d2h, d18o)
                results.append(
                    {
                        "station": row["station_code"],
                        "type": row["station_type"],
                        "d18O": d18o,
                        "d2H": d2h,
                        "d_excess": d_excess,
                        "is_evaporated": d_excess < 10,
                    }
                )

        # LMWL line points for plotting
        x_range = [
            min(r["d18O"] for r in results) - 1,
            max(r["d18O"] for r in results) + 1,
        ]
        lmwl = [{"x": x, "y": 8.0 * x + 10.0} for x in x_range]  # GMWL
        local_mwl = [
            {"x": x, "y": 7.87 * x + 13.61} for x in x_range  # Local MWL (Ghana) - above GMWL
        ]

        return {"samples": results, "gmwl": lmwl, "lmwl": local_mwl}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/objectives/4/transport")
async def get_transport_modeling():
    """Objective 4: FloPy Transport Modeling"""
    try:
        if not check_flopy_available():
            return {"available": False, "message": "FloPy not available"}

        # Synthetic vadose output (simulating 18 months of data)
        n_days = 550
        vadose_times = np.arange(n_days, dtype=float)

        # Simulate seasonal recharge pattern
        recharge_base = 0.001
        recharge_amplitude = 0.002
        recharge = recharge_base + recharge_amplitude * np.sin(
            2 * np.pi * vadose_times / 365
        )
        recharge = np.maximum(recharge, 0)

        # Simulate NO3
        no3_base = 20.0
        no3_concentration = no3_base + 80 * np.exp(
            -(((vadose_times - 90) % 365) ** 2) / 500
        )
        no3_concentration += 60 * np.exp(-(((vadose_times - 270) % 365) ** 2) / 500)

        # Run coupling
        result = couple_vadose_saturated(
            vadose_times=vadose_times,
            vadose_recharge_m_day=recharge,
            vadose_concentration=no3_concentration,
            aquifer_length_m=200.0,
            aquifer_thickness_m=15.0,
            hydraulic_k_m_day=5.0,
            porosity=0.25,
            dispersivity_m=10.0,
            denitrification_k_1_day=0.005,
            head_gradient=0.01,
            use_analytical=False,
        )

        # Convert numpy arrays to lists for JSON
        output_times = (
            result.combined_times.tolist() if result.combined_times is not None else []
        )
        output_conc = (
            result.combined_concentration.tolist()
            if result.combined_concentration is not None
            else []
        )

        # Downsample for frontend if too large
        if len(output_times) > 500:
            step = len(output_times) // 500
            output_times = output_times[::step]
            output_conc = output_conc[::step]

        return {
            "available": True,
            "success": result.success,
            "input_times": vadose_times.tolist()[::5],  # Downsample
            "input_conc": no3_concentration.tolist()[::5],
            "output_times": output_times,
            "output_conc": output_conc,
            "metrics": {
                "attenuation": result.attenuation_factor,
                "peak_conc": result.peak_concentration,
                "travel_time": result.mean_travel_time_days,
            },
            "warnings": result.warnings,
        }
    except Exception as e:
        # Fallback to analytical if FloPy fails in production
        return {
            "available": True,
            "success": False,
            "error": str(e),
            "message": "Transport simulation failed",
        }
