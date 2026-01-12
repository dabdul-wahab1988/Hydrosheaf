"""
Analysis Router - Handles hydrogeochemical analysis endpoints
Now integrated with real Hydrosheaf core engine!
"""

from fastapi import APIRouter, HTTPException, BackgroundTasks
from pydantic import BaseModel, Field
from typing import Dict, List, Optional, Any
from enum import Enum
import uuid
from datetime import datetime
import traceback
import sys

from .. import project_store

HYDROSHEAF_AVAILABLE = None
HYDROSHEAF_IMPORT_ERROR = None
Config = None
fit_network_pipeline = None
auto_disable_missing_modules = None
infer_edges_probabilistic = None
HydrosheafEdge = None


def _load_hydrosheaf() -> None:
    """Lazy-load Hydrosheaf to avoid heavy imports during app startup."""
    global HYDROSHEAF_AVAILABLE, HYDROSHEAF_IMPORT_ERROR
    global Config, fit_network_pipeline, auto_disable_missing_modules
    global infer_edges_probabilistic, HydrosheafEdge

    if HYDROSHEAF_AVAILABLE is not None:
        return
    try:
        from hydrosheaf import (
            Config as HydrosheafConfig,
            fit_network_pipeline as hydrosheaf_fit_network_pipeline,
            auto_disable_missing_modules as hydrosheaf_auto_disable_missing_modules,
            infer_edges_probabilistic as hydrosheaf_infer_edges_probabilistic,
        )
        from hydrosheaf.graph.types import Edge as HydrosheafEdgeType

        Config = HydrosheafConfig
        fit_network_pipeline = hydrosheaf_fit_network_pipeline
        auto_disable_missing_modules = hydrosheaf_auto_disable_missing_modules
        infer_edges_probabilistic = hydrosheaf_infer_edges_probabilistic
        HydrosheafEdge = HydrosheafEdgeType
        HYDROSHEAF_AVAILABLE = True
        HYDROSHEAF_IMPORT_ERROR = None
        try:
            with open("analysis_debug.log", "a") as f:
                f.write(f"\n[{datetime.utcnow().isoformat()}] Hydrosheaf loaded successfully.\n")
        except Exception:
            pass
    except Exception as exc:
        HYDROSHEAF_AVAILABLE = False
        HYDROSHEAF_IMPORT_ERROR = str(exc)
        HydrosheafEdge = None
        print(f"WARNING: Hydrosheaf not available: {exc}", file=sys.stderr)
        print("Install with: pip install ../../ (from web/backend directory)", file=sys.stderr)
        try:
            with open("analysis_debug.log", "a") as f:
                f.write(f"\n[{datetime.utcnow().isoformat()}] WARNING: Hydrosheaf not available: {exc}\n")
        except Exception:
            pass

# Import our adapter
try:
    from ..hydrosheaf_adapter import ConfigAdapter, SampleAdapter, ResultAdapter
    ADAPTER_AVAILABLE = True
except ImportError as e:
    ADAPTER_AVAILABLE = False
    print(f"WARNING: Adapter not available: {e}", file=sys.stderr)

router = APIRouter()

# In-memory storage for analysis jobs (in production, use a database)
analysis_jobs: Dict[str, Dict] = {}


class AnalysisType(str, Enum):
    TRANSPORT = "transport"
    REACTION_PATH = "reaction_path"
    NETWORK_INFERENCE = "network_inference"
    NITRATE_SOURCE = "nitrate_source"
    FULL_PIPELINE = "full_pipeline"


class AnalysisConfig(BaseModel):
    """Configuration for analysis"""
    # Core settings
    lasso_penalty: float = Field(default=0.1, ge=0, le=10)

    # Feature flags - Core modules
    enable_phreeqc: bool = False
    enable_isotopes: bool = True
    enable_uncertainty: bool = False  # Disabled by default (resource intensive)

    # Uncertainty settings
    uncertainty_method: str = Field(default="bootstrap", pattern="^(none|bootstrap|bayesian|monte_carlo)$")
    bootstrap_iterations: int = Field(default=100, ge=10, le=10000)
    bayesian_samples: int = Field(default=5000, ge=100, le=50000)
    bayesian_chains: int = Field(default=4, ge=1, le=8)

    # Nitrate source discrimination
    enable_nitrate_source: bool = False
    nitrate_isotope_mixing: bool = True

    # Temporal analysis
    enable_temporal: bool = False
    temporal_window_days: int = Field(default=365, ge=30, le=3650)
    residence_time_method: str = Field(default="cross_correlation", pattern="^(gradient|cross_correlation|bayesian_lag|ttd|tracer_decay)$")

    # 3D network / layered aquifer
    enable_3d_network: bool = False
    vertical_flow_enabled: bool = True
    vertical_anisotropy: float = Field(default=0.1, ge=0.001, le=1.0)
    enable_layer_system: bool = False

    # Reactive transport validation
    enable_reactive_transport: bool = False
    rt_simulator: str = Field(default="phreeqc_kinetic", pattern="^(phreeqc_kinetic|mt3dms)$")
    rt_time_steps: int = Field(default=100, ge=10, le=1000)

    # Vadose zone modeling (for nitrate)
    enable_vadose_zone: bool = False

    # CoDA (Compositional Data Analysis)
    enable_coda: bool = False

    # Edge inference settings
    edge_radius_km: float = Field(default=10.0, ge=0.1, le=100)
    edge_max_neighbors: int = Field(default=3, ge=1, le=10)
    edge_p_min: float = Field(default=0.75, ge=0.0, le=1.0)
    edge_head_inference: str = Field(default="heuristic", pattern="^(heuristic|bayesian|bayesian_mcmc)$")

    # Gibbs diagram analysis
    enable_gibbs: bool = True
    gibbs_weight: float = Field(default=0.5, ge=0, le=2.0)

    # Ion exchange
    enable_exchange: bool = True


class AnalysisRequest(BaseModel):
    """Request model for starting an analysis"""
    name: str = Field(..., min_length=1, max_length=100)
    analysis_type: AnalysisType
    samples: List[Dict[str, Any]]
    edges: Optional[List[Dict[str, Any]]] = None
    config: Optional[AnalysisConfig] = None
    project_id: Optional[str] = None
    dataset_id: Optional[str] = None
    dataset_name: Optional[str] = None
    dataset_description: Optional[str] = None


class AnalysisStatus(str, Enum):
    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"


def run_analysis_task(job_id: str, request: AnalysisRequest):
    """Background task to run REAL Hydrosheaf analysis"""
    try:
        _load_hydrosheaf()

        analysis_jobs[job_id]["status"] = AnalysisStatus.RUNNING
        analysis_jobs[job_id]["started_at"] = datetime.utcnow().isoformat()

        # Check if Hydrosheaf is available
        if not HYDROSHEAF_AVAILABLE:
            raise RuntimeError(
                f"Hydrosheaf engine not installed on backend. "
                f"Error: {HYDROSHEAF_IMPORT_ERROR if 'HYDROSHEAF_IMPORT_ERROR' in globals() else 'Unknown'}. "
                f"Install with: pip install ../../ from web/backend directory"
            )

        if not ADAPTER_AVAILABLE:
            raise RuntimeError("Hydrosheaf adapter module not available")

        print(f"[Job {job_id}] Starting Hydrosheaf analysis: {request.name}")

        # Convert frontend config to Hydrosheaf Config
        frontend_config = request.config.model_dump() if request.config else {}
        config = ConfigAdapter.frontend_to_hydrosheaf(frontend_config)
        print(f"[Job {job_id}] Config created: lambda_l1={config.lambda_l1}, phreeqc={config.phreeqc_enabled}")

        # Convert frontend samples to Hydrosheaf format
        hydrosheaf_samples = SampleAdapter.frontend_to_hydrosheaf(request.samples)
        print(f"[Job {job_id}] Converted {len(hydrosheaf_samples)} samples to Hydrosheaf format")

        # Auto-disable features if required data is missing
        config = auto_disable_missing_modules(hydrosheaf_samples, config)
        print(f"[Job {job_id}] Auto-disable: phreeqc={config.phreeqc_enabled}, isotope={config.isotope_enabled}")

        # Build or infer edges
        if request.edges and len(request.edges) > 0:
            # Use provided edges
            print(f"[Job {job_id}] Using {len(request.edges)} provided edges")
            edges = []
            for e in request.edges:
                u = e.get('source', e.get('u', ''))
                v = e.get('target', e.get('v', ''))
                weight = e.get('weight', 1.0)
                edge = HydrosheafEdge(
                    edge_id=f"{u}->{v}",
                    u=u,
                    v=v,
                    attrs={"weight": weight},
                )
                edges.append(edge)
            edge_source = "provided"
        else:
            # Infer edges from spatial data if coordinates available
            print(f"[Job {job_id}] No edges provided, attempting to infer from spatial data...")

            # Check if samples have coordinates
            has_coords = any(
                ('x' in s and 'y' in s) or ('latitude' in s and 'longitude' in s)
                for s in hydrosheaf_samples
            )

            if has_coords:
                edges = infer_edges_probabilistic(
                    hydrosheaf_samples,
                    radius_km=config.edge_radius_km,
                    max_neighbors=config.edge_max_neighbors,
                    p_min=config.edge_p_min,
                )
                print(f"[Job {job_id}] Inferred {len(edges)} edges from coordinates")
                edge_source = "inferred"
            else:
                # Create simple sequential edges if no coordinates
                print(f"[Job {job_id}] No coordinates found, creating sequential edges")
                edges = []
                for i in range(len(hydrosheaf_samples) - 1):
                    source_id = hydrosheaf_samples[i].get('site_id', f'sample_{i}')
                    target_id = hydrosheaf_samples[i + 1].get('site_id', f'sample_{i+1}')
                    edges.append(
                        HydrosheafEdge(
                            edge_id=f"{source_id}->{target_id}",
                            u=source_id,
                            v=target_id,
                            attrs={"weight": 1.0},
                        )
                    )
                edge_source = "sequential"

        if not edges:
            raise ValueError(
                "No edges to analyze. Either provide edges explicitly or include "
                "spatial coordinates (x,y or latitude,longitude) in your samples for automatic edge inference."
            )

        if request.project_id:
            edges_for_storage = []
            for edge in edges:
                attrs = edge.attrs or {}
                edges_for_storage.append(
                    {
                        "u": edge.u,
                        "v": edge.v,
                        "weight": attrs.get("weight"),
                        "edge_id": edge.edge_id,
                    }
                )
            project_store.add_input(
                request.project_id,
                {
                    "analysis_job_id": job_id,
                    "analysis_name": request.name,
                    "analysis_type": request.analysis_type.value,
                    "dataset_id": request.dataset_id,
                    "dataset_name": request.dataset_name,
                    "dataset_description": request.dataset_description,
                    "sample_count": len(request.samples),
                    "samples": request.samples,
                    "edges": edges_for_storage,
                    "edge_source": edge_source,
                    "config": frontend_config,
                },
            )

        print(f"[Job {job_id}] Running fit_network_pipeline with {len(edges)} edges...")

        # Run the REAL Hydrosheaf analysis pipeline!
        edge_results, extras = fit_network_pipeline(
            hydrosheaf_samples,
            edges,
            config,
            auto_disable_missing=True,
        )

        print(f"[Job {job_id}] Analysis complete! Got {len(edge_results)} edge results")

        # Convert results to frontend format
        frontend_results = ResultAdapter.hydrosheaf_to_frontend(edge_results, extras)

        # Add metadata
        frontend_results['metadata'] = {
            'hydrosheaf_version': '0.3.0',
            'analysis_engine': 'Hydrosheaf Core',
            'mock_data': False,
            'edges_analyzed': len(edge_results),
            'samples_analyzed': len(hydrosheaf_samples),
        }

        # Store results
        analysis_jobs[job_id]["status"] = AnalysisStatus.COMPLETED
        analysis_jobs[job_id]["completed_at"] = datetime.utcnow().isoformat()
        analysis_jobs[job_id]["results"] = frontend_results

        print(f"[Job {job_id}] SUCCESS - Results saved")

    except Exception as e:
        error_msg = str(e)
        error_trace = traceback.format_exc()

        analysis_jobs[job_id]["status"] = AnalysisStatus.FAILED
        analysis_jobs[job_id]["error"] = error_msg
        analysis_jobs[job_id]["traceback"] = error_trace

        print(f"[Job {job_id}] FAILED: {error_msg}", file=sys.stderr)
        print(error_trace, file=sys.stderr)
        try:
            with open("analysis_debug.log", "a") as f:
                f.write(f"\n[{datetime.utcnow().isoformat()}] Job {job_id} FAILED:\n")
                f.write(f"{error_msg}\n")
                f.write(f"{error_trace}\n")
        except Exception:
            pass


@router.post("/run", response_model=Dict[str, Any])
async def start_analysis(request: AnalysisRequest, background_tasks: BackgroundTasks):
    """Start a new analysis job using REAL Hydrosheaf engine"""
    if request.project_id and not project_store.get_project(request.project_id):
        raise HTTPException(status_code=404, detail="Project not found")

    _load_hydrosheaf()

    job_id = str(uuid.uuid4())

    analysis_jobs[job_id] = {
        "job_id": job_id,
        "name": request.name,
        "analysis_type": request.analysis_type,
        "status": AnalysisStatus.PENDING,
        "created_at": datetime.utcnow().isoformat(),
        "config": request.config.model_dump() if request.config else {},
        "project_id": request.project_id,
        "dataset_id": request.dataset_id,
    }

    if request.project_id:
        project_store.add_analysis_job(request.project_id, job_id)

    background_tasks.add_task(run_analysis_task, job_id, request)

    return {
        "job_id": job_id,
        "status": AnalysisStatus.PENDING,
        "message": "Analysis job created successfully (using REAL Hydrosheaf engine)",
        "hydrosheaf_available": HYDROSHEAF_AVAILABLE,
        "adapter_available": ADAPTER_AVAILABLE,
    }


@router.get("/status/{job_id}")
async def get_analysis_status(job_id: str):
    """Get the status of an analysis job"""
    if job_id not in analysis_jobs:
        raise HTTPException(status_code=404, detail="Analysis job not found")

    job = analysis_jobs[job_id].copy()

    # Don't expose full traceback to frontend unless debugging
    if 'traceback' in job:
        job['has_traceback'] = True
        # Only send first few lines of error
        lines = job['traceback'].split('\n')
        job['error_preview'] = '\n'.join(lines[:10])
        del job['traceback']

    return job


@router.get("/results/{job_id}")
async def get_analysis_results(job_id: str):
    """Get the results of a completed analysis"""
    if job_id not in analysis_jobs:
        raise HTTPException(status_code=404, detail="Analysis job not found")

    job = analysis_jobs[job_id]

    if job["status"] != AnalysisStatus.COMPLETED:
        raise HTTPException(
            status_code=400,
            detail=f"Analysis not completed. Current status: {job['status']}"
        )

    return job.get("results", {})


@router.get("/jobs")
async def list_analysis_jobs():
    """List all analysis jobs"""
    jobs = []
    for job in analysis_jobs.values():
        job_copy = job.copy()
        # Remove traceback from list view
        if 'traceback' in job_copy:
            del job_copy['traceback']
        jobs.append(job_copy)
    return jobs


@router.delete("/jobs/{job_id}")
async def delete_analysis_job(job_id: str):
    """Delete an analysis job"""
    if job_id not in analysis_jobs:
        raise HTTPException(status_code=404, detail="Analysis job not found")

    del analysis_jobs[job_id]
    return {"message": "Job deleted successfully"}


@router.get("/health")
async def health_check():
    """Check if Hydrosheaf engine is available"""
    _load_hydrosheaf()
    return {
        "hydrosheaf_available": HYDROSHEAF_AVAILABLE,
        "adapter_available": ADAPTER_AVAILABLE,
        "status": "ready" if (HYDROSHEAF_AVAILABLE and ADAPTER_AVAILABLE) else "degraded",
        "message": "Real Hydrosheaf engine integrated!" if HYDROSHEAF_AVAILABLE else "Hydrosheaf not installed",
    }
