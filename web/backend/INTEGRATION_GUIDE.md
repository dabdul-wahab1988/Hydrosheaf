# Backend Integration Guide

This guide shows how to integrate the actual Hydrosheaf core engine into the web backend.

## Step 1: Install Hydrosheaf in Backend Environment

```bash
cd web/backend

# Add to requirements.txt:
# hydrosheaf (from parent directory)

# Or install directly:
pip install ../../  # Install from repo root
```

## Step 2: Create Helper Module for Data Transformation

Create `web/backend/app/hydrosheaf_adapter.py`:

```python
"""
Adapter layer between web API and Hydrosheaf core engine.
Handles data transformation and configuration mapping.
"""

from typing import Dict, List, Any, Optional
from hydrosheaf import Config, DEFAULT_ION_ORDER
from hydrosheaf.data.units import mgL_to_mmolL


class ConfigAdapter:
    """Convert frontend config to Hydrosheaf Config"""

    ION_MOLECULAR_WEIGHTS = {
        'Ca': 40.078,
        'Mg': 24.305,
        'Na': 22.990,
        'K': 39.098,
        'HCO3': 61.017,
        'SO4': 96.064,
        'Cl': 35.453,
        'NO3': 62.005,
        'F': 18.998,
        'Fe': 55.845,
        'PO4': 94.971,
    }

    @staticmethod
    def frontend_to_hydrosheaf(frontend_config: Dict[str, Any]) -> Config:
        """Convert frontend config to Hydrosheaf Config object"""
        return Config(
            # Penalty settings
            lambda_l1=frontend_config.get('lasso_penalty', 0.1),

            # Feature flags
            phreeqc_enabled=frontend_config.get('enable_phreeqc', False),
            isotope_enabled=frontend_config.get('enable_isotopes', True),

            # Optional advanced settings (if frontend adds them later)
            ion_order=frontend_config.get('ion_order', DEFAULT_ION_ORDER.copy()),
            weights=frontend_config.get('weights', [1.0] * 10),
            active_minerals=frontend_config.get('active_minerals', None),

            # Detection limit handling
            detection_limit_policy=frontend_config.get('detection_limit_policy', 'half'),
            missing_policy=frontend_config.get('missing_policy', 'skip'),
        )


class SampleAdapter:
    """Convert frontend sample format to Hydrosheaf format"""

    @staticmethod
    def frontend_to_hydrosheaf(frontend_samples: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
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
            if 'location_id' in sample:
                hydrosheaf_sample['site_id'] = sample['location_id']
            elif 'sample_id' in sample:
                hydrosheaf_sample['site_id'] = sample['sample_id']

            # Optional: preserve original sample_id
            if 'sample_id' in sample:
                hydrosheaf_sample['sample_id'] = sample['sample_id']

            # Convert major ions from mg/L to mmol/L
            ions_to_convert = ['ca', 'mg', 'na', 'k', 'hco3', 'so4', 'cl', 'no3', 'f', 'fe', 'po4']

            for ion_lower in ions_to_convert:
                if ion_lower in sample and sample[ion_lower] is not None:
                    ion_capital = ion_lower.capitalize() if ion_lower != 'hco3' else 'HCO3'
                    if ion_lower == 'po4':
                        ion_capital = 'PO4'

                    mg_l_value = sample[ion_lower]
                    mw = ConfigAdapter.ION_MOLECULAR_WEIGHTS.get(ion_capital)

                    if mw:
                        mmol_l_value = mgL_to_mmolL(mg_l_value, mw)
                        hydrosheaf_sample[ion_capital] = mmol_l_value

            # Direct copy for non-concentration fields
            direct_copy_fields = ['ph', 'temperature', 'ec', 'tds', 'date']
            for field in direct_copy_fields:
                if field in sample:
                    hydrosheaf_sample[field] = sample[field]

            # Isotope fields (already in per-mil, no conversion needed)
            if 'd18o' in sample:
                hydrosheaf_sample['18O'] = sample['d18o']
            if 'd2h' in sample:
                hydrosheaf_sample['2H'] = sample['d2h']
            if 'd15n' in sample:
                hydrosheaf_sample['15N_NO3'] = sample['d15n']
            if 'd18o_no3' in sample:
                hydrosheaf_sample['18O_NO3'] = sample['d18o_no3']

            # Spatial coordinates
            for field in ['x', 'y', 'z', 'latitude', 'longitude', 'elevation']:
                if field in sample:
                    hydrosheaf_sample[field] = sample[field]

            hydrosheaf_samples.append(hydrosheaf_sample)

        return hydrosheaf_samples


class ResultAdapter:
    """Convert Hydrosheaf results to frontend format"""

    @staticmethod
    def hydrosheaf_to_frontend(edge_results, extras: Optional[Dict] = None) -> Dict[str, Any]:
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
                'summary': {'total_samples': 0, 'total_edges': 0},
                'edges': [],
                'transport_model': {},
                'reactions': [],
            }

        # Aggregate statistics
        transport_models = [r.transport_model for r in edge_results]
        dominant_process = max(set(transport_models), key=transport_models.count)

        # Average gamma and f across edges
        avg_gamma = sum(r.gamma or 1.0 for r in edge_results) / len(edge_results)
        avg_f = sum(r.f or 0.0 for r in edge_results) / len(edge_results)

        # Collect all unique reactions across network
        all_reactions = {}
        for result in edge_results:
            if result.reaction_labels and result.reaction_extents:
                for label, extent in zip(result.reaction_labels, result.reaction_extents):
                    if abs(extent) > 1e-6:  # Only non-zero reactions
                        if label not in all_reactions:
                            all_reactions[label] = []
                        all_reactions[label].append(extent)

        # Average reaction extents
        reactions_list = [
            {
                'mineral': label,
                'rate_mmol_L': sum(extents) / len(extents),
                'direction': 'dissolution' if sum(extents) > 0 else 'precipitation',
                'occurrences': len(extents),
            }
            for label, extents in all_reactions.items()
        ]

        # Sort by absolute rate
        reactions_list.sort(key=lambda x: abs(x['rate_mmol_L']), reverse=True)

        # Build result structure
        frontend_result = {
            'analysis_type': 'full_pipeline',
            'summary': {
                'total_samples': len(set(r.source for r in edge_results) | set(r.target for r in edge_results)),
                'total_edges': len(edge_results),
            },
            'transport_model': {
                'dominant_process': dominant_process,
                'evaporation_fraction': (avg_gamma - 1.0) if dominant_process == 'evap' else 0.0,
                'mixing_fraction': avg_f if dominant_process == 'mix' else 0.0,
                'average_gamma': avg_gamma,
            },
            'reactions': reactions_list,
            'edges': [
                {
                    'edge_id': r.edge_id,
                    'source': r.source,
                    'target': r.target,
                    'transport_model': r.transport_model,
                    'gamma': r.gamma,
                    'f': r.f,
                    'residual_norm': r.residual_norm,
                    'reactions': [
                        {'label': label, 'extent': extent}
                        for label, extent in zip(r.reaction_labels or [], r.reaction_extents or [])
                        if abs(extent) > 1e-6
                    ]
                }
                for r in edge_results
            ]
        }

        # Add temporal results if available
        if extras and 'temporal_results' in extras:
            frontend_result['temporal'] = {
                'residence_times': [
                    {
                        'edge_id': tr.edge_id,
                        'residence_time_days': tr.residence_time_days,
                        'method': tr.residence_time_method,
                    }
                    for tr in extras['temporal_results']
                ]
            }

        return frontend_result
```

## Step 3: Update Analysis Router

Replace `web/backend/app/routers/analysis.py` with integrated version:

```python
"""
Analysis Router - Handles hydrogeochemical analysis endpoints
"""

from fastapi import APIRouter, HTTPException, BackgroundTasks
from pydantic import BaseModel, Field
from typing import Dict, List, Optional, Any
from enum import Enum
import uuid
from datetime import datetime
import traceback

# Import Hydrosheaf
try:
    from hydrosheaf import (
        Config,
        fit_network_pipeline,
        auto_disable_missing_modules,
        infer_edges_probabilistic,
    )
    from hydrosheaf.graph.types import Edge
    HYDROSHEAF_AVAILABLE = True
except ImportError:
    HYDROSHEAF_AVAILABLE = False
    print("WARNING: Hydrosheaf not installed. Using mock mode.")

from ..hydrosheaf_adapter import ConfigAdapter, SampleAdapter, ResultAdapter

router = APIRouter()

# In-memory storage for analysis jobs
analysis_jobs: Dict[str, Dict] = {}


class AnalysisType(str, Enum):
    TRANSPORT = "transport"
    REACTION_PATH = "reaction_path"
    NETWORK_INFERENCE = "network_inference"
    NITRATE_SOURCE = "nitrate_source"
    FULL_PIPELINE = "full_pipeline"


class AnalysisConfig(BaseModel):
    """Configuration for analysis"""
    lasso_penalty: float = Field(default=0.1, ge=0, le=10)
    enable_phreeqc: bool = False
    enable_isotopes: bool = True
    enable_uncertainty: bool = False  # Resource-intensive, off by default
    bootstrap_iterations: int = Field(default=100, ge=10, le=1000)


class AnalysisRequest(BaseModel):
    """Request model for starting an analysis"""
    name: str = Field(..., min_length=1, max_length=100)
    analysis_type: AnalysisType
    samples: List[Dict[str, Any]]
    edges: Optional[List[Dict[str, Any]]] = None
    config: Optional[AnalysisConfig] = None


class AnalysisStatus(str, Enum):
    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"


def run_analysis_task(job_id: str, request: AnalysisRequest):
    """Background task to run actual Hydrosheaf analysis"""
    try:
        analysis_jobs[job_id]["status"] = AnalysisStatus.RUNNING
        analysis_jobs[job_id]["started_at"] = datetime.utcnow().isoformat()

        if not HYDROSHEAF_AVAILABLE:
            raise RuntimeError("Hydrosheaf engine not installed on backend")

        # Convert frontend config to Hydrosheaf Config
        frontend_config = request.config.model_dump() if request.config else {}
        config = ConfigAdapter.frontend_to_hydrosheaf(frontend_config)

        # Convert frontend samples to Hydrosheaf format
        hydrosheaf_samples = SampleAdapter.frontend_to_hydrosheaf(request.samples)

        # Auto-disable features if required data is missing
        config = auto_disable_missing_modules(hydrosheaf_samples, config)

        # Build or infer edges
        if request.edges and len(request.edges) > 0:
            # Use provided edges
            edges = [
                Edge(
                    u=e.get('source', e.get('u')),
                    v=e.get('target', e.get('v')),
                    weight=e.get('weight', 1.0)
                )
                for e in request.edges
            ]
        else:
            # Infer edges from spatial data
            edges = infer_edges_probabilistic(
                hydrosheaf_samples,
                config,
                radius_km=config.edge_radius_km,
                max_neighbors=config.edge_max_neighbors,
            )

        # Run the analysis pipeline
        edge_results, extras = fit_network_pipeline(
            hydrosheaf_samples,
            edges,
            config,
            auto_disable_missing=True,
        )

        # Convert results to frontend format
        frontend_results = ResultAdapter.hydrosheaf_to_frontend(edge_results, extras)

        # Store results
        analysis_jobs[job_id]["status"] = AnalysisStatus.COMPLETED
        analysis_jobs[job_id]["completed_at"] = datetime.utcnow().isoformat()
        analysis_jobs[job_id]["results"] = frontend_results

    except Exception as e:
        analysis_jobs[job_id]["status"] = AnalysisStatus.FAILED
        analysis_jobs[job_id]["error"] = str(e)
        analysis_jobs[job_id]["traceback"] = traceback.format_exc()
        print(f"Analysis failed for job {job_id}: {e}")
        print(traceback.format_exc())


@router.post("/run", response_model=Dict[str, Any])
async def start_analysis(request: AnalysisRequest, background_tasks: BackgroundTasks):
    """Start a new analysis job using Hydrosheaf engine"""
    job_id = str(uuid.uuid4())

    analysis_jobs[job_id] = {
        "job_id": job_id,
        "name": request.name,
        "analysis_type": request.analysis_type,
        "status": AnalysisStatus.PENDING,
        "created_at": datetime.utcnow().isoformat(),
        "config": request.config.model_dump() if request.config else {},
    }

    background_tasks.add_task(run_analysis_task, job_id, request)

    return {
        "job_id": job_id,
        "status": AnalysisStatus.PENDING,
        "message": "Analysis job created successfully",
        "hydrosheaf_available": HYDROSHEAF_AVAILABLE,
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
        # Only send first few lines
        job['error_preview'] = job['traceback'].split('\n')[:5]
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
    return list(analysis_jobs.values())


@router.delete("/jobs/{job_id}")
async def delete_analysis_job(job_id: str):
    """Delete an analysis job"""
    if job_id not in analysis_jobs:
        raise HTTPException(status_code=404, detail="Analysis job not found")

    del analysis_jobs[job_id]
    return {"message": "Job deleted successfully"}
```

## Step 4: Test the Integration

Create `web/backend/test_integration.py`:

```python
"""
Test script to verify Hydrosheaf integration
"""

import requests
import time
import json

API_BASE = "http://localhost:8000/api"


def test_full_workflow():
    """Test complete analysis workflow"""

    # 1. Upload sample data
    test_samples = [
        {
            "sample_id": "W1",
            "location_id": "Well_1",
            "ca": 80.0,
            "mg": 30.0,
            "na": 40.0,
            "k": 5.0,
            "hco3": 250.0,
            "so4": 70.0,
            "cl": 50.0,
            "no3": 15.0,
            "f": 1.0,
            "fe": 0.1,
            "ph": 7.2,
        },
        {
            "sample_id": "W2",
            "location_id": "Well_2",
            "ca": 95.0,
            "mg": 35.0,
            "na": 50.0,
            "k": 6.0,
            "hco3": 280.0,
            "so4": 85.0,
            "cl": 65.0,
            "no3": 22.0,
            "f": 1.2,
            "fe": 0.15,
            "ph": 7.0,
        },
    ]

    print("Uploading samples...")
    upload_res = requests.post(
        f"{API_BASE}/samples/upload",
        json={"name": "Test Dataset", "samples": test_samples}
    )
    print(f"Upload response: {upload_res.status_code}")
    dataset = upload_res.json()
    print(json.dumps(dataset, indent=2))

    # 2. Start analysis
    print("\nStarting analysis...")
    analysis_req = {
        "name": "Integration Test",
        "analysis_type": "full_pipeline",
        "samples": test_samples,
        "config": {
            "lasso_penalty": 0.1,
            "enable_phreeqc": False,
            "enable_isotopes": False,
        }
    }

    analysis_res = requests.post(
        f"{API_BASE}/analysis/run",
        json=analysis_req
    )
    print(f"Analysis start response: {analysis_res.status_code}")
    job = analysis_res.json()
    job_id = job['job_id']
    print(f"Job ID: {job_id}")
    print(f"Hydrosheaf available: {job.get('hydrosheaf_available')}")

    # 3. Poll for completion
    print("\nPolling for completion...")
    for i in range(30):  # 30 seconds max
        time.sleep(1)
        status_res = requests.get(f"{API_BASE}/analysis/status/{job_id}")
        status = status_res.json()
        print(f"  [{i+1}s] Status: {status['status']}")

        if status['status'] == 'completed':
            print("\n✓ Analysis completed!")
            break
        elif status['status'] == 'failed':
            print(f"\n✗ Analysis failed: {status.get('error')}")
            return

    # 4. Get results
    print("\nFetching results...")
    results_res = requests.get(f"{API_BASE}/analysis/results/{job_id}")
    results = results_res.json()

    print("\nResults:")
    print(json.dumps(results, indent=2))

    # Verify structure
    assert 'transport_model' in results
    assert 'reactions' in results
    assert 'edges' in results

    print("\n✓ Integration test passed!")


if __name__ == "__main__":
    test_full_workflow()
```

## Step 5: Run the Test

```bash
# Terminal 1: Start backend
cd web/backend
uvicorn app.main:app --reload

# Terminal 2: Run test
python test_integration.py
```

## Expected Improvements

After this integration:

✅ **Real Analysis**: Actual Hydrosheaf engine performs calculations
✅ **Accurate Results**: True transport modeling and reaction fitting
✅ **LASSO Sparsity**: Proper sparse reaction selection with L1 penalty
✅ **Auto-Disable**: Gracefully handles missing data (no PHREEQC, no isotopes, etc.)
✅ **Error Handling**: Proper exceptions and traceback logging

## Next Steps

1. Add more comprehensive error messages for user feedback
2. Implement progress tracking (currently just PENDING → RUNNING → COMPLETED)
3. Add support for temporal analysis, uncertainty quantification
4. Implement result caching to avoid re-running identical analyses
5. Add data validation before passing to Hydrosheaf engine

## Troubleshooting

**ImportError: No module named 'hydrosheaf'**
- Install Hydrosheaf in the backend venv: `pip install ../../`

**Results don't match CLI**
- Check data transformation in SampleAdapter
- Verify Config mapping in ConfigAdapter
- Compare edge construction (manual vs inferred)

**Analysis takes too long**
- Reduce number of samples for testing
- Disable PHREEQC if not needed
- Lower bootstrap iterations for uncertainty

**PHREEQC errors**
- Ensure PHREEQC is installed if `enable_phreeqc=True`
- Set `phreeqc_mode` in Config to match your installation
- Check database path in Config
