"""
Samples Router - Handles water sample data endpoints
"""

from fastapi import APIRouter, HTTPException, UploadFile, File
from pydantic import BaseModel, Field
from typing import Dict, List, Optional, Any
import uuid
import json

router = APIRouter()

# In-memory storage for sample datasets
sample_datasets: Dict[str, Dict] = {}


class WaterSample(BaseModel):
    """A single water sample with chemical analysis"""
    sample_id: str
    location_id: str
    date: Optional[str] = None
    
    # Major ions (mg/L)
    ca: Optional[float] = Field(None, description="Calcium")
    mg: Optional[float] = Field(None, description="Magnesium")
    na: Optional[float] = Field(None, description="Sodium")
    k: Optional[float] = Field(None, description="Potassium")
    hco3: Optional[float] = Field(None, description="Bicarbonate")
    so4: Optional[float] = Field(None, description="Sulfate")
    cl: Optional[float] = Field(None, description="Chloride")
    no3: Optional[float] = Field(None, description="Nitrate")
    
    # Isotopes
    d18o: Optional[float] = Field(None, description="δ18O (‰)")
    d2h: Optional[float] = Field(None, description="δ2H (‰)")
    d15n: Optional[float] = Field(None, description="δ15N-NO3 (‰)")
    d18o_no3: Optional[float] = Field(None, description="δ18O-NO3 (‰)")
    
    # Field parameters
    ph: Optional[float] = None
    ec: Optional[float] = Field(None, description="Electrical Conductivity (μS/cm)")
    temperature: Optional[float] = Field(None, description="Temperature (°C)")
    
    # Additional properties
    tds: Optional[float] = Field(None, description="Total Dissolved Solids (mg/L)")


class SampleDataset(BaseModel):
    """A collection of water samples"""
    name: str
    description: Optional[str] = None
    samples: List[WaterSample]


@router.post("/upload")
async def upload_samples(dataset: SampleDataset):
    """Upload a new sample dataset"""
    dataset_id = str(uuid.uuid4())[:8]
    
    sample_datasets[dataset_id] = {
        "id": dataset_id,
        "name": dataset.name,
        "description": dataset.description,
        "samples": [s.model_dump() for s in dataset.samples],
        "sample_count": len(dataset.samples),
    }
    
    return {
        "dataset_id": dataset_id,
        "name": dataset.name,
        "sample_count": len(dataset.samples),
        "message": "Dataset uploaded successfully",
    }


@router.post("/upload-file")
async def upload_samples_file(file: UploadFile = File(...)):
    """Upload samples from a JSON file"""
    try:
        content = await file.read()
        data = json.loads(content)
        
        dataset_id = str(uuid.uuid4())[:8]
        
        if isinstance(data, list):
            samples = data
            name = file.filename or "Uploaded Dataset"
        elif isinstance(data, dict):
            samples = data.get("samples", [])
            name = data.get("name", file.filename or "Uploaded Dataset")
        else:
            raise ValueError("Invalid data format")
        
        sample_datasets[dataset_id] = {
            "id": dataset_id,
            "name": name,
            "description": f"Uploaded from {file.filename}",
            "samples": samples,
            "sample_count": len(samples),
        }
        
        return {
            "dataset_id": dataset_id,
            "name": name,
            "sample_count": len(samples),
            "message": "File uploaded successfully",
        }
    except json.JSONDecodeError:
        raise HTTPException(status_code=400, detail="Invalid JSON file")
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


@router.get("/datasets")
async def list_datasets():
    """List all sample datasets"""
    return [
        {
            "id": ds["id"],
            "name": ds["name"],
            "sample_count": ds["sample_count"],
        }
        for ds in sample_datasets.values()
    ]


@router.get("/datasets/{dataset_id}")
async def get_dataset(dataset_id: str):
    """Get a specific dataset"""
    if dataset_id not in sample_datasets:
        raise HTTPException(status_code=404, detail="Dataset not found")
    
    return sample_datasets[dataset_id]


@router.get("/datasets/{dataset_id}/summary")
async def get_dataset_summary(dataset_id: str):
    """Get statistical summary of a dataset"""
    if dataset_id not in sample_datasets:
        raise HTTPException(status_code=404, detail="Dataset not found")
    
    dataset = sample_datasets[dataset_id]
    samples = dataset["samples"]
    
    # Calculate summary statistics
    numeric_fields = ["ca", "mg", "na", "k", "hco3", "so4", "cl", "no3", "ph", "ec", "d18o", "d2h"]
    summary = {}
    
    for field in numeric_fields:
        values = [s.get(field) for s in samples if s.get(field) is not None]
        if values:
            summary[field] = {
                "count": len(values),
                "min": min(values),
                "max": max(values),
                "mean": sum(values) / len(values),
            }
    
    return {
        "dataset_id": dataset_id,
        "name": dataset["name"],
        "total_samples": len(samples),
        "statistics": summary,
    }


@router.delete("/datasets/{dataset_id}")
async def delete_dataset(dataset_id: str):
    """Delete a dataset"""
    if dataset_id not in sample_datasets:
        raise HTTPException(status_code=404, detail="Dataset not found")
    
    del sample_datasets[dataset_id]
    return {"message": "Dataset deleted successfully"}


# Preload demo data
demo_samples = [
    {
        "sample_id": "S001",
        "location_id": "Well_A",
        "date": "2024-01-15",
        "ca": 85.2,
        "mg": 32.1,
        "na": 45.6,
        "k": 5.2,
        "hco3": 245.0,
        "so4": 78.3,
        "cl": 52.1,
        "no3": 12.5,
        "d18o": -5.2,
        "d2h": -35.1,
        "ph": 7.4,
        "ec": 620,
    },
    {
        "sample_id": "S002",
        "location_id": "Well_B",
        "date": "2024-01-15",
        "ca": 92.1,
        "mg": 28.5,
        "na": 52.3,
        "k": 4.8,
        "hco3": 268.0,
        "so4": 85.2,
        "cl": 61.8,
        "no3": 18.2,
        "d18o": -4.8,
        "d2h": -32.5,
        "ph": 7.2,
        "ec": 685,
    },
    {
        "sample_id": "S003",
        "location_id": "Well_C",
        "date": "2024-01-16",
        "ca": 78.5,
        "mg": 35.8,
        "na": 41.2,
        "k": 6.1,
        "hco3": 232.0,
        "so4": 72.1,
        "cl": 48.5,
        "no3": 8.7,
        "d18o": -5.5,
        "d2h": -37.2,
        "ph": 7.5,
        "ec": 598,
    },
]

sample_datasets["demo"] = {
    "id": "demo",
    "name": "Demo Groundwater Samples",
    "description": "Example dataset for demonstration purposes",
    "samples": demo_samples,
    "sample_count": len(demo_samples),
}
