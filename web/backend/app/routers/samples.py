"""
Samples Router - Handles water sample data endpoints
"""

from fastapi import APIRouter, HTTPException, UploadFile, File
from pydantic import BaseModel, Field
from typing import Dict, List, Optional, Any
import uuid
import json
import csv
import io

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


def parse_csv_samples(content: bytes, filename: str) -> List[Dict]:
    """Parse CSV content into sample dictionaries"""
    # Decode content
    try:
        text = content.decode('utf-8')
    except UnicodeDecodeError:
        text = content.decode('latin-1')

    # Parse CSV
    reader = csv.DictReader(io.StringIO(text))
    samples = []

    # Column name mapping (case-insensitive, common variations)
    column_mapping = {
        'sample_id': ['sample_id', 'sampleid', 'sample', 'id', 'name'],
        'location_id': ['location_id', 'locationid', 'location', 'well', 'well_id', 'site', 'site_id'],
        'date': ['date', 'sample_date', 'collection_date', 'datetime'],
        'ca': ['ca', 'calcium', 'ca2+', 'ca_mgl'],
        'mg': ['mg', 'magnesium', 'mg2+', 'mg_mgl'],
        'na': ['na', 'sodium', 'na+', 'na_mgl'],
        'k': ['k', 'potassium', 'k+', 'k_mgl'],
        'hco3': ['hco3', 'bicarbonate', 'hco3-', 'hco3_mgl', 'alkalinity'],
        'so4': ['so4', 'sulfate', 'sulphate', 'so42-', 'so4_mgl'],
        'cl': ['cl', 'chloride', 'cl-', 'cl_mgl'],
        'no3': ['no3', 'nitrate', 'no3-', 'no3_mgl'],
        'f': ['f', 'fluoride', 'f-', 'f_mgl'],
        'fe': ['fe', 'iron', 'fe2+', 'fe_mgl'],
        'po4': ['po4', 'phosphate', 'po43-', 'po4_mgl'],
        'd18o': ['d18o', 'delta18o', 'o18', 'δ18o', 'oxygen18'],
        'd2h': ['d2h', 'delta2h', 'deuterium', 'δ2h', 'h2'],
        'd15n': ['d15n', 'delta15n', 'n15', 'δ15n', 'nitrogen15'],
        'd18o_no3': ['d18o_no3', 'delta18o_no3', 'o18_no3'],
        'ph': ['ph'],
        'ec': ['ec', 'conductivity', 'electrical_conductivity', 'spc'],
        'temperature': ['temperature', 'temp', 't'],
        'tds': ['tds', 'total_dissolved_solids'],
        'x': ['x', 'longitude', 'lon', 'easting', 'x_coord'],
        'y': ['y', 'latitude', 'lat', 'northing', 'y_coord'],
    }

    # Build reverse mapping from actual column names to standard names
    def find_standard_name(col_name: str) -> Optional[str]:
        col_lower = col_name.lower().strip()
        for standard, variations in column_mapping.items():
            if col_lower in variations:
                return standard
        return None

    for row_num, row in enumerate(reader, start=1):
        sample = {}
        for col_name, value in row.items():
            if value is None or value.strip() == '':
                continue

            standard_name = find_standard_name(col_name)
            if standard_name:
                # Try to convert numeric values
                if standard_name not in ['sample_id', 'location_id', 'date']:
                    try:
                        sample[standard_name] = float(value)
                    except ValueError:
                        sample[standard_name] = value
                else:
                    sample[standard_name] = value
            else:
                # Keep unknown columns with original names (lowercase)
                try:
                    sample[col_name.lower()] = float(value)
                except ValueError:
                    sample[col_name.lower()] = value

        # Generate sample_id if missing
        if 'sample_id' not in sample:
            sample['sample_id'] = f"S{row_num:03d}"

        # Generate location_id if missing
        if 'location_id' not in sample:
            sample['location_id'] = sample['sample_id']

        samples.append(sample)

    return samples


@router.post("/upload-file")
async def upload_samples_file(file: UploadFile = File(...)):
    """Upload samples from a JSON or CSV file"""
    content = await file.read()
    filename = file.filename or "uploaded_file"
    file_ext = filename.lower().split('.')[-1] if '.' in filename else ''

    dataset_id = str(uuid.uuid4())[:8]
    samples = []
    name = filename.rsplit('.', 1)[0] if '.' in filename else filename
    file_type = ""

    # Determine file type and parse accordingly
    if file_ext == 'csv':
        # Parse CSV file
        try:
            samples = parse_csv_samples(content, filename)
            file_type = "CSV"
        except Exception as e:
            raise HTTPException(status_code=400, detail=f"CSV parsing error: {str(e)}")

    elif file_ext == 'json':
        # Parse JSON file
        try:
            data = json.loads(content)
            if isinstance(data, list):
                samples = data
            elif isinstance(data, dict):
                samples = data.get("samples", [])
                name = data.get("name", name)
            else:
                raise HTTPException(status_code=400, detail="Invalid JSON data format - expected array or object with 'samples' key")
            file_type = "JSON"
        except json.JSONDecodeError as e:
            raise HTTPException(status_code=400, detail=f"Invalid JSON format: {str(e)}")

    else:
        # Try to auto-detect format - try JSON first, then CSV
        try:
            data = json.loads(content)
            if isinstance(data, list):
                samples = data
            elif isinstance(data, dict):
                samples = data.get("samples", [])
                name = data.get("name", name)
            else:
                raise ValueError("Not valid JSON structure")
            file_type = "JSON"
        except (json.JSONDecodeError, ValueError):
            # Try CSV
            try:
                samples = parse_csv_samples(content, filename)
                file_type = "CSV"
            except Exception as e:
                raise HTTPException(status_code=400, detail=f"Could not parse file as JSON or CSV: {str(e)}")

    if not samples:
        raise HTTPException(status_code=400, detail="No samples found in the file. Please check the file format.")

    sample_datasets[dataset_id] = {
        "id": dataset_id,
        "name": name,
        "description": f"Uploaded from {filename} ({file_type})",
        "samples": samples,
        "sample_count": len(samples),
    }

    return {
        "dataset_id": dataset_id,
        "name": name,
        "sample_count": len(samples),
        "file_type": file_type,
        "message": f"{file_type} file uploaded successfully",
    }


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


@router.get("/datasets/{dataset_id}/capabilities")
async def get_dataset_capabilities(dataset_id: str):
    """
    Analyze dataset and return available analysis capabilities.
    This helps the frontend show only relevant analysis options.
    """
    if dataset_id not in sample_datasets:
        raise HTTPException(status_code=404, detail="Dataset not found")

    dataset = sample_datasets[dataset_id]
    samples = dataset["samples"]

    # Define field categories
    ion_fields = ['ca', 'mg', 'na', 'k', 'hco3', 'so4', 'cl', 'no3', 'f', 'fe', 'po4']
    isotope_fields = ['d18o', 'd2h']
    nitrate_isotope_fields = ['d15n', 'd18o_no3']
    spatial_fields = ['x', 'y', 'longitude', 'latitude', 'lon', 'lat']
    depth_fields = ['z', 'screen_depth', 'well_depth', 'elevation']
    temporal_fields = ['date', 'datetime', 'sample_date']
    field_param_fields = ['ph', 'ec', 'temperature', 'tds']

    # Count samples with each field
    def count_field(field_list):
        for field in field_list:
            count = sum(1 for s in samples if s.get(field) is not None and s.get(field) != '')
            if count > 0:
                return count, field
        return 0, None

    def get_present_fields(field_list):
        present = []
        for field in field_list:
            count = sum(1 for s in samples if s.get(field) is not None and s.get(field) != '')
            if count > 0:
                present.append(field)
        return present

    # Detect available fields
    available_ions = get_present_fields(ion_fields)
    available_isotopes = get_present_fields(isotope_fields)
    available_nitrate_isotopes = get_present_fields(nitrate_isotope_fields)
    available_spatial = get_present_fields(spatial_fields)
    available_depth = get_present_fields(depth_fields)
    available_temporal = get_present_fields(temporal_fields)
    available_field_params = get_present_fields(field_param_fields)

    # Minimum requirements for each analysis
    has_major_ions = len(set(available_ions) & {'ca', 'mg', 'na', 'hco3', 'cl', 'so4'}) >= 4
    has_coordinates = len(available_spatial) >= 2
    has_isotopes = len(available_isotopes) >= 2
    has_nitrate_isotopes = len(available_nitrate_isotopes) >= 2
    has_nitrate = 'no3' in available_ions
    has_temporal = len(available_temporal) > 0
    has_depth = len(available_depth) > 0
    has_ph = 'ph' in available_field_params

    # Determine available analyses
    available_analyses = {
        "transport": has_major_ions,
        "reaction_path": has_major_ions,
        "network_inference": has_coordinates,
        "isotope_analysis": has_isotopes,
        "phreeqc": has_major_ions and has_ph,
        "nitrate_source": has_nitrate and has_nitrate_isotopes,
        "temporal": has_temporal and has_major_ions,
        "network_3d": has_coordinates and has_depth,
        "gibbs": has_major_ions and ('tds' in available_field_params or 'ec' in available_field_params),
        "exchange": has_major_ions,
        "coda": has_major_ions,
        "vadose_zone": has_nitrate and has_depth,
        "reactive_transport": has_major_ions,
        "uncertainty": has_major_ions,
    }

    # Generate warnings for unavailable analyses
    warnings = []
    if not has_major_ions:
        warnings.append("Major ion analysis requires at least Ca, Mg, Na, HCO3, Cl, SO4")
    if not has_coordinates:
        warnings.append("Network inference requires spatial coordinates (x/y or lat/lon)")
    if not has_isotopes:
        warnings.append("Isotope analysis requires d18o and d2h fields")
    if not has_nitrate_isotopes and has_nitrate:
        warnings.append("Nitrate source discrimination requires d15n and d18o_no3 isotope fields")
    if not has_ph and has_major_ions:
        warnings.append("PHREEQC constraints require pH measurements")
    if not has_depth and has_coordinates:
        warnings.append("3D network analysis requires depth/elevation data")
    if not has_temporal:
        warnings.append("Temporal analysis requires date field in samples")

    return {
        "dataset_id": dataset_id,
        "sample_count": len(samples),
        "available_fields": {
            "ions": available_ions,
            "isotopes": available_isotopes,
            "nitrate_isotopes": available_nitrate_isotopes,
            "spatial": available_spatial,
            "depth": available_depth,
            "temporal": available_temporal,
            "field_params": available_field_params,
        },
        "available_analyses": available_analyses,
        "warnings": warnings,
        "recommendations": {
            "can_run_basic": has_major_ions,
            "can_run_spatial": has_coordinates,
            "can_run_isotope": has_isotopes,
            "can_run_full_pipeline": has_major_ions and has_coordinates,
        }
    }


@router.delete("/datasets/{dataset_id}")
async def delete_dataset(dataset_id: str):
    """Delete a dataset"""
    if dataset_id not in sample_datasets:
        raise HTTPException(status_code=404, detail="Dataset not found")

    del sample_datasets[dataset_id]
    return {"message": "Dataset deleted successfully"}


# Preload demo data with complete fields including coordinates
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
        "temperature": 22.5,
        "tds": 450,
        "x": 0.0,
        "y": 0.0,
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
        "temperature": 23.1,
        "tds": 512,
        "x": 2.5,
        "y": 1.2,
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
        "temperature": 21.8,
        "tds": 425,
        "x": 5.0,
        "y": 3.5,
    },
]

sample_datasets["demo"] = {
    "id": "demo",
    "name": "Demo Groundwater Samples",
    "description": "Example dataset for demonstration purposes",
    "samples": demo_samples,
    "sample_count": len(demo_samples),
}
