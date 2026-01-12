# Hydrosheaf Web Backend

FastAPI backend for the Hydrosheaf web application, now **fully integrated** with the real Hydrosheaf core engine!

## What Changed?

🎉 **The critical issue has been fixed!** The backend now uses the **actual Hydrosheaf geochemical modeling engine** instead of mock data.

### Before (Mock Data)
```python
# Old code - FAKE results
results = {
    "transport_model": {"dominant_process": "mixing", ...},
    "reactions": [
        {"mineral": "Calcite", "rate_mmol_L": 0.23, ...}  # Hardcoded!
    ]
}
```

### After (Real Hydrosheaf)
```python
# New code - REAL analysis!
from hydrosheaf import fit_network_pipeline, Config

edge_results, extras = fit_network_pipeline(
    hydrosheaf_samples,
    edges,
    config,
    auto_disable_missing=True,
)
# Real transport modeling, LASSO reaction fitting, etc.
```

## Setup Instructions

### 1. Install Dependencies

First, install the basic backend requirements:

```bash
cd web/backend

# Create virtual environment (recommended)
python -m venv venv

# Activate it
# Windows:
venv\Scripts\activate
# macOS/Linux:
source venv/bin/activate

# Install base requirements
pip install -r requirements.txt
```

### 2. Install Hydrosheaf Core Engine

**CRITICAL**: You must install the Hydrosheaf package for the backend to work with real data.

```bash
# Still in web/backend directory
pip install ../../

# This installs Hydrosheaf from the repository root
# It includes all dependencies: scipy, pymc, arviz, PyYAML, etc.
```

**Optional** - Install with additional features:
```bash
pip install ../../[plot]      # Add matplotlib for plotting
pip install ../../[phreeqc]   # Add PHREEQC integration
pip install ../../[modpath]   # Add MODPATH support
```

### 3. Verify Installation

```bash
# Check if Hydrosheaf is available
python -c "import hydrosheaf; print(f'Hydrosheaf v{hydrosheaf.__version__} installed!')"

# Should output:
# Hydrosheaf v0.1.0 installed!
```

### 4. Start the Backend Server

```bash
uvicorn app.main:app --reload --port 8000
```

The API will be available at:
- **Main API**: http://localhost:8000
- **Interactive Docs**: http://localhost:8000/api/docs
- **Health Check**: http://localhost:8000/api/health

### 5. Run Integration Tests

In a **separate terminal**, verify everything works:

```bash
cd web/backend
python test_integration.py
```

This will:
1. Check health endpoints
2. Run a full analysis with real Hydrosheaf
3. Test network flow inference
4. Verify results are not mock data

## API Overview

### Analysis Endpoints (`/api/analysis/`)

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/run` | POST | Start real Hydrosheaf analysis |
| `/status/{job_id}` | GET | Check analysis status |
| `/results/{job_id}` | GET | Get analysis results |
| `/jobs` | GET | List all analysis jobs |
| `/jobs/{job_id}` | DELETE | Delete a job |
| `/health` | GET | Check Hydrosheaf availability |

### Network Endpoints (`/api/network/`)

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/create` | POST | Create flow network |
| `/{id}/infer-flow` | POST | Infer flow directions (uses real Hydrosheaf!) |
| `/health` | GET | Check network inference status |

### Sample Endpoints (`/api/samples/`)

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/upload` | POST | Upload sample dataset |
| `/upload-file` | POST | Upload from JSON file |
| `/datasets` | GET | List datasets |
| `/datasets/{id}` | GET | Get dataset details |

### Project Endpoints (`/api/projects/`)

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/create` | POST | Create project |
| `/{id}/download` | GET | Download report |

## Architecture

```
web/backend/
├── app/
│   ├── main.py                  # FastAPI app entry point
│   ├── hydrosheaf_adapter.py    # NEW: Data transformation layer
│   └── routers/
│       ├── analysis.py          # UPDATED: Real Hydrosheaf integration
│       ├── network.py           # UPDATED: Probabilistic inference
│       ├── samples.py           # Sample data management
│       └── projects.py          # Project management
├── requirements.txt             # Dependencies
├── test_integration.py          # NEW: Integration tests
└── README.md                    # This file
```

### Key Components

**`hydrosheaf_adapter.py`** - Transforms data between web API and Hydrosheaf formats:
- `ConfigAdapter`: Maps frontend config → Hydrosheaf `Config`
- `SampleAdapter`: Converts samples (mg/L → mmol/L, field names)
- `ResultAdapter`: Converts Hydrosheaf `EdgeResult` → frontend JSON

**Updated `analysis.py`**:
- Imports real `fit_network_pipeline` from Hydrosheaf
- Converts frontend data using adapters
- Runs actual geochemical analysis
- Returns real results (not mock data!)

**Updated `network.py`**:
- Uses `infer_edges_probabilistic` for Bayesian flow inference
- Falls back to simple gradient method if Hydrosheaf unavailable
- Real graph theory algorithms

## Data Flow

### Frontend → Backend → Hydrosheaf

1. **Frontend sends** (mg/L, lowercase fields):
```json
{
  "sample_id": "W1",
  "location_id": "Well_1",
  "ca": 80.0,
  "mg": 30.0,
  ...
}
```

2. **Adapter converts** (mmol/L, capitalized fields):
```python
{
  "site_id": "Well_1",
  "Ca": 1.996,  # Converted from mg/L
  "Mg": 1.234,
  ...
}
```

3. **Hydrosheaf processes**:
```python
config = Config(lambda_l1=0.1, phreeqc_enabled=False, ...)
edge_results, extras = fit_network_pipeline(samples, edges, config)
```

4. **Adapter converts back**:
```json
{
  "transport_model": {"dominant_process": "evap", ...},
  "reactions": [{"mineral": "Calcite", "rate_mmol_L": 0.234, ...}],
  "metadata": {"mock_data": false, ...}
}
```

## Configuration Mapping

| Frontend Parameter | Hydrosheaf Config | Notes |
|-------------------|-------------------|-------|
| `lasso_penalty` | `lambda_l1` | L1 regularization strength |
| `enable_phreeqc` | `phreeqc_enabled` | Thermodynamic constraints |
| `enable_isotopes` | `isotope_enabled` | δ18O/δ2H analysis |
| `bootstrap_iterations` | N/A | Function parameter, not config |

## Troubleshooting

### "Hydrosheaf not available" error

**Problem**: Backend can't import Hydrosheaf.

**Solution**:
```bash
cd web/backend
pip install ../../
```

Verify with:
```bash
python -c "import hydrosheaf"
```

### Analysis fails immediately

**Check**:
1. Sample data has required fields (Ca, Mg, Na, HCO3, Cl, SO4, NO3, F, Fe, PO4)
2. pH values are present if PHREEQC is enabled
3. Units are correct (mg/L in frontend, will be converted)

**Debug**:
- Check backend console logs for detailed error messages
- Review `error_preview` in job status response
- Use `/api/analysis/health` to verify Hydrosheaf is loaded

### Results look wrong

**Verify**:
1. Check `metadata.mock_data` is `false` in results
2. Compare with CLI results for same input data
3. Ensure edge inference is working (coordinates provided or edges specified)

### PHREEQC errors

If you enable PHREEQC but don't have it installed:
- Set `enable_phreeqc: false` in analysis config, OR
- Install PHREEQC: `pip install ../../[phreeqc]`

## Testing

### Manual API Testing

Visit http://localhost:8000/api/docs for interactive API documentation where you can:
- Upload samples
- Start analyses
- Check results
- All without writing code!

### Automated Testing

```bash
# Run full integration test
python test_integration.py

# Test specific endpoint with curl
curl http://localhost:8000/api/health
curl http://localhost:8000/api/analysis/health
```

### Test with Demo Data

The frontend includes demo samples. Start both backend and frontend:

**Terminal 1 - Backend:**
```bash
cd web/backend
uvicorn app.main:app --reload
```

**Terminal 2 - Frontend:**
```bash
cd web/frontend
npm run dev
```

Then:
1. Open http://localhost:5173
2. Go to "Analysis" page
3. Select "Demo Groundwater Samples"
4. Click "Start Analysis"
5. Wait for real Hydrosheaf results!

## Performance Notes

**Analysis Speed**:
- Small networks (< 5 samples): ~1-2 seconds
- Medium networks (5-20 samples): ~3-10 seconds
- Large networks (> 20 samples): ~10-60 seconds

**Optimization Tips**:
- Disable PHREEQC if not needed (faster)
- Disable uncertainty quantification for quick tests
- Reduce bootstrap iterations (100 instead of 1000)
- Use specific edges instead of full inference

## Development

### Adding New Features

To expose more Hydrosheaf features:

1. Import from Hydrosheaf:
```python
from hydrosheaf.temporal import fit_temporal_edges
```

2. Add endpoint:
```python
@router.post("/run-temporal")
async def run_temporal_analysis(...):
    # Implementation
```

3. Update adapters if needed for data transformation

4. Add tests to `test_integration.py`

### Code Style

- Follow FastAPI best practices
- Use type hints for all functions
- Add docstrings for API endpoints
- Handle errors gracefully (don't expose internal tracebacks to users)

## License

MIT License (same as Hydrosheaf)

## Authors

Backend integration by:
- Dickson Abdul-Wahab
- Ebenezer Aquisman Asare
- Abdul Rashid Dickson

Powered by **Hydrosheaf** - Sheaf-Theoretic Methods in Groundwater Hydrogeochemistry
