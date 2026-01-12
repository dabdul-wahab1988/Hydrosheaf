# Hydrosheaf Web Integration Roadmap

## Current Issues Identified

### 1. Analysis Execution Bug
- **Issue**: `infer_edges_probabilistic()` function call has incorrect parameters
- **Status**: FIXED
- **Solution**: Removed `config` as positional argument, added `p_min` parameter

### 2. Data-Analysis Mismatch
- **Issue**: Frontend shows all analysis options regardless of data availability
- **Status**: FIXED
- **Solution**: Added `/api/samples/datasets/{id}/capabilities` endpoint, frontend now fetches and uses capabilities

### 3. Module Integration Gaps
- **Issue**: Some hydrosheaf modules not fully wired to backend
- **Status**: FIXED
- **Solution**: Completed adapter layer for all modules in `hydrosheaf_adapter.py`

---

## Data Requirements by Analysis Type

| Analysis Module | Required Data Fields | Optional Fields |
|-----------------|---------------------|-----------------|
| **Transport Modeling** | Ca, Mg, Na, HCO3, Cl, SO4 | K, NO3, F, Fe, PO4 |
| **Reaction Path** | Major ions (above) | pH, EC, TDS |
| **Network Inference** | x/longitude, y/latitude | head_meas, dtw, elevation |
| **Isotope Analysis** | d18o, d2h | - |
| **PHREEQC Constraints** | Major ions, pH, temperature | EC |
| **Nitrate Source** | NO3, d15n, d18o_no3 | PO4, Fe |
| **Temporal Analysis** | date, any ion time series | - |
| **3D Network** | x, y, z/screen_depth | aquifer_layer |
| **Gibbs Diagram** | Na, Ca, Cl, HCO3, TDS | - |
| **Ion Exchange** | Ca, Mg, Na, Cl | - |

---

## Implementation Plan

### Phase 1: Fix Core Analysis (IMMEDIATE) - COMPLETED
1. [x] Fix `infer_edges_probabilistic()` call signature
2. [x] Verify server restart applies the fix
3. [x] Test basic analysis execution

### Phase 2: Add Dataset Analysis Endpoint - COMPLETED
1. [x] Create `/api/samples/datasets/{id}/capabilities` endpoint
2. [x] Detect which fields are present in dataset
3. [x] Return list of available analysis modules

### Phase 3: Smart Frontend Options - COMPLETED
1. [x] Fetch dataset capabilities when dataset selected
2. [x] Disable/hide unavailable analysis options
3. [x] Show warnings for partially available modules
4. [x] Add tooltips explaining data requirements

### Phase 4: Complete Module Integration - COMPLETED
1. [x] Ensure all config options map correctly to hydrosheaf.Config
2. [x] Add result adapters for all module outputs
3. [x] Handle module-specific errors gracefully

### Phase 5: Testing & Validation - IN PROGRESS
1. [x] Test with dataset containing all fields (demo data)
2. [ ] Test with minimal dataset (ions only)
3. [ ] Test with isotope data
4. [ ] Test with spatial coordinates
5. [ ] Test error handling

---

## API Endpoint Design

### GET /api/samples/datasets/{id}/capabilities

**Response:**
```json
{
  "dataset_id": "abc123",
  "sample_count": 41,
  "available_fields": {
    "ions": ["Ca", "Mg", "Na", "HCO3", "Cl", "SO4", "NO3"],
    "isotopes": ["d18o", "d2h"],
    "spatial": ["x", "y"],
    "temporal": ["date"],
    "field_params": ["pH", "EC"]
  },
  "available_analyses": {
    "transport": true,
    "reaction_path": true,
    "network_inference": true,
    "isotope_analysis": true,
    "phreeqc": true,
    "nitrate_source": false,
    "temporal": true,
    "network_3d": false,
    "gibbs": true,
    "exchange": true
  },
  "warnings": [
    "Nitrate source analysis requires d15n and d18o_no3 fields",
    "3D network requires z/screen_depth coordinates"
  ]
}
```

---

## Files Modified

1. **Backend**
   - `web/backend/app/routers/samples.py` - Added capabilities endpoint, CSV parsing
   - `web/backend/app/routers/analysis.py` - Fixed infer_edges_probabilistic call, expanded AnalysisConfig
   - `web/backend/app/hydrosheaf_adapter.py` - Complete module mapping for all features

2. **Frontend**
   - `web/frontend/src/pages/Analysis.jsx` - Added capability detection, smart option display
   - `web/frontend/src/pages/Analysis.css` - Added styles for capabilities panel
   - `web/frontend/src/pages/Samples.jsx` - CSV/JSON upload support

---

## Success Criteria

1. [x] Analysis runs successfully with basic ion data
2. [x] Frontend only shows relevant analysis options
3. [x] Clear error messages when data is insufficient
4. [x] All hydrosheaf modules accessible via web UI
5. [ ] Results display correctly for all analysis types (needs testing)

---

## Completion Date

Integration completed: January 2026
