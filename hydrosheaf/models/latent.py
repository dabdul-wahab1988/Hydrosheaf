"""Latent Endmember Identification (Virtual Nodes)."""

import math
from typing import List, Mapping

try:
    import numpy as np
except ImportError:
    np = None

from ..log import get_logger

logger = get_logger(__name__)

def identify_latent_endmembers(
    samples: List[Mapping[str, object]],
    ion_order: List[str],
    n_endmembers: int = 2
) -> List[Mapping[str, object]]:
    """
    Identify potential latent endmembers using Principal Component Analysis.
    
    Returns 'Virtual' samples representing extreme compositions in the dataset.
    """
    if np is None:
        logger.warning("Numpy not available, skipping latent endmember identification.")
        return []

    # 1. Extract Data Matrix
    data_rows = []
    valid_samples = []
    
    for s in samples:
        row = []
        is_valid = True
        for ion in ion_order:
            val = s.get(ion)
            if val is None or not isinstance(val, (int, float)) or math.isnan(val):
                is_valid = False
                break
            row.append(float(val))
        
        if is_valid:
            data_rows.append(row)
            valid_samples.append(s)
            
    if len(data_rows) < n_endmembers + 2:
        return []

    X = np.array(data_rows)
    
    # 2. Log-Ratio Transform (CLR) for compositional validity
    # Replace zeros
    X[X <= 0] = 1e-9
    g = np.exp(np.mean(np.log(X), axis=1, keepdims=True))
    X_clr = np.log(X / g)
    
    # 3. PCA via SVD
    # Center
    mean_clr = np.mean(X_clr, axis=0)
    X_centered = X_clr - mean_clr
    
    U, S, Vt = np.linalg.svd(X_centered, full_matrices=False)
    
    # 4. Identify Extremes along first (n_endmembers - 1) components
    # (n+1 vertices for n-simplex mixing)
    # Heuristic: Pick min and max along PC1, then max along PC2 orthogonal...
    # Simple Heuristic: Points with max projection on PC axes
    
    virtual_samples = []
    
    # We want vertices of the mixing polytope.
    # We construct them artificially by taking mean +/- sigma * PC vector
    # This creates "Pure" endmembers derived from trends.
    
    # Reconstruct in CLR then Inverse CLR
    
    for i in range(n_endmembers):
        # Create a virtual point along PC(i)
        # We go "outwards" from the mean
        # Let's pick the actual sample indices that are extremum
        # This is safer than synthesizing chemistry that might be impossible
        
        pc_scores = U[:, i]
        idx_min = np.argmin(pc_scores)
        idx_max = np.argmax(pc_scores)
        
        # We assume the dataset contains "near-endmembers"
        # Or we can project outwards.
        # Let's project outwards 10% beyond the max/min sample to define a "Source"
        
        vec = Vt[i, :] # Eigenvector
        
        # Max source
        score_max = np.max(pc_scores) * S[i]
        virtual_clr_max = mean_clr + vec * (score_max * 1.1)
        
        # Min source
        score_min = np.min(pc_scores) * S[i]
        virtual_clr_min = mean_clr + vec * (score_min * 1.1)
        
        # Inverse CLR to get concentration profile (sum constrained to 1, need to scale)
        def inv_clr(vec):
            e = np.exp(vec)
            return e / np.sum(e)
            
        comp_max = inv_clr(virtual_clr_max)
        comp_min = inv_clr(virtual_clr_min)
        
        # Scale? We don't know the TDS of the virtual source.
        # Heuristic: Use the TDS of the extremum sample
        tds_max = np.sum(data_rows[idx_max])
        tds_min = np.sum(data_rows[idx_min])
        
        conc_max = comp_max * tds_max
        conc_min = comp_min * tds_min
        
        # Add Virtual Sample 1
        vs1 = {"site_id": f"Virtual_Endmember_PC{i+1}_High", "type": "virtual"}
        for k, ion in enumerate(ion_order):
            vs1[ion] = float(conc_max[k])
        virtual_samples.append(vs1)
        
        # Add Virtual Sample 2
        vs2 = {"site_id": f"Virtual_Endmember_PC{i+1}_Low", "type": "virtual"}
        for k, ion in enumerate(ion_order):
            vs2[ion] = float(conc_min[k])
        virtual_samples.append(vs2)

    return virtual_samples
