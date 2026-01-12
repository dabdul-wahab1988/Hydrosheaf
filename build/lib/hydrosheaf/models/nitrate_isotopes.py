"""
Dual Isotope Mixing Model for Nitrate Source Apportionment.
"""

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
from scipy.stats import multivariate_normal

# Default path to endmembers DB
DEFAULT_DB_PATH = Path(__file__).resolve().parent.parent / "databases" / "nitrate_endmembers.json"


@dataclass
class IsotopeSample:
    d15N: float
    d18O: float


@dataclass
class SourceIsotopes:
    name: str
    d15N_mean: float
    d15N_std: float
    d18O_mean: float
    d18O_std: float


def load_isotope_endmembers(path: Path = DEFAULT_DB_PATH) -> List[SourceIsotopes]:
    """Load endmember definitions from JSON."""
    if not path.exists():
        return []
    
    with open(path, "r") as f:
        data = json.load(f)
    
    sources = []
    for name, params in data.get("sources", {}).items():
        sources.append(SourceIsotopes(
            name=name,
            d15N_mean=params["d15N"]["mean"],
            d15N_std=params["d15N"]["std"],
            d18O_mean=params["d18O"]["mean"],
            d18O_std=params["d18O"]["std"]
        ))
    return sources


def compute_isotope_prob(
    sample: IsotopeSample, 
    sources: List[SourceIsotopes],
    prior_probs: Optional[Dict[str, float]] = None
) -> Dict[str, float]:
    """
    Compute posterior probability of each source using a multivariate normal likelihood.
    
    P(Source|Sample) \propto P(Sample|Source) * P(Source)
    
    P(Sample|Source) ~ N( [d15N, d18O], Sigma )
    where Sigma is diagonal matrix of variances (std^2).
    """
    if not sources:
        return {}

    posteriors = {}
    total_likelihood = 0.0
    
    # 1. Calculate Likelihoods
    likelihoods = {}
    
    for src in sources:
        # Covariance matrix (assuming independence between N and O for simplicity, 
        # though denitrification couples them with slope 0.5)
        cov = [[src.d15N_std**2, 0], [0, src.d18O_std**2]]
        mean = [src.d15N_mean, src.d18O_mean]
        
        # P(data | source)
        try:
            lik = multivariate_normal.pdf([sample.d15N, sample.d18O], mean=mean, cov=cov)
        except Exception:
            lik = 0.0
        
        # Apply prior
        prior = 1.0 / len(sources)
        if prior_probs and src.name in prior_probs:
            prior = prior_probs[src.name]
            
        unnormalized = lik * prior
        likelihoods[src.name] = unnormalized
        total_likelihood += unnormalized
        
    # 2. Normalize
    if total_likelihood <= 0:
        # Fallback to uniform if sample is extremely far from all sources
        uniform_p = 1.0 / len(sources)
        return {src.name: uniform_p for src in sources}
        
    for name, val in likelihoods.items():
        posteriors[name] = val / total_likelihood
        
    return posteriors
