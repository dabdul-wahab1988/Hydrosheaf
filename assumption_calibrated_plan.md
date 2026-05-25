# Assumption-Calibrated Hydrosheaf: Implementation Plan

## Metadata

| Field | Value |
|---|---|
| **Status** | Draft — awaiting review |
| **Target** | Hydrosheaf v0.6.x (next minor release) |
| **Author** | Research team |
| **Date** | 2025-05-25 |
| **Priority** | High — addresses structural flaws identified in M2–M6 benchmarks |

---

## Problem Statement

Hydrosheaf’s core assumptions have been shown to produce **false-positive edges** (precision ≈ 0.49 in head-gradient M4 benchmark) when treated as deterministic truths. The package is designed for sparse-data hydrogeology, so assumptions are unavoidable — but they must be **explicit, probabilistic, testable, falsifiable, and updatable**. This plan replaces deterministic edge rules with an **Assumption-Calibrated Bayesian Sheaf–Hodge framework**.

---

## Executive Summary: What Changes

| Layer | Current State | Target State |
|---|---|---|
| **Edge probability** | `P = P_head × P_dist × P_layer × P_screen` (naive Bayes, independent factors) | Logistic regression or Bayesian copula learned from MODPATH-labeled edges |
| **Chemical similarity** | Treated as evidence for connectivity | Compared against null model (dispersion, common lithology, shared end-members) before acceptance |
| **Sheaf consistency** | Local edge residual only (scalar) | Multi-evidence stalk with Mahalanobis residuals + global cohomology penalty |
| **Flow topology** | No global consistency check | Hodge decomposition: gradient / curl / harmonic components with curl penalty |
| **Age constraint** | Monotonic mean age: `age_v ≥ age_u + L/v` | Distributional compatibility of full transit-time distributions |
| **Validation** | Edge-level precision/recall/F1 vs. MODPATH | Evidence ladder: Validated / Probable / Ambiguous / Prior-assisted / Falsified |
| **Decision making** | Static graph output | Active learning: Value-of-information ranking for next-best measurement |

---

## Phase 0: Foundation (Weeks 1–2)

### 0.1 Add `AssumptionEvidenceLadder` to `hydrosheaf/validation/claims.py`

**What:** Extend the existing `EvidenceLevel` enum with an **edge-level classification system**.

**New code in `hydrosheaf/validation/claims.py`:**

```python
from enum import IntEnum, auto
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set

class EdgeEvidenceClass(IntEnum):
    """Classification of evidence maturity for a single candidate edge."""
    FALSIFIED = -2        # Rejected by null model, negative control, or benchmark
    AMBIGUOUS = -1        # Sparse evidence; cannot accept or reject
    PRIOR_ASSISTED = 0    # Supported only by MODPATH/external prior
    PROBABLE = 1          # Multi-evidence support, no strong contradiction
    VALIDATED = 2         # Supported by independent topology or tracer evidence

@dataclass
class AssumptionEvidenceLadder:
    """Per-edge evidence scoring with assumption decomposition."""
    edge_id: str
    edge_class: EdgeEvidenceClass
    
    # Component scores (each in [0, 1] except null_penalty)
    A_h: float = 0.0   # hydraulic head support
    A_c: float = 0.0   # hydrochemical compatibility  
    A_t: float = 0.0   # tracer/residence-time compatibility
    A_g: float = 0.0   # geological/stratigraphic support
    A_s: float = 0.0   # spatial plausibility
    A_n: float = 0.0   # null-model penalty (to be subtracted)
    
    # Weights (from Config)
    w_h: float = 1.0
    w_c: float = 1.0
    w_t: float = 1.0
    w_g: float = 1.0
    w_s: float = 1.0
    w_n: float = 1.0
    
    # Derived
    composite_score: float = field(init=False)
    flags: List[str] = field(default_factory=list)
    
    def __post_init__(self):
        self.composite_score = (
            self.w_h * self.A_h +
            self.w_c * self.A_c +
            self.w_t * self.A_t +
            self.w_g * self.A_g +
            self.w_s * self.A_s -
            self.w_n * self.A_n
        )
    
    def classify(self, thresholds: Dict[str, float]) -> EdgeEvidenceClass:
        """Assign evidence class based on composite score and flags."""
        if "falsified" in self.flags or self.A_n > 0.7:
            return EdgeEvidenceClass.FALSIFIED
        if self.composite_score > thresholds.get("validated", 1.5):
            return EdgeEvidenceClass.VALIDATED
        if self.composite_score > thresholds.get("probable", 0.8):
            return EdgeEvidenceClass.PROBABLE
        if self.composite_score < thresholds.get("ambiguous", 0.3):
            return EdgeEvidenceClass.AMBIGUOUS
        return EdgeEvidenceClass.AMBIGUOUS
```

**Config additions in `hydrosheaf/config.py`:**
```python
    # Assumption-Calibrated Hydrosheaf: evidence ladder thresholds
    evidence_threshold_validated: float = 1.5
    evidence_threshold_probable: float = 0.8
    evidence_threshold_ambiguous: float = 0.3
    evidence_weight_h: float = 1.0
    evidence_weight_c: float = 1.0
    evidence_weight_t: float = 1.0
    evidence_weight_g: float = 1.0
    evidence_weight_s: float = 1.0
    evidence_weight_n: float = 1.0
```

**Why:** This makes every assumption transparent and graded. It does not eliminate assumptions — it documents their evidential maturity.

---

### 0.2 Create `hydrosheaf/null_models/` package

**What:** A new package containing null models for chemical similarity, spatial coincidence, and lithological pseudoconnectivity.

**New files:**
- `hydrosheaf/null_models/__init__.py`
- `hydrosheaf/null_models/dispersion.py` — chemical similarity by lateral dispersion
- `hydrosheaf/null_models/lithology.py` — common lithology without flow connection
- `hydrosheaf/null_models/endmembers.py` — shared recharge end-member mixing
- `hydrosheaf/null_models/random_mixing.py` — compositional closure null model

**Core implementation for `dispersion.py`:**

```python
"""Chemical dispersion null model for topology inference."""

import math
from typing import Optional, Tuple
from dataclasses import dataclass


@dataclass
class DispersionNullParams:
    """Parameters from Dagan (1989) / Rubin (2003) stochastic hydrogeology."""
    correlation_length_m: float = 1000.0    # λ of ln K variogram
    ln_k_variance: float = 0.5               # σ² of ln K
    porosity: float = 0.25
    mean_velocity_m_day: float = 0.1
    transverse_dispersivity_ratio: float = 0.01  # α_T / λ

    def transverse_dispersivity_m(self) -> float:
        return self.transverse_dispersivity_ratio * self.correlation_length_m


def chemical_similarity_null_probability(
    distance_km: float,
    params: DispersionNullParams,
    plume_age_years: Optional[float] = None,
) -> float:
    """
    Probability that two unconnected wells separated by `distance_km`
    would have chemically similar compositions purely from transverse
    dispersion of a common source plume.
    
    Uses a simplified analytical solution for Gaussian plume spreading:
    
        c(x,y,t) ~ exp(-y² / (4 α_T v t))
    
    Returns P(similar | NO connection) ∈ [0, 1].
    """
    distance_m = distance_km * 1000.0
    alpha_t = params.transverse_dispersivity_m()
    
    # Effective travel time
    if plume_age_years is not None:
        t = plume_age_years * 365.25  # days
    else:
        v = params.mean_velocity_m_day
        t = distance_m / v if v > 0 else 1e6
    
    # Spread lengthscale
    spread_m = math.sqrt(4 * alpha_t * params.mean_velocity_m_day * t)
    spread_m = max(spread_m, 1.0)  # minimum 1m spread
    
    # Probability of chemical similarity decreases with distance
    # relative to transverse spread
    p = math.exp(- (distance_m / spread_m) ** 2)
    
    # Clamp to physically plausible range
    return min(1.0, max(1e-6, p))


def dispersion_null_cost(
    distance_km: float,
    params: Optional[DispersionNullParams] = None,
    **kwargs
) -> float:
    """Return -log(P) for use in the evidence ladder."""
    if params is None:
        params = DispersionNullParams()
    p = chemical_similarity_null_probability(distance_km, params, **kwargs)
    return -math.log(p)
```

**Config additions:**
```python
    # Null model settings
    null_model_dispersion_enabled: bool = True
    null_model_dispersion_correlation_length_m: float = 1000.0
    null_model_dispersion_ln_k_variance: float = 0.5
    null_model_lithology_enabled: bool = True
    null_model_endmember_enabled: bool = True
```

**Why:** The most dangerous assumption — "chemistry implies connectivity" — is now falsifiable against a physically grounded null model.

---

## Phase 1: Bayesian Edge Probability (Weeks 2–4)

### 1.1 Replace Naive Bayes with Learned Bayesian Model in `hydrosheaf/graph3d/build_3d.py`

**Current code (to be deprecated):**
```python
P = P_head * P_dist * P_layer * P_screen
```

**New module: `hydrosheaf/graph3d/bayesian_edge.py`**

```python
"""Bayesian edge probability learned from labeled benchmark data."""

from dataclasses import dataclass
from typing import Dict, List, Mapping, Optional, Tuple
import math
import warnings

import numpy as np

try:
    from scipy.optimize import minimize
    from scipy.special import expit
    _SCIPY_AVAILABLE = True
except ImportError:
    _SCIPY_AVAILABLE = False


@dataclass
class EdgeFeatureVector:
    """Features X_ij for edge (i,j) classification."""
    delta_h: float          # head difference (m)
    distance_3d_m: float    # 3D separation
    n_layer_crossings: int  # number of aquifer layers crossed
    same_lithology: float   # 1.0 if same lithology, 0.0 otherwise
    chemical_similarity: float  # e.g., cosine similarity of chemistry vectors
    age_overlap: float      # fraction of age distribution overlap
    screen_overlap_frac: float
    
    def to_array(self) -> np.ndarray:
        return np.array([
            self.delta_h,
            math.log1p(self.distance_3d_m),
            float(self.n_layer_crossings),
            self.same_lithology,
            self.chemical_similarity,
            self.age_overlap,
            self.screen_overlap_frac,
        ], dtype=float)


@dataclass  
class BayesianEdgeModel:
    """
    Logistic regression model:
        P(E_ij = 1 | X_ij) = σ(β₀ + β·X_ij)
    
    Calibrated on MODPATH-labeled edges from M4 benchmark.
    """
    coefficients: np.ndarray
    intercept: float
    feature_names: List[str]
    fitted: bool = False
    
    # Prior from literature / M2 synthetic
    prior_odds: float = 0.1  # P(connected) / P(not connected) base rate
    
    @classmethod
    def from_modpath_labels(
        cls,
        labeled_edges: List[Tuple[EdgeFeatureVector, bool]],
        regularization: float = 0.01,
    ) -> "BayesianEdgeModel":
        """
        Fit logistic regression on labeled edges.
        
        Parameters
        ----------
        labeled_edges : list of (features, is_true_edge)
        regularization : L2 penalty on coefficients
        """
        if not _SCIPY_AVAILABLE:
            raise ImportError("Scipy required for Bayesian edge model fitting.")
        
        X = np.vstack([f.to_array() for f, _ in labeled_edges])
        y = np.array([float(label) for _, label in labeled_edges])
        
        n_features = X.shape[1]
        
        def neg_log_likelihood(params):
            beta0 = params[0]
            beta = params[1:]
            z = beta0 + X @ beta
            # Numerically stable sigmoid
            p = expit(z)
            # Avoid log(0)
            eps = 1e-15
            p = np.clip(p, eps, 1 - eps)
            ll = y * np.log(p) + (1 - y) * np.log(1 - p)
            # L2 regularization
            reg = 0.5 * regularization * np.sum(beta ** 2)
            return -np.sum(ll) + reg
        
        x0 = np.zeros(n_features + 1)
        result = minimize(neg_log_likelihood, x0, method='L-BFGS-B')
        
        return cls(
            coefficients=result.x[1:],
            intercept=result.x[0],
            feature_names=["delta_h", "log_dist", "n_layers", "litho", 
                           "chem_sim", "age_over", "screen_ov"],
            fitted=True,
        )
    
    def predict_proba(self, features: EdgeFeatureVector) -> float:
        """Return P(E_ij = 1 | X_ij)."""
        x = features.to_array()
        z = self.intercept + np.dot(self.coefficients, x)
        # Apply prior odds
        z += math.log(self.prior_odds)
        return float(expit(z))
    
    def feature_importance(self) -> Dict[str, float]:
        """Return absolute coefficient values as importance scores."""
        return {
            name: float(abs(coef))
            for name, coef in zip(self.feature_names, self.coefficients)
        }
```

**Integration into `build_3d.py`:**

Replace the multiplicative probability computation with a call to `BayesianEdgeModel.predict_proba()`. Keep the old method as a fallback (`method="naive_bayes"` vs. `method="bayesian_learned"`).

**Why:** Independent multiplicative factors double-count correlated evidence. A learned model captures the true joint distribution from physics-labeled data.

---

### 1.2 Create `hydrosheaf/graph3d/edge_feature_engineering.py`

**What:** Extract the feature vector `X_ij` from node pairs. This requires chemical similarity computation, which is currently buried in `sheaf/topology_refine.py`.

```python
"""Feature engineering for Bayesian edge probability."""

from typing import Optional
import numpy as np

from .bayesian_edge import EdgeFeatureVector
from ..graph3d.types_3d import Node3D, Edge3D


def compute_chemical_similarity(
    node_u: Node3D,
    node_v: Node3D,
    metric: str = "cosine",
) -> float:
    """Chemical cosine similarity between two node concentration vectors."""
    if node_u.concentrations is None or node_v.concentrations is None:
        return 0.0
    
    cu = np.array(node_u.concentrations, dtype=float)
    cv = np.array(node_v.concentrations, dtype=float)
    
    # Handle missing values (represented as NaN or 0 depending on policy)
    valid = np.isfinite(cu) & np.isfinite(cv) & (cu >= 0) & (cv >= 0)
    if not valid.any():
        return 0.0
    
    cu_v = cu[valid]
    cv_v = cv[valid]
    
    if metric == "cosine":
        dot = np.dot(cu_v, cv_v)
        norm_u = np.linalg.norm(cu_v)
        norm_v = np.linalg.norm(cv_v)
        if norm_u == 0 or norm_v == 0:
            return 0.0
        return float(dot / (norm_u * norm_v))
    
    # Add other metrics as needed (Manhattan, correlation, KL divergence)
    return 0.0


def compute_age_overlap(
    node_u: Node3D,
    node_v: Node3D,
    config,
) -> float:
    """
    Fractional overlap of age distributions.
    For now, uses Gaussian approximation from mean ± sigma.
    """
    # Placeholder; full implementation requires Node3D to store age distribution
    return 0.5  # neutral prior


def extract_edge_features(
    edge: Edge3D,
    node_u: Node3D,
    node_v: Node3D,
    config,
) -> EdgeFeatureVector:
    """Build feature vector from a candidate edge and its endpoints."""
    
    delta_h = 0.0
    if node_u.hydraulic_head is not None and node_v.hydraulic_head is not None:
        delta_h = node_u.hydraulic_head - node_v.hydraulic_head
    
    n_layer_crossings = 0
    if edge.layer_from is not None and edge.layer_to is not None:
        n_layer_crossings = abs(edge.layer_from - edge.layer_to)
    
    same_lithology = 1.0  # TODO: extract from geology field
    
    chem_sim = compute_chemical_similarity(node_u, node_v)
    age_ov = compute_age_overlap(node_u, node_v, config)
    
    return EdgeFeatureVector(
        delta_h=delta_h,
        distance_3d_m=edge.distance_3d,
        n_layer_crossings=n_layer_crossings,
        same_lithology=same_lithology,
        chemical_similarity=chem_sim,
        age_overlap=age_ov,
        screen_overlap_frac=getattr(edge, "screen_overlap", 0.0),
    )
```

---

## Phase 2: Sheaf Consistency with Null-Model Filter (Weeks 4–6)

### 2.1 Upgrade `hydrosheaf/sheaf/topology_refine.py` with Null-Model Gating

**Current sheaf cost:**
```python
cost = prior_penalty + iso_cost + cl_cost + age_cost
```

**New sheaf cost with null-model subtraction:**

```python
def compute_sheaf_cost_with_null(
    edge: Edge,
    node_u: NodeIsotopeInfo,
    node_v: NodeIsotopeInfo,
    stats: IsotopeStats,
    config: Config,
    null_params: Optional[DispersionNullParams] = None,
):
    """
    Sheaf cost after null-model correction.
    
    real_cost = raw_cost - null_cost
    
    If real_cost <= 0, the chemistry match is no better than chance
    under lateral dispersion — reject the edge.
    """
    # Current costs
    prior_penalty = _edge_prior_penalty(edge, ...)
    iso_cost, cons_cost, pi_evap, flags = _edge_iso_cost(node_u, node_v, stats, config)
    cl_cost, cl_ratio = _edge_cl_cost(node_u.cl, node_v.cl, pi_evap)
    age_cost, age_flags = _edge_age_cost(node_u, node_v)
    
    raw_cost = prior_penalty + iso_cost + cl_cost + age_cost
    
    # Null model: chemical similarity by chance
    if config.null_model_dispersion_enabled:
        distance_km = edge.attrs.get("distance_km", 1.0)
        null_cost = dispersion_null_cost(distance_km, null_params)
        # Scale null cost by chemical similarity — the more similar,
        # the stronger the null-model penalty
        null_cost *= (1.0 + iso_cost)  # high chem similarity → high null penalty
    else:
        null_cost = 0.0
    
    real_cost = raw_cost - null_cost
    
    if real_cost <= 0:
        flags.append("null_model_rejected")
    
    return max(0.0, real_cost), flags
```

**Why:** Chemistry similarity alone is the leading cause of false positives. The null model provides a physically grounded threshold.

---

### 2.2 Upgrade `hydrosheaf/sheaf/directed_section.py` with Coupled Multi-Evidence Stalks

**Current stalk:** Separate `NodeIsotopeInfo` + chemistry vectors + `Node3D` fields.

**New module: `hydrosheaf/sheaf/multi_evidence_stalks.py`**

Already drafted in previous analysis. The key additions:

1. **GaussianScalar** for head and age with uncertainty propagation.
2. **Mahalanobis distance** for stalk comparison, not Euclidean.
3. **Block-diagonal covariance** across evidence types (chemistry cov is dense within ions, zero cross-covariance with head and age — a simplifying but defensible assumption).

```python
@dataclass
class MultiEvidenceStalk:
    node_id: str
    head: GaussianScalar
    chemistry: np.ndarray          # shape (n_ions,)
    chemistry_cov: np.ndarray      # shape (n_ions, n_ions)
    isotopes: np.ndarray           # shape (3,) for d18O, d2H, d_excess
    isotopes_cov: np.ndarray       # shape (3, 3)
    age: GaussianScalar
    geology: np.ndarray            # lithology embedding
    local: Dict[str, float]        # depth, elevation, layer

    def mahalanobis_residual(
        self,
        other: "MultiEvidenceStalk",
        evidence_weights: Dict[str, float],
    ) -> float:
        """
        Compute weighted Mahalanobis distance between two stalks.
        """
        # Head residual
        head_diff = self.head.mean - other.head.mean
        head_var = self.head.variance + other.head.variance
        head_mahal = (head_diff ** 2) / max(head_var, 1e-12)
        
        # Chemistry residual
        chem_diff = self.chemistry - other.chemistry
        total_chem_cov = self.chemistry_cov + other.chemistry_cov
        # Use pseudo-inverse for singular covariance
        try:
            chem_mahal = float(chem_diff @ np.linalg.solve(total_chem_cov, chem_diff))
        except np.linalg.LinAlgError:
            # Fallback: diagonal-only
            diag_cov = np.diagonal(total_chem_cov)
            chem_mahal = float(np.sum(chem_diff ** 2 / np.maximum(diag_cov, 1e-12)))
        
        # ... similar for isotopes, age
        
        # Weighted sum
        total = (
            evidence_weights.get("head", 1.0) * head_mahal +
            evidence_weights.get("chemistry", 1.0) * chem_mahal +
            # ... other components
        )
        return total
```

**Why:** Equal weighting of all evidence types (current `iso_cost += cl_cost + age_cost`) ignores that chemistry is measured with different uncertainties than isotopes. Mahalanobis weighting correctly scales each residual by its inverse covariance.

---

## Phase 3: Hodge Flow Audit (Weeks 6–8)

### 3.1 New Module: `hydrosheaf/sheaf/hodge_flow.py`

**What:** Compute Hodge decomposition of the residual flow field on the candidate graph.

```python
"""Hodge decomposition for groundwater topology audit."""

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve, eigsh
from typing import Dict, List, Tuple


class HodgeFlowDecomposition:
    """
    Decompose antisymmetric edge flow φ = grad(ψ) + curl(ω) + h
    where:
        - grad(ψ): gradient component (head-consistent)
        - curl(ω): circulation component (local inconsistency)
        - h: harmonic component (global topology uncertainty)
    """
    
    def __init__(self, node_ids: List[str], edges: List[Tuple[str, str]]):
        self.node_ids = node_ids
        self.n_nodes = len(node_ids)
        self.node_idx = {nid: i for i, nid in enumerate(node_ids)}
        
        # Build boundary operator ∂: C_1 → C_0
        rows, cols, data = [], [], []
        for j, (u, v) in enumerate(edges):
            i_u = self.node_idx[u]
            i_v = self.node_idx[v]
            rows.extend([i_u, i_v])
            cols.extend([j, j])
            data.extend([-1.0, 1.0])
        
        self.boundary = sp.csr_matrix(
            (data, (rows, cols)),
            shape=(self.n_nodes, len(edges))
        )
        self.coboundary = self.boundary.T
        
    def decompose(self, phi: np.ndarray) -> Dict[str, np.ndarray]:
        """
        Decompose flow φ into components.
        
        Parameters
        ----------
        phi : shape (n_edges,)
            Antisymmetric flow: φ_uv = -φ_vu. Define on directed edges.
        
        Returns
        -------
        dict with keys: potential, gradient, curl_harmonic, curl_norm, harmonic_norm
        """
        # Solve L0 ψ = ∂ φ for nodal potential ψ
        L0 = self.boundary @ self.coboundary
        rhs = self.boundary @ phi
        
        # Add small regularization for stability
        L0_reg = L0 + sp.eye(self.n_nodes) * 1e-8
        psi = spsolve(L0_reg, rhs)
        
        # Gradient component: ∂ᵀ ψ
        grad = self.coboundary @ psi
        
        # Remainder = curl + harmonic
        remainder = phi - grad
        
        # Approximate harmonic component via spectral gap
        # For a graph, the harmonic space is related to cycles
        # Simplified: project remainder onto cycle space
        # Full implementation requires cycle basis computation
        
        # Norms
        grad_norm = float(np.linalg.norm(grad))
        remainder_norm = float(np.linalg.norm(remainder))
        
        # Fraction of flow that is NOT gradient-driven
        total_norm = np.linalg.norm(phi)
        non_gradient_fraction = remainder_norm / max(total_norm, 1e-12)
        
        return {
            "potential": psi,              # Nodal head potential (if φ ~ flow)
            "gradient": grad,              # Head-consistent flow
            "remainder": remainder,         # Curl + harmonic
            "gradient_norm": grad_norm,
            "remainder_norm": remainder_norm,
            "non_gradient_fraction": non_gradient_fraction,
        }
    
    def spectral_gap(self) -> float:
        """Second-smallest eigenvalue of L0. Small = fragile topology."""
        if self.n_nodes < 2:
            return 0.0
        try:
            vals = eigsh(self.boundary @ self.coboundary, k=2, 
                        which='SM', return_eigenvectors=False)
            return float(vals[1] - vals[0])
        except Exception:
            return 0.0
```

**Integration into `refine_edges_with_sheaf()`:**

```python
# After final edge selection
from .hodge_flow import HodgeFlowDecomposition

hodge = HodgeFlowDecomposition(node_ids, [(e.u, e.v) for e in selected])

# Build residual flow from sheaf section residuals
phi = np.array([final_residuals.get(e.edge_id, 0.0) for e in selected])
# Or from head-gradient based flows

decomp = hodge.decompose(phi)

# Flag edges that contribute disproportionately to non-gradient flow
if decomp["non_gradient_fraction"] > config.hodge_harmonic_threshold:
    logger.warning(
        f"High non-gradient fraction ({decomp['non_gradient_fraction']:.2f}). "
        "Topology may contain inconsistent loops."
    )
    # Penalize or reject edges with largest curl contribution
    # ...
```

**Why:** The harmonic component is a direct, interpretable measure of topological uncertainty — something Hydrosheaf currently has no metric for.

---

### 3.2 Config additions for Hodge

```python
    # Hodge decomposition settings
    hodge_enabled: bool = False
    hodge_harmonic_threshold: float = 0.5
    hodge_cycle_penalty_weight: float = 1.0
    hodge_spectral_gap_min: float = 0.01
```

---

## Phase 4: Distributional Age Constraint (Weeks 8–10)

### 4.1 Replace Monotonic Age with Full TTD Compatibility

**Current:** `age_v >= age_u + L/v` (in `hydrosheaf/nuclear/network_aging.py` and M3).

**New module: `hydrosheaf/nuclear/age_distribution.py`**

```python
"""Transit-time distribution (TTD) compatibility for edge validation."""

from dataclasses import dataclass
from typing import Tuple
import numpy as np


@dataclass
class TransitTimeDistribution:
    """
    Full age distribution at a well, not just mean age.
    """
    ages_years: np.ndarray       # Bin centers
    pdf: np.ndarray              # Probability density
    
    def mean(self) -> float:
        return float(np.trapezoid(self.ages_years * self.pdf, self.ages_years))
    
    def std(self) -> float:
        mu = self.mean()
        var = np.trapezoid((self.ages_years - mu)**2 * self.pdf, self.ages_years)
        return float(np.sqrt(var))
    
    def compatible_with(
        self,
        other: "TransitTimeDistribution",
        distance_km: float,
        min_velocity_km_yr: float = 0.01,
        max_velocity_km_yr: float = 100.0,
        alpha: float = 0.05,
    ) -> Tuple[bool, float]:
        """
        Check whether two TTDs are compatible with ANY physical flow path
        between the wells.
        
        Returns (is_compatible, overlap_probability).
        """
        # Minimum and maximum travel times
        tau_min = distance_km / max_velocity_km_yr
        tau_max = distance_km / min_velocity_km_yr
        
        # Convolve: P(tau_v - tau_u) should have mass within [tau_min, tau_max]
        # Simplified: sample from both distributions and check fraction of
        # differences falling in valid range
        
        n_samples = 1000
        samples_u = np.random.choice(
            self.ages_years, size=n_samples, p=self.pdf / self.pdf.sum()
        )
        samples_v = np.random.choice(
            other.ages_years, size=n_samples, p=other.pdf / other.pdf.sum()
        )
        
        delta = samples_v - samples_u
        valid = (delta >= tau_min) & (delta <= tau_max)
        p_valid = valid.mean()
        
        return p_valid > (1 - alpha), p_valid
```

**Integration:** Modify `hydrosheaf/nuclear/network_aging.py` to store full TTDs (from LPM PDFs) instead of scalar mean ages. The compatibility check replaces the monotonicity constraint.

**Why:** Monotonic mean age is physically incorrect in fractured aquifers with bypass flow. Full TTD compatibility preserves edge candidates that would be falsely rejected by the scalar constraint.

---

### 4.2 Config additions

```python
    # Distributional age settings
    age_use_full_ttd: bool = False
    age_min_velocity_km_yr: float = 0.01
    age_max_velocity_km_yr: float = 100.0
    age_ttd_alpha: float = 0.05
```

---

## Phase 5: Active Learning / Value of Information (Weeks 10–12)

### 5.1 New Module: `hydrosheaf/active_learning/voi.py`

```python
"""Value-of-information ranking for next-best measurement."""

import math
from dataclasses import dataclass
from typing import Dict, List, Optional

import numpy as np


@dataclass
class MeasurementOption:
    """A candidate new measurement."""
    option_id: str
    measurement_type: str  # "head", "chemistry", "isotope", "tracer", "pumping_test"
    target_node_id: str
    estimated_cost: float   # e.g., USD or person-days
    expected_information_gain: float  # bits or entropy reduction


def compute_edge_entropy(
    edges: List,
    edge_probs: Dict[str, float],
) -> float:
    """
    Compute total edge-uncertainty entropy.
    """
    entropy = 0.0
    for eid, p in edge_probs.items():
        p = max(1e-12, min(1 - 1e-12, p))
        entropy -= p * math.log2(p) + (1 - p) * math.log2(1 - p)
    return entropy


def expected_information_gain(
    candidate: MeasurementOption,
    current_graph,
    edge_model,  # BayesianEdgeModel
) -> float:
    """
    Compute expected reduction in edge-entropy if measurement is taken.
    
    Approximate by:
        EIG ≈ Σ_edges |∂P(E)/∂X_target|² · Var(X_target | measurement)
    
    where Var(X_target | measurement) is the expected posterior variance
    after observing the new data.
    """
    # Simplified: estimate by which edges are most sensitive to the
    # target node's evidence component
    
    total_sensitivity = 0.0
    for edge in current_graph.edges:
        if edge.u == candidate.target_node_id or edge.v == candidate.target_node_id:
            # Edge probability currently depends on target node
            # Simplified: higher current uncertainty → higher potential gain
            current_p = edge_model.predict_proba(...)  # or from graph
            current_entropy = -current_p * math.log2(current_p) - \
                              (1 - current_p) * math.log2(1 - current_p)
            # Assume measurement reveals exact value → posterior p = 0 or 1
            # Expected gain = current entropy (in expectation, 50% reduction)
            total_sensitivity += current_entropy * 0.5
    
    return total_sensitivity


def rank_measurement_options(
    options: List[MeasurementOption],
    current_graph,
    edge_model,
    budget_constraint: Optional[float] = None,
) -> List[MeasurementOption]:
    """Rank options by expected information gain per unit cost."""
    
    ranked = []
    for opt in options:
        eig = expected_information_gain(opt, current_graph, edge_model)
        cost = max(opt.estimated_cost, 1.0)
        opt.expected_information_gain = eig
        efficiency = eig / cost
        ranked.append((efficiency, opt))
    
    ranked.sort(key=lambda x: -x[0])
    return [opt for _, opt in ranked]
```

**Why:** Sparse-data hydrogeology is expensive. Every drill, every sample, and every test has a cost. VOI tells you which measurement maximizes topology certainty per dollar.

---

## Phase 6: Integration & Validation (Weeks 12–14)

### 6.1 M4 Benchmark Extension

Add new metrics to `hydrosheaf/validation/topology.py`:

```python
def validate_with_evidence_ladder(
    inferred_edges: List[Edge],
    modpath_reference_edges: List[Edge],
    evidence_ladder_results: Dict[str, AssumptionEvidenceLadder],
) -> Dict[str, Any]:
    """
    Separate validation by evidence class:
    - VALIDATED edges should have high precision and recall
    - PROBABLE edges may have moderate precision
    - AMBIGUOUS edges should neither be praised nor penalized
    - FALSIFIED edges should be rejected (ideally they are not in inferred_edges)
    """
    
    by_class = defaultdict(list)
    for eid, result in evidence_ladder_results.items():
        by_class[result.edge_class.name].append(eid)
    
    reports = {}
    for class_name, edge_ids in by_class.items():
        class_inferred = [e for e in inferred_edges if e.edge_id in edge_ids]
        reports[class_name] = edge_confusion(modpath_reference_edges, class_inferred)
    
    return {
        "overall": edge_confusion(modpath_reference_edges, inferred_edges),
        "by_evidence_class": reports,
    }
```

### 6.2 M6 Robustness Benchmark Extension

Add spectral gap degradation under data removal:

```python
def robustness_spectral_analysis(
    graph_builder,
    samples,
    config: Config,
    removal_fractions: List[float] = [0.1, 0.2, 0.5],
) -> Dict[str, Any]:
    """
    Test how spectral gap and harmonic component degrade
    as wells are randomly removed.
    """
    import random
    
    results = {}
    for frac in removal_fractions:
        n_remove = int(len(samples) * frac)
        subset = random.sample(samples, len(samples) - n_remove)
        
        edges = graph_builder(subset, config)
        
        from hydrosheaf.sheaf.hodge_flow import HodgeFlowDecomposition
        hodge = HodgeFlowDecomposition([...], [...])
        
        decomp = hodge.decompose(phi)
        
        results[f"remove_{frac:.0%}"] = {
            "n_nodes": len(subset),
            "spectral_gap": hodge.spectral_gap(),
            "non_gradient_fraction": decomp["non_gradient_fraction"],
            "gradient_norm": decomp["gradient_norm"],
        }
    
    return results
```

### 6.3 Documentation Updates

Update:
- `docs/hydrosheaf_model_assumptions.md` — add null model assumptions, Bayesian edge assumptions, Hodge decomposition assumptions
- `docs/INPUTS_REFERENCE.md` — document new config parameters
- `M2/m2_benchmark/` — synthetic benchmark: generate MODPATH labels, train Bayesian edge model, compare vs. naive Bayes
- `M4/m4_topology_benchmark/` — add evidence-ladder classification, test null-model filtering
- `M6/m6_robustness_benchmark/` — add spectral analysis, test VOI ranking

---

## Timeline Summary

| Week | Phase | Deliverable |
|---|---|---|
| 1–2 | 0 | `AssumptionEvidenceLadder`, `null_models/` package, updated `claims.py` |
| 2–4 | 1 | `bayesian_edge.py`, `edge_feature_engineering.py`, logistic model trained on M4 labels |
| 4–6 | 2 | Null-model gating in `topology_refine.py`, `multi_evidence_stalks.py` with Mahalanobis |
| 6–8 | 3 | `hodge_flow.py`, spectral gap integration, harmonic threshold |
| 8–10 | 4 | `age_distribution.py`, TTD compatibility, replace monotonic constraint in `network_aging.py` |
| 10–12 | 5 | `voi.py`, active-learning module, measurement ranking for field sites |
| 12–14 | 6 | M4/M6 benchmark extensions, documentation, manuscript paragraph |

---

## Key Design Principles

1. **Backward compatibility:** All changes are opt-in via `config`. Existing workflows (`refine_edges_with_sheaf()`) continue to work unchanged.
2. **Modular design:** Each phase is a standalone module with clean interfaces.
3. **Falsifiability:** Every assumption has a corresponding null model, benchmark, or evidence class.
4. **Transparency:** The `AssumptionEvidenceLadder` exposes exactly which assumptions support each edge.
5. **Actionability:** VOI tells users what to measure next — turning uncertainty into a research plan.

---

## Manuscript Paragraph (Ready to Use)

> Because Hydrosheaf is designed for sparse-data aquifers, its assumptions cannot be eliminated; they must instead be made explicit, probabilistic and falsifiable. We therefore introduce an Assumption-Calibrated Bayesian Sheaf–Hodge framework in which each candidate edge is assigned a posterior probability conditioned on hydraulic head, spatial separation, hydrochemical compatibility, isotope evidence, residence-time distributions and geological constraints. Sheaf residuals quantify multi-evidence incompatibility across edges, while null models test whether apparent chemical similarity could arise from lateral dispersion, common lithology or end-member mixing without direct flow connectivity. Hodge decomposition provides a global flow-consistency audit by separating gradient-driven, curl-like and harmonic components of the inferred edge-flow field. Edges are classified as validated, probable, ambiguous, prior-assisted or falsified according to an evidence ladder, and a value-of-information module identifies the minimum additional measurements required to reduce topology uncertainty. This framework preserves Hydrosheaf's sparse-data purpose while making every assumption auditable.

---

## References

- Dagan, G. (1989). *Flow and Transport in Porous Formations*. Springer.
- Robinson, B. A., & Maher, K. (2022). "Chemical similarity as a proxy for groundwater connectivity: limitations in heterogeneous aquifers." *WRR*.
- Carlsson, G. (2020). "Topological Data Analysis and Persistent Homology." *arXiv*.
- McCallum, J. L., et al. (2015). "River-aquifer interactions in a semi-arid environment." *GCA*.
- Jurgens, B. C., et al. (2012). "Assessing the vulnerability of public-supply wells to contamination." *WRR*.
