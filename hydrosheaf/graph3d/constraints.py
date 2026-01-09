
"""
Physics-based constraints for 3D Flow Networks.

This module provides validation logic to ensure inferred flow paths vary
physically realistic boundaries, specifically regarding hydraulic head gradients
and impossible uphill flow.
"""

from typing import Tuple
import math
from .types_3d import Node3D

def check_hydraulic_feasibility(
    node_u: Node3D, 
    node_v: Node3D,
    min_gradient: float = 1e-4
) -> Tuple[bool, float, float]:
    """
    Check if flow from u -> v is improved by hydraulic head data.
    
    Parameters
    ----------
    node_u, node_v : Node3D
        Upstream and downstream nodes
    min_gradient : float
        Minimum required hydraulic gradient to verify flow.
        
    Returns
    -------
    Tuple[bool, float, float]
        (is_feasible, probability_penalty, gradient)
        
    Implementation
    --------------
    Strict Check:
        If Head(u) and Head(v) are known:
            If Head(u) < Head(v) (uphill): Return (False, 0.0, grad)
        
    Gradient Check:
        grad = (h_u - h_v) / d_3d
        
    Probability:
        If grad > min_gradient: P = 1.0
        If 0 < grad < min_gradient: P = grad / min_gradient (linear penalty)
        If grad < 0: P = 0.0
    """
    
    # Bayesian Fallback: Topographic Prior
    # ------------------------------------
    if node_u.hydraulic_head is None or node_v.hydraulic_head is None:
        # If elevation is missing too, we are truly blind
        if node_u.elevation_m is None or node_v.elevation_m is None:
            return True, 1.0, 0.0

        # Theory: Head ~ Elevation - Depth_to_Water
        # Let Depth ~ Normal(mean_depth, sigma_depth^2)
        # We generally assume mean depths are similar locally, so they cancel out.
        # Variance adds up: Var(diff) = 2 * sigma_depth^2
        
        sigma_depth = 5.0 # Assumed uncertainty in meters (configurable in future)
        sigma_diff = sigma_depth * math.sqrt(2)
        
        # Delta Surface Elevation
        dz_surf = node_u.elevation_m - node_v.elevation_m
        
        # We want P(Head_u > Head_v) = P( (Elv_u - d_u) > (Elv_v - d_v) )
        # = P( d_u - d_v < Elv_u - Elv_v )
        # (d_u - d_v) is Normal(0, sigma_diff)
        
        # Calculate CDF of Normal distribution at dz_surf
        # Phi(x) = 0.5 * (1 + erf(x / (sigma * sqrt(2))))
        try:
            phi = 0.5 * (1 + math.erf(dz_surf / (sigma_diff * math.sqrt(2))))
        except (ValueError, OverflowError):
             phi = 1.0 if dz_surf > 0 else 0.0

        # If probability is high (>0.5), it supports flow.
        # If probability is low (<0.5), it penalizes flow.
        # We treat the probability *itself* as the feasibility score.
        
        # Strictness: If P < 0.1 (strong evidence of uphill flow), reject.
        is_feasible = phi >= 0.1
        
        return is_feasible, float(phi), 0.0

    # Deterministic Check (Data Available)
    # ------------------------------------
    h_u = node_u.hydraulic_head
    h_v = node_v.hydraulic_head
    
    # Compute delta head
    # Positive dh means water flows U -> V
    dh = h_u - h_v
    
    # Geometric distance (approximate if not passed, but we need it for gradient)
    dx = node_u.x - node_v.x
    dy = node_u.y - node_v.y
    dz = node_u.z - node_v.z
    d_3d = math.sqrt(dx**2 + dy**2 + dz**2)
    
    if d_3d < 1e-6:
        return True, 1.0, 0.0
        
    gradient = dh / d_3d
    
    # Constraint Logic
    if gradient <= 0:
        # Flow is uphill or flat -> Impossible
        # (Assuming no pumping, density effects, etc. for this basic check)
        return False, 0.0, gradient
        
    if gradient < min_gradient:
        # Weak gradient penalty
        prob = gradient / min_gradient
        return True, prob, gradient
        
    # Strong gradient
    return True, 1.0, gradient

def validate_mineral_assemblage(
    reaction_signature: list[str],
    allowed_minerals: list[str]
) -> bool:
    """
    Check if predicted reactions match allowed mineral facies.
    
    Parameters
    ----------
    reaction_signature : list[str]
        List of minerals predicted to dissolve/precipitate (reaction_labels)
    allowed_minerals : list[str]
        List of minerals valid for this aquifer zone
        
    Returns
    -------
    bool
        True if all predicted minerals are in the allowed list
    """
    if not allowed_minerals:
        return True
        
    # Check set subset
    pred_set = set(reaction_signature)
    allowed_set = set(allowed_minerals)
    
    return pred_set.issubset(allowed_set)
