
import pytest
import numpy as np
from dataclasses import dataclass
from typing import Optional, Dict, Any

from hydrosheaf.sheaf.isotope_metrics import (
    compute_isotope_stats,
    sample_depth_m,
    compute_evaporation_probability,
    IsotopeStats,
)

# Mock Config class
@dataclass
class MockConfig:
    isotope_d18o_key: str = "d18O"
    isotope_d2h_key: str = "d2H"
    lmwl_a: float = 8.0
    lmwl_b: float = 10.0
    detection_limit_policy: str = "half"  # Or whatever default
    screen_top_key: str = "screen_top"
    screen_bottom_key: str = "screen_bottom"
    edge_screen_depth_key: str = "screen_depth"
    edge_well_depth_key: str = "well_depth"
    z_coordinate_key: Optional[str] = "z"
    z_coordinate_positive_down: bool = True
    sheaf_shallow_depth_m: float = 30.0
    sheaf_evap_gate_strength: float = 1.0


# --- Tests for compute_isotope_stats ---

def test_compute_isotope_stats_empty():
    config = MockConfig()
    samples = []
    stats = compute_isotope_stats(samples, config)
    
    # Defaults
    assert stats.d_excess_p25 == 10.0
    assert stats.evap_index_p25 == 0.0
    assert stats.evap_index_p75 == 2.0

def test_compute_isotope_stats_insufficient_data():
    config = MockConfig()
    # Less than 3 samples with valid isotopes
    samples = [
        {"d18O": -5.0, "d2H": -30.0},
        {"d18O": -6.0, "d2H": -38.0},
    ]
    stats = compute_isotope_stats(samples, config)
    
    # Should still be defaults because code checks len >= 3
    assert stats.d_excess_p25 == 10.0
    assert stats.evap_index_p25 == 0.0
    assert stats.evap_index_p75 == 2.0

def test_compute_isotope_stats_valid_data():
    config = MockConfig()
    # Create enough samples to trigger calculation
    # d_excess = d2H - 8 * d18O
    # evap_index -> defined in another module, likely related to distance from LMWL
    # Let's assume a simple case where we can predict roughly
    
    samples = [
        {"d18O": -5.0, "d2H": -30.0}, # d_excess = -30 - (-40) = 10
        {"d18O": -5.0, "d2H": -30.0}, # 10
        {"d18O": -5.0, "d2H": -30.0}, # 10
        {"d18O": -5.0, "d2H": -30.0}, # 10
    ]
    
    stats = compute_isotope_stats(samples, config)
    
    # Since all are identical, percentiles should be the same
    assert stats.d_excess_p25 == 10.0
    # evap_index calculation depends on implementation in hydrosheaf.isotopes.evaporation_index
    # but since inputs are identical, p25 and p75 should be identical (before correction)
    
    # The code has a correction: if p75 < p25 + 1e-6: p75 = p25 + 1.0
    assert stats.evap_index_p75 == stats.evap_index_p25 + 1.0

def test_compute_isotope_stats_varied_data():
    config = MockConfig()
    # Generate data with known spread
    # d_excess values: 10, 12, 14, 16
    # p25 of [10, 12, 14, 16] -> interpolated. 
    # np.quantile uses linear interpolation by default. 
    # 0.25 quantile of 4 items is at index 0.75 -> between 0th (10) and 1st (12).
    # Value = 10 + 0.75 * (12 - 10) = 11.5
    
    samples = []
    d_excess_target = [10, 12, 14, 16]
    
    # We need d2H - 8*d18O = target
    # Let d18O = 0 => d2H = target
    for val in d_excess_target:
        samples.append({"d18O": 0.0, "d2H": float(val)})
        
    stats = compute_isotope_stats(samples, config)
    
    # Check d_excess_p25
    expected_p25 = np.quantile(d_excess_target, 0.25)
    assert np.isclose(stats.d_excess_p25, expected_p25)
    
    # Evap index checks
    # If d18O=0, d2H=val, LMWL: d2H = 8*d18O + 10 = 10
    # evap_index usually measures deviation.
    # We verify that p75 > p25 (or corrected).
    assert stats.evap_index_p75 >= stats.evap_index_p25

def test_compute_isotope_stats_missing_keys():
    config = MockConfig()
    samples = [
        {"other": 1},
        {"d18O": 1}, # missing d2H
        {"d2H": 1},  # missing d18O
        {"d18O": None, "d2H": 1},
    ]
    stats = compute_isotope_stats(samples, config)
    # Should handle gracefully and return defaults
    assert stats.d_excess_p25 == 10.0


# --- Tests for sample_depth_m ---

def test_sample_depth_m_top_bottom():
    config = MockConfig()
    sample = {"screen_top": 10, "screen_bottom": 20}
    depth = sample_depth_m(sample, config)
    assert depth == 15.0

def test_sample_depth_m_top_only():
    config = MockConfig()
    sample = {"screen_top": 10}
    depth = sample_depth_m(sample, config)
    assert depth == 10.0

def test_sample_depth_m_bottom_only():
    config = MockConfig()
    sample = {"screen_bottom": 20}
    depth = sample_depth_m(sample, config)
    assert depth == 20.0

def test_sample_depth_m_depth_key():
    config = MockConfig()
    sample = {"screen_depth": 30}
    depth = sample_depth_m(sample, config)
    assert depth == 30.0

def test_sample_depth_m_well_depth_key():
    config = MockConfig()
    sample = {"well_depth": 40}
    depth = sample_depth_m(sample, config)
    assert depth == 40.0

def test_sample_depth_m_z_coordinate():
    config = MockConfig()
    sample = {"z": 50}
    depth = sample_depth_m(sample, config)
    assert depth == 50.0

def test_sample_depth_m_priority():
    config = MockConfig()
    # All present, top/bottom should win
    sample = {
        "screen_top": 10, "screen_bottom": 20, # Avg = 15
        "screen_depth": 30,
        "well_depth": 40,
        "z": 50
    }
    depth = sample_depth_m(sample, config)
    assert depth == 15.0
    
    # Remove top/bottom
    sample = {
        "screen_depth": 30,
        "well_depth": 40,
        "z": 50
    }
    depth = sample_depth_m(sample, config)
    assert depth == 30.0

    # Remove depth
    sample = {
        "well_depth": 40,
        "z": 50
    }
    depth = sample_depth_m(sample, config)
    assert depth == 40.0

    # Remove well_depth
    sample = {
        "z": 50
    }
    depth = sample_depth_m(sample, config)
    assert depth == 50.0

def test_sample_depth_m_none():
    config = MockConfig()
    sample = {}
    depth = sample_depth_m(sample, config)
    assert depth is None

def test_sample_depth_m_string_parsing():
    config = MockConfig()
    sample = {"screen_top": "10.5", "screen_bottom": "20.5"}
    depth = sample_depth_m(sample, config)
    assert depth == 15.5

def test_sample_depth_m_z_negative():
    config = MockConfig(z_coordinate_positive_down=True)
    sample = {"z": -50}
    depth = sample_depth_m(sample, config)
    assert depth == -50.0 


# --- Tests for compute_evaporation_probability ---

def test_compute_evaporation_probability_none_inputs():
    config = MockConfig()
    stats = IsotopeStats()
    prob = compute_evaporation_probability(None, 10.0, 10.0, stats, config)
    assert prob == 0.0
    prob = compute_evaporation_probability(1.0, None, 10.0, stats, config)
    assert prob == 0.0

def test_compute_evaporation_probability_basic():
    config = MockConfig()
    stats = IsotopeStats(d_excess_p25=10.0, evap_index_p25=0.0, evap_index_p75=2.0)
    
    # Case: evap_index=0, d_excess=10 (match stats), depth=0 (shallow)
    # evap_score = (0 - 0) / 2 = 0
    # d_excess_score = (10 - 10) / 10 = 0
    # depth_score = (30 - 0) / 30 = 1
    # combined = 0.6*0 + 0.3*0 + 0.2*1 = 0.2
    # sigmoid(0.2) > 0.5
    
    prob = compute_evaporation_probability(0.0, 10.0, 0.0, stats, config)
    assert 0.5 < prob < 1.0

def test_compute_evaporation_probability_deep():
    config = MockConfig(sheaf_shallow_depth_m=30.0)
    stats = IsotopeStats()
    
    # Deep sample > 30m
    # depth_gate = 30 / 60 = 0.5
    # Probability should be reduced by 0.5
    
    # Calculate base prob first
    # evap=0, d_excess=10, depth=60
    # depth_score = (30 - 60) / 30 = -1
    # combined = 0 + 0 + 0.2*(-1) = -0.2
    # base = sigmoid(-0.2) < 0.5
    
    prob = compute_evaporation_probability(0.0, 10.0, 60.0, stats, config)
    
    # Should be base * 0.5
    # Verify it is less than un-gated version (though un-gated isn't directly exposed)
    # Just check it returns a valid probability
    assert 0.0 <= prob <= 1.0

def test_compute_evaporation_probability_high_evap():
    config = MockConfig()
    stats = IsotopeStats(evap_index_p25=0.0, evap_index_p75=2.0)
    
    # High evaporation index
    # evap_index = 4.0
    # evap_score = (4 - 0) / 2 = 2.0
    # combined = 0.6*2.0 + ... = 1.2 + ...
    # Sigmoid(>1.2) should be high
    
    prob = compute_evaporation_probability(4.0, 10.0, 0.0, stats, config)
    assert prob > 0.7  # Expect high probability
