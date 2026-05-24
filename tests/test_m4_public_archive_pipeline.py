"""Unit tests for the M4 public archive pipeline."""

from __future__ import annotations

import os
from pathlib import Path
import pandas as pd
import pytest

from hydrosheaf.validation.modpath_archive import (
    scan_modpath_archive,
    standardise_endpoint_columns,
    standardise_pathline_columns,
    build_node_mapping,
    build_reference_edges,
)

PROJECT_ROOT = Path(__file__).resolve().parents[1]
CONFIG_DIR = PROJECT_ROOT / "M4" / "m4_topology_benchmark" / "configs"
RESULTS_DIR = PROJECT_ROOT / "M4" / "m4_topology_benchmark" / "results"


def test_configs_exist_and_scan():
    """Verify that YAML configurations scan correctly and validate against the schema."""
    configs = ["savage.yaml", "great_miami.yaml", "long_island.yaml"]
    for config_name in configs:
        config_path = CONFIG_DIR / config_name
        assert config_path.exists(), f"Missing config file: {config_name}"
        
        config = scan_modpath_archive(str(config_path))
        assert isinstance(config, dict)
        assert "archive_name" in config
        assert "archive_id" in config
        assert "validation_tier" in config
        assert "local_archive_root" in config
        assert "allowed_claim" in config
        assert "required_guardrail" in config


def test_standardise_endpoint_columns():
    """Test that endpoint columns are standardised correctly even if input columns are missing."""
    df = pd.DataFrame([{"particle_id": 1, "x0": 100.0, "y0": 200.0, "time": 10.0}])
    std_df = standardise_endpoint_columns(df, {})
    
    assert "particle_id" in std_df.columns
    assert "initial_cell" in std_df.columns
    assert std_df.loc[0, "initial_cell"] is None
    assert std_df.loc[0, "x0"] == 100.0


def test_standardise_pathline_columns():
    """Test that pathline columns are standardised correctly."""
    df = pd.DataFrame([{"particle_id": 1, "x": 10.0, "y": 20.0, "time": 5.0}])
    std_df = standardise_pathline_columns(df, {})
    
    assert "particle_id" in std_df.columns
    assert "cell" in std_df.columns
    assert std_df.loc[0, "cell"] is None
    assert std_df.loc[0, "x"] == 10.0


def test_build_node_mapping_dummy():
    """Verify that node mapping extracts unique nodes from endpoints and pathlines."""
    endpoints = pd.DataFrame([
        {"initial_cell": 101, "x0": 1.0, "y0": 2.0, "z0": 3.0, "final_cell": 102, "x": 4.0, "y": 5.0, "z": 6.0}
    ])
    pathlines = pd.DataFrame([
        {"cell": 103, "x": 7.0, "y": 8.0, "z": 9.0}
    ])
    
    node_map = build_node_mapping(endpoints, pathlines, {})
    assert len(node_map) == 3
    node_ids = set(node_map["node_id"])
    assert node_ids == {"cell_101", "cell_102", "cell_103"}


def test_build_reference_edges_dummy():
    """Verify that reference edges are extracted from initial/final cell connectivity."""
    endpoints = pd.DataFrame([
        {"initial_cell": 101, "final_cell": 102},
        {"initial_cell": 101, "final_cell": 102},  # Duplicate, should be ignored
        {"initial_cell": 102, "final_cell": 103},
    ])
    node_map = pd.DataFrame()
    
    edges = build_reference_edges(endpoints, pd.DataFrame(), node_map, {})
    assert len(edges) == 2
    edge_ids = set(edges["edge_id"])
    assert edge_ids == {"cell_101->cell_102", "cell_102->cell_103"}
    assert "particle_count" in edges.columns
    assert "allowed_claim" in edges.columns
    assert "required_guardrail" in edges.columns


def test_pipeline_outputs_generated():
    """Verify that the orchestrated pipeline created all output CSVs and JSONs."""
    expected_files = [
        "external_modpath_archive_summary.csv",
        "external_modpath_source_manifest.json",
        "tier_1_savage_archive_summary.csv",
        "tier_2_great_miami_archive_summary.csv",
        "tier_3_long_island_archive_summary.csv",
    ]
    for filename in expected_files:
        path = RESULTS_DIR / filename
        assert path.exists(), f"Missing output file: {filename}"
