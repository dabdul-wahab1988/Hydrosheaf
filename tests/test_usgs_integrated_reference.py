from __future__ import annotations

from M7.m7_nonuniqueness_benchmark.scripts.run_usgs_integrated_reference import (
    DEFAULT_AGE_ROOT,
    DEFAULT_M4_RESULTS,
    inventory_aiken,
    load_modpath_panels,
    load_usgs_age_panel,
)


def test_usgs_age_panel_preserves_model_reference_and_screen_fields() -> None:
    nodes, observations, audit = load_usgs_age_panel(DEFAULT_AGE_ROOT)
    assert audit["status"] == "COMPLETE"
    assert len(nodes) == 1279
    assert len(observations) == 1279
    assert {"node_id", "lat", "lon", "sample_date"}.issubset(nodes.columns)
    assert observations["reported_age_years"].notna().any()


def test_m4_modpath_panels_remain_separate() -> None:
    nodes, edges, travel, audit = load_modpath_panels(DEFAULT_M4_RESULTS)
    assert audit["panels"]["tier_1_savage"]["n_edges"] == 174
    assert audit["panels"]["tier_2_great_miami"]["n_edges"] == 68
    assert len(edges[edges["panel_id"] == "tier_1_savage"]) == 174
    assert len(edges[edges["panel_id"] == "tier_2_great_miami"]) == 68
    assert len(travel) == len(edges)
    assert set(nodes["panel_id"]) == {"tier_1_savage", "tier_2_great_miami"}


def test_aiken_is_mapping_deferred_without_an_extracted_archive() -> None:
    report = inventory_aiken(None)
    assert report["status"] == "NOT_SUPPLIED"
    assert report["mapping_status"] == "DEFERRED"
