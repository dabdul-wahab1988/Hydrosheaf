"""Regression checks for the bounded UER Objective 5 audit path.

These checks exercise the live runner rather than trusting a previously saved
CSV.  They specifically protect the two corrections made after the first
Objective 5 audit: isotope missingness is based on measured values, and CBE
failures are included in the admitted-section gate.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

import scripts.analysis.run_uer_objective5_full as runner


ROOT = Path(__file__).resolve().parents[1]


def test_uer_objective5_cost_manifest_is_dated_and_traceable() -> None:
    """Every frontier cost must have a source, date and explicit quote status."""
    manifest = runner.load_cost_manifest()

    assert manifest["manifest_id"] == "UER-O5-COST-PUBLIC-LIST-2026-09-15"
    assert manifest["target_currency"] == "USD"
    assert manifest["fx_basis"]["observation_date"] == "2026-09-14"
    assert manifest["fx_basis"]["cad_to_usd"] == 0.7190

    options = {row["option_id"]: row for row in manifest["candidates"]}
    assert set(options) == {"SF6", "CFC12", "3H_resample", "14C", "3H_3He"}
    assert options["SF6"]["usd_per_sample"] == 340.0
    assert options["CFC12"]["usd_per_sample"] == 340.0
    assert options["3H_resample"]["usd_per_sample"] == 359.50
    assert options["14C"]["usd_per_sample"] == 413.43
    assert options["3H_3He"]["quote_status"] == "stale_public_list_requires_reconfirmation"
    for row in options.values():
        assert row["source_url"]
        assert row["price_date"]
        assert row["quote_status"]
        assert row["matrix"]
        assert row["shipping"]
        assert row["taxes"]


def test_uer_objective5_cost_manifest_preserves_unpriced_field_scope() -> None:
    """The analytical frontier must not silently turn into a field budget."""
    manifest = runner.load_cost_manifest()
    unpriced = set(manifest["unpriced_components"])

    assert "Ghana field labor and technician time" in unpriced
    assert "outbound and inbound international courier/freight" in unpriced
    assert "insurance, customs brokerage, Ghana import/export charges, duties and taxes" in unpriced

    campaign = manifest["campaign_cost_model"]
    assert campaign["status"] == "quote_required"
    assert campaign["included_in_current_frontier"] is False
    assert len(campaign["components"]) >= 10
    assert "shared_transport" in campaign["full_cost_formula"]


def test_uer_objective5_audited_missingness_and_qc_gate(tmp_path: Path, monkeypatch) -> None:
    """Nitrate denominators and CBE admission must reflect observed data."""
    table_dir = tmp_path / "tables"
    table_dir.mkdir()
    monkeypatch.setattr(runner, "TAB_DIR", table_dir)

    data = pd.read_csv(runner.DATA_CSV)
    edges = pd.read_csv(runner.EDGES_CSV)
    audit, _, denominators, gating = runner.run_sheaf_and_tripartite_boundary(data, edges)

    assert len(audit) == 1600
    assert denominators["nitrate_dual_tested"] == 116
    assert denominators["nitrate_pass"] == 116
    assert denominators["nitrate_fail"] == 0
    assert denominators["nitrate_untested"] == 1484
    assert int(audit["nitrate_conflict"].sum()) == 0

    assert denominators["cbe_fail"] == 66
    assert int(audit["cbe_conflict"].sum()) == 66
    assert int(audit["is_conflict"].sum()) == 168
    assert int((~audit["cbe_ok"] & ~audit["is_conflict"]).sum()) == 0

    # Total energy has a different edge denominator.  The same-edge metric is
    # the relevant diagnostic and must be retained in the comparison table.
    assert gating.loc[1, "dirichlet_energy"] < gating.loc[0, "dirichlet_energy"]
    assert gating.loc[1, "coherent_edge_energy"] > gating.loc[0, "coherent_edge_energy"]
    assert "edge_set_comparison" in gating.columns
