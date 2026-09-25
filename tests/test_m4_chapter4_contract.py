from __future__ import annotations

import pandas as pd

from scripts.analysis.generate_m4_chapter4_contract import (
    MODE_ORDER,
    _load_contract,
    generate_table_4_7,
)


def test_m4_chapter4_contract_is_truth_blind_and_separates_modes(tmp_path):
    summary, manifest = _load_contract()

    assert list(summary["mode"]) == MODE_ORDER
    assert set(summary["n_nodes"]) == {153}
    assert set(summary["n_candidate_edges"]) == {23256}
    assert set(summary["n_reference_edges"]) == {174}
    assert summary["truth_blind_inference"].all()
    assert summary["reference_edges_used_only_after_inference"].all()
    assert not summary["calibration_transfer_established"].any()
    assert manifest["candidate_contract"]["all_directed_pairs_default"] is True
    assert manifest["candidate_contract"]["reference_edges_in_candidate_generation"] is False
    assert manifest["candidate_contract"]["reference_edges_in_feature_construction"] is False
    assert manifest["candidate_contract"]["reference_edges_in_threshold_selection"] is False

    csv_path, md_path = generate_table_4_7(tmp_path)
    table = pd.read_csv(csv_path)
    assert len(table) == 2
    assert table["Candidate directed pairs"].tolist() == [23256, 23256]
    assert table["Reference edges (evaluation only)"].tolist() == [174, 174]
    assert "Primary limited-data benchmark" in table["Role"].tolist()
    assert "MODPATH" in md_path.read_text(encoding="utf-8")
    assert "legacy benchmark contract" in md_path.read_text(encoding="utf-8")

    docx_csv = tmp_path / "table_4_7_m4_truth_blind_topology_contract_docx.csv"
    assert docx_csv.exists()
    assert len(pd.read_csv(docx_csv).columns) == 5
