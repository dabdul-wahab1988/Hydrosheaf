from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from M7.m7_nonuniqueness_benchmark.scripts.run_integrated_truth_blind import (
    audit_blind_tables,
)


def test_blind_table_audit_accepts_observations_and_rejects_truth(tmp_path: Path) -> None:
    case = tmp_path / "cases" / "locked_test_1"
    case.mkdir(parents=True)
    pd.DataFrame([{"sample_id": "a", "Ca": 1.0}]).to_csv(
        case / "blind_observations.csv", index=False
    )
    report = audit_blind_tables(tmp_path)
    assert report["status"] == "PASS"
    assert report["n_blind_tables"] == 1

    pd.DataFrame([{"sample_id": "a", "true_age_years": 4.0}]).to_csv(
        case / "blind_observations.csv", index=False
    )
    with pytest.raises(ValueError, match="truth fields"):
        audit_blind_tables(tmp_path)
