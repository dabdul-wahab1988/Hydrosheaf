"""Regression tests for the Table 4.13 held-out compatibility aggregation."""

from __future__ import annotations

import numpy as np
import pandas as pd

from scripts.analysis import generate_m8_ttd_tables_figures as m8


def test_table_4_13_counts_compatible_and_incompatible_rows(tmp_path, monkeypatch):
    """The compatibility denominator must not be pre-filtered to COMPATIBLE."""

    source = tmp_path / "m3_development.csv"
    rows = [
        {
            "held_out_tracer": "3H",
            "held_out_status": "COMPATIBLE",
            "held_out_observed": 1.0,
            "held_out_sigma": 0.1,
            "prediction_lower": 0.5,
            "prediction_upper": 1.5,
            "prediction_width": 1.0,
            "held_out_compatible": True,
        },
        {
            "held_out_tracer": "3H",
            "held_out_status": "INCOMPATIBLE",
            "held_out_observed": 3.0,
            "held_out_sigma": 0.1,
            "prediction_lower": 0.5,
            "prediction_upper": 1.5,
            "prediction_width": 1.0,
            "held_out_compatible": False,
        },
        {
            "held_out_tracer": "3H",
            "held_out_status": "ABSTAIN",
            "held_out_observed": 2.0,
            "held_out_sigma": 0.1,
            "prediction_lower": np.nan,
            "prediction_upper": np.nan,
            "prediction_width": np.nan,
            "held_out_compatible": np.nan,
        },
        {
            "held_out_tracer": "3H",
            "held_out_status": np.nan,
            "held_out_observed": np.nan,
            "held_out_sigma": np.nan,
            "prediction_lower": np.nan,
            "prediction_upper": np.nan,
            "prediction_width": np.nan,
            "held_out_compatible": np.nan,
        },
    ]
    pd.DataFrame(rows).to_csv(source, index=False)

    monkeypatch.setattr(m8, "M3_DEV_CSV", source)
    monkeypatch.setattr(m8, "TABLES_DIR", tmp_path / "tables")
    m8.TABLES_DIR.mkdir()

    m8.generate_table_4_13()

    result = pd.read_csv(
        m8.TABLES_DIR / "table_4_13_heldout_tracer_predictive_compatibility.csv"
    )
    row = result.iloc[0]
    assert row["Evaluated Sites (n)"] == 2
    assert row["Empirical Compatibility (%)"] == "50.0% (1/2)"
    assert row["Model Compatibility Gate"] == "WARNING"

    markdown = (
        m8.TABLES_DIR / "table_4_13_heldout_tracer_predictive_compatibility.md"
    ).read_text(encoding="utf-8")
    assert "3H: 2 evaluable (1 compatible; 1 incompatible), 1 abstain" in markdown
    assert "1 of 2 evaluable predictions were compatible (50.0%)" in markdown
