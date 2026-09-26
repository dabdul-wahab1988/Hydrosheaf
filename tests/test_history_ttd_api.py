"""Public import contracts for the recharge-history-aware TTD API."""

from __future__ import annotations

import pytest


def _history_ttd_module():
    return pytest.importorskip(
        "hydrosheaf.nuclear.history_ttd",
        reason=(
            "history-aware TTD API depends on worker-A module "
            "hydrosheaf.nuclear.history_ttd"
        ),
    )


def test_history_ttd_is_public_from_api_package_and_nuclear() -> None:
    history_ttd = _history_ttd_module()

    from hydrosheaf import fit_history_ttd as package_fit_history_ttd
    from hydrosheaf.api import fit_history_ttd as api_fit_history_ttd
    from hydrosheaf.nuclear import (
        HistoryTracerObservation as nuclear_history_observation,
        HistoryTTDResult as nuclear_history_result,
        fit_history_ttd as nuclear_fit_history_ttd,
    )

    assert api_fit_history_ttd is not None
    assert package_fit_history_ttd is api_fit_history_ttd
    assert nuclear_fit_history_ttd is history_ttd.fit_history_ttd
    assert nuclear_history_observation is history_ttd.HistoryTracerObservation
    assert nuclear_history_result is history_ttd.HistoryTTDResult

    import hydrosheaf.nuclear as nuclear

    assert {
        "HistoryTracerObservation",
        "HistoryTTDResult",
        "fit_history_ttd",
    }.issubset(nuclear.__all__)


def test_history_ttd_public_wrapper_runs_minimal_contract() -> None:
    history_ttd = _history_ttd_module()

    from hydrosheaf import fit_history_ttd
    from hydrosheaf.nuclear.input_history import InputHistory
    from hydrosheaf.nuclear.history_ttd import HistoryTracerObservation

    observation = HistoryTracerObservation(
        tracer="d18O",
        value=0.0,
        sigma=0.1,
    )
    result = fit_history_ttd(
        observations=[observation],
        sample_year=2024.0,
        age_grid_years=[0.0, 1.0, 2.0],
        source_histories={
            "d18O": InputHistory(
                dates=[2022.0, 2023.0, 2024.0],
                values=[0.0, 0.0, 0.0],
            ),
        },
    )

    assert isinstance(result, history_ttd.HistoryTTDResult)
    assert result.provenance["claim_scope"] == "conditional_inference"
    assert result.provenance["field_validation_status"] == "not_performed"


def test_importing_package_does_not_eagerly_import_history_stack() -> None:
    import subprocess
    import sys

    subprocess.run(
        [
            sys.executable,
            "-c",
            (
                "import sys, hydrosheaf; "
                "assert 'fit_history_ttd' in hydrosheaf.__all__; "
                "assert 'hydrosheaf.nuclear.history_ttd' not in sys.modules"
            ),
        ],
        check=True,
    )
