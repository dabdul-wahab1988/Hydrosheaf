from __future__ import annotations

from pathlib import Path
import sys

import pandas as pd
import pytest


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "M6" / "m6_field_transfer_benchmark" / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

import m6_common as m6  # noqa: E402
from run_m6_q1 import ANALYSIS_STEPS, R_STEPS  # noqa: E402
from run_m6_refined_cross_sectional import SOURCES  # noqa: E402


def test_legacy_field_loaders_are_disabled() -> None:
    with pytest.raises(RuntimeError, match="old M6 field-transfer experiment is disabled"):
        m6.load_all()
    with pytest.raises(RuntimeError, match="legacy seasonal field branch is retired"):
        m6.load_northern_ghana()


def test_q1_runner_contains_every_authoritative_stage() -> None:
    assert ANALYSIS_STEPS == ("run_m6_refined_cross_sectional.py",)
    assert R_STEPS == ()
    assert set(SOURCES) == {"central_region", "upper_east_region"}


def test_cbe_helper_abstains_on_incomplete_chemistry() -> None:
    assert pd.isna(m6.charge_balance_error({"Ca": 1.0, "Mg": 2.0}))
