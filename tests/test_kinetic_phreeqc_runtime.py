from __future__ import annotations

import numpy as np
import pytest

from hydrosheaf.config import Config
from hydrosheaf.reactive_transport import KineticParameters
from hydrosheaf.reactive_transport.kinetic_phreeqc import (
    build_kinetic_block,
    run_phreeqc_kinetic,
)


def test_phreeqpython_executes_the_pending_kinetic_step():
    pytest.importorskip("phreeqpython")
    config = Config()
    block = build_kinetic_block(
        reaction_labels=["calcite"],
        extents=[1.0],
        residence_time_days=10.0,
        kinetic_params={
            "calcite": KineticParameters(
                reaction_name="calcite",
                rate_constant=1e-10,
                surface_area=0.1,
            )
        },
    )
    result = run_phreeqc_kinetic(
        initial_solution={
            "pH": 6.5,
            "Ca": 0.05,
            "Mg": 0.02,
            "Na": 0.10,
            "K": 0.01,
            "Cl": 0.10,
            "SO4": 0.01,
            "NO3": 0.01,
            "HCO3": 0.10,
            "F": 0.0,
            "Fe": 0.0,
            "PO4": 0.0,
        },
        kinetics_block=block,
        residence_time_days=10.0,
        config=config,
    )
    assert result["success"], result["error_message"]
    values = np.asarray(result["final_composition"], dtype=float)
    assert values.size == len(config.ion_order)
    assert np.all(np.isfinite(values))
