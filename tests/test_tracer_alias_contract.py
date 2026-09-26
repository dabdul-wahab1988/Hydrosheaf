from datetime import datetime
from pathlib import Path

import numpy as np
import pytest

from hydrosheaf.nuclear.tracer_inputs import (
    load_tracer_histories_csv,
    normalize_tracer_key,
)
from hydrosheaf.nuclear.tracer_registry import (
    build_default_tracer_registry,
    canonicalize_tracer_alias,
)
from hydrosheaf.temporal import TemporalNode, TimeSeriesSample
from hydrosheaf.temporal.residence_time import (
    _extract_tracer_series,
    _parse_tracer_candidates,
)
from hydrosheaf.temporal.time_series import load_time_series_csv


@pytest.mark.parametrize(
    ("alias", "canonical"),
    [
        ("d18O", "d18O"),
        ("18O", "d18O"),
        ("D18O", "d18O"),
        ("delta18O", "d18O"),
        ("Δ18O", "d18O"),
        ("δ18O", "d18O"),
        ("d2H", "d2H"),
        ("2H", "d2H"),
        ("D2H", "d2H"),
        ("delta2H", "d2H"),
        ("Δ2H", "d2H"),
        ("δ2H", "d2H"),
        ("H-3", "3H"),
        ("Tritium", "3H"),
        ("C-14", "14C"),
        ("Radiocarbon", "14C"),
        ("SF_6", "SF6"),
        ("CFC-12", "CFC12"),
        ("Krypton-85", "85Kr"),
    ],
)
def test_canonical_alias_function_and_nuclear_wrapper(alias: str, canonical: str):
    assert canonicalize_tracer_alias(alias) == canonical
    assert normalize_tracer_key(alias) == canonical


def test_registry_exposes_canonical_ids_and_aliases():
    registry = build_default_tracer_registry()

    for canonical in ("d18O", "d2H", "3H", "SF6", "14C"):
        spec = registry[canonical]
        assert spec.metadata["canonical_id"] == canonical
        assert canonical in spec.aliases
        assert tuple(spec.metadata["aliases"]) == spec.aliases
        for alias in spec.aliases:
            assert canonicalize_tracer_alias(alias) == canonical


def test_tracer_history_csv_canonicalizes_long_and_wide_aliases(tmp_path: Path):
    long_path = tmp_path / "long.csv"
    long_path.write_text(
        "\n".join(
            [
                "tracer,year,value,sigma",
                "δ18O,2000,-8.0,0.1",
                "delta2H,2000,-55.0,0.8",
                "H-3,2000,5.0,0.2",
                "Radiocarbon,2000,85.0,1.0",
                "SF_6,2000,4.0,0.1",
            ]
        ),
        encoding="utf-8",
    )
    long_histories = load_tracer_histories_csv(long_path)
    assert set(long_histories) == {"d18O", "d2H", "3H", "14C", "SF6"}

    wide_path = tmp_path / "wide.csv"
    wide_path.write_text(
        "\n".join(
            [
                "year,Δ18O,δ2H,Krypton-85",
                "2000,-8.0,-55.0,12.0",
                "2010,-7.5,-52.0,10.0",
            ]
        ),
        encoding="utf-8",
    )
    wide_histories = load_tracer_histories_csv(wide_path)
    assert set(wide_histories) == {"d18O", "d2H", "85Kr"}
    assert wide_histories["d18O"].interpolate([2005.0])[0][0] == pytest.approx(-7.75)


def test_time_series_csv_canonicalizes_stable_and_nuclear_aliases(tmp_path: Path):
    path = tmp_path / "temporal.csv"
    path.write_text(
        "\n".join(
            [
                "node_id,sample_id,timestamp,Cl,δ18O,delta2H,SF6,H-3,Radiocarbon,Krypton-85",
                "A,a1,2020-01-02,10,-8,-55,4,5,85,12",
                "A,a2,2020-01-01,11,-7,-52,5,6,86,11",
            ]
        ),
        encoding="utf-8",
    )

    nodes = load_time_series_csv(path, ion_order=["Cl"], node_id_column="node_id")

    assert list(nodes) == ["A"]
    samples = nodes["A"].samples
    assert [sample.sample_id for sample in samples] == ["a2", "a1"]
    assert samples[0].isotopes == {
        "d18O": -7.0,
        "d2H": -52.0,
        "SF6": 5.0,
        "3H": 6.0,
        "14C": 86.0,
        "85Kr": 11.0,
    }


def test_residence_time_candidate_and_series_extraction_canonicalize_aliases():
    assert _parse_tracer_candidates("Cl, δ18O, 2H, H-3, Radiocarbon, SF_6") == [
        "Cl",
        "d18O",
        "d2H",
        "3H",
        "14C",
        "SF6",
    ]

    node = TemporalNode(
        node_id="A",
        samples=[
            TimeSeriesSample(
                sample_id="a1",
                node_id="A",
                timestamp=datetime(2020, 1, 1),
                concentrations=[1.0],
                isotopes={"δ18O": -8.0, "Δ2H": -55.0, "H-3": 5.0},
            ),
            TimeSeriesSample(
                sample_id="a2",
                node_id="A",
                timestamp=datetime(2020, 1, 2),
                concentrations=[2.0],
                isotopes={"δ18O": -7.0, "Δ2H": -52.0, "H-3": 6.0},
            ),
        ],
    )

    np.testing.assert_allclose(_extract_tracer_series(node, "18O", None), [-8.0, -7.0])
    np.testing.assert_allclose(_extract_tracer_series(node, "d2H", None), [-55.0, -52.0])
    np.testing.assert_allclose(_extract_tracer_series(node, "tritium", None), [5.0, 6.0])
    np.testing.assert_allclose(
        _extract_tracer_series(node, "d18O", ["δ18O"]),
        [1.0, 2.0],
    )
