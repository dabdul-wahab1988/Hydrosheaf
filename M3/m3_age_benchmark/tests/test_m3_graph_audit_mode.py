"""Tests for the opt-in graph compatibility audit in the M3 identified-TTD run.

The audit is Version 1: compatibility/obstruction only.  These tests pin the
two properties that matter scientifically -- the default path is untouched, and
supplying edges never narrows a local identified set.
"""

from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
import pytest
import yaml

BENCHMARK_DIR = Path(__file__).resolve().parents[1]
SCRIPT_DIR = BENCHMARK_DIR / "scripts"
REPO_ROOT = BENCHMARK_DIR.parents[1]
for _path in (str(REPO_ROOT), str(SCRIPT_DIR)):
    if _path not in sys.path:
        sys.path.insert(0, _path)

import run_m3_identified_ttd_benchmark as runner  # noqa: E402

from hydrosheaf.nuclear.ttd_identified import (  # noqa: E402
    AgeFunctional,
    TracerConstraint,
    solve_ttd_identified_set,
)

# A site/tracer fold known to yield a feasible local report inside the first
# rows of the national benchmark dataset.
E2E_SITE = "BISCPAS1-15"
E2E_TRACER = "SF6"
E2E_MAX_ROWS = 16


def _fixed_report(mass):
    """Two-bin local report pinned to ``mass`` (mirrors tests/test_ttd_graph.py)."""
    return solve_ttd_identified_set(
        [0.0, 10.0],
        [TracerConstraint("bin-0", [1.0, 0.0], mass[0], 1.0)],
        [AgeFunctional("young", [1.0, 0.0], 1.0e-10)],
        sigma_multiplier=0.0,
    )


def _edge_document(**overrides):
    edge = {
        "edge_id": "e1",
        "source_site_id": "A",
        "target_site_id": "B",
        "transport_generator": "identity",
        "provenance_tier": "independent",
        "source": "synthetic test fixture",
    }
    edge.update(overrides)
    return {"schema_version": "1.0", "edges": [edge]}


def _write(path: Path, document) -> Path:
    if path.suffix == ".json":
        path.write_text(json.dumps(document), encoding="utf-8")
    else:
        path.write_text(yaml.safe_dump(document), encoding="utf-8")
    return path


# --------------------------------------------------------------------------
# Edge-file loading and validation
# --------------------------------------------------------------------------


@pytest.mark.parametrize("suffix", [".yaml", ".json"])
def test_load_graph_edges_accepts_yaml_and_json(tmp_path, suffix):
    path = _write(tmp_path / f"edges{suffix}", _edge_document())
    edges = runner.load_graph_edges(path, 2)

    assert len(edges) == 1
    edge = edges[0]
    assert (edge.edge_id, edge.source_site_id, edge.target_site_id) == ("e1", "A", "B")
    assert np.allclose(np.asarray(edge.transport.matrix), np.eye(2))
    assert edge.transport.provenance_tier == "independent"


def test_load_graph_edges_accepts_explicit_matrix(tmp_path):
    document = _edge_document(
        transport_generator=None, transport_matrix=[[1.0, 0.5], [0.0, 0.5]]
    )
    path = _write(tmp_path / "edges.yaml", document)
    edges = runner.load_graph_edges(path, 2)

    assert np.allclose(np.asarray(edges[0].transport.matrix), [[1.0, 0.5], [0.0, 0.5]])


def test_load_graph_edges_rejects_columns_not_summing_to_one(tmp_path):
    """Validation is delegated to MassTransportMap's own constructor invariants."""
    document = _edge_document(
        transport_generator=None, transport_matrix=[[0.5, 0.0], [0.0, 0.5]]
    )
    path = _write(tmp_path / "edges.yaml", document)
    with pytest.raises(ValueError, match="columns must sum to one"):
        runner.load_graph_edges(path, 2)


def test_load_graph_edges_rejects_negative_entries(tmp_path):
    document = _edge_document(
        transport_generator=None, transport_matrix=[[1.1, 0.0], [-0.1, 1.0]]
    )
    path = _write(tmp_path / "edges.yaml", document)
    with pytest.raises(ValueError, match="non-negative"):
        runner.load_graph_edges(path, 2)


def test_load_graph_edges_rejects_wrong_shape(tmp_path):
    document = _edge_document(
        transport_generator=None, transport_matrix=[[1.0, 0.0], [0.0, 1.0]]
    )
    path = _write(tmp_path / "edges.yaml", document)
    with pytest.raises(ValueError, match="n_target_bins, n_source_bins"):
        runner.load_graph_edges(path, 3)


def test_load_graph_edges_rejects_ambiguous_and_unknown_transport(tmp_path):
    both = _edge_document(transport_matrix=[[1.0, 0.0], [0.0, 1.0]])
    with pytest.raises(ValueError, match="exactly one"):
        runner.load_graph_edges(_write(tmp_path / "both.yaml", both), 2)

    neither = _edge_document(transport_generator=None)
    with pytest.raises(ValueError, match="exactly one"):
        runner.load_graph_edges(_write(tmp_path / "neither.yaml", neither), 2)

    unknown = _edge_document(transport_generator="upstream_mixing")
    with pytest.raises(ValueError, match="Unsupported transport_generator"):
        runner.load_graph_edges(_write(tmp_path / "unknown.yaml", unknown), 2)


def test_load_graph_edges_rejects_missing_fields_and_duplicates(tmp_path):
    missing = _edge_document(provenance_tier="")
    with pytest.raises(ValueError, match="missing"):
        runner.load_graph_edges(_write(tmp_path / "missing.yaml", missing), 2)

    duplicate = _edge_document()
    duplicate["edges"] = [duplicate["edges"][0], dict(duplicate["edges"][0])]
    with pytest.raises(ValueError, match="Duplicate edge_id"):
        runner.load_graph_edges(_write(tmp_path / "dupe.yaml", duplicate), 2)

    empty = {"schema_version": "1.0", "edges": []}
    with pytest.raises(ValueError, match="no edges"):
        runner.load_graph_edges(_write(tmp_path / "empty.yaml", empty), 2)


# --------------------------------------------------------------------------
# Audit behaviour over frozen local reports
# --------------------------------------------------------------------------


def test_identity_map_on_a_site_against_itself_is_compatible(tmp_path):
    path = _write(
        tmp_path / "edges.yaml",
        _edge_document(edge_id="self", source_site_id="A", target_site_id="A"),
    )
    edges = runner.load_graph_edges(path, 2)
    report = _fixed_report([0.25, 0.75])

    records = runner.run_graph_compatibility_audits(
        edges, {("A", "SF6"): report}, {}, ["SF6"]
    )

    assert len(records) == 1
    assert records[0]["status"] == "COMPATIBLE"
    assert records[0]["minimum_l1_gap"] == pytest.approx(0.0, abs=1e-7)
    assert records[0]["tightening_authorized"] is False
    assert records[0]["graph_mode"] == "compatibility_audit_only"


def test_missing_endpoint_report_yields_abstain_with_reason(tmp_path):
    path = _write(tmp_path / "edges.yaml", _edge_document())
    edges = runner.load_graph_edges(path, 2)
    reports = {("A", "SF6"): _fixed_report([0.25, 0.75])}
    fold_reasons = {("B", "SF6"): "NOT_ELIGIBLE:missing_sample_year"}

    records = runner.run_graph_compatibility_audits(
        edges, reports, fold_reasons, ["SF6"]
    )

    assert len(records) == 1, "abstentions must be preserved, never dropped"
    record = records[0]
    assert record["status"] == "ABSTAIN"
    assert "endpoint_report_unavailable" in record["reason"]
    assert "target=B(NOT_ELIGIBLE:missing_sample_year)" in record["reason"]
    assert np.isnan(record["minimum_l1_gap"])
    assert record["tightening_authorized"] is False


def test_absent_site_abstains_rather_than_being_skipped(tmp_path):
    path = _write(tmp_path / "edges.yaml", _edge_document())
    edges = runner.load_graph_edges(path, 2)

    records = runner.run_graph_compatibility_audits(edges, {}, {}, ["SF6", "3H"])

    assert len(records) == 2
    assert {record["status"] for record in records} == {"ABSTAIN"}
    assert all("site_not_in_run" in record["reason"] for record in records)


def test_incompatible_edge_is_reported_without_tightening(tmp_path):
    path = _write(tmp_path / "edges.yaml", _edge_document())
    edges = runner.load_graph_edges(path, 2)
    source = _fixed_report([1.0, 0.0])
    target = _fixed_report([0.0, 1.0])
    source_bound_before = source.bounds["young"].lower
    target_bound_before = target.bounds["young"].lower

    records = runner.run_graph_compatibility_audits(
        edges, {("A", "SF6"): source, ("B", "SF6"): target}, {}, ["SF6"]
    )

    assert records[0]["status"] == "INCOMPATIBLE"
    assert records[0]["minimum_l1_gap"] == pytest.approx(2.0)
    assert records[0]["tightening_authorized"] is False
    assert source.bounds["young"].lower == source_bound_before
    assert target.bounds["young"].lower == target_bound_before


def test_report_sink_is_write_only_and_does_not_change_the_record():
    """The graph side channel must not perturb the local fold record."""
    import pandas as pd

    import run_m3_usgs_benchmark as usgs

    protocol = runner.load_protocol()
    data = usgs.load_benchmark_dataset(sources="national").head(E2E_MAX_ROWS)
    row = dict(data.iloc[E2E_MAX_ROWS - 3])

    without = runner.evaluate_held_out_row_safely(row, E2E_TRACER, protocol)
    sink: dict = {}
    with_sink = runner.evaluate_held_out_row_safely(
        row, E2E_TRACER, protocol, report_sink=sink
    )

    assert pd.Series(without).equals(pd.Series(with_sink))
    assert without["graph_mode"] == "disabled_local_only"


# --------------------------------------------------------------------------
# End-to-end: the default path is untouched and edges never alter local bounds
# --------------------------------------------------------------------------


@pytest.fixture(scope="module")
def default_run(tmp_path_factory):
    """One shared baseline run with no --graph-edges supplied."""
    directory = tmp_path_factory.mktemp("default_run")
    output = directory / "run.csv"
    exit_code = runner.main(
        [
            "--max-rows",
            str(E2E_MAX_ROWS),
            "--withhold-tracer",
            E2E_TRACER,
            "--output",
            str(output),
        ]
    )
    assert exit_code == 0
    manifest = json.loads(
        output.with_name(f"{output.stem}_manifest.json").read_text(encoding="utf-8")
    )
    return output, manifest


def test_default_path_is_local_only_and_writes_no_edge_output(default_run):
    output, manifest = default_run

    assert manifest["graph_mode"] == "disabled_local_only"
    assert "graph_audit" not in manifest
    assert not output.with_name(f"{output.stem}_graph_edges.csv").exists()

    import pandas as pd

    results = pd.read_csv(output)
    assert set(results["graph_mode"].unique()) == {"disabled_local_only"}
    assert not [column for column in results.columns if column.startswith("edge_")]


def test_graph_output_requires_graph_edges(tmp_path):
    with pytest.raises(SystemExit):
        runner.main(
            [
                "--max-rows",
                "1",
                "--output",
                str(tmp_path / "run.csv"),
                "--graph-output",
                str(tmp_path / "edges.csv"),
            ]
        )


def test_supplying_edges_leaves_every_local_bound_unchanged(default_run, tmp_path):
    """The whole point of Version 1: an edge audits, it never tightens."""
    baseline_output, baseline_manifest = default_run
    edge_file = _write(
        tmp_path / "edges.yaml",
        _edge_document(
            edge_id="self-identity",
            source_site_id=E2E_SITE,
            target_site_id=E2E_SITE,
        ),
    )
    output = tmp_path / "run.csv"
    exit_code = runner.main(
        [
            "--max-rows",
            str(E2E_MAX_ROWS),
            "--withhold-tracer",
            E2E_TRACER,
            "--output",
            str(output),
            "--graph-edges",
            str(edge_file),
        ]
    )
    assert exit_code == 0

    # Local per-site fold output must be byte-for-byte the baseline.
    assert output.read_bytes() == baseline_output.read_bytes()

    manifest = json.loads(
        output.with_name(f"{output.stem}_manifest.json").read_text(encoding="utf-8")
    )
    assert manifest["graph_mode"] == "compatibility_audit_only"
    assert manifest["output_sha256"] == baseline_manifest["output_sha256"]
    assert manifest["summary"] == baseline_manifest["summary"]
    assert manifest["graph_audit"]["tightening_authorized"] is False
    assert manifest["graph_audit"]["conditioning_performed"] is False
    assert manifest["graph_audit"]["n_declared_edges"] == 1


def test_edge_output_schema_and_real_site_self_edge_is_compatible(
    default_run, tmp_path
):
    edge_file = _write(
        tmp_path / "edges.yaml",
        _edge_document(
            edge_id="self-identity",
            source_site_id=E2E_SITE,
            target_site_id=E2E_SITE,
        ),
    )
    output = tmp_path / "run.csv"
    graph_output = tmp_path / "custom_edges.csv"
    assert (
        runner.main(
            [
                "--max-rows",
                str(E2E_MAX_ROWS),
                "--withhold-tracer",
                E2E_TRACER,
                "--output",
                str(output),
                "--graph-edges",
                str(edge_file),
                "--graph-output",
                str(graph_output),
            ]
        )
        == 0
    )

    import pandas as pd

    edges = pd.read_csv(graph_output)
    required = {
        "edge_id",
        "source_site_id",
        "target_site_id",
        "status",
        "minimum_l1_gap",
        "maximum_bin_gap",
        "provenance_tier",
        "transport_source",
        "tightening_authorized",
        "message",
    }
    assert required.issubset(edges.columns)
    assert len(edges) == 1
    assert not edges["tightening_authorized"].any()

    record = edges.iloc[0]
    assert record["source_site_id"] == E2E_SITE
    assert record["status"] == "COMPATIBLE"
    assert record["minimum_l1_gap"] == pytest.approx(0.0, abs=1e-6)
