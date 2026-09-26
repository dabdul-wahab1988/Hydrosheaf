"""Regression tests for the output-reproduction contract."""

from __future__ import annotations

import struct
import zlib

from hydrosheaf.reproducibility import (
    DEFAULT_EPOCH,
    apply_deterministic_environment,
    compare_trees,
    reproducible_now_iso,
)
from scripts.reproduce_outputs import _replace_png_text_chunks, _scoped_manifest


def _png_chunk(kind: bytes, payload: bytes) -> bytes:
    return (
        struct.pack(">I", len(payload))
        + kind
        + payload
        + struct.pack(">I", zlib.crc32(kind + payload) & 0xFFFFFFFF)
    )


def test_reproducible_environment_freezes_timestamp(monkeypatch):
    monkeypatch.delenv("HYDROSHEAF_REPRODUCIBLE", raising=False)
    environment = apply_deterministic_environment(seed=31415, epoch=DEFAULT_EPOCH)

    assert environment["HYDROSHEAF_REPRODUCIBLE_SEED"] == "31415"
    assert reproducible_now_iso() == DEFAULT_EPOCH


def test_png_metadata_alignment_preserves_pixel_payload():
    signature = b"\x89PNG\r\n\x1a\n"
    ihdr = _png_chunk(b"IHDR", b"header")
    old_text = _png_chunk(b"tEXt", b"Software\x00old renderer")
    idat = _png_chunk(b"IDAT", b"pixel bytes")
    end = _png_chunk(b"IEND", b"")
    candidate = signature + ihdr + old_text + idat + end

    aligned = _replace_png_text_chunks(candidate, [b"Software\x00historical renderer"])

    assert b"Software\x00historical renderer" in aligned
    assert b"Software\x00old renderer" not in aligned
    assert idat in aligned


def test_scoped_manifest_excludes_known_unresolved_families():
    manifest = {
        "schema": "test",
        "files": [
            {"path": "chapter4/tables/a.csv", "bytes": 1},
            {"path": "ttd_graph_virtual_benchmark_v1/programme/a.json", "bytes": 2},
            {"path": "objective3_closure_2026-09-15/a.csv", "bytes": 3},
        ],
    }

    scoped = _scoped_manifest(
        manifest,
        include_families={"chapter4", "ttd_graph_virtual_benchmark_v1"},
        exclude_families={"ttd_graph_virtual_benchmark_v1/programme"},
    )

    assert [item["path"] for item in scoped["files"]] == ["chapter4/tables/a.csv"]


def test_identical_manifests_are_exact():
    manifest = {"files": [{"path": "a.csv", "bytes": 1, "sha256": "x"}]}
    comparison = compare_trees(manifest, manifest)

    assert comparison["exact"] is True
    assert comparison["changed"] == []
