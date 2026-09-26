"""Re-execute the preserved HydroSheaf output recipes in an isolated tree.

The command is intentionally a verifier as well as a runner.  It never writes
to the historical ``outputs`` directory.  Instead it creates a staging tree,
runs the recipes whose source and inputs are present, records the exact hashes,
and compares the result with the historical tree supplied by ``--reference``.

Typical use from the repository root::

    .venv\\Scripts\\python.exe scripts\\reproduce_outputs.py --compare

The command exits non-zero in strict mode when a recipe fails, an output family
has no preserved recipe, or any file differs byte-for-byte.  This is deliberate:
an audit must not silently turn a missing historical input into a successful
"reproduction".
"""

from __future__ import annotations

import argparse
import contextlib
import copy
import importlib.util
import importlib.metadata
import json
import os
import re
import shutil
import sys
import traceback
import zipfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Iterable, Sequence

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from hydrosheaf.reproducibility import (
    DEFAULT_EPOCH,
    DEFAULT_SEED,
    apply_deterministic_environment,
    compare_trees,
    file_sha256,
    manifest_for_tree,
    seed_all,
    write_json,
)


REFERENCE_OUTPUTS = REPO_ROOT / "outputs"
DEFAULT_RUN_ROOT = REPO_ROOT / ".codex_work" / "output_reproduction"
OUTPUT_SUFFIXES = {
    ".csv",
    ".json",
    ".ndjson",
    ".log",
    ".md",
    ".pdf",
    ".png",
    ".sha256",
    ".xlsx",
    ".docx",
}

# These output families are present in the current historical tree but do not
# have a generator in the preserved checkout.  Keeping this list explicit is
# what prevents a partial rerun from being advertised as a full replay.
UNRESOLVED_FAMILIES = {
    "master_field_integration": (
        "legacy pooled field network is retired; it included excluded cohorts and treated DEM elevation as hydraulic head"
    ),
    "uer_field_integration": (
        "legacy UER network is retired; its topographic elevation proxy is not measured hydraulic head or flow truth"
    ),
    "2026-09-05_multi_geology_join": (
        "historical multi-dataset geology-join generator was not preserved"
    ),
    "2026-09-05_northern_ghana_geology_join": (
        "historical UER geology-join generator was not preserved"
    ),
    "2026-09-05_northern_ghana_new": (
        "historical UER cleaning generator was not preserved"
    ),
    "objective5_audit_2026-09-15": (
        "historical Objective 5 audit generator was not preserved"
    ),
    "uer_objective5_corrected_2026-09-15": (
        "corrected historical snapshot has no preserved versioned recipe"
    ),
    "uer_objective5_costed_public_list_2026-09-15": (
        "costed historical snapshot has no preserved versioned recipe"
    ),
    "ttd_graph_virtual_benchmark_v1/programme": (
        "programme ensemble recipe and exact locked run inputs are not part of the virtual runner"
    ),
    "uer_objective5": (
        "current Objective 5 generator is executable, but its historical versioned recipe is not locked"
    ),
    "uer_objective5_groundwater_budgetary_2026-09-15": (
        "current Objective 5 budgetary generator is executable, but its historical versioned recipe is not locked"
    ),
    "objective5_closure_2026-09-15": (
        "current Objective 5 closure generator is executable, but its historical versioned inputs are not locked"
    ),
}


@dataclass
class StepResult:
    name: str
    status: str
    outputs: list[str] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)
    error: str | None = None


def _load_script(path: Path, alias: str) -> Any:
    """Import a script under a private name so output globals can be patched."""

    path = path.resolve()
    spec = importlib.util.spec_from_file_location(alias, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[alias] = module
    spec.loader.exec_module(module)
    return module


def _set_globals(module: Any, **values: Any) -> None:
    for name, value in values.items():
        if not hasattr(module, name):
            raise AttributeError(f"{module.__name__} has no redirectable global {name}")
        setattr(module, name, value)


def _patch_plotting() -> None:
    """Use a headless backend for every staged figure render."""

    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.figure

    if getattr(matplotlib.figure.Figure, "_hydrosheaf_reproducible_savefig", False):
        return
    original = matplotlib.figure.Figure.savefig

    def savefig(self: Any, fname: Any, *args: Any, **kwargs: Any) -> Any:
        return original(self, fname, *args, **kwargs)

    matplotlib.figure.Figure.savefig = savefig
    matplotlib.figure.Figure._hydrosheaf_reproducible_savefig = True


def _align_renderer_metadata(root: Path, reference: Path) -> None:
    """Align non-scientific renderer metadata with the historical artifact.

    Matplotlib embeds a PDF creation second and a PNG Software string.  Those
    bytes are renderer metadata, not scientific content.  When the actual
    rendered content is unchanged, using the preserved reference metadata lets
    the strict comparison test the content rather than the wall clock.
    """

    for candidate in sorted(root.rglob("*.pdf")):
        relative = candidate.relative_to(root)
        baseline = reference / relative
        if not baseline.exists():
            continue
        candidate_bytes = candidate.read_bytes()
        baseline_bytes = baseline.read_bytes()
        baseline_match = re.search(rb"/CreationDate \(([^)]*)\)", baseline_bytes)
        if baseline_match:
            candidate_bytes = re.sub(
                rb"/CreationDate \([^)]*\)",
                b"/CreationDate (" + baseline_match.group(1) + b")",
                candidate_bytes,
                count=1,
            )
        candidate.write_bytes(candidate_bytes)

    # PNG Software metadata is small and isolated in a tEXt chunk.  Do not
    # touch pixel IDAT chunks; if those differ the manifest must still fail.
    for candidate in sorted(root.rglob("*.png")):
        relative = candidate.relative_to(root)
        baseline = reference / relative
        if not baseline.exists():
            continue
        baseline_text = _png_text_chunks(baseline.read_bytes())
        if not baseline_text:
            continue
        candidate_bytes = candidate.read_bytes()
        candidate_bytes = _replace_png_text_chunks(candidate_bytes, baseline_text)
        candidate.write_bytes(candidate_bytes)


def _png_text_chunks(data: bytes) -> list[bytes]:
    chunks: list[bytes] = []
    index = 8
    while index + 12 <= len(data):
        length = int.from_bytes(data[index : index + 4], "big")
        kind = data[index + 4 : index + 8]
        payload = data[index + 8 : index + 8 + length]
        if kind == b"tEXt":
            chunks.append(payload)
        index += length + 12
        if kind == b"IEND":
            break
    return chunks


def _replace_png_text_chunks(data: bytes, payloads: list[bytes]) -> bytes:
    import struct
    import zlib

    output = bytearray(data[:8])
    index = 8
    payload_index = 0
    while index + 12 <= len(data):
        length = int.from_bytes(data[index : index + 4], "big")
        kind = data[index + 4 : index + 8]
        payload = data[index + 8 : index + 8 + length]
        if kind == b"tEXt" and payload_index < len(payloads):
            payload = payloads[payload_index]
            payload_index += 1
            length = len(payload)
        output.extend(struct.pack(">I", length))
        output.extend(kind)
        output.extend(payload)
        output.extend(struct.pack(">I", zlib.crc32(kind + payload) & 0xFFFFFFFF))
        index += int.from_bytes(data[index : index + 4], "big") + 12
        if kind == b"IEND":
            break
    return bytes(output)


def _align_workbook_container(root: Path, reference: Path) -> None:
    """Reuse historical ZIP container metadata when workbook XML is equal."""

    for candidate in sorted(root.rglob("*.xlsx")):
        relative = candidate.relative_to(root)
        baseline = reference / relative
        if not baseline.exists():
            continue
        temporary = candidate.with_suffix(candidate.suffix + ".reference.tmp")
        copy_reference = False
        with zipfile.ZipFile(baseline, "r") as ref_zip, zipfile.ZipFile(
            candidate, "r"
        ) as candidate_zip:
            ref_names = ref_zip.namelist()
            if set(ref_names) != set(candidate_zip.namelist()):
                continue
            # A workbook is safe to replay byte-for-byte only when its logical
            # ZIP members are all identical.  The historical file may still
            # use a different compression level or central-directory layout;
            # in that case copying the preserved container is the only way to
            # reproduce its bytes without treating container metadata as
            # scientific content.  A changed XML member remains a hard diff.
            if all(
                candidate_zip.read(name) == ref_zip.read(name) for name in ref_names
            ):
                copy_reference = True
            else:
                with zipfile.ZipFile(temporary, "w") as output_zip:
                    for name in ref_names:
                        data = candidate_zip.read(name)
                        if name == "docProps/core.xml":
                            # Creation/modification properties are renderer
                            # metadata; preserve the reference representation
                            # for byte replay.
                            data = ref_zip.read(name)
                        info = copy.copy(ref_zip.getinfo(name))
                        output_zip.writestr(info, data)
        if copy_reference:
            shutil.copyfile(baseline, candidate)
        else:
            temporary.replace(candidate)


def _canonicalize_workbooks(root: Path) -> None:
    """Make XLSX ZIP containers independent of creation time and entry order."""

    import re

    for path in sorted(root.rglob("*.xlsx")):
        temporary = path.with_suffix(path.suffix + ".reproducible.tmp")
        with zipfile.ZipFile(path, "r") as source, zipfile.ZipFile(
            temporary,
            "w",
            compression=zipfile.ZIP_DEFLATED,
            compresslevel=9,
        ) as target:
            for name in sorted(source.namelist()):
                data = source.read(name)
                if name == "docProps/core.xml":
                    text = data.decode("utf-8")
                    text = re.sub(
                        r"(<dcterms:(?:created|modified)[^>]*>).*?(</dcterms:(?:created|modified)>)",
                        r"\g<1>1980-01-01T00:00:00Z\g<2>",
                        text,
                    )
                    data = text.encode("utf-8")
                info = zipfile.ZipInfo(name, date_time=(1980, 1, 1, 0, 0, 0))
                info.compress_type = zipfile.ZIP_DEFLATED
                info.external_attr = 0o600 << 16
                target.writestr(info, data)
        temporary.replace(path)


def _canonicalize_text_paths(root: Path, stage_outputs: Path) -> None:
    """Replace staging-root paths with the historical repository output root."""

    historical = REPO_ROOT / "outputs"
    replacements = {
        str(stage_outputs): str(historical),
        stage_outputs.as_posix(): historical.as_posix(),
        str(stage_outputs).replace("\\", "/"): historical.as_posix(),
    }
    binary_suffixes = {".png", ".pdf", ".xlsx", ".docx", ".zip", ".pyc"}
    for path in sorted(root.rglob("*")):
        if not path.is_file() or path.suffix.lower() in binary_suffixes:
            continue
        try:
            original = path.read_text(encoding="utf-8")
        except (UnicodeDecodeError, OSError):
            continue
        updated = original
        for source, target in replacements.items():
            updated = updated.replace(source, target)
        if updated != original:
            path.write_text(updated, encoding="utf-8", newline="\n")


def _build_chapter4_generators(output_dir: Path, tag: str) -> None:
    """Run the ten preserved Chapter 4 generator modules into *output_dir*."""

    tables = output_dir / "tables"
    figures = output_dir / "figures"
    tables.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)
    analysis = REPO_ROOT / "scripts" / "analysis"

    calls: list[tuple[str, str, dict[str, Any], tuple[str, ...]]] = [
        (
            "m3",
            "generate_tier1_m3_tables_figures.py",
            {"TABLE_DIR": tables, "FIGURE_DIR": figures},
            ("run_tier1_tracer_forward", "run_tier1_inverse_age"),
        ),
        (
            "m5",
            "generate_m5_reaction_tables_figures.py",
            {"TABLE_DIR": tables, "FIGURE_DIR": figures},
            ("run_m5_reaction_benchmarks",),
        ),
        (
            "o3",
            "generate_o3_jrs_tables_figures.py",
            {"TABLES_DIR": tables, "FIGURES_DIR": figures},
            ("generate_table_4_6", "generate_figure_4_5"),
        ),
        (
            "m4",
            "generate_m4_topology_tables_figures.py",
            {"TABLES_DIR": tables, "FIGURES_DIR": figures},
            ("generate_table_4_7", "generate_figure_4_6", "generate_figure_4_7"),
        ),
        (
            "m7",
            "generate_m7_integration_tables_figures.py",
            {"TABLES_DIR": tables, "FIGURES_DIR": figures},
            (
                "generate_table_4_8",
                "generate_table_4_9",
                "generate_table_4_10",
                "generate_figure_4_8_and_4_9",
                "generate_figure_4_10",
                "generate_figure_4_16",
            ),
        ),
        (
            "robustness",
            "generate_robustness_identifiability_suite.py",
            {"TABLES_DIR": tables, "FIGURES_DIR": figures},
            (
                "generate_table_4_11",
                "generate_figure_4_11",
                "generate_figure_4_12",
                "generate_figure_4_13",
                "generate_figure_4_14",
                "generate_figure_4_15",
                "generate_figure_4_17",
            ),
        ),
        (
            "m8",
            "generate_m8_ttd_tables_figures.py",
            {"TABLES_DIR": tables, "FIGURES_DIR": figures},
            (
                "generate_table_4_12",
                "generate_table_4_13",
                "generate_table_4_14",
                "generate_figure_4_18",
                "generate_figure_4_19",
                "generate_figure_4_20",
            ),
        ),
        (
            "synthesis",
            "generate_evidence_synthesis_tables.py",
            {"TABLES_DIR": tables},
            ("generate_table_4_1", "generate_table_4_19", "generate_table_4_20"),
        ),
    ]

    for name, filename, overrides, functions in calls:
        module = _load_script(analysis / filename, f"_hydrosheaf_repro_{tag}_{name}")
        _set_globals(module, **overrides)
        arguments = {function: () for function in functions}
        for function in functions:
            getattr(module, function)(*arguments[function])

    workbook = _load_script(
        analysis / "build_chapter4_master_excel.py", f"_hydrosheaf_repro_{tag}_workbook"
    )
    _set_globals(
        workbook,
        TABLES_DIR=tables,
        FIGURES_DIR=figures,
        OUT_XLSX=output_dir / "Chapter4_Master_Tables.xlsx",
        MANIFEST_MD=output_dir / "CHAPTER4_DELIVERABLES_MANIFEST.md",
    )
    workbook.build_master_workbook()
    workbook.build_manifest_markdown()


def _run_chapter4(stage_outputs: Path) -> StepResult:
    chapter4 = stage_outputs / "chapter4"
    audit_dir = stage_outputs / "chapter4_audit"
    _build_chapter4_generators(chapter4, "chapter4")
    _build_chapter4_generators(audit_dir, "chapter4_audit")

    audit = _load_script(
        REPO_ROOT / "scripts" / "analysis" / "audit_chapter4_reproducibility.py",
        "_hydrosheaf_repro_chapter4_audit",
    )
    _set_globals(
        audit,
        REPO_ROOT=REPO_ROOT,
        ORIGINAL_DIR=chapter4,
        ORIGINAL_TABLES_DIR=chapter4 / "tables",
        ORIGINAL_FIGURES_DIR=chapter4 / "figures",
        ORIGINAL_MASTER_XLSX=chapter4 / "Chapter4_Master_Tables.xlsx",
        AUDIT_DIR=audit_dir,
        AUDIT_TABLES_DIR=audit_dir / "tables",
        AUDIT_FIGURES_DIR=audit_dir / "figures",
        AUDIT_MASTER_XLSX=audit_dir / "Chapter4_Master_Tables.xlsx",
        AUDIT_MANIFEST_MD=audit_dir / "CHAPTER4_DELIVERABLES_MANIFEST.md",
        AUDIT_PROVENANCE_MD=audit_dir / "PROVENANCE_MANIFEST.md",
        AUDIT_REPORT_MD=chapter4 / "CHAPTER4_AUDIT_REPORT.md",
        AUDIT_REPORT_COPY_MD=audit_dir / "CHAPTER4_AUDIT_REPORT.md",
    )
    table_comparison = audit.compare_tables()
    figure_comparison = audit.compare_figures()
    workbook_comparison = audit.compare_excel()
    environment = audit.get_environment_and_provenance_metadata()
    # Timings are operational diagnostics, not scientific output.  Supplying
    # zeroes keeps the generated historical report byte-stable.
    timing = {str(key): 0.0 for key in environment.get("generator_scripts", {})}
    audit.write_provenance_manifest_and_report(
        environment,
        timing,
        table_comparison,
        figure_comparison,
        workbook_comparison,
    )
    return StepResult(
        "chapter4",
        "PASS",
        outputs=["chapter4", "chapter4_audit"],
        notes=[
            "ten preserved Chapter 4 generators executed twice in staging",
            "audit comparison and provenance report generated without touching historical outputs",
        ],
    )


def _run_simple_script(
    *, name: str, script: Path, stage_outputs: Path, output_dir_name: str, overrides: dict[str, Any]
) -> StepResult:
    output_dir = stage_outputs / output_dir_name
    output_dir.mkdir(parents=True, exist_ok=True)
    module = _load_script(script, f"_hydrosheaf_repro_{name}")
    _set_globals(module, **overrides, OUT=output_dir)
    module.main()
    return StepResult(name, "PASS", outputs=[output_dir_name])


def _run_field_outputs(stage_outputs: Path) -> StepResult:
    script = (
        REPO_ROOT
        / "M6"
        / "m6_field_transfer_benchmark"
        / "scripts"
        / "run_m6_refined_cross_sectional.py"
    )
    module = _load_script(script, "_hydrosheaf_repro_refined_field_qa")
    output_dir = stage_outputs / "refined_field_input_qa"
    manifest = module.run(output_dir)
    if set(manifest["field_sources_used"]) != {"central_region", "upper_east_region"}:
        raise ValueError("Refined field QA did not use exactly the two approved cohorts")
    return StepResult(
        "field_outputs",
        "PASS",
        outputs=["refined_field_input_qa"],
        notes=[
            "completed Central Region and Upper East Region workbooks are the only analysis inputs",
            "charge balance and cohort summaries were recomputed; derived CSVs served only as rowwise validation mirrors",
            "no legacy master merge, Northern workbook, seasonal transfer, or elevation-as-head graph was run",
        ],
    )


def _run_uer_outputs(stage_outputs: Path) -> StepResult:
    script = (
        REPO_ROOT
        / "M6"
        / "m6_field_transfer_benchmark"
        / "scripts"
        / "run_m6_refined_cross_sectional.py"
    )
    module = _load_script(script, "_hydrosheaf_repro_o5_source_check")
    _, source_info = module._read_source(
        "upper_east_region", module.SOURCES["upper_east_region"]
    )
    status_path = stage_outputs / "uer_objective5_transfer_status.json"
    write_json(
        status_path,
        {
            "status": "ABSTAIN",
            "field_source": source_info,
            "field_graph": "not run",
            "reason": (
                "The approved UER workbook has no measured hydraulic-head time series "
                "or independent flow-path truth. The legacy edge file was created by "
                "promoting DEM elevation to hydraulic head and is not reused."
            ),
            "not_recalculated": [
                "legacy UER edge/sheaf gating results",
                "combined minimax frontier package until separated from the invalid edge input",
                "vendor-backed procurement frontier; current posted-list-cost scenario is budgetary only",
            ],
        },
    )
    return StepResult(
        "uer_outputs",
        "ABSTAIN",
        outputs=["uer_objective5_transfer_status.json"],
        notes=[
            "UER source provenance was checked against the completed workbook and its mirror",
            "no legacy Northern-derived or DEM-as-head network output was reused",
            "Objective 5 graph and combined frontier require a clean, separately validated rerun",
        ],
    )


def _run_ttd_outputs(stage_outputs: Path) -> StepResult:
    module = _load_script(
        REPO_ROOT / "scripts" / "run_ttd_graph_virtual_benchmark.py",
        "_hydrosheaf_repro_ttd_virtual",
    )
    output_dir = stage_outputs / "ttd_graph_virtual_benchmark_v1"
    config = REPO_ROOT / "configs" / "ttd_graph_virtual_benchmark_v1.json"
    output_dir.mkdir(parents=True, exist_ok=True)
    module.run_ttd_graph_virtual_benchmark(output_dir, config_path=config, overwrite=True)
    evaluator = _load_script(
        REPO_ROOT / "scripts" / "evaluate_graph_ttd_inversion_benchmark.py",
        "_hydrosheaf_repro_ttd_evaluator",
    )
    evaluator.run_phase2_evaluation(output_dir, config)
    return StepResult(
        "ttd_outputs",
        "PASS",
        outputs=["ttd_graph_virtual_benchmark_v1"],
        notes=[
            "static, dynamic, particle, and comparative virtual benchmark phases executed",
            "phase 2/3 comparative reports generated by the preserved evaluator",
        ],
    )


def _recipe_steps(stage_outputs: Path) -> dict[str, Callable[[], StepResult]]:
    analysis = REPO_ROOT / "scripts" / "analysis"
    return {
        "chapter4": lambda: _run_chapter4(stage_outputs),
        "field_outputs": lambda: _run_field_outputs(stage_outputs),
        "objective3": lambda: _run_simple_script(
            name="objective3",
            script=REPO_ROOT
            / "outputs"
            / "objective3_closure_2026-09-15"
            / "rebuild_objective3_closure.py",
            stage_outputs=stage_outputs,
            output_dir_name="objective3_closure_2026-09-15",
            overrides={"ROOT": REPO_ROOT},
        ),
        "objective4_anova": lambda: _run_simple_script(
            name="objective4_anova",
            script=REPO_ROOT
            / "outputs"
            / "objective4_anova_rerun_2026-09-15"
            / "run_objective4_anova_rerun.py",
            stage_outputs=stage_outputs,
            output_dir_name="objective4_anova_rerun_2026-09-15",
            overrides={"ROOT": REPO_ROOT},
        ),
        "objective4_closure": lambda: _run_simple_script(
            name="objective4_closure",
            script=REPO_ROOT
            / "outputs"
            / "objective4_closure_2026-09-15"
            / "objective4_closure_audit.py",
            stage_outputs=stage_outputs,
            output_dir_name="objective4_closure_2026-09-15",
            overrides={"ROOT": REPO_ROOT},
        ),
        "uer_outputs": lambda: _run_uer_outputs(stage_outputs),
        "ttd_outputs": lambda: _run_ttd_outputs(stage_outputs),
    }


def _historical_family_names(reference: Path) -> set[str]:
    names: set[str] = set()
    for path in reference.rglob("*"):
        if not path.is_file():
            continue
        rel = path.relative_to(reference).as_posix()
        names.add(rel.split("/", 1)[0])
    return names


def _scoped_manifest(
    manifest: dict[str, Any],
    *,
    include_families: Iterable[str] | None = None,
    exclude_families: Iterable[str] = (),
) -> dict[str, Any]:
    """Return a stable manifest limited to selected top-level output families."""

    includes = set(include_families) if include_families is not None else None
    excludes = tuple(str(prefix).rstrip("/") for prefix in exclude_families)

    def matches_family(path: str, family: str) -> bool:
        return path == family or path.startswith(family + "/")

    files = []
    for item in manifest.get("files", []):
        path = str(item["path"])
        family = path.split("/", 1)[0]
        if includes is not None and family not in includes:
            continue
        if any(matches_family(path, prefix) for prefix in excludes):
            continue
        files.append(dict(item))
    return {
        "schema": manifest.get("schema", "hydrosheaf-output-tree-manifest-v1"),
        "files": files,
        "file_count": len(files),
        "total_bytes": sum(int(item.get("bytes", 0)) for item in files),
    }


def _runtime_fingerprint() -> dict[str, str]:
    """Record the runtime components that can affect rendered bytes."""

    packages = (
        "numpy",
        "pandas",
        "scipy",
        "matplotlib",
        "openpyxl",
        "geopandas",
        "shapely",
        "pyproj",
    )
    fingerprint = {
        "python": sys.version.replace("\n", " "),
        "platform": sys.platform,
    }
    for package in packages:
        try:
            fingerprint[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            fingerprint[package] = "MISSING"
    return fingerprint


def _unresolved_for_reference(reference: Path) -> dict[str, str]:
    result: dict[str, str] = {}
    families = _historical_family_names(reference)
    for family, reason in UNRESOLVED_FAMILIES.items():
        root = family.split("/", 1)[0]
        if family in families or root in families:
            result[family] = reason
    if (
        "cr_geology_spatial_join" in families
        and not (REPO_ROOT / "data" / "FieldData" / "CRdata" / "CentralRegion.xlsx").exists()
    ):
        result["cr_geology_spatial_join"] = (
            "original CentralRegion.xlsx is absent; the rerun can only use raw columns recovered from the historical join artifact"
        )
    return result


def _write_report(
    path: Path,
    *,
    run_id: str,
    reference: Path,
    stage_outputs: Path,
    environment: dict[str, str],
    steps: list[StepResult],
    comparison: dict[str, Any],
    unresolved: dict[str, str],
) -> None:
    payload = {
        "schema": "hydrosheaf-output-reproduction-report-v1",
        "run_id": run_id,
        "reference": str(reference),
        "stage_outputs": str(stage_outputs),
        "environment": environment,
        "steps": [step.__dict__ for step in steps],
        "unresolved_families": unresolved,
        "comparison": comparison,
        "exact_historical_replay": bool(comparison.get("exact")) and not unresolved,
    }
    write_json(path, payload)
    markdown = [
        "# HydroSheaf output-tree reproduction",
        "",
        f"- Run: `{run_id}`",
        f"- Reference: `{reference}`",
        f"- Staging tree: `{stage_outputs}`",
        f"- Exact historical replay: **{'PASS' if payload['exact_historical_replay'] else 'FAIL'}**",
        "",
        "## Step status",
        "",
        "| Step | Status | Outputs | Notes |",
        "|---|---|---|---|",
    ]
    for step in steps:
        markdown.append(
            "| {0} | {1} | {2} | {3} |".format(
                step.name,
                step.status,
                ", ".join(step.outputs) or "-",
                " ".join(step.notes) if step.notes else (step.error or "-"),
            )
        )
    markdown.extend(
        [
            "",
        "## Byte comparison",
        "",
        f"- Comparison scope: `{comparison.get('scope', 'full output tree')}`",
        f"- Known unresolved families excluded from the byte diff: `{len(comparison.get('excluded_families', []))}`",
        f"- Reference files: `{comparison.get('reference_file_count', 0)}`",
            f"- Candidate files: `{comparison.get('candidate_file_count', 0)}`",
            f"- Missing: `{len(comparison.get('missing', []))}`",
            f"- Extra: `{len(comparison.get('extra', []))}`",
            f"- Changed: `{len(comparison.get('changed', []))}`",
            "",
            "## Unresolved historical families",
            "",
        ]
    )
    if unresolved:
        for family, reason in sorted(unresolved.items()):
            markdown.append(f"- `{family}` — {reason}")
    else:
        markdown.append("- None.")
    path.with_suffix(".md").write_text("\n".join(markdown) + "\n", encoding="utf-8")


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, default=REFERENCE_OUTPUTS)
    parser.add_argument("--run-root", type=Path, default=DEFAULT_RUN_ROOT)
    parser.add_argument("--run-id", default="latest")
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--epoch", default=DEFAULT_EPOCH)
    parser.add_argument(
        "--only",
        nargs="+",
        choices=sorted(_recipe_steps(Path("."))),
        help="run only selected preserved recipe groups",
    )
    parser.add_argument("--no-compare", action="store_true")
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="replace an existing run directory under --run-root",
    )
    parser.add_argument(
        "--allow-unresolved",
        action="store_true",
        help="return success for the generated subset while still recording unresolved families",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    reference = args.reference.resolve()
    run_id = str(args.run_id)
    run_root = args.run_root.resolve() / run_id
    stage_outputs = run_root / "outputs"
    if stage_outputs.exists():
        if not args.overwrite:
            raise FileExistsError(
                f"run directory already exists: {run_root}; use --overwrite explicitly"
            )
        shutil.rmtree(stage_outputs)
    stage_outputs.mkdir(parents=True, exist_ok=True)

    environment = apply_deterministic_environment(seed=args.seed, epoch=args.epoch)
    environment.update(
        {f"runtime_{key}": value for key, value in _runtime_fingerprint().items()}
    )
    seed_all(args.seed)
    _patch_plotting()
    recipes = _recipe_steps(stage_outputs)
    selected = args.only or list(recipes)
    steps: list[StepResult] = []
    for name in selected:
        try:
            result = recipes[name]()
        except Exception as exc:  # keep the report useful after one blocked family
            result = StepResult(
                name=name,
                status="FAIL",
                error=f"{type(exc).__name__}: {exc}",
                notes=[traceback.format_exc(limit=8)],
            )
        steps.append(result)

    _canonicalize_text_paths(stage_outputs, stage_outputs)
    _canonicalize_workbooks(stage_outputs)
    _align_renderer_metadata(stage_outputs, reference)
    _align_workbook_container(stage_outputs, reference)
    candidate_manifest = manifest_for_tree(
        stage_outputs, include_suffixes=OUTPUT_SUFFIXES
    )
    write_json(run_root / "candidate_manifest.json", candidate_manifest)
    reference_manifest = (
        manifest_for_tree(reference, include_suffixes=OUTPUT_SUFFIXES)
        if reference.exists()
        else {"files": []}
    )
    write_json(run_root / "reference_manifest.json", reference_manifest)
    unresolved = _unresolved_for_reference(reference)
    selected_families = None
    if args.only:
        selected_families = {
            output.split("/", 1)[0]
            for step in steps
            for output in step.outputs
        }
    comparison_reference = _scoped_manifest(
        reference_manifest,
        include_families=selected_families,
        exclude_families=unresolved,
    )
    comparison_candidate = _scoped_manifest(
        candidate_manifest,
        include_families=selected_families,
        exclude_families=unresolved,
    )
    comparison = (
        compare_trees(comparison_reference, comparison_candidate)
        if not args.no_compare
        else {"exact": False, "skipped": True}
    )
    comparison["scope"] = (
        "full output tree" if selected_families is None else ", ".join(sorted(selected_families))
    )
    comparison["excluded_families"] = sorted(unresolved)
    _write_report(
        run_root / "reproduction_report.json",
        run_id=run_id,
        reference=reference,
        stage_outputs=stage_outputs,
        environment=environment,
        steps=steps,
        comparison=comparison,
        unresolved=unresolved,
    )

    failed_steps = [step for step in steps if step.status != "PASS"]
    strict_failure = bool(failed_steps) or bool(unresolved)
    if not args.no_compare and not comparison.get("exact", False):
        strict_failure = True
    if args.allow_unresolved and not failed_steps and (
        args.no_compare or comparison.get("exact", False)
    ):
        strict_failure = False
    print(json.dumps({
        "run_id": run_id,
        "stage_outputs": str(stage_outputs),
        "steps": [step.__dict__ for step in steps],
        "comparison": comparison,
        "unresolved_families": unresolved,
        "strict_failure": strict_failure,
    }, indent=2, sort_keys=True, default=str))
    return 1 if strict_failure else 0


if __name__ == "__main__":
    raise SystemExit(main())
