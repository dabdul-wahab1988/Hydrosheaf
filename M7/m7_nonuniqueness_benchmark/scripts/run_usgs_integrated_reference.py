"""Prepare a provenance-preserving USGS model-reference benchmark package.

This command intentionally combines *panels*, not labels.  The USGS national
age release is a model-derived residence-time panel, while the M4 MODPATH
archives are model-derived topology/travel-time panels.  They have no shared
well-to-cell crosswalk in this repository, so the output keeps them separate
and refuses to emit an integrated field-accuracy score.

The optional ``--aiken-root`` argument points to a user-supplied Aiken
MODFLOW-NWT/MODPATH5 release (the supplied ZIP-backed release is supported).
The command parses its compact workbooks and validated text/binary MODPATH
panels without extracting the multi-gigabyte archives.  It never guesses a
well-to-cell crosswalk and never converts model-derived pathways into direct
well--well truth.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
from typing import Any, Iterable
import zipfile

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[3]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from hydrosheaf.validation.integrated_benchmark import (  # noqa: E402
    audit_reference_panel,
    build_run_manifest,
    source_manifest_entry,
    write_json,
)
from hydrosheaf.validation.aiken_reference import (  # noqa: E402
    PUBLICATION_DOI as AIKEN_PUBLICATION_DOI,
    REFERENCE_TYPE as AIKEN_REFERENCE_TYPE,
    SOURCE_DOI as AIKEN_SOURCE_DOI,
    SOURCE_URL as AIKEN_SOURCE_URL,
    AikenReference,
    load_aiken_reference,
)


DEFAULT_AGE_ROOT = (
    REPO_ROOT
    / "M2"
    / "m2_benchmark"
    / "external"
    / "usgs_age"
    / "input"
    / "DataForNationalGroundwaterAge_1_1"
)
DEFAULT_M4_RESULTS = REPO_ROOT / "M4" / "m4_topology_benchmark" / "results"
DEFAULT_OUTPUT = REPO_ROOT / ".codex_work" / "runs" / "RUN-USGS-INTEGRATED-REFERENCE-20260908-01"

AGE_TABLES = {
    "sites": "Table_1_Sites.txt",
    "ages": "Table_2_Ages.txt",
    "tracers": "Table_3_Tracers.txt",
}
MODPATH_TABLES = {
    "tier_1_savage": "tier_1_savage_edge_agreement.csv",
    "tier_2_great_miami": "tier_2_great_miami_edge_agreement.csv",
    "tier_3_long_island": "tier_3_long_island_edge_agreement.csv",
}

USGS_AGE_DOI = "10.5066/P9W7T0DN"
USGS_AGE_URL = "https://www.usgs.gov/data/data-distribution-groundwater-age-aquifers-used-public-supply-continental-united-states-2004"


def _read_tsv(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(
        path,
        sep="\t",
        na_values=["na", "NA", ""],
        encoding="latin-1",
        low_memory=False,
    )


def _read_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path, low_memory=False)


def _finite_numeric(frame: pd.DataFrame, columns: Iterable[str]) -> pd.DataFrame:
    result = frame.copy()
    for column in columns:
        if column in result:
            result[column] = pd.to_numeric(result[column], errors="coerce")
    return result


def load_usgs_age_panel(age_root: Path) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    """Load the three compact USGS tables and retain model-output semantics."""

    tables = {name: _read_tsv(age_root / filename) for name, filename in AGE_TABLES.items()}
    missing_files = [
        filename for name, filename in AGE_TABLES.items() if tables[name].empty
    ]
    if any(frame.empty for frame in tables.values()):
        return (
            pd.DataFrame(columns=["node_id"]),
            pd.DataFrame(columns=["node_id"]),
            {
                "status": "MISSING" if len(missing_files) == len(AGE_TABLES) else "PARTIAL",
                "missing_files": missing_files,
                "duplicate_counts": {},
                "n_rows": 0,
            },
        )

    sites = tables["sites"].drop_duplicates("SampleID", keep="first").copy()
    ages = tables["ages"].drop_duplicates("SampleID", keep="first").copy()
    tracers = tables["tracers"].drop_duplicates("SampleID", keep="first").copy()
    merged = sites.merge(ages, on="SampleID", how="inner", suffixes=("", "_age"))
    merged = merged.merge(tracers, on="SampleID", how="left", suffixes=("", "_tracer"))
    merged = _finite_numeric(
        merged,
        (
            "LatDD83",
            "LongDD83",
            "TopOfScrn_m",
            "Depth_m",
            "ScrnLg_m",
            "MidPt_m",
            "Rpt_TotAge_yrs",
            "Rept_TotAge_Err_yrs",
            "Rpt_Probability",
            "3H_TU",
            "3H_err_TU",
            "3He_trit_TU",
            "3He_trit_err_TU",
            "SF6_pptv",
            "SF6_err_pptv",
            "CFC-11_pptv",
            "CFC-11_err_pptv",
            "CFC-12_pptv",
            "CFC-12_err_pptv",
            "14C_pmC",
            "14C_err_pmC",
        ),
    )
    merged["node_id"] = "USGSAGE:" + merged["SampleID"].astype(str).str.strip()
    merged["reference_panel"] = "usgs_national_age"
    merged["reference_type"] = "calibrated_model_reference"
    merged["age_reference_kind"] = "USGS reported LPM model output"
    # Canonical names are written alongside the native USGS names so the
    # generic contract can audit coordinate/time coverage without guessing.
    merged["lat"] = merged["LatDD83"]
    merged["lon"] = merged["LongDD83"]
    merged["sample_date"] = merged["SampleDate"]
    merged["screen_top_m"] = merged["TopOfScrn_m"]
    merged["screen_bottom_m"] = merged["Depth_m"]
    merged["screen_length_m"] = merged["ScrnLg_m"]
    merged["reported_age_years"] = merged["Rpt_TotAge_yrs"]
    merged["reported_age_sigma_years"] = merged["Rept_TotAge_Err_yrs"]
    merged["tracer_count"] = merged.get("LPM_TracersMod", "").fillna("").map(
        lambda value: len([item for item in str(value).split(",") if item.strip()])
    )

    keep = [
        "node_id",
        "reference_panel",
        "reference_type",
        "SampleID",
        "StudyUnit",
        "StudyArea",
        "AqGroup",
        "State",
        "LatDD83",
        "LongDD83",
        "TopOfScrn_m",
        "Depth_m",
        "ScrnLg_m",
        "MidPt_m",
        "SampleDate",
        "LPM_Name",
        "LPM_TracersMod",
        "AgeCat",
        "FracAnthropocene",
        "FracHolocene",
        "FracPleistocene",
        "Rpt_TotAge_yrs",
        "Rept_TotAge_Err_yrs",
        "Rpt_ChiSquare",
        "Rpt_Probability",
        "age_reference_kind",
        "lat",
        "lon",
        "sample_date",
        "screen_top_m",
        "screen_bottom_m",
        "screen_length_m",
        "reported_age_years",
        "reported_age_sigma_years",
        "tracer_count",
        "3H_TU",
        "3H_err_TU",
        "3He_trit_TU",
        "3He_trit_err_TU",
        "SF6_pptv",
        "SF6_err_pptv",
        "CFC-11_pptv",
        "CFC-11_err_pptv",
        "CFC-12_pptv",
        "CFC-12_err_pptv",
        "14C_pmC",
        "14C_err_pmC",
    ]
    for column in keep:
        if column not in merged:
            merged[column] = pd.NA
    observations = merged[keep].copy()
    nodes = observations[
        [
            "node_id",
            "reference_panel",
            "reference_type",
            "SampleID",
            "StudyUnit",
            "AqGroup",
            "State",
            "LatDD83",
            "LongDD83",
            "lat",
            "lon",
            "TopOfScrn_m",
            "Depth_m",
            "ScrnLg_m",
            "MidPt_m",
            "SampleDate",
            "sample_date",
            "age_reference_kind",
        ]
    ].copy()
    duplicates = {
        name: int(len(tables[name]) - tables[name]["SampleID"].nunique())
        for name in tables
        if "SampleID" in tables[name]
    }
    audit = {
        "status": "COMPLETE",
        "missing_files": missing_files,
        "duplicate_counts": duplicates,
        "n_rows": int(len(observations)),
        "n_age_values": int(observations["Rpt_TotAge_yrs"].notna().sum()),
        "n_coordinate_pairs": int(
            (observations["LatDD83"].notna() & observations["LongDD83"].notna()).sum()
        ),
    }
    return nodes, observations, audit


def load_modpath_panels(results_root: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    """Load the already-ingested compact M4 edge tables, one archive at a time."""

    node_rows: list[dict[str, Any]] = []
    edge_frames: list[pd.DataFrame] = []
    travel_frames: list[pd.DataFrame] = []
    panel_audits: dict[str, Any] = {}
    summary_path = results_root / "external_modpath_archive_summary.csv"
    summary = _read_csv(summary_path)
    summary_by_tier = {
        str(row.get("validation_tier")): row.to_dict()
        for _, row in summary.iterrows()
    }
    for tier, filename in MODPATH_TABLES.items():
        frame = _read_csv(results_root / filename)
        meta = summary_by_tier.get(tier, {})
        if frame.empty:
            panel_audits[tier] = {
                "status": "MISSING",
                "file": filename,
                "n_edges": 0,
                "source_doi": meta.get("source_doi"),
            }
            continue
        frame = frame.copy()
        frame["panel_id"] = tier
        frame["reference_type"] = "calibrated_model_reference"
        frame["u"] = frame["source_node"].astype(str)
        frame["v"] = frame["target_node"].astype(str)
        frame["edge_id"] = frame["u"] + "->" + frame["v"]
        frame["travel_time_days"] = pd.to_numeric(frame.get("travel_time_mean"), errors="coerce")
        frame["travel_time_p10_days"] = pd.to_numeric(frame.get("travel_time_p10"), errors="coerce")
        frame["travel_time_p90_days"] = pd.to_numeric(frame.get("travel_time_p90"), errors="coerce")
        frame["model_name"] = meta.get("archive_name", tier)
        frame["source_doi"] = meta.get("source_doi")
        frame["source_url"] = meta.get("source_url")
        frame["edge_reference_kind"] = "MODPATH endpoint source-receptor pair"
        edge_frames.append(
            frame[
                [
                    "panel_id",
                    "reference_type",
                    "edge_id",
                    "u",
                    "v",
                    "endpoint_particle_count",
                    "pathline_particle_count",
                    "classification",
                    "direction_agrees",
                    "source_receptor_overlap",
                    "travel_time_days",
                    "travel_time_p10_days",
                    "travel_time_p90_days",
                    "model_name",
                    "source_doi",
                    "source_url",
                    "edge_reference_kind",
                ]
            ]
        )
        for node_id in sorted(set(frame["u"]) | set(frame["v"])):
            node_rows.append(
                {
                    "node_id": node_id,
                    "panel_id": tier,
                    "reference_type": "calibrated_model_reference",
                    "model_name": meta.get("archive_name", tier),
                    "source_doi": meta.get("source_doi"),
                    "node_reference_kind": "MODPATH grid-cell endpoint projection",
                }
            )
        travel_frames.append(
            frame[
                [
                    "panel_id",
                    "edge_id",
                    "u",
                    "v",
                    "travel_time_days",
                    "travel_time_p10_days",
                    "travel_time_p90_days",
                    "source_doi",
                ]
            ]
        )
        panel_audits[tier] = {
            "status": "COMPLETE",
            "file": filename,
            "n_edges": int(len(frame)),
            "n_nodes": int(len(set(frame["u"]) | set(frame["v"]))),
            "source_doi": meta.get("source_doi"),
            "source_url": meta.get("source_url"),
            "model_name": meta.get("archive_name", tier),
        }
    nodes = pd.DataFrame(node_rows, columns=["node_id", "panel_id", "reference_type", "model_name", "source_doi", "node_reference_kind"])
    edges = pd.concat(edge_frames, ignore_index=True) if edge_frames else pd.DataFrame()
    travel = pd.concat(travel_frames, ignore_index=True) if travel_frames else pd.DataFrame()
    return nodes, edges, travel, {"summary_path": str(summary_path), "panels": panel_audits}


def _aiken_archive_metadata(path: Path) -> dict[str, Any]:
    """Read ZIP central-directory metadata without hashing or extracting it.

    The Aiken release contains multi-gigabyte archives.  A normal benchmark
    run must not silently stream those archives merely to create an inventory;
    archive/member hashes remain an explicit provenance operation in
    ``inventory_aiken_source``.
    """

    result: dict[str, Any] = {
        "name": path.name,
        "path": str(path),
        "size_bytes": path.stat().st_size,
        "hash_status": "not_requested",
    }
    try:
        with zipfile.ZipFile(path) as handle:
            infos = handle.infolist()
        files = [item for item in infos if not item.is_dir()]
        result.update(
            {
                "status": "READABLE",
                "member_count": len(files),
                "directory_count": len(infos) - len(files),
                "uncompressed_member_bytes": int(sum(item.file_size for item in files)),
                "extension_counts": {
                    suffix: sum(item.filename.lower().endswith(suffix) for item in files)
                    for suffix in (".xlsx", ".loc", ".pth", ".ept", ".sum", ".nam")
                },
            }
        )
    except (OSError, zipfile.BadZipFile) as exc:
        result.update({"status": "UNREADABLE", "error": str(exc)})
    return result


def inventory_aiken(
    root: Path | None,
    *,
    reference: AikenReference | None = None,
) -> dict[str, Any]:
    """Inventory and, when supplied, attest the parsed Aiken reference.

    This is deliberately a *lightweight* inventory: it reports archive sizes
    and central-directory counts but does not hash multi-gigabyte archives.
    Pass ``reference`` after ``load_aiken_reference`` to record the parsed
    panel counts and evidence boundary.  Full archive/member hashing remains
    an explicit call to :func:`inventory_aiken_source`.
    """

    expected = [
        "ancillary/geochemistry/Table 10 Split into 7 tables revised 1 12 21.xlsx",
        "ancillary/geochemistry/Tables 2 and 8 thru 15 v 8 18 2020.xlsx",
        "readme.txt",
        "model.zip",
        "output.zip",
        "ancillary.zip",
    ]
    if root is None:
        return {
            "status": "NOT_SUPPLIED",
            "expected_files": expected,
            "mapping_status": "DEFERRED",
            "reason": "Supply the local Aiken release root before parsing its panel-separated reference data.",
        }
    root = root.resolve()
    if not root.exists():
        return {
            "status": "MISSING",
            "root": str(root),
            "expected_files": [],
            "mapping_status": "DEFERRED",
            "reason": "The supplied Aiken root does not exist; no download or inference is attempted.",
        }

    # Resolve both extracted members and members retained in the release ZIPs.
    archive_paths: list[Path]
    if root.is_file() and root.suffix.lower() == ".zip":
        archive_paths = [root]
    else:
        archive_paths = sorted(root.glob("*.zip")) if root.is_dir() else []
    archive_names: dict[str, tuple[Path, int]] = {}
    for archive in archive_paths:
        try:
            with zipfile.ZipFile(archive) as handle:
                for info in handle.infolist():
                    if not info.is_dir():
                        archive_names[info.filename.replace("\\", "/").lower()] = (archive, int(info.file_size))
        except (OSError, zipfile.BadZipFile):
            continue

    found: list[dict[str, Any]] = []
    for relative in expected:
        path = root / relative if root.is_dir() else None
        direct = path.is_file() if path is not None else False
        member_match = next(
            (
                (name, archive, size)
                for name, (archive, size) in archive_names.items()
                if name == relative.lower() or name.endswith("/" + relative.lower())
            ),
            None,
        )
        item: dict[str, Any] = {
            "relative_path": relative,
            "exists": bool(direct or member_match),
            "size_bytes": path.stat().st_size if direct and path is not None else (member_match[2] if member_match else None),
            "location": "extracted_file" if direct else ("zip_member" if member_match else None),
        }
        if member_match:
            item["archive"] = str(member_match[1])
            item["member"] = member_match[0]
        found.append(item)
    archives = [_aiken_archive_metadata(path) for path in archive_paths]
    report: dict[str, Any] = {
        "status": "PARSED" if reference is not None else ("READY_FOR_PARSING" if any(item["exists"] for item in found) else "MISSING"),
        "root": str(root),
        "expected_files": found,
        "archives": archives,
        "mapping_status": "MODEL_REFERENCE_ONLY" if reference is not None else "DEFERRED",
        "source_doi": AIKEN_SOURCE_DOI,
        "publication_doi": AIKEN_PUBLICATION_DOI,
        "source_url": AIKEN_SOURCE_URL,
        "hash_policy": "Archive/member SHA-256 not requested by the normal run; use inventory_aiken_source explicitly.",
        "reason": (
            "Aiken panels were parsed as calibrated-model references. They are not independent "
            "well-to-well adjacency, age, flow, or reaction truth."
            if reference is not None
            else "No Aiken columns are mapped until the release root is supplied to the parser."
        ),
    }
    if reference is not None:
        report["reference_audit"] = reference.audit
        report["readme"] = reference.readme
        report["modelgeoref"] = reference.modelgeoref
        report["parsed_counts"] = {
            "well_metadata": int(len(reference.well_metadata)),
            "field_samples": int(len(reference.field_samples)),
            "cfc_ages": int(len(reference.cfc_ages)),
            "pathways": int(len(reference.pathways)),
            "chemistry": int(len(reference.chemistry)),
            "model_locations": int(len(reference.model_locations)),
            "model_references": int(len(reference.model_references)),
            "modpath5_endpoints": int(len(getattr(reference, "endpoints", pd.DataFrame()))),
            "modpath5_pathlines": int(len(getattr(reference, "pathlines", pd.DataFrame()))),
            "modpath5_binary_diagnostics": int(len(getattr(reference, "binary_parse_diagnostics", pd.DataFrame()))),
        }
    return report


def _write_csv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, index=False)


def _write_aiken_reference_outputs(
    reference: AikenReference,
    output: Path,
) -> dict[str, str]:
    """Persist Aiken's canonical panels without collapsing evidence types."""

    tables: dict[str, pd.DataFrame] = {
        "aiken_well_metadata": reference.well_metadata,
        "aiken_field_samples": reference.field_samples,
        "aiken_cfc_ages": reference.cfc_ages,
        "aiken_pathways": reference.pathways,
        "aiken_chemistry": reference.chemistry,
        "aiken_model_locations": reference.model_locations,
        "aiken_model_references": reference.model_references,
        "aiken_model_name_files": reference.model_name_files,
    }
    # Newer adapters may expose validated binary MODPATH panels.  Keep the
    # runner backward-compatible with the table-only adapter while writing
    # those panels whenever they are present.
    optional_tables = {
        "aiken_modpath_pathlines": ("modpath_pathlines", "pathlines"),
        "aiken_modpath_endpoints": ("modpath_endpoints", "endpoints"),
        "aiken_modpath_binary_diagnostics": ("binary_parse_diagnostics", "run_audits"),
    }
    for output_name, attributes in optional_tables.items():
        frame = next(
            (
                getattr(reference, attribute)
                for attribute in attributes
                if isinstance(getattr(reference, attribute, None), pd.DataFrame)
            ),
            None,
        )
        if frame is not None:
            tables[output_name] = frame
    generated: dict[str, str] = {}
    for name, frame in tables.items():
        filename = f"{name}.csv"
        _write_csv(frame, output / filename)
        generated[name] = filename
    write_json(output / "aiken_reference_manifest.json", reference.to_manifest())
    generated["aiken_reference_manifest"] = "aiken_reference_manifest.json"
    return generated


def _aiken_panel_audits(reference: AikenReference) -> list[Any]:
    """Build separate audits for age, transport, and chemistry panels."""

    metadata = {
        "required_components": ("nodes", "observations"),
        "model_name": "USGS Aiken County calibrated MODFLOW-NWT/MODPATH5 reference",
        "source_doi": AIKEN_SOURCE_DOI,
        "independent_age_truth": False,
        "independent_direct_adjacency_truth": False,
        "independent_reaction_truth": False,
        "transport_time_available": True,
        "screen_intervals_available": False,
    }
    return [
        audit_reference_panel(
            "aiken_cfc_apparent_age",
            AIKEN_REFERENCE_TYPE,
            nodes=reference.well_metadata,
            observations=reference.cfc_ages,
            metadata=metadata,
        ),
        audit_reference_panel(
            "aiken_modpath_pathway_reference",
            AIKEN_REFERENCE_TYPE,
            nodes=reference.well_metadata,
            observations=reference.pathways,
            metadata=metadata,
        ),
        audit_reference_panel(
            "aiken_chemistry_reference",
            AIKEN_REFERENCE_TYPE,
            nodes=reference.well_metadata,
            observations=reference.chemistry,
            metadata=metadata,
        ),
        audit_reference_panel(
            "aiken_modpath5_binary_transport",
            AIKEN_REFERENCE_TYPE,
            nodes=reference.well_metadata,
            observations=getattr(reference, "endpoints", pd.DataFrame()),
            metadata=metadata,
        ),
    ]


def run(
    *,
    output: Path,
    age_root: Path,
    m4_results: Path,
    aiken_root: Path | None = None,
) -> dict[str, Any]:
    output = output.resolve()
    if output.exists() and any(output.iterdir()):
        raise FileExistsError(f"Refusing to overwrite non-empty run directory: {output}")
    output.mkdir(parents=True, exist_ok=True)

    age_nodes, age_observations, age_audit = load_usgs_age_panel(age_root.resolve())
    mod_nodes, mod_edges, travel, mod_audit = load_modpath_panels(m4_results.resolve())
    _write_csv(age_nodes, output / "usgs_age_nodes.csv")
    _write_csv(age_observations, output / "usgs_age_observations.csv")
    _write_csv(mod_nodes, output / "modpath_reference_nodes.csv")
    _write_csv(mod_edges, output / "modpath_reference_edges.csv")
    _write_csv(travel, output / "modpath_reference_travel_times.csv")

    age_panel = audit_reference_panel(
        "usgs_national_age",
        "calibrated_model_reference",
        nodes=age_nodes,
        observations=age_observations,
        metadata={
            "required_components": ("nodes", "observations"),
            "model_name": "USGS reported lumped-parameter-model ages",
            "source_doi": USGS_AGE_DOI,
            "independent_age_truth": False,
            "screen_intervals_available": bool(
                age_nodes.get("TopOfScrn_m", pd.Series(dtype=float)).notna().any()
            ),
        },
    )
    model_panels = []
    summary_by_tier = mod_audit["panels"]
    for tier, item in summary_by_tier.items():
        panel_edges = mod_edges[mod_edges["panel_id"] == tier] if not mod_edges.empty else pd.DataFrame()
        panel_nodes = mod_nodes[mod_nodes["panel_id"] == tier] if not mod_nodes.empty else pd.DataFrame()
        model_panels.append(
            audit_reference_panel(
                tier,
                "calibrated_model_reference",
                nodes=panel_nodes,
                edges=panel_edges,
                metadata={
                    "model_name": item.get("model_name", tier),
                    "source_doi": item.get("source_doi"),
                    "transport_time_available": True,
                    "independent_direct_adjacency_truth": False,
                },
            )
        )
    aiken_reference: AikenReference | None = None
    aiken_panels: list[Any] = []
    aiken_outputs: dict[str, str] = {}
    if aiken_root is not None:
        # Parse the local release directly from its extracted files/ZIP
        # members.  No download, extraction, binary decoding, or inferred
        # crosswalk is hidden inside the benchmark run.
        aiken_reference = load_aiken_reference(aiken_root.resolve())
        aiken_outputs = _write_aiken_reference_outputs(aiken_reference, output)
        aiken_panels = _aiken_panel_audits(aiken_reference)

    audits = [age_panel, *model_panels, *aiken_panels]

    sources = [
        source_manifest_entry(
            age_root / filename,
            source_id=f"usgs_age_{name}",
            role="USGS public-supply groundwater-age input table",
            doi=USGS_AGE_DOI,
            url=USGS_AGE_URL,
            license_name="CC0 / USGS data release",
        )
        for name, filename in AGE_TABLES.items()
    ]
    sources.extend(
        source_manifest_entry(
            m4_results / filename,
            source_id=f"m4_{tier}",
            role="Existing M4 compact MODPATH reference-edge table",
            doi=summary_by_tier.get(tier, {}).get("source_doi"),
            url=summary_by_tier.get(tier, {}).get("source_url"),
            license_name="USGS public data release / repository-derived projection",
        )
        for tier, filename in MODPATH_TABLES.items()
    )
    if aiken_root is not None:
        sources.append(
            source_manifest_entry(
                aiken_root,
                source_id="usgs_aiken_extracted_root",
                role="USGS Aiken calibrated-model reference release root",
                doi=AIKEN_SOURCE_DOI,
                url=AIKEN_SOURCE_URL,
                license_name="CC0 / USGS data release",
            )
        )
        # Hash only small standalone metadata files.  The multi-gigabyte ZIP
        # archives remain size-inventoried; full hashes are explicit via the
        # adapter's inventory_aiken_source() function.
        for metadata_name in ("readme.txt", "modelgeoref.txt", "sir2022-5036.xml"):
            metadata_path = aiken_root.resolve() / metadata_name
            if metadata_path.is_file():
                sources.append(
                    source_manifest_entry(
                        metadata_path,
                        source_id=f"usgs_aiken_{metadata_name.replace('.', '_')}",
                        role=f"Aiken release metadata: {metadata_name}",
                        doi=AIKEN_SOURCE_DOI,
                        url=AIKEN_SOURCE_URL,
                        license_name="CC0 / USGS data release",
                    )
                )

    aiken = inventory_aiken(aiken_root, reference=aiken_reference)
    write_json(output / "aiken_inventory.json", aiken)
    write_json(
        output / "coverage.json",
        {
            "schema": "usgs-integrated-reference-coverage-v1",
            "usgs_age": age_audit,
            "modpath": mod_audit,
            "panels": [audit.to_dict() for audit in audits],
            "aiken": aiken,
            "aiken_reference_panels": [audit.to_dict() for audit in aiken_panels],
        },
    )
    manifest = build_run_manifest(
        run_id=output.name,
        protocol="usgs-integrated-reference-v1",
        reference_audits=audits,
        sources=sources,
        claim_boundary=(
            "Point 1 is a panel-separated model-reference benchmark. It supports "
            "USGS age/LPM and MODPATH topology reproduction diagnostics. No common "
            "well-to-cell crosswalk or independent field labels are present, so an "
            "integrated field-accuracy score is prohibited."
        ),
        extra={
            "age_panel": {"n_nodes": int(len(age_nodes)), "n_observations": int(len(age_observations))},
            "modpath_panels": {tier: int(len(mod_edges[mod_edges["panel_id"] == tier])) if not mod_edges.empty else 0 for tier in MODPATH_TABLES},
            "aiken_inventory": aiken,
            "aiken_reference": aiken_reference.audit if aiken_reference is not None else None,
            "source_doi_registry": {
                "usgs_age": USGS_AGE_DOI,
                "usgs_aiken": "10.5066/P9U0GHLU",
                "savage_modpath": "10.5066/F7J102FK",
                "great_miami_modpath": "10.5066/P9X4C9R6",
                "long_island_modpath": "10.5066/P97VFXZ4",
            },
            "generated_files": {
                "usgs_age_nodes": "usgs_age_nodes.csv",
                "usgs_age_observations": "usgs_age_observations.csv",
                "modpath_reference_nodes": "modpath_reference_nodes.csv",
                "modpath_reference_edges": "modpath_reference_edges.csv",
                "modpath_reference_travel_times": "modpath_reference_travel_times.csv",
                "coverage": "coverage.json",
                "aiken_inventory": "aiken_inventory.json",
                **aiken_outputs,
            },
            "created_by": "run_usgs_integrated_reference.py",
        },
    )
    write_json(output / "manifest.json", manifest)
    (output / "README.md").write_text(
        "# USGS integrated reference package\n\n"
        "This run keeps the national USGS model-derived age panel and the M4 "
        "MODPATH model-derived topology panels separate. They are not joined: "
        "no verified well-to-model-cell crosswalk, coeval chemistry panel, or "
        "independent direct-adjacency labels are available. Use the CSVs for "
        "component/model-reference diagnostics only. When --aiken-root is "
        "supplied, its well, chemistry, CFC, pathway, and validated text/binary "
        "panels are written separately; Aiken remains calibrated-model evidence "
        "and is never scored as independent direct-adjacency or reaction truth.\n",
        encoding="utf-8",
    )
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--age-root", type=Path, default=DEFAULT_AGE_ROOT)
    parser.add_argument("--m4-results", type=Path, default=DEFAULT_M4_RESULTS)
    parser.add_argument(
        "--aiken-root",
        type=Path,
        default=None,
        help="optional local Aiken release root (ZIP-backed or extracted); no download is performed",
    )
    args = parser.parse_args()
    manifest = run(
        output=args.output,
        age_root=args.age_root,
        m4_results=args.m4_results,
        aiken_root=args.aiken_root,
    )
    print(json.dumps(manifest, indent=2, ensure_ascii=False, default=str))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())


__all__ = [
    "DEFAULT_AGE_ROOT",
    "DEFAULT_M4_RESULTS",
    "DEFAULT_OUTPUT",
    "inventory_aiken",
    "load_modpath_panels",
    "load_usgs_age_panel",
    "main",
    "run",
]
