from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
import time
import zipfile
from dataclasses import asdict
from pathlib import Path
from typing import Iterable
from urllib.request import urlopen

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[3]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from hydrosheaf.api import fit_network_pipeline
from hydrosheaf.config import default_config
from hydrosheaf.graph.types import Edge

BENCHMARK_ROOT = Path(__file__).resolve().parents[1]
EXTERNAL_ROOT = BENCHMARK_ROOT / "external" / "usgs_public_chem"
INPUT_DIR = EXTERNAL_ROOT / "input"
RESULT_DIR = EXTERNAL_ROOT / "results"
FIGURE_DIR = BENCHMARK_ROOT / "figures"

SOURCE_DOI = "10.5066/P9J7I9DH"
SOURCE_PAGE = (
    "https://www.usgs.gov/data/select-groundwater-quality-and-quality-control-data-"
    "national-water-quality-assessment-project"
)
SCIENCEBASE_ITEM = "606f0b21d34ef99870188821"
SCIENCEBASE_JSON = f"https://www.sciencebase.gov/catalog/item/{SCIENCEBASE_ITEM}?format=json"
DATA_ZIP_NAME = "NWQP_GW_QW_Data_Release_v4.zip"
DATA_ZIP_URL = (
    "https://www.sciencebase.gov/catalog/file/get/606f0b21d34ef99870188821"
    "?f=__disk__71%2F47%2F72%2F714772ae14299ae14ceebb5058b2dafeb04edb0e"
)
METADATA_URL = (
    "https://www.sciencebase.gov/catalog/file/get/606f0b21d34ef99870188821"
    "?f=__disk__08%2F2d%2Fcd%2F082dcd4ef8aa31014a8c5c8af12511d1448faf87"
)
REVISION_URL = (
    "https://www.sciencebase.gov/catalog/file/get/606f0b21d34ef99870188821"
    "?f=__disk__88%2Fbe%2Fde%2F88bede51df45f1d1105d4eec9cd5ca67b485965a"
)

MOLECULAR_WEIGHTS = {
    "Ca": 40.08,
    "Mg": 24.31,
    "Na": 22.99,
    "K": 39.10,
    "HCO3": 61.02,
    "Cl": 35.45,
    "SO4": 96.06,
    "F": 19.00,
    "NO3": 62.00,
}

KEY_COLUMNS = [
    "Network_type",
    "Network_name",
    "NAWQA_well_identification_number",
    "State",
    "Sample_date",
]


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download_if_missing(url: str, path: Path) -> None:
    if path.exists() and path.stat().st_size > 0:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with urlopen(url, timeout=120) as response:
        path.write_bytes(response.read())


def ensure_source_files() -> Path:
    INPUT_DIR.mkdir(parents=True, exist_ok=True)
    download_if_missing(DATA_ZIP_URL, INPUT_DIR / DATA_ZIP_NAME)
    download_if_missing(METADATA_URL, INPUT_DIR / "NWQP_GW_QW_Metadata_V4.xml")
    download_if_missing(REVISION_URL, INPUT_DIR / "NWQP_GW_QW_RevisionHistory_v4.txt")

    extract_root = INPUT_DIR / "NWQP_GW_QW_Data_Release_v4"
    if not (extract_root / "NWQP_GW_QW_Data_Release_v4").exists():
        with zipfile.ZipFile(INPUT_DIR / DATA_ZIP_NAME) as archive:
            archive.extractall(extract_root)
    return extract_root / "NWQP_GW_QW_Data_Release_v4"


def read_table(data_dir: Path, name: str) -> pd.DataFrame:
    return pd.read_csv(data_dir / name, sep="\t", dtype=str, na_values=["--", "NC", "na", "NA", ""])


def numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def parse_degrees_minutes(value: object) -> float:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return math.nan
    text = str(value).strip()
    if not text:
        return math.nan
    sign = -1 if text.startswith("-") else 1
    cleaned = text.replace("-", " ").replace("+", " ")
    parts = [part for part in cleaned.split() if part]
    try:
        if len(parts) == 1:
            return float(text)
        degrees = abs(float(parts[0]))
        minutes = float(parts[1])
        return sign * (degrees + minutes / 60.0)
    except ValueError:
        return math.nan


def first_numeric(df: pd.DataFrame, columns: Iterable[str]) -> pd.Series:
    out = pd.Series(np.nan, index=df.index, dtype=float)
    for column in columns:
        if column in df.columns:
            out = out.fillna(numeric(df[column]))
    return out


def build_public_samples(data_dir: Path, network_name: str | None, max_sites: int) -> tuple[pd.DataFrame, dict]:
    sites = read_table(data_dir, "Table_1_site_list_v4.txt")
    indicators = read_table(data_dir, "Table_3_qw_indicators_v4.txt")
    nutrients = read_table(data_dir, "Table_4_nutrients_v4.txt")
    majors = read_table(data_dir, "Table_5_major_ions_v4.txt")

    merge_keys = KEY_COLUMNS + ["Sample_time"]
    df = sites.merge(majors, on=KEY_COLUMNS, how="inner", suffixes=("", "_major"))
    df = df.merge(indicators, on=merge_keys, how="left", suffixes=("", "_indicator"))
    df = df.merge(nutrients, on=merge_keys, how="left", suffixes=("", "_nutrient"))

    harmonized = pd.DataFrame(
        {
            "site_id": df["NAWQA_well_identification_number"],
            "sample_id": df["NAWQA_well_identification_number"],
            "network_type": df["Network_type"],
            "network_name": df["Network_name"],
            "state": df["State"],
            "sample_date": pd.to_datetime(df["Sample_date"], errors="coerce"),
            "aquifer": df["Principal_and_regional_and_or_other_aquifer_information"],
            "well_depth_m": numeric(df["Well_depth"]),
            "screen_top_m": numeric(df["Depth_to_perforation_top"]),
            "screen_bottom_m": numeric(df["Depth_to_perforation_bottom"]),
            "latitude": df["Latitude"].map(parse_degrees_minutes),
            "longitude": df["Longitude"].map(parse_degrees_minutes),
            "pH": first_numeric(df, ["V_P00400_pH", "V_P00403_pH_wu_lab"]),
            "EC": first_numeric(df, ["V_P00095_Specific_cond_at_25C"]),
            "TDS": first_numeric(df, ["V_P70300_Total_diss_solids_dry_at_180C"]),
            "temp_c": first_numeric(df, ["V_P00010_Temperature_water"]),
            "Ca_mg_l": first_numeric(df, ["V_P00915_Calcium_wf"]),
            "Mg_mg_l": first_numeric(df, ["V_P00925_Magnesium_wf"]),
            "Na_mg_l": first_numeric(df, ["V_P00930_Sodium_wf"]),
            "K_mg_l": first_numeric(df, ["V_P00935_Potassium_wf"]),
            "HCO3_mg_l": first_numeric(
                df,
                ["V_P63786_Bicarbonate_wf_Gran_field", "V_P00453_Bicarbonate_wf_inflect_pt_fld"],
            ),
            "Cl_mg_l": first_numeric(df, ["V_P00940_Chloride_wf"]),
            "SO4_mg_l": first_numeric(df, ["V_P00945_Sulfate_wf"]),
            "F_mg_l": first_numeric(df, ["V_P00950_Fluoride_wf"]),
            "NO3_mg_l": first_numeric(df, ["V_P00631_NO3_and_NO2_wf"]),
        }
    )
    for ion, molecular_weight in MOLECULAR_WEIGHTS.items():
        harmonized[ion] = harmonized[f"{ion}_mg_l"] / molecular_weight

    required = ["Ca", "Mg", "Na", "K", "HCO3", "Cl", "SO4", "F", "pH", "EC", "latitude", "longitude"]
    complete = harmonized.dropna(subset=required).copy()
    if network_name:
        complete = complete[complete["network_name"].str.casefold() == network_name.casefold()].copy()
    else:
        counts = complete.groupby("network_name")["site_id"].count().sort_values(ascending=False)
        if counts.empty:
            raise ValueError("No public USGS network has the required major-ion fields.")
        selected_network = str(counts.index[0])
        complete = complete[complete["network_name"] == selected_network].copy()

    complete = complete.sort_values(["sample_date", "site_id"]).head(max_sites).reset_index(drop=True)
    if len(complete) < 10:
        raise ValueError(
            f"Selected public network has only {len(complete)} complete samples; choose another network."
        )

    source_summary = {
        "source": "USGS NAWQA groundwater-quality public data release",
        "doi": SOURCE_DOI,
        "source_page": SOURCE_PAGE,
        "sciencebase_item": SCIENCEBASE_ITEM,
        "selected_network_name": str(complete["network_name"].iloc[0]),
        "selected_network_type": str(complete["network_type"].iloc[0]),
        "selected_states": sorted(complete["state"].dropna().unique().tolist()),
        "selected_aquifers": sorted(complete["aquifer"].dropna().unique().tolist()),
        "complete_public_samples": int(len(complete)),
        "available_complete_public_samples": int(
            harmonized.dropna(subset=required).shape[0]
        ),
        "data_zip_sha256": file_sha256(INPUT_DIR / DATA_ZIP_NAME),
    }
    return complete, source_summary


def candidate_edges(samples: pd.DataFrame, k_nearest: int) -> list[Edge]:
    coords = samples[["latitude", "longitude", "well_depth_m"]].copy()
    coords["well_depth_m"] = coords["well_depth_m"].fillna(coords["well_depth_m"].median())
    scale = coords.std(numeric_only=True).replace(0, 1.0)
    scaled = (coords - coords.mean(numeric_only=True)) / scale
    edges: list[Edge] = []
    site_ids = samples["site_id"].tolist()
    for i, source in enumerate(site_ids):
        distances = []
        for j, target in enumerate(site_ids):
            if i == j:
                continue
            distance = float(np.linalg.norm(scaled.iloc[i].to_numpy() - scaled.iloc[j].to_numpy()))
            depth_delta = float(coords.iloc[j]["well_depth_m"] - coords.iloc[i]["well_depth_m"])
            distances.append((distance, depth_delta, target))
        for distance, depth_delta, target in sorted(distances)[:k_nearest]:
            edges.append(
                Edge(
                    edge_id=f"{source}->{target}",
                    u=source,
                    v=target,
                    attrs={
                        "distance_scaled": distance,
                        "depth_delta_m": depth_delta,
                        "source_tier": "public_usgs_chemistry_nearest_neighbor",
                    },
                )
            )
    return edges


def write_report(samples: pd.DataFrame, results: pd.DataFrame, source_summary: dict) -> None:
    report_path = RESULT_DIR / "e4b_public_chem_report.md"
    top = results.nsmallest(20, "objective_score")
    summary_cols = ["edge_id", "transport_model", "objective_score", "chemistry_r2"]
    lines = [
        "# E4b Public Groundwater-Chemistry Pilot Report",
        "",
        f"Run timestamp: {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}",
        "",
        f"Source: USGS NAWQA groundwater-quality public data release, DOI `{SOURCE_DOI}`.",
        f"Selected network: `{source_summary['selected_network_name']}`.",
        f"Samples fitted: {len(samples)}.",
        f"Candidate edge fits: {len(results)}.",
        "",
        "## Top 20 Inferred Connectivity",
        "",
        top[summary_cols].to_markdown(index=False),
        "",
        "## Reviewer Interpretation",
        "",
        "This run replaces the local E4 pilot with a fully public groundwater-chemistry dataset. "
        "It validates Hydrosheaf's public-data ingestion, unit harmonization, sparse candidate graph "
        "construction, and chemistry-fit behavior under real field observations. It does not provide "
        "an independently known reaction-path truth graph, so it should be described as public-field "
        "demonstration evidence rather than full process-truth validation.",
        "",
    ]
    report_path.write_text("\n".join(lines), encoding="utf-8")


def update_workplan(source_summary: dict, result_count: int) -> None:
    path = BENCHMARK_ROOT / "tables" / "external_validation_workplan.csv"
    if not path.exists():
        return
    df = pd.read_csv(path)
    row = {
        "validation_tier": "E4b_public_field_chemistry",
        "required_by_m2_section": "3.6 and Figure S8",
        "primary_dataset": "USGS NAWQA groundwater-quality public data release",
        "source_or_doi": SOURCE_DOI,
        "source_url": SOURCE_PAGE,
        "hydrosheaf_task": "Run Hydrosheaf on public major-ion groundwater chemistry under sparse-field conditions",
        "required_outputs": (
            "external/usgs_public_chem/results/public_chem_samples.csv; "
            "external/usgs_public_chem/results/public_chem_edge_results.csv; "
            "figures/figure_s8_e4_sparse_pilot_network.png"
        ),
        "status": "completed_public_demonstration",
        "note": (
            f"Processed {source_summary['complete_public_samples']} public samples from "
            f"{source_summary['selected_network_name']}; generated {result_count} candidate fits. "
            "No independent reaction-path truth graph is available."
        ),
    }
    df = df[df["validation_tier"] != row["validation_tier"]]
    df = pd.concat([df, pd.DataFrame([row])], ignore_index=True)
    df.to_csv(path, index=False)


def run_public_chem(network_name: str | None, max_sites: int, k_nearest: int) -> dict:
    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    data_dir = ensure_source_files()
    samples, source_summary = build_public_samples(data_dir, network_name, max_sites)
    edges = candidate_edges(samples, k_nearest)
    print(f"Generated {len(edges)} candidate edges for {len(samples)} samples.")

    config = default_config()
    config.ion_order = ["Ca", "Mg", "Na", "K", "HCO3", "Cl", "SO4", "F"]
    config.weights = [1.0] * len(config.ion_order)
    # Align conservative weights to the 8-ion order (index 5 is Cl)
    config.conservative_weights = [0.01, 0.01, 0.01, 0.01, 0.01, 1.0, 0.01, 0.01]
    config.active_minerals = ["calcite", "dolomite", "gypsum", "halite"]
    config.latent_endmembers_enabled = False
    config.latent_endmembers_count = 0

    results, _extras = fit_network_pipeline(samples.to_dict("records"), edges, config)
    rows = []
    for result in results:
        row = asdict(result)
        for label, value in zip(result.z_labels, result.z_extents):
            row[f"reaction_{label}"] = value
        for nested in [
            "transport_probabilities",
            "candidate_scores",
            "constraints_active",
            "si_u",
            "si_v",
            "isotope_metrics",
            "gibbs_metrics",
            "qc_flags",
            "reaction_fit",
        ]:
            row.pop(nested, None)
        rows.append(row)

    result_df = pd.DataFrame(rows).sort_values("objective_score")
    samples.to_csv(RESULT_DIR / "public_chem_samples.csv", index=False)
    result_df.to_csv(RESULT_DIR / "public_chem_edge_results.csv", index=False)

    manifest = {
        "run_timestamp_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "source_summary": source_summary,
        "source_files": {
            DATA_ZIP_NAME: {
                "path": str(INPUT_DIR / DATA_ZIP_NAME),
                "sha256": file_sha256(INPUT_DIR / DATA_ZIP_NAME),
                "bytes": (INPUT_DIR / DATA_ZIP_NAME).stat().st_size,
            },
            "NWQP_GW_QW_Metadata_V4.xml": {
                "path": str(INPUT_DIR / "NWQP_GW_QW_Metadata_V4.xml"),
                "sha256": file_sha256(INPUT_DIR / "NWQP_GW_QW_Metadata_V4.xml"),
                "bytes": (INPUT_DIR / "NWQP_GW_QW_Metadata_V4.xml").stat().st_size,
            },
            "NWQP_GW_QW_RevisionHistory_v4.txt": {
                "path": str(INPUT_DIR / "NWQP_GW_QW_RevisionHistory_v4.txt"),
                "sha256": file_sha256(INPUT_DIR / "NWQP_GW_QW_RevisionHistory_v4.txt"),
                "bytes": (INPUT_DIR / "NWQP_GW_QW_RevisionHistory_v4.txt").stat().st_size,
            },
        },
        "hydrosheaf_config": {
            "ion_order": config.ion_order,
            "active_minerals": config.active_minerals,
            "latent_endmembers_enabled": config.latent_endmembers_enabled,
            "candidate_graph": f"{k_nearest}-nearest neighbors in latitude/longitude/depth space",
        },
        "outputs": {
            "samples": str(RESULT_DIR / "public_chem_samples.csv"),
            "edge_results": str(RESULT_DIR / "public_chem_edge_results.csv"),
            "report": str(RESULT_DIR / "e4b_public_chem_report.md"),
        },
    }
    (RESULT_DIR / "e4b_public_chem_source_manifest.json").write_text(
        json.dumps(manifest, indent=2),
        encoding="utf-8",
    )
    write_report(samples, result_df, source_summary)
    update_workplan(source_summary, len(result_df))
    return manifest


from hydrosheaf.log import setup_logging

def main() -> None:
    setup_logging(verbose=True)
    parser = argparse.ArgumentParser(description="Run E4b public USGS groundwater-chemistry validation.")
    parser.add_argument("--network-name", default=None, help="Optional NAWQA network name to fit.")
    parser.add_argument("--max-sites", type=int, default=30, help="Maximum public samples to fit.")
    parser.add_argument("--k-nearest", type=int, default=5, help="Candidate outgoing nearest-neighbor edges.")
    args = parser.parse_args()
    manifest = run_public_chem(args.network_name, args.max_sites, args.k_nearest)
    summary = manifest["source_summary"]
    print(
        "E4b public chemistry run complete: "
        f"{summary['complete_public_samples']} samples from {summary['selected_network_name']}."
    )


if __name__ == "__main__":
    main()
