"""Small, non-extracting tests for the Aiken calibrated-reference adapter."""

from __future__ import annotations

from io import BytesIO
from pathlib import Path
import hashlib
import zipfile

import openpyxl

from hydrosheaf.validation.aiken_reference import (
    CORE_WORKBOOK_NAME,
    VOC_WORKBOOK_NAME,
    inventory_aiken_source,
    inventory_zip,
    load_aiken_reference,
    parse_apparent_age_label,
    parse_qualified_value,
)


def _xlsx(sheets: dict[str, list[list[object]]]) -> bytes:
    workbook = openpyxl.Workbook()
    workbook.remove(workbook.active)
    for name, rows in sheets.items():
        worksheet = workbook.create_sheet(name)
        for row in rows:
            worksheet.append(row)
    stream = BytesIO()
    workbook.save(stream)
    return stream.getvalue()


def _core_workbook() -> bytes:
    return _xlsx(
        {
            "Table 2": [
                ["Well ID", "County number for well, AK-n or LEX-n", "USGS site ID number", "Year Installed", "Altitude of land surface, feet NGVD 29*", "Total depth of completed well, feet below land-surface altitude+", "Pump type", "Aquifer screened by well (open hole if Bedrock)"],
                ["Public-supply wells"],
                ["Office well", "AK-455", 333124081554500, 1971, 458, 285, "Submersible", "Crouch Branch and McQueen Branch"],
                ["Monitoring wells"],
                ["AK-848", "AK-848", 333233081290704, 1985, 295, 124, "Submersible", "Gordon"],
            ],
            "Table 8": [
                ["Well ID", "County number, AK-n or LEX-n", "USGS site ID number", "Sample date", "Sample time", "Sample collection method", "Water temperature (°C)", "Specific Conductance (µs/cm)", "pH", "Dissolved oxygen (mg/L)", "Dissolved oxygen (% saturation at sample temperature)", "Comments"],
                ["Office well", "AK-455", 333124081554500, "2015-07-22", "09:47", "Wellhead", 21.5, 53, 4.46, 8.92, 103.3, ""],
            ],
            "Table 9": [
                ["Well ID", "County number, AK-n or LEX-n", "USGS site ID number", "Sample date", "Sample time", "Sample collection method", "Radium-226 (226Ra)", "Radium-228 (228Ra)", "Radium-226 + Radium-228 (226Ra + 228Ra)", "Radium-228/Radium-226 (228Ra/226Ra)"],
                ["Office well", "AK-455", 333124081554500, "2015-07-22", "09:47", "Wellhead", 0.416, 1.69, 2.106, 4.0625],
            ],
            "Table 11": [
                ["Well ID", "County number, AK-n", "USGS site ID number", "Sample date", "Sample time", "Sample collection location", "dH, permil", "dO, permil", "Aquifer screened by well (open hole if Bedrock)"],
                ["Office well", "AK-455", 333124081554500, "2015-07-22", "09:47", "Wellhead", -25.63, -4.88, "Crouch Branch and McQueen Branch"],
            ],
            "Table 12": [
                ["Well ID", "County number, AK-n", "USGS site ID number", "Sample date", "Sample time", "Sample collection method", "Nitrogen as nitrate, mg/L"],
                ["Office well", "AK-455", 333124081554500, "2016-10-19", "09:26", "Wellhead", "<0.04"],
            ],
            "Table 13": [
                ["Well ID", "County number, AK-n", "USGS site ID number", "Sample date", "Sample time", "Concentration in solution", None, None, "Piston-type flow recharge dates", None, None, "CFCs used for ages", "Assigned CFC apparent groundwater age date", "Percent young water in mixture", None, None],
                [None, None, None, None, None, "(pg/kg)", None, None, None, "(elapsed time, in years, before sample collection)"],
                [None, None, None, None, None, "CFC-11", "CFC-12", "CFC-113", "CFC-11", "CFC-12", "CFC-113", None, None, "From CFC11/12", "From CFC113/12", "From CFC113/11"],
                ["Office well", "AK-455", 333124081554500, "2015-07-22", "09:47", 6151.04, 3722.48, 6.23, "C", "C", 43.6, "CFC-113", "Early 1970's", "NP", "NP", "NP"],
            ],
            "Table 14": [
                ["Well ID", "County number, AK-n or LEX-n", "USGS site ID", "Sample date", "Sample time", "Sample temperature", "Recharge", None, "Methane", "Carbon dioxide", "Nitrogen", "Oxygen", "Argon", "Calculated recharge temperature (°C)"],
                ["Office well", "AK-455", 333124081554500, "2015-07-22", "10:00", 21.5, 458, "15Y0517", 0, 21.65, 17.15, 8.56, 0.60, 18.1],
            ],
            "Table 15": [
                ["Well ID", "County number for well, AK-n or LEX-n", "USGS site ID number", "Groundwater age, elapsed time since recharge, yrs before 2015", "Recharge date", "Groundwater-flow pathway extent from most distal recharge area to the well, ft", "Groundwater flow velocity, estimated, ft/yr", "Comment"],
                ["Office well", "AK-455", 333124081554500, 46, "Early 1970's", 2900, 63, "Low levels of VOCs, MTBE, nitrate"],
                ["Oakwood well", "AK-2715", 333051081350101, "NP", "Early 1960's", "NP", "NP", "Did not reach water table by end of simulation time"],
            ],
        }
    )


def _voc_workbook() -> bytes:
    return _xlsx(
        {
            "Table 10 aromatics": [
                ["Well ID", "County number, AK-n", "USGS site ID number", "Sample date", "Sample time", "Pump Status", "benzene", "toluene"],
                ["Office well", "AK-455", 333124081554500, "2015-07-22", "09:47", "On", "<0.026", 0.05],
            ]
        }
    )


def _make_aiken_source(root: Path, *, extracted: bool = False) -> None:
    root.mkdir(parents=True, exist_ok=True)
    readme = (
        "Publication https://doi.org/10.3133/sir20225036\n"
        "Model Archive Data Release https://doi.org/10.5066/P9U0GHLU\n"
        "Archive created: 2020-06-23\nArchive updated: 2022-05-12\nArchive released: 2023-05-23\n"
        "six MODFLOW-NWT simulations and fourteen MODPATH simulations\n"
        "16 PSWs and 4 monitoring wells; CFC and dissolved-gas concentrations were used to estimate groundwater age and particle tracking."
    ).encode("utf-16")
    georef = b"# Datum NAD 83\nupper_left -82.218802 33.560357\nupper_right -81.547336 34.089942\nlower_right -81.047445 33.643655\nlower_left -81.718494 33.116774\n"
    (root / "readme.txt").write_bytes(readme)
    (root / "modelgeoref.txt").write_bytes(georef)
    if extracted:
        geochemistry = root / "ancillary" / "geochemistry"
        geochemistry.mkdir(parents=True)
        (geochemistry / CORE_WORKBOOK_NAME).write_bytes(_core_workbook())
        (geochemistry / VOC_WORKBOOK_NAME).write_bytes(_voc_workbook())
        model = root / "model" / "AK-455_MP"
        model.mkdir(parents=True)
        (model / "AK-455.loc").write_text("138 177 7 0.1 0.2 0.3 0 0 0 0\n", encoding="ascii")
        (model / "usgs.model.reference").write_text("xul 1\nyul 2\nrotation 3\nlength_units feet\ntime_units days\nmodel MODPATH5\nESRI:102733\n", encoding="ascii")
        return
    ancillary = BytesIO()
    with zipfile.ZipFile(ancillary, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        archive.writestr(f"ancillary/geochemistry/{CORE_WORKBOOK_NAME}", _core_workbook())
        archive.writestr(f"ancillary/geochemistry/{VOC_WORKBOOK_NAME}", _voc_workbook())
    (root / "ancillary.zip").write_bytes(ancillary.getvalue())
    model = BytesIO()
    with zipfile.ZipFile(model, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        archive.writestr("model/AK-455_MP/AK-455.loc", "138 177 7 0.1 0.2 0.3 0 0 0 0\n")
        archive.writestr("model/AK-455_MP/usgs.model.reference", "xul 1\nyul 2\nrotation 3\nlength_units feet\ntime_units days\nmodel MODPATH5\nESRI:102733\n")
        archive.writestr("model/AK-455_MP/AK-455_MP.nam", "list 97 AK-455.sum\n")
    (root / "model.zip").write_bytes(model.getvalue())


def test_qualified_values_preserve_censoring_and_status() -> None:
    assert parse_qualified_value("<0.04") == {
        "raw_value": "<0.04",
        "value_numeric": None,
        "value_lower": None,
        "value_upper": 0.04,
        "qualifier": "<",
        "status": "below_reporting_limit",
        "flags": ["left_censored"],
    }
    assert parse_qualified_value("C")["status"] == "above_equilibrium"
    assert "NP" in parse_qualified_value("NP")["flags"]
    assert parse_qualified_value("3-5")["value_numeric"] == 4.0


def test_apparent_age_label_is_interval_not_a_precise_age() -> None:
    parsed = parse_apparent_age_label("Late 1960's to early 1970's")
    assert parsed["age_parse_status"] == "interval"
    assert parsed["recharge_year_min"] == 1967
    assert parsed["recharge_year_max"] == 1973
    assert "qualitative_recharge_date_interval" in parsed["age_flags"]


def test_inventory_hashes_members_without_extracting(tmp_path: Path) -> None:
    archive = tmp_path / "tiny.zip"
    archive.write_bytes(b"")
    with zipfile.ZipFile(archive, "w") as handle:
        handle.writestr("nested/file.txt", b"aiken")
    report = inventory_zip(archive)
    member = next(item for item in report["members"] if item["name"] == "nested/file.txt")
    assert member["sha256"] == hashlib.sha256(b"aiken").hexdigest()
    assert member["hash_status"] == "complete"


def test_inventory_can_skip_archive_and_member_hashes(tmp_path: Path) -> None:
    archive = tmp_path / "tiny.zip"
    with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED) as handle:
        handle.writestr("nested/file.txt", b"aiken")
    report = inventory_zip(archive, hash_members=False)
    assert report["sha256"] is None
    assert report["members"][0]["sha256"] is None
    assert report["members"][0]["hash_status"] == "not_requested"


def test_load_aiken_zip_backed_source_is_panel_separated(tmp_path: Path) -> None:
    root = tmp_path / "AikenCounty"
    _make_aiken_source(root)
    before = {path.relative_to(root) for path in root.rglob("*")}
    reference = load_aiken_reference(root)
    after = {path.relative_to(root) for path in root.rglob("*")}

    assert before == after
    assert reference.reference_type == "calibrated_model_reference"
    assert reference.readme["publication_doi"] == "10.3133/sir20225036"
    assert reference.readme["data_release_doi"] == "10.5066/P9U0GHLU"
    assert reference.readme["n_modflow_simulations"] == 6
    assert reference.readme["n_modpath_simulations"] == 14
    assert reference.modelgeoref["status"] == "COMPLETE"
    assert set(reference.well_metadata["well_id"]) == {"AK-455", "AK-848"}
    assert reference.field_samples.loc[0, "sample_datetime"] == "2015-07-22T09:47:00"
    assert reference.field_samples.loc[0, "water_temperature_c"] == 21.5
    cfc = reference.cfc_ages.iloc[0]
    assert cfc["cfc_11_status"] == "observed"
    assert cfc["piston_cfc_11_raw"] == "C"
    assert cfc["apparent_age_label"] == "Early 1970's"
    assert cfc["recharge_year_min"] == 1970
    assert cfc["age_likelihood_status"] == "screening_interval_available"
    assert set(reference.pathways["pathway_status"]) == {"reported", "not_possible"}
    assert set(reference.pathways["direct_adjacency_truth_status"]) == {"ABSTAIN"}
    assert "benzene" in set(reference.chemistry["parameter"])
    assert reference.chemistry[reference.chemistry["parameter"] == "benzene"].iloc[0]["measurement_status"] == "below_reporting_limit"
    chemistry = reference.chemistry
    # Identifier/metadata columns stay in the native table but must not be
    # emitted as analytes in the canonical chemistry panel.
    assert not chemistry["parameter"].str.contains(
        "county number|usgs site|recharge$", case=False, na=False, regex=True
    ).any()
    isotope = chemistry[chemistry["parameter"].str.startswith("dH")].iloc[0]
    assert isotope["chemical_family"] == "stable_isotope"
    assert isotope["parameter_role"] == "stable_isotope"
    assert isotope["unit"] == "permil"
    radium_ratio = chemistry[chemistry["parameter_role"] == "radium_ratio"].iloc[0]
    assert radium_ratio["unit"] == "dimensionless"
    dissolved_oxygen = chemistry[chemistry["parameter"].str.contains("Dissolved oxygen", case=False)].iloc[0]
    assert dissolved_oxygen["chemical_family"] == "field_physical_or_other"
    assert dissolved_oxygen["unit"] == "mg/L"
    oxygen_saturation = chemistry[chemistry["parameter"].str.contains("saturation", case=False)].iloc[0]
    assert oxygen_saturation["unit"] == "%"
    cfc_roles = set(chemistry.loc[chemistry["source_sheet"] == "Table 13", "parameter_role"])
    assert cfc_roles == {"cfc_concentration", "piston_elapsed_years", "young_water_fraction"}
    assert set(chemistry.loc[chemistry["parameter_role"] == "piston_elapsed_years", "unit"]) == {"years"}
    assert set(chemistry.loc[chemistry["parameter_role"] == "young_water_fraction", "unit"]) == {"%"}
    assert len(reference.model_locations) == 1
    assert reference.model_locations.loc[0, "direct_adjacency_truth_status"] == "ABSTAIN"
    assert reference.audit["integrated_scoring_allowed"] is False
    assert reference.audit["capabilities"]["independent_direct_adjacency_truth"] is False


def test_load_aiken_extracted_source_has_same_minimal_contract(tmp_path: Path) -> None:
    root = tmp_path / "extracted"
    _make_aiken_source(root, extracted=True)
    reference = load_aiken_reference(root)
    assert len(reference.well_metadata) == 2
    assert len(reference.cfc_ages) == 1
    assert len(reference.model_locations) == 1
    assert reference.model_references.iloc[0]["coordinate_reference_system"] == "ESRI:102733"


def test_full_member_hashing_is_explicit_on_aiken_inventory(tmp_path: Path) -> None:
    root = tmp_path / "AikenCounty"
    _make_aiken_source(root)
    report = inventory_aiken_source(root, hash_members=True)
    assert report["hash_policy"]["member_sha256"] == "complete"
    assert all(
        item["sha256"] is not None
        for archive in report["archives"]
        for item in archive["members"]
        if not item["is_dir"]
    )


def test_aiken_inventory_metadata_mode_does_not_hash_archives(tmp_path: Path) -> None:
    root = tmp_path / "AikenCounty"
    _make_aiken_source(root)
    report = inventory_aiken_source(root, hash_members=False)
    assert report["hash_policy"]["archive_sha256"] == "not_requested"
    assert report["hash_policy"]["member_sha256"] == "not_requested"
    assert all(archive["sha256"] is None for archive in report["archives"])
