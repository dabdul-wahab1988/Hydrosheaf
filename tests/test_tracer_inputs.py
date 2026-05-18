from pathlib import Path

from hydrosheaf.nuclear.joint_lpm import fit_lpm_model, predict_lpm_tracers
from hydrosheaf.nuclear.tracer_inputs import (
    SiteInputContext,
    build_site_tracer_histories,
    dissolved_gas_to_atmospheric_equivalent,
    load_tracer_histories_csv,
    site_input_history_metadata,
    standardize_gas_observations,
)


def test_load_local_tracer_histories_long_csv(tmp_path: Path):
    path = tmp_path / "histories.csv"
    path.write_text(
        "\n".join(
            [
                "tracer,year,value,sigma",
                "SF6,2000,4.0,0.1",
                "SF6,2010,7.0,0.2",
                "CFC12,2000,540,5",
                "CFC12,2010,528,5",
            ]
        ),
        encoding="utf-8",
    )

    histories = load_tracer_histories_csv(path)

    assert set(histories) == {"SF6", "CFC12"}
    assert histories["SF6"].interpolate([2005.0])[0][0] == 5.5


def test_standardize_gas_observations_uses_explicit_solubility_only():
    observations = standardize_gas_observations(
        {
            "sf6_dissolved": 20.0,
            "sf6_solubility_per_pptv": 2.0,
            "sf6_excess_air_fraction": 0.25,
            "sf6_dissolved_sigma": 1.0,
        }
    )

    assert observations["sf6_pptv"] == 8.0
    assert observations["sf6_sigma_pptv"] == 0.4
    assert observations["gas_correction_notes"]
    assert dissolved_gas_to_atmospheric_equivalent(20.0, 2.0, excess_air_fraction=0.25) == 8.0


def test_joint_fit_accepts_local_history_paths(tmp_path: Path):
    history_path = tmp_path / "local_histories.csv"
    history_path.write_text(
        "\n".join(
            [
                "year,SF6",
                "1990,2.0",
                "2000,4.0",
                "2010,6.0",
                "2020,8.0",
            ]
        ),
        encoding="utf-8",
    )
    histories = load_tracer_histories_csv(history_path)
    predicted = predict_lpm_tracers(
        "PFM",
        {"mean_age_years": 10.0},
        2020.0,
        ["SF6"],
        histories=histories,
        max_age_years=50.0,
    )
    fit = fit_lpm_model(
        {"sf6_atm_equiv_pptv": predicted["SF6"], "sf6_sigma_pptv": 0.1},
        sample_year=2020.0,
        model="PFM",
        history_paths=[str(history_path)],
        max_age_years=50.0,
    )

    assert fit.converged
    assert abs(fit.parameters["mean_age_years"] - 10.0) < 1.5


# --- Phase 3: Site-Specific Input Histories ---

def test_northern_latitude_returns_northern_hemisphere_histories():
    ctx = SiteInputContext(site_id="A", sample_year=2020.0, latitude=45.0)
    histories = build_site_tracer_histories(ctx)
    meta = site_input_history_metadata(ctx)
    assert "3H" in histories
    assert "SF6" in histories
    assert meta["input_history_region"] == "northern_hemisphere"
    assert "northern" in meta["input_history_source"]


def test_southern_latitude_records_fallback():
    ctx = SiteInputContext(site_id="B", sample_year=2020.0, latitude=-35.0)
    histories = build_site_tracer_histories(ctx)
    meta = site_input_history_metadata(ctx)
    assert meta["input_history_region"] == "southern_hemisphere"
    assert "fallback" in meta["input_history_source"]


def test_build_site_tracer_histories_returns_non_empty_mapping():
    ctx = SiteInputContext(site_id="C", sample_year=2020.0)
    histories = build_site_tracer_histories(ctx)
    assert len(histories) >= 2
    for key, hist in histories.items():
        assert hasattr(hist, "years")
        assert hasattr(hist, "values")
