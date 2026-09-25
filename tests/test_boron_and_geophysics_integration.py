"""
Comprehensive integration tests for Boron, Diagnostic Tracers, Adaptive Sheaf Stalks,
and Continuous Geophysics conditioning in Hydrosheaf.
"""

import math
import numpy as np
import pandas as pd
import pytest

from hydrosheaf.config import Config
from hydrosheaf.data.units import (
    ChemicalRegistry,
    ChemicalSpecies,
    SPECIES_REGISTRY,
    MOLAR_MASS_G_MOL,
    CHARGE_EQUIV,
    get_species_charge_equivalent,
    get_species_molar_mass,
    get_who_drinking_water_limit,
    mgL_to_mmolL,
    mmolL_to_mgL,
)
from hydrosheaf.data.field import (
    CANONICAL_IONS,
    _base_record,
    _fill_optional_diagnostics,
)
from hydrosheaf.data.minerals import (
    MINERAL_LIBRARY,
    check_strontium_provenance,
    check_cl_br_source,
    get_mineral_stoich,
)
from hydrosheaf.models.reactions import (
    build_reaction_dictionary,
    INDICATOR_IONS,
)
from hydrosheaf.models.nitrate_isotopes import (
    IsotopeSample,
    SourceIsotopes,
    compute_isotope_prob,
    load_endmember_database,
)
from hydrosheaf.graph3d.types_3d import Node3D, Edge3D
from hydrosheaf.graph3d.build_3d import (
    attach_geophysical_observations,
    build_network_3d,
    geophysical_barrier_check,
    compute_geophysical_transit_time,
    load_geophysical_observations,
)
from hydrosheaf.graph.types import Edge
from hydrosheaf.sheaf.directed_section import DirectedEdgeMap, build_edge_maps
from hydrosheaf.sheaf.cohomology import (
    build_coboundary_matrix,
    build_rhs_vector,
    compute_cohomology,
    _is_uniform_sheaf,
    _build_adaptive_stalks,
)
from hydrosheaf.inference.topology_posterior import _edge_prior_probability
from hydrosheaf.inference.edge_fit import fit_edge
from hydrosheaf.inference.network_fit import fit_network, infer_edges
from hydrosheaf.outputs.tables import edge_results_table
from hydrosheaf.calibration.well_active_learning import (
    DEFAULT_MEASUREMENT_COSTS,
    SUPPORTED_MEASUREMENT_TYPES,
    topology_contrast_surrogate,
)
from hydrosheaf.calibration.active_learning import (
    FLAG_TO_MEASUREMENTS,
    _score_evidence_ambiguity,
    _score_geophysical_uncertainty,
)
from hydrosheaf.nitrate_source_v2 import infer_node_posteriors


class TestChemicalRegistryAndTracers:
    """Validate chemical registry, diagnostic tracers, speciation and WHO thresholds."""

    def test_registry_registered_tracers(self):
        """Ensure B, Br, SiO2, Sr, Mn, As, Li, Ba are properly registered."""
        for sp in ["B", "Br", "SiO2", "Sr", "Mn", "As", "Li", "Ba"]:
            assert sp in SPECIES_REGISTRY
            assert sp in MOLAR_MASS_G_MOL
            assert get_species_molar_mass(sp) > 0.0
            assert sp in CANONICAL_IONS

    def test_boron_speciation_ph_dependence(self):
        """Boron is neutral B(OH)3 at circumneutral pH, and B(OH)4^- (charge -1) at high pH."""
        # At pH 7.0 (well below pKa ~ 9.24), predominantly B(OH)3^0 -> charge eq close to 0
        charge_ph7 = get_species_charge_equivalent("B", ph=7.0)
        assert abs(charge_ph7) < 0.05

        # At pH 9.24 (pKa), exactly 50% dissociation -> charge eq close to -0.5
        charge_pka = get_species_charge_equivalent("B", ph=9.24)
        assert abs(charge_pka - (-0.5)) < 0.05

        # At pH 11.0 (well above pKa), predominantly B(OH)4^- -> charge eq close to -1.0
        charge_ph11 = get_species_charge_equivalent("B", ph=11.0)
        assert abs(charge_ph11 - (-1.0)) < 0.05

    def test_arsenic_speciation_redox_dependence(self):
        """Arsenic has charge 0 under reducing conditions (H3AsO3^0) and negative charge under oxic conditions."""
        charge_reducing = get_species_charge_equivalent("As", ph=7.0, do_mg_l=0.1)
        assert abs(charge_reducing) < 0.05

        charge_oxic = get_species_charge_equivalent("As", ph=7.0, do_mg_l=6.0)
        assert charge_oxic < -0.8

    def test_who_guidelines(self):
        """Check WHO thresholds for Boron and Arsenic."""
        assert get_who_drinking_water_limit("B") == 2.4
        assert get_who_drinking_water_limit("As") == 0.01
        assert get_who_drinking_water_limit("Mn") == 0.4
        assert get_who_drinking_water_limit("Ba") == 1.3
        assert get_who_drinking_water_limit("UnknownElement") is None

    def test_live_registry_conversion_for_custom_species(self):
        """Conversions must see species registered after module import."""
        symbol = "Xx_live_test"
        SPECIES_REGISTRY.register(
            ChemicalSpecies(
                symbol=symbol,
                name="Live test species",
                molar_mass_g_mol=5.0,
                base_valence=1.0,
                charge_equiv=1.0,
            )
        )

        assert mgL_to_mmolL(10.0, symbol) == pytest.approx(2.0)
        assert mmolL_to_mgL(2.0, symbol) == pytest.approx(10.0)

    def test_field_boundary_preserves_explicit_sr_ratio(self):
        record = _base_record(
            dataset="test",
            sample_id="S1",
            site_id="S1",
        )
        _fill_optional_diagnostics(record, {"87Sr/86Sr": 0.7085})
        assert record["sr_ratio_87_86"] == pytest.approx(0.7085)


class TestNitrateBoronBayesianForensics:
    """Validate 4D isotopic & trace forensics discriminating sewage vs animal manure."""

    def test_endmembers_database_split(self):
        """Check that nitrate_endmembers.json has Septic_Sewage and Animal_Manure with Boron."""
        db = load_endmember_database()
        assert "Septic_Sewage" in db
        assert "Animal_Manure" in db
        assert "Manure" in db  # backward compatibility alias

        sewage = db["Septic_Sewage"]
        manure = db["Animal_Manure"]

        # Septic sewage typically has higher boron concentration and lower d11B than manure
        assert sewage.d11B_mean is not None and manure.d11B_mean is not None
        assert sewage.d11B_mean < manure.d11B_mean
        assert sewage.ln_B_mean is not None and manure.ln_B_mean is not None
        assert sewage.ln_B_mean > manure.ln_B_mean

    def test_4d_likelihood_discrimination(self):
        """A sample with high boron and low d11B should be attributed to sewage rather than animal manure."""
        # Wastewater sample: d15N=12, d18O=5, d11B=6.0, B=500 ug/L (ln_B ~ 6.21)
        sample_sewage = IsotopeSample(
            sample_id="W1",
            d15N=12.0,
            d18O=5.0,
            d11B=6.0,
            B=500.0,
        )

        db = load_endmember_database()
        sewage_source = db["Septic_Sewage"]
        manure_source = db["Animal_Manure"]

        probs = compute_isotope_prob(sample_sewage, [sewage_source, manure_source])
        p_sewage = probs["Septic_Sewage"]
        p_manure = probs["Animal_Manure"]

        assert p_sewage > 10.0 * p_manure, f"p_sewage={p_sewage}, p_manure={p_manure}"

    def test_production_node_posterior_uses_boron_subsources(self):
        """The public node inference path must expose sewage/manure posteriors."""
        frame = pd.DataFrame(
            [
                {
                    "site_id": "W1",
                    "NO3": 1.0,
                    "Cl": 1.0,
                    "Ca": 1.0,
                    "Mg": 1.0,
                    "Na": 1.0,
                    "K": 1.0,
                    "HCO3": 1.0,
                    "SO4": 1.0,
                    "d15N": 12.0,
                    "d18O_NO3": 5.0,
                    "B_ug_L": 500.0,
                    "d11B": 6.0,
                }
            ]
        ).set_index("site_id", drop=False)
        cfg = Config()
        cfg.nitrate_isotope_include_subsources = True
        cfg.nitrate_isotope_boron_enabled = True

        result = infer_node_posteriors(frame, [], config=cfg)["W1"]

        assert result.boron_used is True
        assert result.d11b_used is True
        assert result.p_septic_sewage is not None
        assert result.p_animal_manure is not None
        assert result.p_septic_sewage > result.p_animal_manure
        assert "Septic_Sewage" in (result.source_fractions or {})

    def test_network_fit_propagates_boron_source_result_fields(self):
        """Network fitting must copy node posterior fields onto EdgeResult."""
        ions = {
            "Ca": 1.0,
            "Mg": 1.0,
            "Na": 1.0,
            "K": 1.0,
            "HCO3": 1.0,
            "Cl": 1.0,
            "SO4": 1.0,
            "NO3": 1.0,
            "F": 1.0,
            "Fe": 1.0,
            "PO4": 1.0,
        }
        samples = []
        for site_id, d15, boron, d11b in (
            ("A", 8.0, 100.0, 20.0),
            ("B", 12.0, 500.0, 6.0),
        ):
            samples.append(
                {
                    "site_id": site_id,
                    "EC": 100.0,
                    "TDS": 100.0,
                    "pH": 7.0,
                    "d15N": d15,
                    "d18O_NO3": 5.0,
                    "B_ug_L": boron,
                    "d11B": d11b,
                    **ions,
                }
            )
        cfg = Config(
            phreeqc_enabled=False,
            nitrate_source_enabled=True,
            transport_models_enabled=["evap"],
        )

        results = fit_network(samples, [Edge(edge_id="A->B", u="A", v="B")], cfg)

        assert len(results) == 1
        assert results[0].nitrate_source_boron_used is True
        assert results[0].nitrate_source_d11b_used is True
        assert results[0].nitrate_source_p_septic_sewage is not None
        assert results[0].nitrate_source_p_animal_manure is not None


class TestMineralsAndRedoxLadder:
    """Validate silicate stoichiometry, redox ladder ordering, and isotope diagnostics."""

    def test_silicate_stoichiometry_sio2(self):
        """Serpentine, talc, biotite, chlorite must include SiO2."""
        for mineral in ["serpentine", "talc", "biotite", "chlorite"]:
            assert mineral in MINERAL_LIBRARY
            assert "SiO2" in MINERAL_LIBRARY[mineral]
            assert MINERAL_LIBRARY[mineral]["SiO2"] > 0

    def test_manganese_reduction_reaction(self):
        """Manganese reduction stoichiometry should be present in reaction definitions."""
        assert "manganese_reduction" in MINERAL_LIBRARY
        rxn = MINERAL_LIBRARY["manganese_reduction"]
        assert rxn["Mn"] == 1.0
        assert rxn["HCO3"] == 1.0

        # Verify build_reaction_dictionary includes manganese_reduction when Mn and HCO3 are measured
        cfg = Config()
        cfg.measured_ions = ["Ca", "Mg", "Na", "K", "HCO3", "Cl", "SO4", "NO3", "Mn"]
        matrix, labels, mineral_mask, penalty_scales = build_reaction_dictionary(
            cfg, sample={"DO": 0.2, "NO3": 0.05}
        )
        assert "manganese_reduction" in labels

    def test_strontium_provenance_check(self):
        """Test Strontium isotope provenance invariant."""
        # Same source within tolerance
        inv, penalty = check_strontium_provenance(0.7085, 0.70852, tolerance=0.0002)
        assert inv is True
        assert penalty == 1.0

        # Different aquifer provenance
        inv_diff, penalty_diff = check_strontium_provenance(0.7085, 0.7120, tolerance=0.0002)
        assert inv_diff is False
        assert penalty_diff < 0.5

    def test_cl_br_source_diagnostic(self):
        """Test Cl/Br molar ratio classification."""
        # Halite dissolution typically Cl/Br > 1000
        assert check_cl_br_source(2000.0) == "halite_dissolution"
        # Sewage / agricultural wastewater typically 300 - 800
        assert check_cl_br_source(500.0) == "domestic_or_animal_wastewater"
        # Atmospheric recharge typically 50 - 150
        assert check_cl_br_source(100.0) == "precipitation_or_recharge"


class TestAdaptiveSheafStalks:
    """Validate dimension-varying cellular sheaf stalks and intersection restriction maps."""

    def test_adaptive_stalk_setup(self):
        """Build edge maps with mismatched species sets and verify intersection stalks."""
        edge_ab = Edge(edge_id="A->B", u="A", v="B")
        edge_bc = Edge(edge_id="B->C", u="B", v="C")

        # Node A has 3 ions, Node B has 5 ions, Node C has 3 ions
        # Edge A->B measures Ca, Mg, Cl
        em_ab = DirectedEdgeMap(
            edge=edge_ab,
            alpha=1.0,
            offset=[0.0, 0.0, 0.0],
            weight=1.0,
            objective=0.01,
            transport_model="mix",
            endmember_id=None,
            residual_norm=0.01,
            species=["Ca", "Mg", "Cl"],
        )
        # Edge B->C measures Ca, Cl, B
        em_bc = DirectedEdgeMap(
            edge=edge_bc,
            alpha=0.9,
            offset=[0.0, 0.0, 0.0],
            weight=1.0,
            objective=0.02,
            transport_model="mix",
            endmember_id=None,
            residual_norm=0.02,
            species=["Ca", "Cl", "B"],
        )

        edge_maps = [em_ab, em_bc]
        assert not _is_uniform_sheaf(edge_maps, 3)

        node_species, col_offsets, total_cols = _build_adaptive_stalks(edge_maps)
        assert "A" in node_species and "B" in node_species and "C" in node_species
        # Node B should contain all species involved in its incident edges: Ca, Mg, Cl, B
        assert set(node_species["B"]) == {"Ca", "Mg", "Cl", "B"}

        # Coboundary matrix D and rhs should build without errors
        D = build_coboundary_matrix(edge_maps)
        rhs = build_rhs_vector(edge_maps)
        assert D.shape[0] == 6  # 3 rows for AB + 3 rows for BC
        assert D.shape[1] == total_cols
        assert len(rhs) == 6

        # Cohomology computation must complete cleanly
        res = compute_cohomology(edge_maps)
        assert "obstruction_energy" in res
        assert "h0_dim" in res
        assert res["obstruction_energy"] >= 0.0

    def test_production_edge_builder_preserves_species_labels(self):
        """Normal edge construction must reach the adaptive sheaf contract."""
        cfg = Config()
        edge = Edge(edge_id="A->B", u="A", v="B")
        maps = build_edge_maps(
            [edge],
            {
                "A": {"Ca": 1.0, "Mg": 2.0, "Cl": 3.0, "B": 0.1},
                "B": {"Ca": 1.1, "Cl": 2.8, "B": 0.2},
            },
            cfg,
        )

        assert maps["A->B"].species == ["Ca", "Cl", "B"]
        assert len(maps["A->B"].offset) == 3

    def test_expanded_ion_order_uses_species_aligned_weights(self):
        """Expanded tracer panels must receive deterministic aligned weights."""
        cfg = Config(ion_order=["Ca", "Cl", "Br", "B", "Sr"])
        cfg.validate()

        assert len(cfg.weights) == 5
        assert len(cfg.conservative_weights) == 5
        assert cfg.get_conservative_weights()[1] == pytest.approx(1.0)
        assert cfg.get_conservative_weights()[2] == pytest.approx(2.0)
        assert cfg.get_conservative_weights()[3] == pytest.approx(0.01)

    def test_edge_result_attaches_cl_br_and_sr_diagnostics(self):
        """Fitted edge results must retain tracer diagnostics for downstream evidence."""
        cfg = Config(
            ion_order=["Ca", "HCO3", "Cl", "Br"],
            transport_models_enabled=["evap"],
            phreeqc_enabled=False,
            cl_br_ratio_enabled=True,
            sr_provenance_enabled=True,
        )
        result = fit_edge(
            [1.0, 1.0, 500.0, 1.0],
            [1.0, 1.0, 1000.0, 2.0],
            cfg,
            edge_id="A->B",
            u="A",
            v="B",
            obs_u={"sr_ratio_87_86": 0.7085},
            obs_v={"sr_ratio_87_86": 0.7120},
        )

        assert result.cl_br_source_u == "domestic_or_animal_wastewater"
        assert result.cl_br_source_v == "domestic_or_animal_wastewater"
        assert result.cl_br_metrics["ratio_cl_br_v"] == pytest.approx(500.0)
        assert result.sr_provenance_invariant is False
        assert result.sr_provenance_score_penalty > 0.0
        row = edge_results_table([result])[0]
        assert "cl_br_metrics" in row
        assert row["sr_provenance_invariant"] is False


class TestGeophysics3DAndConditioning:
    """Validate 3D geophysics conditioning, barrier checks, and transit times."""

    def test_geophysical_barrier_check(self):
        """Test ERT structural barrier and dyke penalties."""
        cfg = Config()
        cfg.geophysics_enabled = True
        cfg.geophysics_barrier_enabled = True
        cfg.geophysics_barrier_dyke_penalty = 0.05
        cfg.geophysics_sigma_struct = 50.0

        node_u = Node3D(
            node_id="N1", x=100.0, y=100.0, z=10.0, elevation_m=50.0,
            attrs={"apparent_resistivity_ohm_m": 40.0}
        )
        node_v_ok = Node3D(
            node_id="N2", x=200.0, y=200.0, z=10.0, elevation_m=50.0,
            attrs={"apparent_resistivity_ohm_m": 55.0}
        )
        node_v_barrier = Node3D(
            node_id="N3", x=300.0, y=300.0, z=10.0, elevation_m=50.0,
            attrs={"apparent_resistivity_ohm_m": 350.0, "geophysical_barrier": True}
        )

        # OK edge
        is_ok, p_ok = geophysical_barrier_check(node_u, node_v_ok, cfg)
        assert is_ok is True
        assert p_ok > 0.9

        # Barrier edge
        is_barrier, p_barrier = geophysical_barrier_check(node_u, node_v_barrier, cfg)
        assert is_barrier is False
        assert p_barrier < 0.1

    def test_bedrock_elevation_conversion_uses_positive_down_depth(self):
        """Bedrock elevations must be converted before comparing node depth."""
        cfg = Config()
        cfg.geophysics_enabled = True
        cfg.geophysics_barrier_enabled = True

        node_u = Node3D(
            node_id="N1", x=0.0, y=0.0, z=20.0, elevation_m=100.0,
            attrs={"bedrock_elevation_m": 35.0},
        )
        node_v = Node3D(
            node_id="N2", x=1.0, y=1.0, z=20.0, elevation_m=100.0,
            attrs={"bedrock_elevation_m": 35.0},
        )
        is_ok, p_ok = geophysical_barrier_check(node_u, node_v, cfg)
        assert is_ok is True
        assert p_ok == pytest.approx(1.0)

        node_v_deep = Node3D(
            node_id="N3", x=1.0, y=1.0, z=80.0, elevation_m=100.0,
            attrs={"bedrock_elevation_m": 35.0},
        )
        is_deep, p_deep = geophysical_barrier_check(node_u, node_v_deep, cfg)
        assert is_deep is False
        assert p_deep < p_ok

    def test_snmr_hydraulic_conductivity_and_travel_time(self):
        """Test SDR formula K = C_SDR * phi^4 * (T2*)^2 and advective transit time."""
        cfg = Config()
        cfg.geophysics_enabled = True
        cfg.geophysics_snmr_k_enabled = True
        cfg.geophysics_snmr_csdr = 1.0e-9

        # phi = 0.25, T2* = 120 ms
        node_u = Node3D(
            node_id="N1", x=0.0, y=0.0, z=10.0, elevation_m=100.0,
            hydraulic_head=100.0,
            attrs={"phi_sNMR": 0.25, "t2_star_ms": 120.0}
        )
        node_v = Node3D(
            node_id="N2", x=500.0, y=0.0, z=10.0, elevation_m=100.0,
            hydraulic_head=98.0,
            attrs={"phi_sNMR": 0.25, "t2_star_ms": 120.0}
        )

        d_3d = 500.0
        delta_h = 2.0

        k_m_day, tau_years = compute_geophysical_transit_time(node_u, node_v, d_3d, delta_h, cfg)

        assert k_m_day is not None and k_m_day > 0.0
        assert tau_years is not None and tau_years > 0.0

        # Verify manual SDR calculation:
        # K_m_s = 1e-9 * (0.25^4) * (120^2) = 1e-9 * 0.00390625 * 14400 = 5.625e-8 m/s
        # K_m_day = 5.625e-8 * 86400 = 0.00486 m/day
        expected_k_day = 1e-9 * (0.25**4) * (120.0**2) * 86400.0
        assert abs(k_m_day - expected_k_day) < 1e-6

    def test_explicit_geophysical_k_without_porosity_does_not_make_up_transit_time(self):
        """A direct K observation still needs porosity for an advective tau."""
        cfg = Config()
        cfg.geophysics_enabled = True
        cfg.geophysics_snmr_k_enabled = True
        node_u = Node3D(
            node_id="N1", x=0.0, y=0.0, z=10.0, elevation_m=100.0,
            attrs={"k_geophys_m_day": 0.5},
        )
        node_v = Node3D(
            node_id="N2", x=100.0, y=0.0, z=10.0, elevation_m=100.0,
            attrs={"k_geophys_m_day": 0.5},
        )

        k_m_day, tau_years = compute_geophysical_transit_time(
            node_u, node_v, 100.0, 1.0, cfg
        )
        assert k_m_day == pytest.approx(0.5)
        assert tau_years is None

    def test_3d_builder_preserves_sample_geophysics_and_no_hydraulic_fallback(self):
        """The normal 3D builder must retain sNMR observations at node/edge level."""
        cfg = Config()
        cfg.geophysics_enabled = True
        cfg.geophysics_snmr_k_enabled = True
        cfg.edge_p_min = 0.1
        cfg.edge_radius_km = 5.0

        network = build_network_3d(
            [
                {
                    "site_id": "N1",
                    "x": 0.0,
                    "y": 0.0,
                    "screen_depth": 10.0,
                    "elevation": 100.0,
                    "head_meas": 100.0,
                    "phi_sNMR": 0.25,
                    "t2_star_ms": 120.0,
                },
                {
                    "site_id": "N2",
                    "x": 100.0,
                    "y": 0.0,
                    "screen_depth": 10.0,
                    "elevation": 100.0,
                    "head_meas": 99.0,
                    "phi_sNMR": 0.25,
                    "t2_star_ms": 120.0,
                },
            ],
            cfg,
            use_haversine=False,
        )

        assert network.nodes["N1"].attrs["phi_sNMR"] == pytest.approx(0.25)
        assert network.edges
        assert any(edge.k_geophys_m_day is not None for edge in network.edges)

        no_fallback = Node3D(
            node_id="N3", x=0.0, y=0.0, z=10.0, elevation_m=100.0,
            attrs={"effective_porosity": 0.25},
        )
        no_fallback_v = Node3D(
            node_id="N4", x=100.0, y=0.0, z=10.0, elevation_m=100.0,
            attrs={"effective_porosity": 0.25},
        )
        k_missing, tau_missing = compute_geophysical_transit_time(
            no_fallback, no_fallback_v, 100.0, 1.0, cfg
        )
        assert k_missing is None
        assert tau_missing is None

    def test_network_fit_edge_adapter_preserves_geophysical_contract(self):
        """The ordinary infer_edges API must not discard 3-D geophysical fields."""
        cfg = Config()
        cfg.network_3d_enabled = True
        cfg.geophysics_enabled = True
        cfg.geophysics_snmr_k_enabled = True
        cfg.edge_p_min = 0.05
        cfg.edge_radius_km = 5.0
        edges = infer_edges(
            [
                {
                    "site_id": "A",
                    "x": 0.0,
                    "y": 0.0,
                    "screen_depth": 10.0,
                    "elevation": 100.0,
                    "head_meas": 100.0,
                    "phi_sNMR": 0.25,
                    "t2_star_ms": 120.0,
                    "formation_conductivity": 100.0,
                    "geophysical_k_cv": 0.6,
                },
                {
                    "site_id": "B",
                    "x": 0.001,
                    "y": 0.0,
                    "screen_depth": 10.0,
                    "elevation": 100.0,
                    "head_meas": 99.0,
                    "phi_sNMR": 0.25,
                    "t2_star_ms": 120.0,
                    "formation_conductivity": 500.0,
                    "geophysical_k_cv": 0.7,
                },
            ],
            config=cfg,
        )
        assert len(edges) == 1
        attrs = edges[0].attrs
        assert attrs["physics_source"] == "geophysics_snmr"
        assert attrs["k_geophys_m_day"] == pytest.approx(0.00486)
        assert attrs["formation_conductivity_u"] == pytest.approx(100.0)
        assert attrs["formation_conductivity_v"] == pytest.approx(500.0)
        assert attrs["geophysical_k_cv"] == pytest.approx(0.7)

    def test_geophysical_table_adapter_attaches_node_observations(self, tmp_path):
        """Configured JSON/CSV-style observations must reach Node3D attributes."""
        source = tmp_path / "ert_observations.json"
        source.write_text(
            '{"observations": [{"node_id": "N1", '
            '"apparent_resistivity_ohm_m": 350, '
            '"bedrock_elevation_m": 35.0}]}',
            encoding="utf-8",
        )
        observations = load_geophysical_observations(source)
        nodes = attach_geophysical_observations(
            [Node3D(node_id="N1", x=0.0, y=0.0, z=20.0, elevation_m=100.0)],
            observations,
        )

        assert nodes[0].attrs["apparent_resistivity_ohm_m"] == 350.0
        assert nodes[0].attrs["bedrock_elevation_m"] == 35.0

    def test_topology_posterior_prior_geophysics(self):
        """Test geophysics barrier modulation of edge prior in topology posterior."""
        cfg = Config()
        cfg.geophysics_enabled = True
        cfg.geophysics_barrier_dyke_penalty = 0.02
        cfg.geophysics_sigma_struct = 50.0

        # Normal edge
        edge_normal = Edge(
            u="A", v="B", edge_id="A->B",
            attrs={"p_uv": 0.8, "apparent_resistivity_u": 50.0, "apparent_resistivity_v": 60.0}
        )
        p_normal = _edge_prior_probability(edge_normal, cfg)
        assert abs(p_normal - 0.8) < 1e-4

        # Edge crossing a detected dyke / barrier
        edge_dyke = Edge(
            u="A", v="C", edge_id="A->C",
            attrs={"p_uv": 0.8, "geophysical_barrier": True}
        )
        p_dyke = _edge_prior_probability(edge_dyke, cfg)
        assert abs(p_dyke - (0.8 * 0.02)) < 1e-4

    def test_topology_prior_consumes_fluid_ec_alias_and_sr_invariant(self):
        """Fluid-EC and Sr endpoint diagnostics must modulate topology priors."""
        ec_cfg = Config()
        ec_cfg.geophysics_enabled = True
        ec_cfg.geophysics_fluid_ec_enabled = True
        ec_cfg.geophysics_fluid_ec_weight = 0.5
        ec_edge = Edge(
            u="A", v="B", edge_id="A->B",
            attrs={
                "p_uv": 0.8,
                "formation_conductivity_u": 100.0,
                "formation_conductivity_v": 500.0,
            },
        )
        assert _edge_prior_probability(ec_edge, ec_cfg) == pytest.approx(0.4)

        sr_cfg = Config()
        sr_cfg.sr_provenance_enabled = True
        sr_cfg.sr_provenance_tolerance = 0.0002
        sr_edge = Edge(
            u="A", v="C", edge_id="A->C",
            attrs={
                "p_uv": 0.8,
                "sr_ratio_u": 0.7085,
                "sr_ratio_v": 0.7120,
            },
        )
        assert _edge_prior_probability(sr_edge, sr_cfg) < 0.8

    def test_active_learning_campaign_geophysics_and_boron(self):
        """Test active learning measurement options and costs for boron and geophysics."""
        assert "trace_elements_boron" in DEFAULT_MEASUREMENT_COSTS
        assert "nitrate_boron_isotopes" in DEFAULT_MEASUREMENT_COSTS
        assert "geophysics_ert_survey" in DEFAULT_MEASUREMENT_COSTS
        assert "geophysics_snmr" in DEFAULT_MEASUREMENT_COSTS

        assert "geophysical_barrier_ambiguity" in FLAG_TO_MEASUREMENTS
        assert "nitrate_source_ambiguity" in FLAG_TO_MEASUREMENTS
        assert "salinity_source_ambiguity" in FLAG_TO_MEASUREMENTS

        # Test surrogate prediction
        context_ert = {
            "well_id": "W1",
            "measurement_type": "geophysics_ert_survey",
            "active_edge_ids": ["W1->W2"],
            "candidate_edges": [Edge(u="W1", v="W2", edge_id="W1->W2")],
        }
        mean_ert, sd_ert = topology_contrast_surrogate(context_ert)
        assert mean_ert == 45.0
        assert sd_ert == 5.0

    def test_diagnostic_flags_reach_active_learning_scoring(self):
        edge = Edge(
            u="W1",
            v="W2",
            edge_id="W1->W2",
            attrs={
                "cl_br_metrics": {
                    "source_class_u": "precipitation_or_recharge",
                    "source_class_v": "domestic_or_animal_wastewater",
                },
                "sr_provenance_invariant": False,
            },
        )
        score, reason = _score_evidence_ambiguity(edge)
        assert score >= 0.6
        assert "salinity_source_ambiguity" in reason

    def test_geophysical_uncertainty_reaches_active_learning_utility(self):
        edge = Edge(
            u="W1",
            v="W2",
            edge_id="W1->W2",
            attrs={
                "geophysics_enabled": True,
                "k_geophys_m_day": 0.1,
                "k_geophys_std_m_day": 0.08,
            },
        )
        score = _score_geophysical_uncertainty(edge)
        assert score == pytest.approx(0.8)
        score_evidence, reason = _score_evidence_ambiguity(edge)
        assert score_evidence >= 0.6
        assert "geophysical_uncertainty" in reason
