"""Exact small-model checks of conditional hydrochemical ACMD certificates."""

from __future__ import annotations

from dataclasses import replace
import hashlib
import math

import numpy as np
import pytest

from hydrosheaf import ACMDStatus, ChemACMD, ChemicalAction
from hydrosheaf.config import Config
from hydrosheaf.reactive_transport import (
    build_isotope_mass_action,
    compile_geochemical_polytope,
    reaction_stoichiometry_from_config,
)


def _confounded_sulfate_system():
    # Gypsum and anhydrite have identical aqueous Ca/SO4 stoichiometry.
    return compile_geochemical_polytope(
        major_ions_u={"Ca": 0.0, "SO4": 0.0},
        major_ions_v={"Ca": 1.0, "SO4": 1.0},
        reaction_stoichiometry={
            "gypsum": {"Ca": 1.0, "SO4": 1.0},
            "anhydrite": {"Ca": 1.0, "SO4": 1.0},
        },
        extent_bounds={"gypsum": (0.0, 1.1), "anhydrite": (0.0, 1.1)},
        ion_error_bounds={"Ca": 0.0, "SO4": 0.0},
    )


def _loop(actions, tolerance=0.05):
    return ChemACMD.create_loop(
        _confounded_sulfate_system(),
        actions,
        {"extent:gypsum": 1.0},
        tolerance,
        target_name="gypsum_extent",
    )


def test_stoichiometric_ambiguity_and_isotope_mass_certification():
    polytope = _confounded_sulfate_system()
    isotope = build_isotope_mass_action(
        polytope,
        action_id="sulfur_isotope",
        target_well="W1",
        isotope_system="34S",
        element="SO4",
        isotope_atoms_per_species=1.0,
        # Synthetic, site-declared fractions; these are not universal signatures.
        source_atom_fractions={"extent:gypsum": 0.04, "extent:anhydrite": 0.06},
        error_bound=0.0002,
        lab_cost=150.0,
    )
    redundant = ChemicalAction(
        action_id="repeat_sulfate",
        measurement_type="SO4",
        target_well="W1",
        response={"extent:gypsum": 1.0, "extent:anhydrite": 1.0},
        error_bound=0.001,
        lab_cost=50.0,
    )
    loop = ChemACMD.create_loop(
        polytope, [redundant, isotope], {"extent:gypsum": 1.0},
        0.05, target_name="gypsum_extent",
    )
    assert loop.status == ACMDStatus.ACTIVE
    assert math.isclose(loop.current_ambiguity, 1.0, abs_tol=1e-7)
    action, diagnostics = loop.recommend_next_action()
    assert action.action_id == "sulfur_isotope"
    assert diagnostics["details"]["repeat_sulfate"]["reduction"] == pytest.approx(0.0)

    # Equal reaction extents produce 0.05 mmol/L heavy-isotope equivalent.
    loop.step(action.action_id, observed_value=0.05)
    assert loop.status == ACMDStatus.CERTIFIED
    assert loop.current_ambiguity == pytest.approx(0.02, abs=1e-7)
    certificate = loop.get_certificate()
    assert len(certificate.certificate_hash) == 64
    assert certificate.verify_hash()
    assert not replace(certificate, certificate_hash="0" * 64).verify_hash()
    assert len(certificate.metadata["problem_hash"]) == 64
    assert certificate.certificate_hash != certificate.metadata["problem_hash"]
    assert all(char in "0123456789abcdef" for char in certificate.certificate_hash)
    assert certificate.steps[0].observed_value == pytest.approx(0.05)
    assert hashlib.sha256(b"unrelated").hexdigest() != certificate.certificate_hash
    same_design = ChemACMD.create_loop(
        polytope, [redundant, isotope], {"extent:gypsum": 1.0},
        0.05, target_name="gypsum_extent",
    )
    same_design.step("sulfur_isotope", observed_value=0.052)
    assert same_design.get_certificate().certificate_hash != certificate.certificate_hash


def test_impossibility_witness_satisfies_mass_balance_and_candidate():
    redundant = ChemicalAction(
        action_id="repeat_sulfate",
        measurement_type="SO4",
        target_well="W1",
        response={"extent:gypsum": 1.0, "extent:anhydrite": 1.0},
        error_bound=0.001,
        lab_cost=50.0,
    )
    loop = _loop([redundant], tolerance=0.02)
    action, diagnostics = loop.recommend_next_action()
    assert action is None
    assert diagnostics["status"] == "IMPOSSIBILITY_WITNESS"
    certificate = loop.get_certificate()
    assert certificate.status == ACMDStatus.IMPOSSIBILITY_WITNESS
    low, high = certificate.lower_witness, certificate.upper_witness
    assert high[0] - low[0] > 0.02
    for state in (low, high):
        assert np.isclose(state[0] + state[1], 1.0)
        assert np.isclose(state[2], 1.0)
        assert np.all(state >= -1e-8)
    assert abs((high[0] + high[1]) - (low[0] + low[1])) <= 0.002


def test_inconsistent_initial_and_new_evidence_are_distinct():
    impossible = compile_geochemical_polytope(
        {"Ca": 0.0, "SO4": 0.0},
        {"Ca": 1.0, "SO4": 0.2},
        {"gypsum": {"Ca": 1.0, "SO4": 1.0}},
        {"gypsum": (0.0, 2.0)},
        {"Ca": 0.0, "SO4": 0.0},
    )
    initial = ChemACMD.create_loop(
        impossible, [], {"extent:gypsum": 1.0}, 0.1,
        target_name="gypsum_extent",
    )
    assert initial.status == ACMDStatus.INCONSISTENT_EVIDENCE

    isotope = ChemicalAction(
        action_id="isotope", measurement_type="isotope_mass_34S",
        target_well="W1", response={"extent:gypsum": 0.04},
        error_bound=0.0001, lab_cost=10.0,
    )
    loop = _loop([isotope])
    loop.step("isotope", observed_value=0.2)
    assert loop.status == ACMDStatus.INCONSISTENT_EVIDENCE
    assert math.isnan(loop.current_ambiguity)
    assert loop.get_certificate().status == ACMDStatus.INCONSISTENT_EVIDENCE


def test_mixing_simplex_and_absolute_exchange_cap():
    mixture = compile_geochemical_polytope(
        {"Ca": 0.0}, {"Ca": 0.25},
        {"inert": {}}, {"inert": (0.0, 0.0)}, {"Ca": 0.0},
        endmembers={"recharge_B": {"Ca": 1.0}},
    )
    mix_loop = ChemACMD.create_loop(
        mixture, [], {"mix:recharge_B": 1.0}, 0.01,
        target_name="recharge_B_fraction", target_units="fraction",
    )
    assert mix_loop.status == ACMDStatus.ALREADY_RESOLVED
    assert mix_loop.current_ambiguity == pytest.approx(0.0)

    capped = compile_geochemical_polytope(
        {"Ca": 0.0}, {"Ca": 0.0},
        {"exchange_A": {}, "exchange_B": {}},
        {"exchange_A": (-1.0, 1.0), "exchange_B": (-1.0, 1.0)},
        {"Ca": 0.0},
        absolute_extent_cap=(("exchange_A", "exchange_B"), 0.5),
    )
    cap_loop = ChemACMD.create_loop(
        capped, [], {"extent:exchange_A": 1.0}, 0.1,
        target_name="exchange_A_extent",
    )
    assert cap_loop.current_ambiguity == pytest.approx(1.0)


def test_colocated_well_mobilization_discount_and_model_guards():
    coarse = ChemicalAction(
        action_id="coarse", measurement_type="sulfate_mass", target_well="W1",
        response={"extent:gypsum": 1.0}, error_bound=0.1,
        lab_cost=1.0, field_travel_cost=20.0,
    )
    precise = ChemicalAction(
        action_id="precise", measurement_type="sulfate_mass", target_well="W1",
        response={"extent:gypsum": 1.0}, error_bound=0.001,
        lab_cost=100.0, field_travel_cost=20.0,
    )
    precise_other_well = ChemicalAction(
        action_id="precise_W2", measurement_type="sulfate_mass", target_well="W2",
        response={"extent:gypsum": 1.0}, error_bound=0.001,
        lab_cost=100.0, field_travel_cost=20.0,
    )
    loop = _loop([coarse, precise, precise_other_well], tolerance=0.01)
    chosen, _ = loop.recommend_next_action()
    assert chosen.action_id == "coarse"
    loop.step("coarse", observed_value=0.5)
    assert loop.status == ACMDStatus.ACTIVE
    assert loop.cumulative_cost == pytest.approx(21.0)
    assert loop.candidates["precise"].incremental_cost(loop.visited_targets) == pytest.approx(100.0)
    assert loop.candidates["precise_W2"].incremental_cost(loop.visited_targets) == pytest.approx(120.0)
    next_action, _ = loop.recommend_next_action()
    assert next_action.action_id == "precise"

    with pytest.raises(ValueError, match="explicit response"):
        ChemicalAction(
            action_id="bad", measurement_type="34S", target_well="W1",
            response={}, error_bound=0.01, lab_cost=1.0,
        )
    with pytest.raises(ValueError, match="unknown states"):
        ChemACMD.create_loop(
            _confounded_sulfate_system(), [ChemicalAction(
                action_id="bad", measurement_type="34S", target_well="W1",
                response={"extent:unknown": 1.0}, error_bound=0.01, lab_cost=1.0,
            )], {"extent:gypsum": 1.0}, 0.1, target_name="gypsum_extent",
        )

    guarded = _loop([coarse], tolerance=0.01)
    with pytest.raises(ValueError, match="observed_value must be finite"):
        guarded.step("coarse", float("nan"))
    with pytest.raises(ValueError, match="different LP intervals"):
        guarded.step("coarse", 0.5, error_bound=0.1, sigma=0.1)
    assert "coarse" in guarded.candidates
    assert guarded.cumulative_cost == 0.0


def test_isotope_helper_requires_complete_signatures_and_refuses_sinks():
    polytope = _confounded_sulfate_system()
    with pytest.raises(ValueError, match="exactly these contributing states"):
        build_isotope_mass_action(
            polytope, action_id="bad", target_well="W1", isotope_system="34S",
            element="SO4", isotope_atoms_per_species=1.0,
            source_atom_fractions={"extent:gypsum": 0.04},
            error_bound=0.001, lab_cost=1.0,
        )
    sinking = compile_geochemical_polytope(
        {"SO4": 1.0}, {"SO4": 0.5},
        {"sulfate_reduction": {"SO4": -1.0}},
        {"sulfate_reduction": (0.0, 1.0)}, {"SO4": 0.0},
    )
    with pytest.raises(ValueError, match="removal or fractionation"):
        build_isotope_mass_action(
            sinking, action_id="bad", target_well="W1", isotope_system="34S",
            element="SO4", isotope_atoms_per_species=1.0,
            source_atom_fractions={}, error_bound=0.001,
            lab_cost=1.0,
        )


def test_existing_reaction_dictionary_orientation_is_preserved():
    config = Config(active_minerals=["gypsum"])
    reactions = reaction_stoichiometry_from_config(config)
    assert reactions["gypsum"]["Ca"] == pytest.approx(1.0)
    assert reactions["gypsum"]["SO4"] == pytest.approx(1.0)


def test_isotope_mass_response_counts_atoms_in_carrier_species():
    nitrate = compile_geochemical_polytope(
        {"NO3": 0.0}, {"NO3": 1.0},
        {"nitrate_source": {"NO3": 1.0}},
        {"nitrate_source": (0.0, 1.0)}, {"NO3": 0.0},
    )
    oxygen = build_isotope_mass_action(
        nitrate, action_id="18O", target_well="W1", isotope_system="18O",
        element="NO3", isotope_atoms_per_species=3.0,
        source_atom_fractions={"extent:nitrate_source": 0.01},
        error_bound=0.001, lab_cost=10.0,
    )
    assert oxygen.response["extent:nitrate_source"] == pytest.approx(0.03)
