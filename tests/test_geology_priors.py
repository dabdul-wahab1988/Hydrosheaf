from __future__ import annotations

import pytest

from hydrosheaf.config import Config
from hydrosheaf.models.geology_context import (
    add_boundary_confidence,
    compare_geology_contexts,
    normalise_geology_context,
)
from hydrosheaf.models.geology_priors import (
    apply_geology_penalty_scales,
    geology_reaction_prior,
    load_geology_family_multipliers,
)
from hydrosheaf.models.reactions import build_reaction_dictionary


def test_field_geology_aliases_use_declared_precedence() -> None:
    context = normalise_geology_context(
        {
            "geology_join_status": "MATCHED",
            "geology_stratigraphic_unit": "  S1  ",
            "geology_symbol": "SYM",
            "geology_code_1000": "C1000",
            "geology_legend_text": "Basement",
        }
    )

    assert context.status == "matched"
    assert context.unit == "S1"

    fallback = normalise_geology_context(
        {
            "geology_join_status": "MATCHED",
            "geology_stratigraphic_unit": None,
            "geology_symbol": "SYM",
            "geology_code_1000": "C1000",
        }
    )
    assert fallback.unit == "SYM"


def test_no_match_and_unknown_contexts_are_neutral() -> None:
    profile = {"S1": {"carbonate": 0.25, "silicate": 4.0}}
    for context in (
        {"geology_join_status": "NO_MATCH", "geology_symbol": "stale"},
        {"geology_join_status": "NOT_JOINED", "geology_symbol": "stale"},
        None,
    ):
        prior = geology_reaction_prior(
            "calcite",
            context,
            family_multipliers=profile,
        )
        assert prior.is_neutral
        assert prior.multiplier == pytest.approx(1.0)
        assert prior.penalty_scale == pytest.approx(1.0)


def test_explicit_matched_context_applies_soft_family_multiplier() -> None:
    profile = {"S1": {"carbonate": 0.25, "silicate": 4.0}}
    scales, priors = apply_geology_penalty_scales(
        ["calcite", "albite", "unknown_reaction"],
        [1.0, 2.0, 3.0],
        {"geology_join_status": "MATCHED", "geology_symbol": "S1"},
        family_multipliers=profile,
    )

    assert scales == pytest.approx([4.0, 0.5, 3.0])
    assert priors["calcite"].status == "matched_prior"
    assert priors["albite"].status == "matched_prior"
    assert priors["unknown_reaction"].status == "neutral_unmapped"


def test_geology_pair_comparison_does_not_assign_unknown_units() -> None:
    result = compare_geology_contexts(
        {"geology_join_status": "MATCHED", "geology_symbol": "S1"},
        {"geology_join_status": "NO_MATCH", "geology_symbol": "stale"},
    )

    assert result.relation == "unknown"
    assert result.status == "ambiguous"
    assert result.same_unit is None
    assert result.neutral_prior_required is True


def test_reaction_dictionary_keeps_tuple_shape_and_honours_explicit_context() -> None:
    config = Config(
        active_minerals=["calcite", "albite"],
        exchange_enabled=False,
        use_thermodynamic_logic_gates=False,
        reaction_processes_enabled=[],
    )
    legacy = build_reaction_dictionary(config)
    neutral = build_reaction_dictionary(
        config,
        geology_context={"geology_join_status": "NO_MATCH", "geology_symbol": "S1"},
        geology_family_multipliers={"S1": {"carbonate": 0.25, "silicate": 4.0}},
    )
    matched = build_reaction_dictionary(
        config,
        geology_context={"geology_join_status": "MATCHED", "geology_symbol": "S1"},
        geology_family_multipliers={"S1": {"carbonate": 0.25, "silicate": 4.0}},
    )

    assert len(legacy) == 4
    assert len(neutral) == 4
    assert len(matched) == 4
    _, labels, _, legacy_scales = legacy
    _, neutral_labels, _, neutral_scales = neutral
    _, matched_labels, _, matched_scales = matched
    assert labels == neutral_labels == matched_labels
    calcite = labels.index("calcite")
    albite = labels.index("albite")
    assert neutral_scales[calcite] == pytest.approx(1.0)
    assert neutral_scales[albite] == pytest.approx(1.0)
    assert matched_scales[calcite] == pytest.approx(4.0)
    assert matched_scales[albite] == pytest.approx(0.25)
    # The old omitted-context call retains the historical Config bias.
    assert legacy_scales[calcite] == pytest.approx(50.0)
    assert legacy_scales[albite] == pytest.approx(0.5)


def test_user_geology_dictionary_is_loaded_and_preferred_key_is_honoured(tmp_path) -> None:
    dictionary = tmp_path / "geology_dictionary.json"
    dictionary.write_text(
        '{"family_multipliers": {"tmbt": {"multipliers": {"silicate": 2.0, "carbonate": 0.5}}}}',
        encoding="utf-8",
    )
    profile = load_geology_family_multipliers(dictionary)
    assert profile["tmbt"]["silicate"] == pytest.approx(2.0)
    prior = geology_reaction_prior(
        "calcite",
        {
            "geology_join_status": "MATCHED",
            "geology_symbol": "tmbt",
            "geology_stratigraphic_unit": "a different descriptive name",
        },
        family_multipliers=profile,
        preferred_key="geology_symbol",
    )
    assert prior.status == "matched_prior"
    assert prior.multiplier == pytest.approx(0.5)

    config = Config(
        active_minerals=["calcite", "albite"],
        exchange_enabled=False,
        use_thermodynamic_logic_gates=False,
        reaction_processes_enabled=[],
        geology_dictionary_path=str(dictionary),
        geology_prior_strength=1.0,
    )
    _, labels, _, scales = build_reaction_dictionary(
        config,
        sample={
            "geology_join_status": "MATCHED",
            "geology_symbol": "tmbt",
        },
    )
    assert scales[labels.index("calcite")] == pytest.approx(2.0)
    assert scales[labels.index("albite")] == pytest.approx(0.5)


def test_boundary_confidence_softens_but_does_not_create_a_geology_match() -> None:
    context = add_boundary_confidence(
        {
            "geology_join_status": "MATCHED",
            "geology_symbol": "S1",
            "geology_boundary_distance_m": 0.0,
        },
        sigma_m=500.0,
    )
    assert context["geology_confidence"] == pytest.approx(0.0)
    prior = geology_reaction_prior(
        "calcite",
        context,
        family_multipliers={"S1": {"carbonate": 0.25}},
        strength=1.0,
        preferred_key="geology_symbol",
    )
    assert prior.status == "matched_prior"
    assert prior.multiplier == pytest.approx(1.0)
