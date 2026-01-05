"""PHREEQC input templates."""

from typing import Mapping


def build_solution_block(sample: Mapping[str, object], solution_id: int, temp_default_c: float) -> str:
    temp_c = sample.get("temp_c")
    if temp_c is None:
        temp_c = temp_default_c
    lines = [
        f"SOLUTION {solution_id}",
        f"temp      {temp_c}",
        f"pH        {sample.get('pH')}",
        "units     mg/L",
        f"Ca        {sample.get('Ca')}",
        f"Mg        {sample.get('Mg')}",
        f"Na        {sample.get('Na')}",
        f"Cl        {sample.get('Cl')}",
        f"S(6)      {sample.get('SO4')}   as SO4",
        f"N(5)      {sample.get('NO3')}   as NO3",
        f"F         {sample.get('F')}",
        f"Alkalinity {sample.get('HCO3')} as HCO3",
        "END",
    ]
    return "\n".join(lines)


def build_selected_output_block() -> str:
    return """
SELECTED_OUTPUT
    -file                 selected_output.csv
    -reset                false
    -high_precision       true
    -user_punch           true
    -pH                   true
    -ionic_strength       true
    -charge_balance       true
    -totals               Ca Mg Na Cl S(6) N(5) F Alkalinity
    -saturation_indices   Calcite Dolomite Gypsum Anhydrite Halite Fluorite

USER_PUNCH
    -headings  sample_id  pH  I  charge_balance  SI_Calcite  SI_Dolomite  SI_Gypsum  SI_Anhydrite  SI_Halite  SI_Fluorite
10  PUNCH      SOLN,     -LA("H+"),  MU,  CB,
20  PUNCH      SI("Calcite"), SI("Dolomite"), SI("Gypsum"), SI("Anhydrite"), SI("Halite"), SI("Fluorite")
""".strip()
