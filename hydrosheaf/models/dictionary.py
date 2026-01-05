"""Reaction dictionary helpers."""

from typing import List, Tuple

from ..config import Config
from .reactions import build_reaction_dictionary as _build


def build_reaction_dictionary(config: Config) -> Tuple[List[List[float]], List[str], List[bool]]:
    return _build(config)
