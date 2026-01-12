"""Graph dataclasses."""

from dataclasses import dataclass, field
from typing import Dict


@dataclass
class Node:
    node_id: str
    attrs: Dict[str, object] = field(default_factory=dict)


@dataclass
class Edge:
    edge_id: str
    u: str
    v: str
    attrs: Dict[str, object] = field(default_factory=dict)
