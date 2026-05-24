"""Utility module for plotting configuration and metadata handling."""

import json
from dataclasses import dataclass, field, asdict
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

import matplotlib.pyplot as plt
from ..log import get_logger

logger = get_logger(__name__)

@dataclass
class PlotConfig:
    """Configuration for plot styling and output."""
    style: str = "seaborn-v0_8-whitegrid"
    palette: str = "colorblind"
    font_scale: float = 1.0
    figsize: Tuple[float, float] = (10.0, 6.0)
    dpi: int = 300
    file_format: str = "png"
    transparent: bool = False
    
    # Metadata context
    context: Dict[str, Any] = field(default_factory=dict)

    def apply(self):
        """Apply style settings to matplotlib context."""
        try:
            plt.style.use(self.style)
        except OSError:
            # Fallback if style not found
            plt.style.use("default")
            
        import seaborn as sns
        try:
            sns.set_palette(self.palette)
            sns.set_context("paper", font_scale=self.font_scale)
        except ImportError:
            pass

def save_with_metadata(
    fig: plt.Figure,
    path: str,
    config: PlotConfig,
    extra_metadata: Optional[Dict[str, Any]] = None
) -> None:
    """Save figure with sidecar JSON metadata."""
    if not path:
        return

    path_obj = Path(path)
    
    # Ensure directory exists
    path_obj.parent.mkdir(parents=True, exist_ok=True)
    
    # Save image
    fig.savefig(
        path,
        dpi=config.dpi,
        bbox_inches="tight",
        transparent=config.transparent,
        format=config.file_format or path_obj.suffix.lstrip(".")
    )
    logger.info(f"Saved plot to {path}")

    # Save metadata
    meta_path = path_obj.with_suffix(".json")
    metadata = {
        "timestamp": datetime.now().isoformat(),
        "file_path": str(path_obj),
        "config": asdict(config),
        "extra_metadata": extra_metadata or {}
    }
    
    try:
        with open(meta_path, "w", encoding="utf-8") as f:
            json.dump(metadata, f, indent=2)
    except IOError as e:
        logger.warning(f"Failed to save plot metadata to {meta_path}: {e}")
