"""Repository-relative paths for the public analysis code."""

from __future__ import annotations

import os
import tempfile
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]


def get_data_root() -> Path:
    """Return the authorized CHARLS data root configured by the researcher."""
    configured = os.environ.get("CHARLS_DATA_ROOT", "").strip()
    if configured:
        return Path(configured).expanduser().resolve()
    return REPOSITORY_ROOT / "data" / "raw"


def get_derived_dir() -> Path:
    """Return the local, version-control-excluded derived-data directory."""
    return REPOSITORY_ROOT / "data" / "derived"


def get_results_dir() -> Path:
    """Return the local analysis-output directory."""
    return REPOSITORY_ROOT / "results"


def get_temp_dir() -> Path:
    """Return a temporary directory used to unpack the harmonized data file."""
    return Path(tempfile.gettempdir()) / "charls_hcap_analysis"
