#!/usr/bin/env python3
"""Build the local participant-level input for HCAP measurement modelling."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np

from project_paths import get_derived_dir, get_results_dir


SCRIPT = Path(__file__).resolve().parent / "11_construct_hcap_indicators.py"
OUTPUT = get_derived_dir() / "b1_hcap_fiml_input_2026-06-19.csv"
MANIFEST = get_results_dir() / "manifests" / "hcap_model_input.json"

INDICATORS = [
    "test_short_learning_trial1",
    "test_short_delayed",
    "test_ten_learning_trial1",
    "test_ten_learning_trial2",
    "test_ten_learning_trial3",
    "test_ten_delayed",
    "test_recognition",
    "test_serial7",
    "test_animal_fluency",
    "test_language",
    "test_visuospatial",
    "test_orientation",
]


def load_indicator_module():
    """Load the indicator-construction script without duplicating its logic."""
    specification = importlib.util.spec_from_file_location(
        "hcap_indicators",
        SCRIPT,
    )
    if specification is None or specification.loader is None:
        raise RuntimeError(f"Cannot import indicator module: {SCRIPT}")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


def main() -> None:
    indicators = load_indicator_module()
    data = indicators.build_scores()
    data = data.loc[data["age60"]].copy()

    standardized_columns = []
    scaling = {}
    for indicator in INDICATORS:
        values = data[indicator]
        mean = float(values.mean())
        standard_deviation = float(values.std())
        if not np.isfinite(standard_deviation) or standard_deviation <= 0:
            raise ValueError(f"Invalid SD for {indicator}: {standard_deviation}")
        standardized_name = f"z_{indicator}"
        data[standardized_name] = (values - mean) / standard_deviation
        standardized_columns.append(standardized_name)
        scaling[indicator] = {
            "mean": mean,
            "sd": standard_deviation,
            "n_observed": int(values.notna().sum()),
            "missing_pct": round(float(values.isna().mean() * 100), 3),
        }

    keep = [
        "ID",
        "r4agey",
        "female",
        "edu3",
        "rural",
        "vision_problem",
        "hearing_problem",
        "sensory_impairment_count",
        "dual_sensory_impairment",
        "r4wtrespb",
        "r4adlab_c",
        "r4iadla",
        "r4imrc",
        "r4dlrc",
        "r4ser7",
        "r4orient",
        "r4draw",
        "iqcode_mean",
        "blessed_score",
    ] + standardized_columns
    output = data[keep].copy()
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(OUTPUT, index=False)

    manifest = {
        "purpose": "Local participant-level input for HCAP measurement modelling.",
        "n_age_ge_60": int(len(output)),
        "n_with_any_indicator": int(
            output[standardized_columns].notna().any(axis=1).sum()
        ),
        "n_complete_indicators": int(
            output[standardized_columns].notna().all(axis=1).sum()
        ),
        "indicator_scaling": scaling,
        "output": str(OUTPUT),
        "data_policy": "Participant-level output remains outside version control.",
    }
    MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    MANIFEST.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
