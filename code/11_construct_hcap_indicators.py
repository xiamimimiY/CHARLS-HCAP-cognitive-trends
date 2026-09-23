#!/usr/bin/env python3
"""Construct the 12 scored CHARLS-HCAP indicators used in the study."""

from __future__ import annotations

import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

from project_paths import get_data_root, get_temp_dir


DATA_ROOT = get_data_root()
TMP_HARMONIZED = get_temp_dir() / "H_CHARLS_D_Data.dta"

COGNITION_FILE = DATA_ROOT / "2018 CHARLS Wave 4" / "Cognition.dta"
INFORMANT_FILE = DATA_ROOT / "2018 CHARLS Wave 4" / "Insider.dta"
HEALTH_FILE = (
    DATA_ROOT
    / "2018 CHARLS Wave 4"
    / "Health_Status_and_Functioning.dta"
)
HARMONIZED_ZIP = DATA_ROOT / "Harmonized CHARLS" / "H_CHARLS_D_Data.zip"

TARGET_WORDS = {
    "dc050_w4",
    "dc052_w4",
    "dc053_w4",
    "dc055_w4",
    "dc058_w4",
    "dc059_w4",
    "dc061_w4",
    "dc064_w4",
    "dc066_w4",
    "dc067_w4",
}


def ensure_harmonized() -> Path:
    """Extract the harmonized file to a temporary location when needed."""
    if TMP_HARMONIZED.exists():
        return TMP_HARMONIZED
    TMP_HARMONIZED.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(HARMONIZED_ZIP) as archive:
        archive.extract("H_CHARLS_D_Data.dta", TMP_HARMONIZED.parent)
    return TMP_HARMONIZED


def count_correct_one(
    frame: pd.DataFrame,
    columns: list[str],
    minimum: int = 1,
) -> pd.Series:
    """Count responses coded correct while enforcing minimum availability."""
    values = frame[columns].copy()
    valid = values.isin([1, 5])
    score = values.eq(1).sum(axis=1).astype(float)
    score[valid.sum(axis=1).lt(minimum)] = np.nan
    return score


def score_multiselect(
    frame: pd.DataFrame,
    correct_columns: list[str],
    administration_columns: list[str],
    invalid_columns: list[str] | None = None,
) -> pd.Series:
    """Score expanded multiple-response fields without treating blanks as errors."""
    values = frame[correct_columns]
    administered = frame[administration_columns].notna().any(axis=1)
    score = values.fillna(0).ne(0).sum(axis=1).astype(float)
    invalid = pd.Series(False, index=frame.index)
    for column in invalid_columns or []:
        invalid |= frame[column].fillna(0).ne(0)
    score[~administered | invalid] = np.nan
    return score


def score_serial7(frame: pd.DataFrame) -> pd.Series:
    """Score the five serial-seven subtraction responses."""
    targets = [93, 86, 79, 72, 65]
    result_columns = [f"dc014_w4_{index}_1" for index in range(1, 6)]
    reason_columns = [f"dc014_w4_{index}" for index in range(1, 6)]
    score = pd.Series(0.0, index=frame.index)
    observed = pd.Series(0, index=frame.index)
    for result, reason, target in zip(
        result_columns,
        reason_columns,
        targets,
    ):
        available = frame[result].notna() | frame[reason].notna()
        observed += available.astype(int)
        score += frame[result].eq(target).astype(float)
    score[observed.eq(0)] = np.nan
    return score


def score_recognition(frame: pd.DataFrame) -> pd.Series:
    """Score the 20 recognition fields using target and distractor coding."""
    columns = [f"dc{index:03d}_w4" for index in range(48, 68)]
    valid = frame[columns].isin([1, 5])
    correct = pd.DataFrame(False, index=frame.index, columns=columns)
    for column in columns:
        correct[column] = (
            frame[column].eq(1)
            if column in TARGET_WORDS
            else frame[column].eq(5)
        )
    score = correct.sum(axis=1).astype(float)
    score[valid.sum(axis=1).lt(10)] = np.nan
    return score


def score_informant(frame: pd.DataFrame) -> pd.DataFrame:
    """Construct the informant measures used in missingness sensitivities."""
    output = frame.copy()
    iq_columns = [f"dd{index:03d}_w4" for index in range(12, 38)]
    iq_values = output[iq_columns].where(
        output[iq_columns].isin([1, 2, 3, 4, 5])
    )
    output["iqcode_mean"] = iq_values.mean(axis=1)
    output.loc[iq_values.notna().sum(axis=1).lt(13), "iqcode_mean"] = np.nan

    blessed_columns = [f"dd{index:03d}_w4" for index in range(42, 48)]
    blessed_values = output[blessed_columns].where(
        output[blessed_columns].isin([1, 2])
    )
    output["blessed_score"] = blessed_values.eq(2).sum(axis=1).astype(float)
    output.loc[
        blessed_values.notna().sum(axis=1).lt(4),
        "blessed_score",
    ] = np.nan
    return output


def build_scores() -> pd.DataFrame:
    """Read authorized inputs and construct the finalized HCAP indicators."""
    cognition = pd.read_stata(COGNITION_FILE, convert_categoricals=False)
    informant = score_informant(
        pd.read_stata(INFORMANT_FILE, convert_categoricals=False)
    )
    harmonized = pd.read_stata(
        ensure_harmonized(),
        columns=[
            "ID",
            "r4agey",
            "ragender",
            "raeduc_c",
            "h4rural",
            "r4wtrespb",
            "r4adlab_c",
            "r4iadla",
            "r4imrc",
            "r4dlrc",
            "r4ser7",
            "r4orient",
            "r4draw",
        ],
        convert_categoricals=False,
    )
    health = pd.read_stata(
        HEALTH_FILE,
        columns=["ID", "da005_3_", "da005_4_", "da006_w4_3_", "da006_w4_4_"],
        convert_categoricals=False,
    )
    frame = (
        cognition.merge(informant, on="ID", how="outer", suffixes=("", "_inf"))
        .merge(harmonized, on="ID", how="left")
        .merge(health, on="ID", how="left")
    )

    orientation_columns = [
        "dc001_w4",
        "dc002_w4",
        "dc003_w4",
        "dc005_w4",
        "dc006_w4",
        "dc007_w4",
        "dc008_w4",
        "dc009_w4",
        "dc010_w4",
        "dc012_w4",
    ]
    frame["test_orientation"] = count_correct_one(
        frame,
        orientation_columns,
        minimum=5,
    )

    for repetition in range(1, 6):
        correct = [f"dc013_w4_{repetition}_s{index}" for index in [1, 2, 3]]
        administration = correct + [
            f"dc013_w4_{repetition}_s4",
            f"dc013_w4_{repetition}_s97",
        ]
        frame[f"test_short_learning_trial{repetition}"] = score_multiselect(
            frame,
            correct_columns=correct,
            administration_columns=administration,
            invalid_columns=[f"dc013_w4_{repetition}_s97"],
        )

    short_delayed = [f"dc015_w4_s{index}" for index in [1, 2, 3]]
    frame["test_short_delayed"] = score_multiselect(
        frame,
        correct_columns=short_delayed,
        administration_columns=short_delayed + ["dc015_w4_s4", "dc015_w4_s97"],
        invalid_columns=["dc015_w4_s97"],
    )

    for question in [28, 29, 30]:
        correct = [f"dc{question:03d}_w4_s{index}" for index in range(1, 11)]
        administration = correct + [
            f"dc{question:03d}_w4_s11",
            f"dc{question:03d}_w4_s12",
        ]
        frame[f"test_ten_learning_trial{question - 27}"] = score_multiselect(
            frame,
            correct_columns=correct,
            administration_columns=administration,
            invalid_columns=[f"dc{question:03d}_w4_s12"],
        )

    ten_delayed = [f"dc047_w4_s{index}" for index in range(1, 11)]
    frame["test_ten_delayed"] = score_multiselect(
        frame,
        correct_columns=ten_delayed,
        administration_columns=ten_delayed + ["dc047_w4_s11", "dc047_w4_s12"],
        invalid_columns=["dc047_w4_s12"],
    )
    frame["test_recognition"] = score_recognition(frame)
    frame["test_serial7"] = score_serial7(frame)

    animal_columns = ["dc033_w4", "dc035_w4", "dc037_w4", "dc039_w4"]
    animal_values = frame[animal_columns].mask(frame[animal_columns].gt(100))
    frame["test_animal_fluency"] = animal_values.sum(axis=1, min_count=2)

    language_columns = [
        "dc016_w4",
        "dc017_w4",
        "dc018_w4",
        "dc019_w4",
        "dc023_w4",
        "dc025_w4",
        "dc026_w4_revised",
        "dc027_w4",
        "dc042_w4",
        "dc043_w4",
    ]
    frame["test_language"] = count_correct_one(
        frame,
        language_columns,
        minimum=5,
    )

    visuospatial_columns = ["dc020_w4", "dc021_w4", "dc022_w4", "dc024_w4"]
    frame["test_visuospatial"] = count_correct_one(
        frame,
        visuospatial_columns,
        minimum=2,
    )

    frame["age60"] = frame["r4agey"].ge(60)
    frame["female"] = frame["ragender"].eq(2).astype(float)
    frame["edu3"] = np.select(
        [
            frame["raeduc_c"].eq(1),
            frame["raeduc_c"].isin([2, 3, 4]),
            frame["raeduc_c"].ge(5) & frame["raeduc_c"].le(10),
        ],
        [1, 2, 3],
        default=np.nan,
    )
    frame["rural"] = frame["h4rural"].where(frame["h4rural"].isin([0, 1]))
    frame["vision_problem"] = (
        frame[["da005_3_", "da006_w4_3_"]].eq(1).any(axis=1).astype(float)
    )
    frame["hearing_problem"] = (
        frame[["da005_4_", "da006_w4_4_"]].eq(1).any(axis=1).astype(float)
    )
    frame["sensory_impairment_count"] = (
        frame["vision_problem"] + frame["hearing_problem"]
    )
    frame["dual_sensory_impairment"] = (
        frame["sensory_impairment_count"].eq(2).astype(float)
    )
    return frame
