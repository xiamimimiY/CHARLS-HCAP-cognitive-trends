#!/usr/bin/env python3
"""Build the local longitudinal input for calibrated-cognition inference.

Raw and harmonized CHARLS files are read only. The participant-level output
is written under the version-control-excluded data/derived directory.
"""

from __future__ import annotations

import hashlib
import json
import zipfile
from pathlib import Path

from project_paths import (
    REPOSITORY_ROOT,
    get_data_root,
    get_derived_dir,
    get_results_dir,
    get_temp_dir,
)

import numpy as np
import pandas as pd


DATA_ROOT = get_data_root()
TMP_HARMONIZED = get_temp_dir() / "H_CHARLS_D_Data.dta"
OUT = get_derived_dir() / "b1_frozen_transport_mi_input_2026-06-28.csv"
MANIFEST = get_results_dir() / "manifests" / "longitudinal_input.json"

WAVES = {1: 2011, 2: 2013, 3: 2015, 4: 2018}
COGNITION = ["imrc", "dlrc", "ser7", "orient", "draw"]


def ensure_harmonized() -> Path:
    if TMP_HARMONIZED.exists():
        return TMP_HARMONIZED
    TMP_HARMONIZED.parent.mkdir(parents=True, exist_ok=True)
    harmonized_zip = DATA_ROOT / "Harmonized CHARLS" / "H_CHARLS_D_Data.zip"
    with zipfile.ZipFile(harmonized_zip) as zf:
        zf.extract("H_CHARLS_D_Data.dta", TMP_HARMONIZED.parent)
    return TMP_HARMONIZED


def pseudonym(value: object) -> str:
    return hashlib.sha256(f"B1-frozen-transport::{value}".encode("utf-8")).hexdigest()[:20]


def main() -> None:
    columns = [
        "ID",
        "householdID",
        "communityID",
        "ragender",
        "raeduc_c",
    ]
    for wave in WAVES:
        columns.extend(
            [
                f"r{wave}agey",
                f"r{wave}iwstat",
                f"r{wave}wtrespb",
                f"r{wave}hukou",
                f"h{wave}rural",
                f"r{wave}adlab_c",
                f"r{wave}iadla" if wave >= 2 else "r1adla_c",
                f"hh{wave}cperc",
            ]
        )
        columns.extend(f"r{wave}{stem}" for stem in COGNITION)
    columns = list(dict.fromkeys(columns))

    harmonized = pd.read_stata(
        ensure_harmonized(),
        columns=columns,
        convert_categoricals=False,
    )

    records = []
    counts = []
    for wave, year in WAVES.items():
        iadl_name = f"r{wave}iadla" if wave >= 2 else "r1adla_c"
        frame = pd.DataFrame(
            {
                "ID": harmonized["ID"].astype(str),
                "pid": harmonized["ID"].map(pseudonym),
                "household_id": harmonized["householdID"].astype(str),
                "community_id": harmonized["communityID"].astype(str),
                "wave": wave,
                "year": year,
                "age": harmonized[f"r{wave}agey"],
                "interview_status": harmonized[f"r{wave}iwstat"],
                "female": harmonized["ragender"].eq(2).astype(int),
                "education3": np.select(
                    [
                        harmonized["raeduc_c"].eq(1),
                        harmonized["raeduc_c"].isin([2, 3, 4]),
                        harmonized["raeduc_c"].ge(5)
                        & harmonized["raeduc_c"].le(10),
                    ],
                    [1, 2, 3],
                    default=np.nan,
                ),
                "rural": harmonized[f"h{wave}rural"],
                "hukou": harmonized[f"r{wave}hukou"],
                "weight": harmonized[f"r{wave}wtrespb"],
                "adl": harmonized[f"r{wave}adlab_c"],
                "iadl": harmonized[iadl_name],
                "consumption_per_capita": harmonized[f"hh{wave}cperc"],
            }
        )
        for stem in COGNITION:
            frame[stem] = harmonized[f"r{wave}{stem}"]

        eligible = frame["age"].ge(60) & frame["interview_status"].eq(1)
        frame = frame.loc[eligible].copy()
        frame["log_consumption"] = np.log(
            frame["consumption_per_capita"].where(
                frame["consumption_per_capita"].gt(0)
            )
        )
        records.append(frame)
        counts.append(
            {
                "wave": wave,
                "year": year,
                "n_eligible": int(len(frame)),
                "n_communities": int(frame["community_id"].nunique(dropna=True)),
                "n_positive_weight": int(frame["weight"].gt(0).sum()),
                **{
                    f"{stem}_missing_pct": round(float(frame[stem].isna().mean() * 100), 3)
                    for stem in COGNITION
                },
            }
        )

    output = pd.concat(records, ignore_index=True)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(OUT, index=False)

    manifest = {
        "purpose": "Ignored individual-level input for frozen-anchor longitudinal MI and survey inference.",
        "output": str(OUT.relative_to(REPOSITORY_ROOT)),
        "rows": int(len(output)),
        "participants": int(output["pid"].nunique()),
        "wave_counts": counts,
        "raw_data_policy": "Raw files are read only; data/derived and results are excluded from version control.",
    }
    MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    MANIFEST.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
