"""Quality-control rules. Each rule looks at one plate and returns findings.

A finding is a pointer for a human reviewer: "look at this well / group /
plate, here is why". Rules never delete or change data.

All thresholds live in QCConfig. They are demonstration heuristics that can
be changed, not accepted industry standards.
"""

from dataclasses import dataclass

import numpy as np
import pandas as pd

from assay_qc.metrics import plate_z_prime, summarize_replicates
from assay_qc.schema import BLANK, VEHICLE

FINDING_COLUMNS = [
    "finding_id", "plate_id", "level", "well", "compound", "concentration_uM",
    "rule", "observed", "threshold", "explanation", "suggested_check",
]


@dataclass(frozen=True)
class QCConfig:
    # The first two defaults are the medians of the thresholds chosen on tuning
    # plates in scripts/run_evaluation.py. Z' keeps a conventional 0.5.
    control_outlier_z: float = 7.0    # robust z-score above which a vehicle/blank well is flagged
    replicate_cv_pct: float = 17.0    # replicate CV (%) above which a group is flagged
    cv_floor_multiple: float = 3.0    # CV is only judged if group mean >= this x blank mean
    z_prime_min: float = 0.5          # plates with Z' below this are flagged


def _finding(plate_id, level, rule, observed, threshold, explanation, suggested_check,
             well=None, compound=None, concentration_uM=None) -> dict:
    return {
        "plate_id": plate_id, "level": level, "well": well, "compound": compound,
        "concentration_uM": concentration_uM, "rule": rule, "observed": observed,
        "threshold": threshold, "explanation": explanation, "suggested_check": suggested_check,
    }


def check_missing_signal(plate: pd.DataFrame) -> list[dict]:
    """R1: wells without a measurement. A structural problem, no statistics."""
    plate_id = plate["plate_id"].iloc[0]
    return [
        _finding(
            plate_id, "well", "missing_signal", observed="missing", threshold="-",
            explanation=f"Well {row.well} ({row.sample_type}) has no signal value.",
            suggested_check="Check the reader export and whether the well was dispensed.",
            well=row.well,
        )
        for row in plate[plate["signal"].isna()].itertuples()
    ]


def check_control_outliers(plate: pd.DataFrame, config: QCConfig) -> list[dict]:
    """R2: a vehicle or blank well far from the other wells of the same type.

    Uses the median and the MAD (median absolute deviation) instead of mean
    and SD, because one extreme well shifts the mean and inflates the SD,
    which would hide the very outlier we are looking for.

        robust z = |signal - median| / (1.4826 * MAD)

    The factor 1.4826 makes the MAD comparable to an SD for normal data.
    """
    plate_id = plate["plate_id"].iloc[0]
    findings = []

    for control in (VEHICLE, BLANK):
        wells = plate[(plate["sample_type"] == control) & plate["signal"].notna()]
        if len(wells) < 4:
            continue  # too few wells to say what is typical
        median = wells["signal"].median()
        mad = (wells["signal"] - median).abs().median()
        if mad == 0:
            continue
        robust_z = (wells["signal"] - median).abs() / (1.4826 * mad)

        for index in wells.index[robust_z > config.control_outlier_z]:
            findings.append(_finding(
                plate_id, "well", "control_outlier",
                observed=round(float(robust_z[index]), 1), threshold=config.control_outlier_z,
                explanation=(
                    f"{control} well {wells.at[index, 'well']} has signal "
                    f"{wells.at[index, 'signal']:.0f}; the median of the {control} wells is {median:.0f}."
                ),
                suggested_check=(
                    "Inspect the well (dispensing, bubbles, splash from a neighbour). "
                    "Normalisation and Z' of this plate use this well."
                ),
                well=wells.at[index, "well"],
            ))
    return findings


def check_replicate_cv(plate: pd.DataFrame, config: QCConfig) -> list[dict]:
    """R3: replicate groups whose raw-signal CV is above the threshold.

    Groups whose mean signal is close to the blank (background) level are
    not judged, because CV becomes unstable near the background.
    """
    plate_id = plate["plate_id"].iloc[0]
    blank_mean = plate.loc[plate["sample_type"] == BLANK, "signal"].mean()
    floor = 0.0 if np.isnan(blank_mean) else config.cv_floor_multiple * blank_mean

    summary = summarize_replicates(plate, min_mean_signal=floor)
    flagged = summary[summary["cv_pct"] > config.replicate_cv_pct]

    return [
        _finding(
            plate_id, "group", "replicate_cv",
            observed=round(float(row.cv_pct), 1), threshold=config.replicate_cv_pct,
            explanation=(
                f"Replicates of {row.compound} at {row.concentration_uM:g} uM disagree: "
                f"CV {row.cv_pct:.1f} % over {row.n} wells."
            ),
            suggested_check=(
                "Compare the replicate wells. High CV can be a pipetting problem, but also "
                "appears on the steep part of a dose-response curve."
            ),
            compound=row.compound, concentration_uM=row.concentration_uM,
        )
        for row in flagged.itertuples()
    ]


def check_z_prime(plate: pd.DataFrame, config: QCConfig) -> list[dict]:
    """R4: plates where the cell signal is poorly separated from the background (low Z')."""
    plate_id = plate["plate_id"].iloc[0]
    z = plate_z_prime(plate)

    if np.isnan(z):
        return [_finding(
            plate_id, "plate", "z_prime", observed="not computable", threshold=config.z_prime_min,
            explanation="Z' could not be computed: fewer than two usable vehicle or blank wells.",
            suggested_check="Check that vehicle controls and blanks are present and labelled.",
        )]
    if z < config.z_prime_min:
        return [_finding(
            plate_id, "plate", "z_prime", observed=round(z, 2), threshold=config.z_prime_min,
            explanation=f"Vehicle signal is poorly separated from the blank: Z' = {z:.2f}.",
            suggested_check=(
                "Look at the vehicle and blank wells: are the vehicle wells noisy, is the "
                "background high, is a single outlier well responsible?"
            ),
        )]
    return []


def run_qc(data: pd.DataFrame, config: QCConfig = QCConfig()) -> pd.DataFrame:
    """Run all rules on every plate and return one table of findings."""
    findings = []
    for _, plate in data.groupby("plate_id", sort=True):
        findings += check_missing_signal(plate)
        findings += check_control_outliers(plate, config)
        findings += check_replicate_cv(plate, config)
        findings += check_z_prime(plate, config)

    table = pd.DataFrame(findings, columns=FINDING_COLUMNS[1:])
    table.insert(0, "finding_id", [f"F{n:03d}" for n in range(1, len(table) + 1)])
    return table
