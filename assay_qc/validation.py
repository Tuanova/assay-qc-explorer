"""Checks that an uploaded table has the expected structure.

Validation answers "can we analyse this file at all?". It reports problems as
readable messages instead of letting the application crash. A missing signal
value is NOT a validation error: the file is still usable and the QC rules
report the well.
"""

import pandas as pd

from assay_qc.schema import BLANK, COLUMNS, COMPOUND, SAMPLE_TYPES, VEHICLE

WELL_PATTERN = r"^[A-H]([1-9]|1[0-2])$"  # A1 ... H12, no leading zero


def validate_data(data: pd.DataFrame) -> tuple[list[str], list[str]]:
    """Return (errors, warnings).

    Errors mean the file cannot be analysed. Warnings mean it can, but some
    plates will not be normalised.
    """
    errors, warnings = [], []

    missing_columns = [column for column in COLUMNS if column not in data.columns]
    if missing_columns:
        return [f"Missing required columns: {', '.join(missing_columns)}."], warnings
    if data.empty:
        return ["The file has no rows."], warnings

    wells = data["well"].astype(str)
    bad_wells = sorted(set(wells[~wells.str.match(WELL_PATTERN)]))
    if bad_wells:
        errors.append(f"Invalid well addresses (expected A1-H12): {', '.join(bad_wells[:10])}.")

    duplicates = data[data.duplicated(["plate_id", "well"], keep=False)]
    if not duplicates.empty:
        examples = sorted({f"{row.plate_id}/{row.well}" for row in duplicates.itertuples()})
        errors.append(f"The same well appears more than once on a plate: {', '.join(examples[:10])}.")

    unknown_types = sorted(set(data["sample_type"].astype(str)) - set(SAMPLE_TYPES))
    if unknown_types:
        errors.append(
            f"Unknown sample_type values: {', '.join(unknown_types)}. "
            f"Allowed: {', '.join(SAMPLE_TYPES)}."
        )

    for column in ("signal", "concentration_uM"):
        as_number = pd.to_numeric(data[column], errors="coerce")
        not_numeric = data[column].notna() & as_number.isna()
        if not_numeric.any():
            errors.append(f"Column '{column}' contains {int(not_numeric.sum())} non-numeric value(s).")

    compounds = data[data["sample_type"] == COMPOUND]
    if compounds["compound"].isna().any() or compounds["concentration_uM"].isna().any():
        errors.append("Some 'compound' wells have no compound name or no concentration.")

    for plate_id, plate in data.groupby("plate_id"):
        for required in (VEHICLE, BLANK):
            if required not in set(plate["sample_type"]):
                warnings.append(f"Plate {plate_id} has no '{required}' wells and cannot be normalised.")

    return errors, warnings


def clean_types(data: pd.DataFrame) -> pd.DataFrame:
    """Return a copy with consistent column types. Call after validation passed."""
    result = data[COLUMNS].copy()
    result["plate_id"] = result["plate_id"].astype(str)
    result["well"] = result["well"].astype(str)
    result["signal"] = pd.to_numeric(result["signal"])
    result["concentration_uM"] = pd.to_numeric(result["concentration_uM"])
    return result
