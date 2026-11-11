"""Data schema and 96-well plate layout.

Everything else in the project relies on the names defined here, so this is
the single place to look up what a column or a sample type means.
"""

import pandas as pd

# --- Columns of the input table (one row = one well) -----------------------
COLUMNS = [
    "plate_id",          # text, identifies one physical plate
    "well",              # "A1" ... "H12"
    "sample_type",       # one of SAMPLE_TYPES
    "compound",          # compound name; empty for controls and blanks
    "concentration_uM",  # compound concentration in micromolar; empty for controls
    "replicate",         # 1-based replicate number within its group
    "signal",            # raw luminescence (RLU); empty = measurement missing
]

# --- Sample types ----------------------------------------------------------
# We avoid "positive/negative control" because their meaning depends on the assay.
VEHICLE = "vehicle_control"  # cells without compound (solvent only) -> defines 100 % viability
BLANK = "blank"              # medium + reagent, no cells -> background signal
COMPOUND = "compound"        # cells + compound (test compound or reference)
SAMPLE_TYPES = (VEHICLE, BLANK, COMPOUND)

# --- Plate geometry --------------------------------------------------------
ROWS = "ABCDEFGH"
N_COLS = 12

# --- Demo experiment -------------------------------------------------------
# STS = staurosporine, a well-known cytotoxic compound used here as a reference.
# It is tested like any other compound; it does NOT define 0 % viability.
COMPOUNDS = ("STS", "Compound_A", "Compound_B")
N_REPLICATES = 3
# Highest concentration of each compound (row A). Every following row is a
# 3-fold dilution, giving an 8-point series. STS is potent, so it starts lower.
TOP_CONCENTRATION_UM = {"STS": 10.0, "Compound_A": 100.0, "Compound_B": 100.0}
DILUTION_FACTOR = 3


def concentration_uM(compound: str, row_index: int) -> float:
    """Concentration of `compound` in a given row (0 = row A = highest)."""
    return round(TOP_CONCENTRATION_UM[compound] / DILUTION_FACTOR**row_index, 4)


def well_name(row: str, col: int) -> str:
    """Return a well address such as 'A1' or 'H12'."""
    return f"{row}{col}"


def build_layout() -> pd.DataFrame:
    """Return the 96 wells of the demo plate, without plate_id and signal.

    Layout:
      column 1 and 12  vehicle controls (16), on both sides of the plate
      columns 2-10     3 compounds x 3 replicate columns, rows = concentrations
      column 11        blanks (8)
    """
    wells = []
    vehicle_count = 0

    for row_index, row in enumerate(ROWS):
        for col in range(1, N_COLS + 1):
            compound, concentration = None, None

            if col in (1, 12):
                sample_type = VEHICLE
                vehicle_count += 1
                replicate = vehicle_count
            elif col == 11:
                sample_type = BLANK
                replicate = row_index + 1
            else:
                sample_type = COMPOUND
                compound = COMPOUNDS[(col - 2) // N_REPLICATES]
                replicate = (col - 2) % N_REPLICATES + 1
                concentration = concentration_uM(compound, row_index)

            wells.append(
                {
                    "well": well_name(row, col),
                    "sample_type": sample_type,
                    "compound": compound,
                    "concentration_uM": concentration,
                    "replicate": replicate,
                }
            )

    return pd.DataFrame(wells)
