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
VEHICLE = "vehicle_control"  # cells + solvent only -> defines 100 % viability
KILL = "kill_control"        # cells + cytotoxic reference -> defines 0 % viability
BLANK = "blank"              # medium + reagent, no cells -> background, QC display only
COMPOUND = "compound"        # cells + test compound
SAMPLE_TYPES = (VEHICLE, KILL, BLANK, COMPOUND)

# --- Plate geometry --------------------------------------------------------
ROWS = "ABCDEFGH"
N_COLS = 12

# --- Demo experiment -------------------------------------------------------
COMPOUNDS = ("Compound_A", "Compound_B", "Compound_C")
N_REPLICATES = 3
# 8-point, 3-fold dilution series starting at 100 uM (row A = highest).
CONCENTRATIONS_UM = tuple(round(100 / 3**i, 4) for i in range(len(ROWS)))


def well_name(row: str, col: int) -> str:
    """Return a well address such as 'A1' or 'H12'."""
    return f"{row}{col}"


def build_layout() -> pd.DataFrame:
    """Return the 96 wells of the demo plate, without plate_id and signal.

    Layout:
      column 1 and 12  vehicle and kill controls, alternating, so that both
                       control types sit on both sides of the plate (8 + 8)
      columns 2-10     3 compounds x 3 replicate columns, rows = concentrations
      column 11        blanks (8)
    """
    wells = []
    control_count = {VEHICLE: 0, KILL: 0}

    for row_index, row in enumerate(ROWS):
        for col in range(1, N_COLS + 1):
            compound, concentration = None, None

            if col in (1, 12):
                # Column 1 starts with vehicle, column 12 starts with kill.
                even_row = row_index % 2 == 0
                starts_with_vehicle = col == 1
                sample_type = VEHICLE if even_row == starts_with_vehicle else KILL
                control_count[sample_type] += 1
                replicate = control_count[sample_type]
            elif col == 11:
                sample_type = BLANK
                replicate = row_index + 1
            else:
                sample_type = COMPOUND
                compound = COMPOUNDS[(col - 2) // N_REPLICATES]
                replicate = (col - 2) % N_REPLICATES + 1
                concentration = CONCENTRATIONS_UM[row_index]

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
