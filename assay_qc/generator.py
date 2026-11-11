"""Synthetic data generator for a luminescence cell-viability assay.

The generator returns two tables:
  * data  - what the application sees (columns from schema.COLUMNS)
  * truth - which problems were injected on purpose (never shown to the rules)

All numbers here are invented to look plausible. They are not taken from any
real experiment.
"""

from dataclasses import dataclass

import numpy as np
import pandas as pd

from assay_qc.schema import BLANK, COLUMNS, COMPOUND, VEHICLE, build_layout

TRUTH_COLUMNS = [
    "plate_id", "level", "issue_type", "well", "compound", "concentration_uM", "detail",
]


@dataclass(frozen=True)
class CompoundProfile:
    """Dose-response shape of a simulated compound (Hill curve)."""

    ic50_uM: float  # concentration giving a half-maximal effect
    hill: float     # steepness of the curve
    bottom: float   # viability fraction left at very high concentration


PROFILES = {
    "STS": CompoundProfile(ic50_uM=0.05, hill=1.0, bottom=0.03),           # reference, potent
    "Compound_A": CompoundProfile(ic50_uM=5.0, hill=1.2, bottom=0.0),      # active
    "Compound_B": CompoundProfile(ic50_uM=10_000.0, hill=1.0, bottom=0.0), # inactive
}

BACKGROUND_NOISE_RLU = 50  # additive reader noise, same for every well


def viability_fraction(concentration_uM, profile: CompoundProfile):
    """Expected viability (1 = like vehicle control, 0 = no living cells)."""
    effect = 1 / (1 + (concentration_uM / profile.ic50_uM) ** profile.hill)
    return profile.bottom + (1 - profile.bottom) * effect


def generate_plate(plate_id: str, rng: np.random.Generator, poor_separation: str | None = None) -> pd.DataFrame:
    """Simulate one plate without well-level problems.

    poor_separation: None, "noisy_vehicle" (vehicle wells scatter a lot) or
    "high_background" (blank signal is high and variable, e.g. contaminated
    medium), so the window between cells and background is poor.
    """
    plate = build_layout()

    # Plate-to-plate differences: overall brightness, background, noise level.
    vehicle_level = 100_000 * rng.lognormal(mean=0, sigma=0.10)
    blank_level = rng.uniform(300, 600)
    noise_cv = rng.uniform(0.04, 0.09)
    vehicle_cv, blank_cv = noise_cv, 0.05

    if poor_separation == "noisy_vehicle":
        vehicle_cv = rng.uniform(0.22, 0.32)
    elif poor_separation == "high_background":
        blank_level = vehicle_level * rng.uniform(0.40, 0.60)
        blank_cv = 0.15
    elif poor_separation is not None:
        raise ValueError(f"Unknown poor_separation mode: {poor_separation}")

    # Expected (noise-free) signal of every well.
    expected = np.zeros(len(plate))
    cv = np.full(len(plate), noise_cv)
    for i, well in plate.iterrows():
        if well["sample_type"] == VEHICLE:
            expected[i], cv[i] = vehicle_level, vehicle_cv
        elif well["sample_type"] == BLANK:
            expected[i], cv[i] = blank_level, blank_cv
        else:
            fraction = viability_fraction(well["concentration_uM"], PROFILES[well["compound"]])
            expected[i] = blank_level + fraction * (vehicle_level - blank_level)

    # Noise: a part proportional to the signal plus a constant background part.
    proportional = expected * cv * rng.standard_normal(len(plate))
    background = BACKGROUND_NOISE_RLU * rng.standard_normal(len(plate))
    signal = np.clip(expected + proportional + background, 0, None).round()

    plate.insert(0, "plate_id", plate_id)
    plate["signal"] = signal
    return plate[COLUMNS]


def _perturbation_factor(rng: np.random.Generator) -> float:
    """Factor applied to one well: clearly too low or clearly too high."""
    if rng.random() < 0.5:
        return float(rng.uniform(0.45, 0.80))
    return float(rng.uniform(1.25, 1.90))


def inject_well_issues(
    plate: pd.DataFrame,
    rng: np.random.Generator,
    n_perturbed_groups: int = 0,
    n_missing: int = 0,
    control_outlier: bool = False,
) -> list[dict]:
    """Modify `plate` in place and return the ground-truth rows.

    A replicate group gets at most one problem, so the truth stays unambiguous.
    """
    plate_id = plate["plate_id"].iloc[0]
    truth = []

    compound_wells = plate[plate["sample_type"] == COMPOUND]
    groups = list(compound_wells.groupby(["compound", "concentration_uM"]).groups.items())
    chosen = rng.choice(len(groups), size=n_perturbed_groups + n_missing, replace=False)

    for position, group_number in enumerate(chosen):
        (compound, concentration), well_indices = groups[group_number]
        index = rng.choice(well_indices)
        row = {
            "plate_id": plate_id,
            "well": plate.at[index, "well"],
            "compound": compound,
            "concentration_uM": concentration,
        }
        if position < n_perturbed_groups:
            factor = _perturbation_factor(rng)
            plate.at[index, "signal"] = round(plate.at[index, "signal"] * factor)
            truth.append({**row, "level": "group", "issue_type": "replicate_perturbation",
                          "detail": f"signal x {factor:.2f}"})
        else:
            plate.at[index, "signal"] = np.nan
            truth.append({**row, "level": "well", "issue_type": "missing_signal", "detail": ""})

    if control_outlier:
        control = str(rng.choice([VEHICLE, BLANK]))
        index = rng.choice(plate.index[plate["sample_type"] == control])
        # A vehicle well that lost signal, or a blank well with extra signal
        # (for example a splash from a neighbouring well).
        factor = float(rng.uniform(0.40, 0.70) if control == VEHICLE else rng.uniform(3.0, 6.0))
        plate.at[index, "signal"] = round(plate.at[index, "signal"] * factor)
        truth.append({"plate_id": plate_id, "level": "well", "issue_type": "control_outlier",
                      "well": plate.at[index, "well"], "compound": None, "concentration_uM": None,
                      "detail": f"{control} signal x {factor:.2f}"})

    return truth


def generate_dataset(
    n_plates: int = 10,
    seed: int = 0,
    clean_fraction: float = 0.3,
    poor_separation_prob: float = 0.15,
    plate_prefix: str = "P",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Generate several independent plates. Returns (data, truth).

    The same seed always gives the same dataset. About `clean_fraction` of the
    plates get no injected problem at all.
    """
    rng = np.random.default_rng(seed)
    plates, truth = [], []

    for number in range(1, n_plates + 1):
        plate_id = f"{plate_prefix}{number:03d}"
        is_clean = rng.random() < clean_fraction

        poor_separation = None
        if not is_clean and rng.random() < poor_separation_prob:
            poor_separation = str(rng.choice(["noisy_vehicle", "high_background"]))

        plate = generate_plate(plate_id, rng, poor_separation)

        if poor_separation is not None:
            truth.append({"plate_id": plate_id, "level": "plate",
                          "issue_type": "poor_control_separation", "well": None,
                          "compound": None, "concentration_uM": None, "detail": poor_separation})
        if not is_clean:
            truth += inject_well_issues(
                plate,
                rng,
                n_perturbed_groups=int(rng.integers(1, 4)),
                n_missing=int(rng.integers(0, 3)),
                # Keep plate-level and control-well problems on separate plates.
                control_outlier=poor_separation is None and rng.random() < 0.5,
            )
        plates.append(plate)

    data = pd.concat(plates, ignore_index=True)
    return data, pd.DataFrame(truth, columns=TRUTH_COLUMNS)
