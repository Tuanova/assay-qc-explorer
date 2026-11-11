"""Calculations on a single plate: control statistics, normalisation, Z'-factor
and replicate summaries.

All functions take the rows of ONE plate. Missing signals (NaN) are ignored in
the statistics; reporting them is the job of the QC rules, not of this module.
"""

import numpy as np
import pandas as pd

from assay_qc.schema import COMPOUND, KILL, VEHICLE


class NormalizationError(ValueError):
    """Raised when a plate cannot be normalised (e.g. a control is missing)."""


def control_stats(plate: pd.DataFrame) -> pd.DataFrame:
    """Return n, mean, SD and CV (%) of the signal for each sample type."""
    stats = plate.groupby("sample_type")["signal"].agg(
        n="count", mean="mean", sd="std"
    )
    stats["cv_pct"] = 100 * stats["sd"] / stats["mean"]
    return stats


def percent_viability(signal, mean_vehicle: float, mean_kill: float):
    """Convert raw signal to percent viability.

        viability % = 100 * (signal - mean_kill) / (mean_vehicle - mean_kill)

    The vehicle control mean maps to 100 %, the kill control mean to 0 %.
    Results are NOT clipped to 0-100 %: values outside that range are
    information about the data and must stay visible.
    """
    return 100 * (signal - mean_kill) / (mean_vehicle - mean_kill)


def normalize_plate(plate: pd.DataFrame, min_window_fraction: float = 0.1) -> pd.DataFrame:
    """Return a copy of the plate with a 'viability_pct' column.

    Raises NormalizationError when a control is missing or when the window
    between the controls (vehicle mean - kill mean) is too small to divide by.
    'Too small' means below `min_window_fraction` of the vehicle mean; this is
    a configurable demo heuristic, not an accepted standard.
    """
    means = plate.groupby("sample_type")["signal"].mean()
    for control in (VEHICLE, KILL):
        if control not in means or np.isnan(means[control]):
            raise NormalizationError(f"No usable '{control}' wells on this plate.")

    window = means[VEHICLE] - means[KILL]
    if window <= min_window_fraction * means[VEHICLE]:
        raise NormalizationError(
            f"Control window too small to normalise: vehicle mean {means[VEHICLE]:.0f}, "
            f"kill mean {means[KILL]:.0f}."
        )

    result = plate.copy()
    result["viability_pct"] = percent_viability(
        result["signal"], means[VEHICLE], means[KILL]
    )
    return result


def z_prime(vehicle_signals, kill_signals) -> float:
    """Z'-factor: how well the two controls are separated on a plate.

        Z' = 1 - 3 * (SD_vehicle + SD_kill) / |mean_vehicle - mean_kill|

    1 is the theoretical maximum; values fall as the controls get noisier or
    closer together. Returns NaN when it cannot be computed (fewer than two
    wells in a control, or identical means). With only a few control wells
    the SDs are uncertain, so Z' is an estimate, not an exact property.
    """
    vehicle = pd.Series(vehicle_signals, dtype=float).dropna()
    kill = pd.Series(kill_signals, dtype=float).dropna()
    if len(vehicle) < 2 or len(kill) < 2:
        return float("nan")

    separation = abs(vehicle.mean() - kill.mean())
    if separation == 0:
        return float("nan")
    return float(1 - 3 * (vehicle.std() + kill.std()) / separation)


def plate_z_prime(plate: pd.DataFrame) -> float:
    """Z'-factor of a plate, taken from its vehicle and kill control wells."""
    return z_prime(
        plate.loc[plate["sample_type"] == VEHICLE, "signal"],
        plate.loc[plate["sample_type"] == KILL, "signal"],
    )


def summarize_replicates(plate: pd.DataFrame, min_mean_signal: float = 0.0) -> pd.DataFrame:
    """Summarise each replicate group (compound x concentration) of a plate.

    Returns one row per group with the number of valid wells, mean, SD and
    CV (%) of the RAW signal. If the plate has a 'viability_pct' column, its
    mean and SD are added too.

    CV is computed on the raw signal, not on percent viability, because
    viability can be close to zero and dividing by a near-zero mean is
    meaningless. CV is left empty (NaN) when a group has fewer than two valid
    wells or when its mean signal is below `min_mean_signal`.
    """
    compounds = plate[plate["sample_type"] == COMPOUND]
    groups = compounds.groupby(["compound", "concentration_uM"])

    summary = groups["signal"].agg(n="count", mean_signal="mean", sd_signal="std")
    summary["cv_pct"] = 100 * summary["sd_signal"] / summary["mean_signal"]
    summary.loc[summary["mean_signal"] < min_mean_signal, "cv_pct"] = np.nan

    if "viability_pct" in plate.columns:
        viability = groups["viability_pct"].agg(mean_viability="mean", sd_viability="std")
        summary = summary.join(viability)

    return summary.reset_index()
