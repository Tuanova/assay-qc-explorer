"""Calculations on a single plate: control statistics, normalisation, Z'-factor
and replicate summaries.

All functions take the rows of ONE plate. Missing signals (NaN) are ignored in
the statistics; reporting them is the job of the QC rules, not of this module.
"""

import numpy as np
import pandas as pd

from assay_qc.schema import BLANK, COMPOUND, VEHICLE


class NormalizationError(ValueError):
    """Raised when a plate cannot be normalised (e.g. the blank is missing)."""


def control_stats(plate: pd.DataFrame) -> pd.DataFrame:
    """Return n, mean, SD and CV (%) of the signal for each sample type."""
    stats = plate.groupby("sample_type")["signal"].agg(
        n="count", mean="mean", sd="std"
    )
    stats["cv_pct"] = 100 * stats["sd"] / stats["mean"]
    return stats


def percent_viability(signal, mean_vehicle: float, mean_blank: float):
    """Convert raw signal to percent viability in two steps.

        net signal  = signal - mean_blank                     (background correction)
        viability % = 100 * net signal / (mean_vehicle - mean_blank)   (normalisation)

    The blank mean is the background of the plate (medium + reagent, no cells).
    The vehicle control mean maps to 100 %; a well with no signal above the
    background maps to 0 %.
    Results are NOT clipped to 0-100 %: values outside that range are
    information about the data and must stay visible.
    """
    net_signal = signal - mean_blank
    return 100 * net_signal / (mean_vehicle - mean_blank)


def normalize_plate(plate: pd.DataFrame, min_window_fraction: float = 0.1) -> pd.DataFrame:
    """Return a copy of the plate with 'net_signal' and 'viability_pct' columns.

    Raises NormalizationError when the vehicle control or the blank is missing,
    or when the window between them (vehicle mean - blank mean) is too small
    to divide by. 'Too small' means below `min_window_fraction` of the vehicle
    mean; this is a configurable demo heuristic, not an accepted standard.
    """
    means = plate.groupby("sample_type")["signal"].mean()
    for required in (VEHICLE, BLANK):
        if required not in means or np.isnan(means[required]):
            raise NormalizationError(f"No usable '{required}' wells on this plate.")

    window = means[VEHICLE] - means[BLANK]
    if window <= min_window_fraction * means[VEHICLE]:
        raise NormalizationError(
            f"Signal window too small to normalise: vehicle mean {means[VEHICLE]:.0f}, "
            f"blank mean {means[BLANK]:.0f}."
        )

    result = plate.copy()
    result["net_signal"] = result["signal"] - means[BLANK]
    result["viability_pct"] = percent_viability(
        result["signal"], means[VEHICLE], means[BLANK]
    )
    return result


def z_prime(vehicle_signals, blank_signals) -> float:
    """Z'-factor: how well the cell signal is separated from the background.

        Z' = 1 - 3 * (SD_vehicle + SD_blank) / |mean_vehicle - mean_blank|

    1 is the theoretical maximum; values fall as the wells get noisier or the
    cell signal gets closer to the background. Here Z' is computed from
    vehicle and blank wells, so it describes the signal-to-background window
    of the plate, not the response to a cytotoxic control.
    Returns NaN when it cannot be computed (fewer than two wells in a group,
    or identical means). With only a few wells the SDs are uncertain, so Z'
    is an estimate, not an exact property.
    """
    vehicle = pd.Series(vehicle_signals, dtype=float).dropna()
    blank = pd.Series(blank_signals, dtype=float).dropna()
    if len(vehicle) < 2 or len(blank) < 2:
        return float("nan")

    separation = abs(vehicle.mean() - blank.mean())
    if separation == 0:
        return float("nan")
    return float(1 - 3 * (vehicle.std() + blank.std()) / separation)


def plate_z_prime(plate: pd.DataFrame) -> float:
    """Z'-factor of a plate, taken from its vehicle and blank wells."""
    return z_prime(
        plate.loc[plate["sample_type"] == VEHICLE, "signal"],
        plate.loc[plate["sample_type"] == BLANK, "signal"],
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
