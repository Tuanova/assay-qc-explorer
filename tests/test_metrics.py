import numpy as np
import pandas as pd
import pytest

from assay_qc.metrics import (
    NormalizationError,
    normalize_plate,
    summarize_replicates,
    z_prime,
)
from assay_qc.schema import BLANK, COMPOUND, VEHICLE


def make_plate(vehicle, blank, compound_signals):
    """Tiny hand-made plate: vehicle wells, blanks and one compound group."""
    rows = [{"sample_type": VEHICLE, "signal": s} for s in vehicle]
    rows += [{"sample_type": BLANK, "signal": s} for s in blank]
    rows += [
        {"sample_type": COMPOUND, "compound": "X", "concentration_uM": 1.0, "signal": s}
        for s in compound_signals
    ]
    return pd.DataFrame(rows)


def test_normalization_maps_vehicle_to_100_and_blank_to_0():
    plate = make_plate(vehicle=[900, 1100], blank=[90, 110], compound_signals=[550])
    result = normalize_plate(plate)

    means = result.groupby("sample_type")["viability_pct"].mean()
    assert means[VEHICLE] == pytest.approx(100)
    assert means[BLANK] == pytest.approx(0)
    # 550 is exactly halfway between 100 and 1000.
    assert result["viability_pct"].iloc[-1] == pytest.approx(50)


def test_normalization_does_not_clip_values_outside_0_100():
    plate = make_plate(vehicle=[1000, 1000], blank=[100, 100], compound_signals=[1450, 10])
    viability = normalize_plate(plate)["viability_pct"].tolist()

    assert viability[-2] == pytest.approx(150)
    assert viability[-1] == pytest.approx(-10)


def test_normalization_fails_clearly_without_blank():
    plate = make_plate(vehicle=[1000, 1000], blank=[], compound_signals=[500])
    with pytest.raises(NormalizationError, match="blank"):
        normalize_plate(plate)


def test_normalization_fails_when_signal_is_not_above_background():
    plate = make_plate(vehicle=[1000, 1000], blank=[990, 995], compound_signals=[500])
    with pytest.raises(NormalizationError, match="window"):
        normalize_plate(plate)


def test_z_prime_matches_hand_calculation():
    # SD_vehicle = 10, SD_blank = 2, means 100 and 10 -> 1 - 3 * 12 / 90 = 0.6
    assert z_prime([90, 100, 110], [8, 10, 12]) == pytest.approx(0.6)


def test_z_prime_is_nan_with_a_single_vehicle_well():
    assert np.isnan(z_prime([100], [8, 10, 12]))


def test_replicate_cv_matches_hand_calculation():
    plate = make_plate(vehicle=[1000, 1000], blank=[100, 100], compound_signals=[90, 100, 110])
    group = summarize_replicates(plate).iloc[0]

    assert group["n"] == 3
    assert group["mean_signal"] == pytest.approx(100)
    assert group["cv_pct"] == pytest.approx(10)  # SD 10 / mean 100


def test_replicate_cv_is_empty_when_only_one_well_is_left():
    plate = make_plate(vehicle=[1000, 1000], blank=[100, 100], compound_signals=[100, np.nan, np.nan])
    group = summarize_replicates(plate).iloc[0]

    assert group["n"] == 1
    assert np.isnan(group["cv_pct"])


def test_replicate_cv_is_empty_below_signal_floor():
    plate = make_plate(vehicle=[1000, 1000], blank=[100, 100], compound_signals=[90, 100, 110])
    group = summarize_replicates(plate, min_mean_signal=200).iloc[0]

    assert np.isnan(group["cv_pct"])
