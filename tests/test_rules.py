import numpy as np
import pandas as pd

from assay_qc.rules import FINDING_COLUMNS, QCConfig, run_qc
from assay_qc.schema import BLANK, COMPOUND, VEHICLE

GOOD_VEHICLE = [980, 1000, 1020, 990, 1010, 1000, 995, 1005]
GOOD_BLANK = [98, 100, 102, 99, 101, 100, 97, 103]


def make_plate(vehicle=GOOD_VEHICLE, blank=GOOD_BLANK, compound_signals=(500, 505, 495)):
    """Hand-made plate: 8 vehicle wells, 8 blanks and one compound group."""
    rows = [{"sample_type": VEHICLE, "signal": s} for s in vehicle]
    rows += [{"sample_type": BLANK, "signal": s} for s in blank]
    rows += [
        {"sample_type": COMPOUND, "compound": "X", "concentration_uM": 1.0, "signal": s}
        for s in compound_signals
    ]
    plate = pd.DataFrame(rows)
    plate.insert(0, "plate_id", "T1")
    plate.insert(1, "well", [f"A{i + 1}" for i in range(len(plate))])
    return plate


def rules_found(plate, config=QCConfig()):
    return run_qc(plate, config)["rule"].tolist()


def test_good_plate_has_no_findings():
    findings = run_qc(make_plate())
    assert findings.empty
    assert list(findings.columns) == FINDING_COLUMNS


def test_missing_signal_is_reported_with_its_well():
    plate = make_plate(compound_signals=(500, np.nan, 495))
    findings = run_qc(plate)

    assert findings["rule"].tolist() == ["missing_signal"]
    assert findings["well"].iloc[0] == plate.loc[plate["signal"].isna(), "well"].iloc[0]


def test_vehicle_outlier_is_flagged():
    vehicle = GOOD_VEHICLE[:-1] + [600]
    assert "control_outlier" in rules_found(make_plate(vehicle=vehicle))


def test_blank_outlier_is_flagged():
    blank = GOOD_BLANK[:-1] + [400]
    assert "control_outlier" in rules_found(make_plate(blank=blank))


def test_replicate_cv_depends_on_threshold():
    plate = make_plate(compound_signals=(450, 500, 550))  # CV = 10 %

    assert "replicate_cv" in rules_found(plate, QCConfig(replicate_cv_pct=5))
    assert "replicate_cv" not in rules_found(plate, QCConfig(replicate_cv_pct=20))


def test_replicate_cv_is_not_judged_near_background():
    plate = make_plate(compound_signals=(90, 120, 150))  # CV 25 %, but mean 120 < 3 x blank 100
    assert "replicate_cv" not in rules_found(plate)


def test_low_z_prime_is_flagged():
    noisy_vehicle = [600, 1400, 800, 1200, 700, 1300, 900, 1100]
    assert "z_prime" in rules_found(make_plate(vehicle=noisy_vehicle))


def test_rules_do_not_change_the_data():
    plate = make_plate(vehicle=GOOD_VEHICLE[:-1] + [600], compound_signals=(500, np.nan, 900))
    before = plate.copy()
    run_qc(plate)
    pd.testing.assert_frame_equal(plate, before)
