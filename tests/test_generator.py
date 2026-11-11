import pandas as pd

from assay_qc.generator import generate_dataset
from assay_qc.schema import BLANK, COLUMNS, COMPOUND, KILL, VEHICLE


def test_every_plate_has_96_unique_wells_and_expected_layout():
    data, _ = generate_dataset(n_plates=3, seed=1)

    assert list(data.columns) == COLUMNS
    for _, plate in data.groupby("plate_id"):
        assert plate["well"].nunique() == 96
        counts = plate["sample_type"].value_counts()
        assert counts[VEHICLE] == 8
        assert counts[KILL] == 8
        assert counts[BLANK] == 8
        assert counts[COMPOUND] == 72


def test_same_seed_gives_same_dataset():
    data_1, truth_1 = generate_dataset(n_plates=4, seed=7)
    data_2, truth_2 = generate_dataset(n_plates=4, seed=7)

    pd.testing.assert_frame_equal(data_1, data_2)
    pd.testing.assert_frame_equal(truth_1, truth_2)


def test_clean_dataset_has_no_truth_rows_and_no_missing_signal():
    data, truth = generate_dataset(n_plates=3, seed=1, clean_fraction=1.0)

    assert truth.empty
    assert data["signal"].notna().all()


def test_injected_missing_signals_are_recorded_in_truth():
    data, truth = generate_dataset(n_plates=10, seed=3, clean_fraction=0.0)

    missing_in_data = set(
        map(tuple, data.loc[data["signal"].isna(), ["plate_id", "well"]].values)
    )
    missing_in_truth = set(
        map(tuple, truth.loc[truth["issue_type"] == "missing_signal", ["plate_id", "well"]].values)
    )
    assert missing_in_data == missing_in_truth
