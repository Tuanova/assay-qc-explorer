import pandas as pd
import pytest

from assay_qc.evaluation import choose_threshold, confusion, evaluate
from assay_qc.generator import generate_dataset
from assay_qc.rules import QCConfig


def test_confusion_counts_and_ratios():
    units = {"a", "b", "c", "d", "e"}
    result = confusion(truth_units={"a", "b"}, flagged_units={"b", "c"}, units=units)

    assert (result["TP"], result["FP"], result["FN"], result["TN"]) == (1, 1, 1, 2)
    assert result["precision"] == pytest.approx(0.5)
    assert result["recall"] == pytest.approx(0.5)
    assert result["accuracy_if_nothing_flagged"] == pytest.approx(0.6)


def test_clean_plates_have_nothing_to_find():
    data, truth = generate_dataset(n_plates=5, seed=1, clean_fraction=1.0)
    result = evaluate(data, truth, QCConfig())

    assert (result["TP"] == 0).all()
    assert (result["FN"] == 0).all()


def test_choose_threshold_takes_most_recall_within_precision_limit():
    sweep_table = pd.DataFrame({
        "threshold": [10, 15, 20, 25],
        "precision": [0.40, 0.92, 0.95, 1.00],
        "recall": [1.00, 0.90, 0.80, 0.60],
    })
    assert choose_threshold(sweep_table, min_precision=0.90) == 15
