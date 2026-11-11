"""Compare findings with the generator's ground truth.

Each rule is evaluated at its own level, with its own unit of counting:

    replicate_cv     unit = replicate group (plate x compound x concentration)
    control_outlier  unit = vehicle or blank well
    z_prime          unit = plate

The three levels are never mixed in one table. `missing_signal` is not
evaluated here: an empty value is found by definition, which is covered by
unit tests instead.
"""

import pandas as pd

from assay_qc.rules import QCConfig, run_qc
from assay_qc.schema import BLANK, COMPOUND, VEHICLE

# rule name -> (issue type in the truth file, columns that identify one unit)
LEVELS = {
    "replicate_cv": ("replicate_perturbation", ["plate_id", "compound", "concentration_uM"]),
    "control_outlier": ("control_outlier", ["plate_id", "well"]),
    "z_prime": ("poor_control_separation", ["plate_id"]),
}


def _keys(table: pd.DataFrame, columns: list[str]) -> set[tuple]:
    """Turn rows into a set of tuples so that truth and findings can be compared."""
    table = table[columns].copy()
    if "concentration_uM" in columns:
        table["concentration_uM"] = table["concentration_uM"].astype(float).round(4)
    return set(map(tuple, table.to_numpy()))


def all_units(data: pd.DataFrame, rule: str) -> set[tuple]:
    """Every unit the rule could flag (needed to count true negatives)."""
    columns = LEVELS[rule][1]
    if rule == "replicate_cv":
        data = data[data["sample_type"] == COMPOUND]
    elif rule == "control_outlier":
        data = data[data["sample_type"].isin([VEHICLE, BLANK])]
    return _keys(data.drop_duplicates(columns), columns)


def confusion(truth_units: set, flagged_units: set, units: set) -> dict:
    """Count TP, FP, FN, TN and derive precision, recall and accuracy.

    precision = TP / (TP + FP)   of everything flagged, how much was a real problem
    recall    = TP / (TP + FN)   of all real problems, how many were flagged
    """
    tp = len(truth_units & flagged_units)
    fp = len(flagged_units - truth_units)
    fn = len(truth_units - flagged_units)
    tn = len(units) - tp - fp - fn
    return {
        "units": len(units), "TP": tp, "FP": fp, "FN": fn, "TN": tn,
        "precision": tp / (tp + fp) if tp + fp else float("nan"),
        "recall": tp / (tp + fn) if tp + fn else float("nan"),
        "accuracy": (tp + tn) / len(units),
        # What you would score by never flagging anything. With rare problems
        # this is already high, which is why accuracy alone says little here.
        "accuracy_if_nothing_flagged": (len(units) - len(truth_units)) / len(units),
    }


def evaluate(data: pd.DataFrame, truth: pd.DataFrame, config: QCConfig) -> pd.DataFrame:
    """Run the rules and return one row of counts per rule."""
    findings = run_qc(data, config)
    rows = []
    for rule, (issue_type, columns) in LEVELS.items():
        truth_units = _keys(truth[truth["issue_type"] == issue_type], columns)
        flagged_units = _keys(findings[findings["rule"] == rule], columns)
        rows.append({"rule": rule, **confusion(truth_units, flagged_units, all_units(data, rule))})
    return pd.DataFrame(rows)


def sweep(data: pd.DataFrame, truth: pd.DataFrame, rule: str, parameter: str, values) -> pd.DataFrame:
    """Evaluate one rule for several values of one threshold."""
    rows = []
    for value in values:
        result = evaluate(data, truth, QCConfig(**{parameter: value}))
        row = result[result["rule"] == rule].iloc[0].to_dict()
        rows.append({"threshold": value, **row})
    return pd.DataFrame(rows)


def choose_threshold(sweep_table: pd.DataFrame, min_precision: float) -> float:
    """Pick the threshold with the most recall among those reaching `min_precision`.

    The criterion says: false alarms cost reviewer time, so at least
    `min_precision` of the flags must be real; within that limit, catch as
    many problems as possible. Ties go to the higher (stricter) threshold.
    """
    good = sweep_table[sweep_table["precision"] >= min_precision]
    if good.empty:
        raise ValueError(f"No threshold reaches precision {min_precision}.")
    best = good[good["recall"] == good["recall"].max()]
    return float(best["threshold"].max())
