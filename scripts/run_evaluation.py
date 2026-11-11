"""Evaluate the QC rules on synthetic plates and save the results.

Run from the repository root:

    python -m scripts.run_evaluation

Protocol
  1. Thresholds are chosen on TUNING plates only.
  2. Results are reported on separate TEST plates, generated with other seeds.
     The test plates are never used to choose anything.
  3. The split is by whole plates, never by wells of the same plate.
  4. Everything is repeated for several seed pairs to show how much the
     numbers move.
"""

from pathlib import Path

import numpy as np
import pandas as pd

from assay_qc.evaluation import choose_threshold, evaluate, sweep
from assay_qc.generator import generate_dataset
from assay_qc.metrics import summarize_replicates
from assay_qc.rules import QCConfig, run_qc

N_PLATES = 40             # per set
N_REPEATS = 5             # seed pairs
MIN_PRECISION = 0.90      # criterion, fixed before looking at any result
CV_GRID = np.arange(5.0, 50.5, 1.0)
OUTLIER_GRID = np.arange(2.5, 12.5, 0.5)
OUT_DIR = Path("results")


def group_examples(data, truth, config, n=5) -> pd.DataFrame:
    """A few missed groups and a few false alarms, for the replicate CV rule."""
    keys = ["plate_id", "compound", "concentration_uM"]
    cv = pd.concat(
        summarize_replicates(plate).assign(plate_id=plate_id)
        for plate_id, plate in data.groupby("plate_id")
    )[keys + ["cv_pct"]]
    findings = run_qc(data, config)
    flagged = findings.loc[findings["rule"] == "replicate_cv", keys].assign(flagged=True)
    injected = truth.loc[truth["issue_type"] == "replicate_perturbation", keys + ["detail"]]

    table = cv.merge(injected, how="left", on=keys).merge(flagged, how="left", on=keys)
    table["flagged"] = table["flagged"].eq(True)
    missed = table[table["detail"].notna() & ~table["flagged"]].assign(kind="missed")
    false_alarm = table[table["detail"].isna() & table["flagged"]].assign(kind="false alarm")
    columns = ["kind", *keys, "cv_pct", "detail"]
    return pd.concat([missed.head(n), false_alarm.head(n)])[columns].round({"cv_pct": 1})


def main() -> None:
    OUT_DIR.mkdir(exist_ok=True)
    summaries = []

    for repeat in range(N_REPEATS):
        tuning, tuning_truth = generate_dataset(N_PLATES, seed=100 + repeat, plate_prefix="TUNE")
        test, test_truth = generate_dataset(N_PLATES, seed=200 + repeat, plate_prefix="TEST")

        cv_sweep = sweep(tuning, tuning_truth, "replicate_cv", "replicate_cv_pct", CV_GRID)
        outlier_sweep = sweep(tuning, tuning_truth, "control_outlier", "control_outlier_z", OUTLIER_GRID)
        config = QCConfig(
            replicate_cv_pct=choose_threshold(cv_sweep, MIN_PRECISION),
            control_outlier_z=choose_threshold(outlier_sweep, MIN_PRECISION),
            # Z' threshold stays at its default: too few plate-level cases to tune.
        )

        result = evaluate(test, test_truth, config)
        result.insert(0, "repeat", repeat + 1)
        result.insert(2, "threshold", [config.replicate_cv_pct, config.control_outlier_z, config.z_prime_min])
        summaries.append(result)

        if repeat == 0:
            # Curves for the app: same thresholds evaluated on tuning and on test plates.
            curves = pd.concat([
                cv_sweep.assign(plates="tuning"),
                sweep(test, test_truth, "replicate_cv", "replicate_cv_pct", CV_GRID).assign(plates="test"),
            ])
            curves.assign(chosen_threshold=config.replicate_cv_pct).to_csv(
                OUT_DIR / "cv_threshold_sweep.csv", index=False)
            group_examples(test, test_truth, config).to_csv(OUT_DIR / "cv_examples.csv", index=False)

    summary = pd.concat(summaries, ignore_index=True)
    summary.to_csv(OUT_DIR / "evaluation_summary.csv", index=False)

    pd.set_option("display.width", 200)
    print(summary.drop(columns=["accuracy", "accuracy_if_nothing_flagged"]).round(2).to_string(index=False))
    print("\nTotals over all repeats (test plates):")
    totals = summary.groupby("rule")[["units", "TP", "FP", "FN"]].sum()
    totals["precision"] = totals["TP"] / (totals["TP"] + totals["FP"])
    totals["recall"] = totals["TP"] / (totals["TP"] + totals["FN"])
    print(totals.round(2).to_string())


if __name__ == "__main__":
    main()
