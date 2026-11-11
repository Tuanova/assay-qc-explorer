# Assay QC Explorer

> **Status: work in progress.** The data model, the synthetic data generator, the
> core calculations, the QC rules and the Streamlit app are in place. The
> evaluation is not written yet, and the rule thresholds have not been tuned.

A small tool for exploring and quality-checking data from a cell-viability
assay run on 96-well plates. It is meant to help a scientist review a plate
quickly: look at the controls, the replicates, and the wells or groups that
deserve a second look.

It supports data review only. It does not make clinical decisions, and it does
not decide on its own that an unusual result is an error.

**All data in this repository are synthetic.** They were produced by the
generator in `assay_qc/generator.py` and do not come from any real experiment.

## Simulated experiment

An ATP-based luminescence viability assay (CellTiter-Glo type) on a 96-well
plate. The signal reflects the amount of ATP in the well. It is an indicator of
viability, not a direct count of living cells.

| Columns | Content |
|---|---|
| 1 and 12 | 16 vehicle controls (cells without compound) |
| 2-10 | 3 compounds x 3 replicates; rows A-H are an 8-point, 3-fold dilution |
| 11 | 8 blanks (medium and reagent, no cells) |

The three compounds are staurosporine (`STS`, a reference cytotoxic compound,
10 uM top concentration), `Compound_A` (active) and `Compound_B` (inactive),
both starting at 100 uM. STS is treated like any other compound: it shows a
dose-response curve and does not define 0 % viability, because its effect
depends on concentration and exposure time.

## Data schema

One row per well.

| Column | Meaning |
|---|---|
| `plate_id` | Identifier of one physical plate |
| `well` | Well address, `A1` to `H12` |
| `sample_type` | `vehicle_control`, `blank` or `compound` |
| `compound` | Compound name; empty for controls and blanks |
| `concentration_uM` | Compound concentration in micromolar; empty for controls and blanks |
| `replicate` | Replicate number within its group, starting at 1 |
| `signal` | Raw luminescence (RLU); empty means the measurement is missing |

Sample types:

- `vehicle_control`: cells with solvent (DMSO) only. Defines 100 % viability.
  The model assumes the same final DMSO concentration in vehicle and compound
  wells, so that differences in viability can be attributed to the compound.
- `blank`: medium and reagent without cells. Defines the background signal.
- `compound`: cells with a compound (test compound or the STS reference).

## Calculations

Per plate, using that plate's own vehicle and blank wells:

```
net signal  = signal - mean_blank                              (background correction)
viability % = 100 * net signal / (mean_vehicle - mean_blank)   (normalisation)
Z'          = 1 - 3 * (SD_vehicle + SD_blank) / |mean_vehicle - mean_blank|
```

- Viability is not clipped to 0-100 %. Values outside that range stay visible.
- A plate with missing vehicle or blank wells, or with a vehicle signal too
  close to the background, is not normalised; the code raises a clear error.
- Replicate CV is computed on the raw signal, because percent viability can be
  close to zero and a CV with a near-zero mean is meaningless.
- Z' is computed from vehicle and blank wells, so it describes the
  signal-to-background window of the plate, not the response to a cytotoxic
  control. It is estimated from a small number of wells, so it is uncertain.
  Any threshold used with it in this project is a configurable demo heuristic.

## QC rules

Each rule returns findings: which plate, well or group to look at, the observed
value, the threshold, a short explanation and a suggested check. Rules never
delete or change data.

| Rule | Unit | Flags |
|---|---|---|
| `missing_signal` | well | a well without a signal value |
| `control_outlier` | vehicle or blank well | robust z-score (median / MAD) above the threshold |
| `replicate_cv` | replicate group | raw-signal CV above the threshold; not judged near the background level |
| `z_prime` | plate | Z' below the threshold, or not computable |

Thresholds are set in `QCConfig` in `assay_qc/rules.py`. They are configurable
demonstration heuristics, not accepted standards, and the current defaults are
placeholders until the evaluation is done.

## Synthetic data and ground truth

`generate_dataset()` creates independent plates with plate-to-plate differences
in brightness and noise, and can inject known problems:

| Issue | Level | What the generator does |
|---|---|---|
| `missing_signal` | well | removes one measurement |
| `replicate_perturbation` | replicate group | multiplies one well of the group by 0.45-0.80 or 1.25-1.90 |
| `control_outlier` | vehicle or blank well | lowers one vehicle well, or raises one blank well |
| `poor_control_separation` | plate | very noisy vehicle wells, or a high and variable background |

The injected problems are saved in a separate truth file. The analysis code
never reads it.

The example in `data/example/` was generated with seed 20, chosen so that the
six plates include clean plates and every issue type. It is a demonstration
set, not an evaluation set.

## Run locally (Windows, PowerShell)

```powershell
git clone https://github.com/Tuanova/assay-qc-explorer.git
cd assay-qc-explorer
py -m venv .venv
.venv\Scripts\Activate.ps1
pip install -r requirements.txt
pytest
python -m scripts.generate_data --n-plates 6 --seed 20
python -m scripts.run_qc data/example/example_plates.csv
streamlit run app.py
```

The last command opens the app in the browser with the built-in example.
Uploaded files are processed in memory and are not saved.

## Repository layout

```
app.py                  Streamlit application (display only)
assay_qc/schema.py      column names, sample types, plate layout
assay_qc/generator.py   synthetic plates and injected problems
assay_qc/metrics.py     normalisation, Z', replicate summary
assay_qc/rules.py       QC rules and the findings table
assay_qc/validation.py  structure checks for uploaded files
scripts/                command-line helpers
tests/                  tests of the calculations and the generator
data/example/           small synthetic example with its truth file
```

## Use of AI

This project is being developed with the help of an AI assistant (Claude).
A precise account of which parts were AI-generated and which were written or
changed by the author will be added here before the project is presented.
