# Assay QC Explorer

> **Status: work in progress.** The data model, the synthetic data generator and
> the core calculations are in place. The QC rules, the Streamlit app and the
> evaluation are not written yet.

A small tool for exploring and quality-checking data from a cell-viability
assay run on 96-well plates. It is meant to help a scientist review a plate
quickly: look at the controls, the replicates, and the wells or groups that
deserve a second look.

It supports data review only. It does not make clinical decisions, and it does
not decide on its own that an unusual result is an error.

**All data in this repository are synthetic.** They were produced by the
generator in `assay_qc/generator.py` and do not come from any real experiment.

## Simulated experiment

A luminescence viability assay (signal proportional to the number of living
cells) on a 96-well plate:

| Columns | Content |
|---|---|
| 1 and 12 | 8 vehicle controls and 8 kill controls, alternating on both sides |
| 2-10 | 3 compounds x 3 replicates; rows A-H are an 8-point, 3-fold dilution from 100 uM |
| 11 | 8 blanks |

## Data schema

One row per well.

| Column | Meaning |
|---|---|
| `plate_id` | Identifier of one physical plate |
| `well` | Well address, `A1` to `H12` |
| `sample_type` | `vehicle_control`, `kill_control`, `blank` or `compound` |
| `compound` | Compound name; empty for controls and blanks |
| `concentration_uM` | Compound concentration in micromolar; empty for controls and blanks |
| `replicate` | Replicate number within its group, starting at 1 |
| `signal` | Raw luminescence (RLU); empty means the measurement is missing |

Sample types:

- `vehicle_control`: cells with solvent only. Defines 100 % viability.
- `kill_control`: cells with a cytotoxic reference. Defines 0 % viability.
- `blank`: medium and reagent without cells. Shows the background; not used in normalisation.
- `compound`: cells with a test compound.

## Calculations

Per plate, using that plate's own controls:

```
viability % = 100 * (signal - mean_kill) / (mean_vehicle - mean_kill)
Z'          = 1 - 3 * (SD_vehicle + SD_kill) / |mean_vehicle - mean_kill|
```

- Viability is not clipped to 0-100 %. Values outside that range stay visible.
- A plate with a missing control, or with controls too close together, is not
  normalised; the code raises a clear error instead.
- Replicate CV is computed on the raw signal, because percent viability can be
  close to zero and a CV with a near-zero mean is meaningless.
- Z' is estimated from 8 wells per control, so it is uncertain. Any threshold
  used with it in this project is a configurable demo heuristic.

## Synthetic data and ground truth

`generate_dataset()` creates independent plates with plate-to-plate differences
in brightness and noise, and can inject known problems:

| Issue | Level | What the generator does |
|---|---|---|
| `missing_signal` | well | removes one measurement |
| `replicate_perturbation` | replicate group | multiplies one well of the group by 0.45-0.80 or 1.25-1.90 |
| `control_outlier` | control well | disturbs one vehicle or kill control well |
| `poor_control_separation` | plate | very noisy controls, or a kill control that is not fully killed |

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
```

## Repository layout

```
assay_qc/schema.py      column names, sample types, plate layout
assay_qc/generator.py   synthetic plates and injected problems
assay_qc/metrics.py     normalisation, Z', replicate summary
scripts/                command-line helpers
tests/                  tests of the calculations and the generator
data/example/           small synthetic example with its truth file
```

## Use of AI

This project is being developed with the help of an AI assistant (Claude).
A precise account of which parts were AI-generated and which were written or
changed by the author will be added here before the project is presented.
