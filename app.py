"""Assay QC Explorer - Streamlit application.

Run locally from the repository root:

    streamlit run app.py

The app only displays what the functions in `assay_qc/` compute. All
calculations and rules live there, so they can be tested without the app.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

from assay_qc.metrics import NormalizationError, control_stats, normalize_plate, plate_z_prime
from assay_qc.rules import QCConfig, run_qc
from assay_qc.schema import COMPOUND, N_COLS, ROWS
from assay_qc.validation import clean_types, validate_data

EXAMPLE_PATH = Path(__file__).parent / "data" / "example" / "example_plates.csv"

st.set_page_config(page_title="Assay QC Explorer", layout="wide")


@st.cache_data
def load_example() -> pd.DataFrame:
    return pd.read_csv(EXAMPLE_PATH)


def plate_heatmap(plate: pd.DataFrame, value: str, flagged_wells: set[str]) -> go.Figure:
    """Draw the plate as an 8 x 12 grid coloured by `value`."""
    grid = np.full((len(ROWS), N_COLS), np.nan)
    hover = np.full((len(ROWS), N_COLS), "", dtype=object)
    for row in plate.itertuples():
        r, c = ROWS.index(row.well[0]), int(row.well[1:]) - 1
        grid[r, c] = getattr(row, value)
        label = row.sample_type if row.sample_type != COMPOUND else f"{row.compound}, {row.concentration_uM:g} uM"
        hover[r, c] = f"{row.well}: {label}"

    columns = [str(c) for c in range(1, N_COLS + 1)]
    figure = go.Figure(go.Heatmap(
        z=grid, x=columns, y=list(ROWS), customdata=hover, colorscale="Viridis",
        hovertemplate="%{customdata}<br>" + value + ": %{z:.1f}<extra></extra>",
        colorbar={"title": value}, xgap=2, ygap=2,
    ))
    if flagged_wells:
        figure.add_trace(go.Scatter(
            x=[well[1:] for well in flagged_wells], y=[well[0] for well in flagged_wells],
            mode="markers", name="has a finding", hoverinfo="skip",
            marker={"symbol": "circle-open", "size": 26, "color": "red", "line": {"width": 3}},
        ))
    figure.update_yaxes(autorange="reversed")
    figure.update_xaxes(side="top", type="category")
    figure.update_layout(height=420, margin={"l": 10, "r": 10, "t": 30, "b": 10}, showlegend=False)
    return figure


def dose_response_plot(normalized: pd.DataFrame) -> go.Figure:
    """Viability against concentration: every replicate as a point, means as a line."""
    figure = go.Figure()
    compounds = normalized[normalized["sample_type"] == COMPOUND]
    for name, wells in compounds.groupby("compound"):
        means = wells.groupby("concentration_uM")["viability_pct"].mean()
        figure.add_trace(go.Scatter(x=means.index, y=means.values, mode="lines", name=name, legendgroup=name))
        figure.add_trace(go.Scatter(
            x=wells["concentration_uM"], y=wells["viability_pct"], mode="markers",
            name=name, legendgroup=name, showlegend=False, text=wells["well"],
            hovertemplate="%{text}: %{y:.1f} %<extra></extra>",
        ))
    figure.add_hline(y=100, line_dash="dot", line_color="grey")
    figure.add_hline(y=0, line_dash="dot", line_color="grey")
    figure.update_xaxes(type="log", title="Concentration (uM, log scale)")
    figure.update_yaxes(title="Viability (% of vehicle control)")
    figure.update_layout(height=420, margin={"l": 10, "r": 10, "t": 30, "b": 10})
    return figure


def wells_with_findings(plate: pd.DataFrame, findings: pd.DataFrame) -> set[str]:
    """Wells named by a finding, plus all wells of a flagged replicate group."""
    wells = set(findings["well"].dropna())
    for finding in findings[findings["level"] == "group"].itertuples():
        in_group = (plate["compound"] == finding.compound) & np.isclose(
            plate["concentration_uM"], finding.concentration_uM
        )
        wells |= set(plate.loc[in_group, "well"])
    return wells


# --- Sidebar: data source and thresholds -----------------------------------
st.sidebar.header("Data")
uploaded = st.sidebar.file_uploader("Upload a CSV file", type="csv")
st.sidebar.download_button(
    "Download example CSV", EXAMPLE_PATH.read_bytes(), "example_plates.csv", "text/csv"
)
st.sidebar.caption("Uploaded files are processed in memory for this session and are not saved.")

st.sidebar.header("QC thresholds")
st.sidebar.caption("Configurable demo heuristics, not accepted standards.")
config = QCConfig(
    replicate_cv_pct=st.sidebar.slider("Replicate CV above (%)", 5.0, 50.0, QCConfig.replicate_cv_pct, 1.0),
    z_prime_min=st.sidebar.slider("Z' below", 0.0, 1.0, QCConfig.z_prime_min, 0.05),
    control_outlier_z=st.sidebar.slider("Control outlier robust z above", 2.0, 10.0, QCConfig.control_outlier_z, 0.5),
)

# --- Header ----------------------------------------------------------------
st.title("Assay QC Explorer")
st.write(
    "Explore plates from a cell-viability assay and see which wells, replicate "
    "groups or plates deserve a second look. The tool points at data to review; "
    "it does not remove or correct anything and does not decide what is an error."
)

if uploaded is None:
    raw = load_example()
    st.info("Showing the built-in example. These data are **synthetic**, with problems injected on purpose.")
else:
    try:
        raw = pd.read_csv(uploaded)
    except Exception as error:  # unreadable file: show the reason instead of crashing
        st.error(f"Could not read the file as CSV: {error}")
        st.stop()

errors, warnings = validate_data(raw)
if errors:
    st.error("The file cannot be analysed:\n\n" + "\n".join(f"- {message}" for message in errors))
    st.stop()
for message in warnings:
    st.warning(message)

data = clean_types(raw)
findings = run_qc(data, config)
# 'observed' and 'threshold' mix numbers and text ("missing"); show them as text.
findings = findings.astype({"observed": str, "threshold": str})

overview_tab, plate_tab, findings_tab, about_tab = st.tabs(["Overview", "Plate", "Findings", "About"])

# --- Overview: one row per plate -------------------------------------------
with overview_tab:
    rows = []
    for plate_id, plate in data.groupby("plate_id"):
        stats = control_stats(plate)
        rows.append({
            "plate_id": plate_id,
            "wells": len(plate),
            "vehicle mean": stats["mean"].get("vehicle_control"),
            "vehicle CV %": stats["cv_pct"].get("vehicle_control"),
            "blank mean": stats["mean"].get("blank"),
            "Z'": plate_z_prime(plate),
            "findings": int((findings["plate_id"] == plate_id).sum()),
        })
    st.subheader("Plates")
    st.dataframe(
        pd.DataFrame(rows), hide_index=True, width="stretch",
        column_config={
            "vehicle mean": st.column_config.NumberColumn(format="%.0f"),
            "vehicle CV %": st.column_config.NumberColumn(format="%.1f"),
            "blank mean": st.column_config.NumberColumn(format="%.0f"),
            "Z'": st.column_config.NumberColumn(format="%.2f"),
        },
    )
    st.caption("Z' here compares vehicle wells with blanks: it describes the signal-to-background window.")

# --- Plate: map, controls, dose-response -----------------------------------
with plate_tab:
    plate_id = st.selectbox("Plate", sorted(data["plate_id"].unique()))
    plate = data[data["plate_id"] == plate_id]
    plate_findings = findings[findings["plate_id"] == plate_id]

    try:
        normalized = normalize_plate(plate)
    except NormalizationError as error:
        normalized = None
        st.warning(f"This plate cannot be normalised: {error}")

    choices = ["signal"] if normalized is None else ["signal", "viability_pct"]
    value = st.radio("Colour by", choices, horizontal=True)

    map_column, curve_column = st.columns(2)
    with map_column:
        st.subheader("Plate map")
        shown = plate if normalized is None else normalized
        st.plotly_chart(plate_heatmap(shown, value, wells_with_findings(plate, plate_findings)), width="stretch")
        st.caption("Red circles mark wells that belong to a finding.")
    with curve_column:
        st.subheader("Dose-response")
        if normalized is None:
            st.write("Not available without normalisation.")
        else:
            st.plotly_chart(dose_response_plot(normalized), width="stretch")
            st.caption("Points are single wells, lines join the replicate means. Values are not clipped to 0-100 %.")

    st.subheader("Controls on this plate")
    st.dataframe(control_stats(plate).round(1), width="stretch")

    st.subheader(f"Findings on this plate ({len(plate_findings)})")
    st.dataframe(plate_findings.drop(columns="plate_id"), hide_index=True, width="stretch")

    with st.expander("All wells of this plate"):
        st.dataframe(shown, hide_index=True, width="stretch")

# --- Findings: everything, with export -------------------------------------
with findings_tab:
    st.subheader(f"All findings ({len(findings)})")
    selected_rules = st.multiselect("Rules", sorted(findings["rule"].unique()), default=sorted(findings["rule"].unique()))
    visible = findings[findings["rule"].isin(selected_rules)]
    st.dataframe(visible, hide_index=True, width="stretch")
    st.download_button("Download findings as CSV", visible.to_csv(index=False), "findings.csv", "text/csv")
    st.caption("A finding is a reason to look, not a verdict. Nothing is removed from the data.")

# --- About ------------------------------------------------------------------
with about_tab:
    st.markdown(
        """
**What the numbers mean**

- `net signal = signal - mean blank` (background correction)
- `viability % = 100 * net signal / (mean vehicle - mean blank)`
- `Z' = 1 - 3 * (SD vehicle + SD blank) / |mean vehicle - mean blank|`
- Replicate CV is computed on the raw signal of the wells sharing compound and concentration.

**Rules**

| Rule | Flags |
|---|---|
| `missing_signal` | a well without a signal value |
| `control_outlier` | a vehicle or blank well far from the others of its type (median / MAD) |
| `replicate_cv` | a replicate group with CV above the threshold |
| `z_prime` | a plate with Z' below the threshold |

**Limitations**

- The example data are synthetic. Results on them say nothing about real assays.
- Thresholds are demonstration heuristics and must be set for each assay.
- A high replicate CV can be a pipetting problem, but also a real effect on the
  steep part of a dose-response curve. The tool cannot tell these apart.
- This is a data-review aid. It makes no clinical or scientific decision.
"""
    )
