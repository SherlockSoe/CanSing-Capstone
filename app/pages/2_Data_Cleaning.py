"""Data Cleaning page: Sankey diagram of how proteins were filtered."""

import json
import sys
from pathlib import Path

import plotly.graph_objects as go
import streamlit as st

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))

CLEANING_SUMMARY_PATH = (
    REPO_ROOT / "data" / "processed" / "cleaning_summary.json"
)

st.set_page_config(page_title="PPI Prediction | Data Cleaning", page_icon="🧬")
st.title("Data Cleaning")

st.markdown(
    "The 6,067 candidate proteins are narrowed down to those with both "
    "sequence data (UniProt) and at least one recorded interaction "
    "(BioGRID) — proteins with no recorded interactions are dropped "
    "before building the interaction matrix used for training."
)

if not CLEANING_SUMMARY_PATH.exists():
    st.info(
        "Not generated yet. Once the notebook's data-cleaning step writes "
        f"its summary to `{CLEANING_SUMMARY_PATH.relative_to(REPO_ROOT)}`, "
        "this page will render the Sankey diagram automatically."
    )
    st.stop()

with open(CLEANING_SUMMARY_PATH) as f:
    summary = json.load(f)

fig = go.Figure(
    go.Sankey(
        node=dict(label=summary["labels"], pad=20, thickness=20),
        link=dict(
            source=summary["source"],
            target=summary["target"],
            value=summary["value"],
        ),
    )
)
fig.update_layout(title="Protein filtering flow", font_size=12)
st.plotly_chart(fig)

n_all = summary["value"][0] + summary["value"][1]
n_kept = summary["value"][0]
col1, col2 = st.columns(2)
col1.metric("Proteins before cleaning", f"{n_all:,}")
col2.metric("Proteins kept", f"{n_kept:,}", f"-{n_all - n_kept:,}")
