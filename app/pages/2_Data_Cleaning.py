"""Data Cleaning page: Sankey diagram of how proteins were filtered."""

import sys
from pathlib import Path

import streamlit as st

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))

CLEANING_SUMMARY_PATH = (
    REPO_ROOT / "data" / "processed" / "cleaning_summary.json"
)

st.set_page_config(page_title="PPI Prediction | Data Cleaning", page_icon="🧬")
st.title("Data Cleaning")

st.markdown(
    "This page will show a Sankey diagram of how the 6,067 candidate "
    "proteins were narrowed down — for example, limiting the allowable "
    "overlap between protein sequences to mitigate data leakage during "
    "training."
)

if not CLEANING_SUMMARY_PATH.exists():
    st.info(
        "Not generated yet. Once the notebook's data-cleaning step writes "
        f"its summary to `{CLEANING_SUMMARY_PATH.relative_to(REPO_ROOT)}`, "
        "this page will render the Sankey diagram automatically."
    )
    st.stop()

st.warning(
    f"Found {CLEANING_SUMMARY_PATH.name}, but Sankey rendering isn't "
    "implemented yet — add it here once the summary schema is finalized."
)
