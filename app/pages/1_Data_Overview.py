"""Data Overview page: real numbers from the raw UniProt/BioGRID data."""

import sys
from pathlib import Path

import streamlit as st

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))

from src.ppi_utils import read_fasta, load_biogrid_interactions  # noqa: E402

st.set_page_config(page_title="PPI Prediction | Data Overview", page_icon="🧬")
st.title("Data Overview")

DATA_DIR = REPO_ROOT / "data" / "raw"
FASTA_PATH = DATA_DIR / "uniprotkb_proteome_UP000002311.fasta"
BIOGRID_PATH = (
    DATA_DIR
    / "BIOGRID-ORGANISM-Saccharomyces_cerevisiae_S288c-5.0.259.mitab.txt"
)

missing = [p for p in (FASTA_PATH, BIOGRID_PATH) if not p.exists()]
if missing:
    st.warning(
        "Raw data not found. Run the download scripts first:\n\n"
        "```bash\n"
        "python scripts/download_biogrid.py\n"
        "python scripts/download_uniprot.py\n"
        "```\n\n"
        f"Missing: {', '.join(p.name for p in missing)}"
    )
    st.stop()


@st.cache_data
def _load_fasta(path):
    return read_fasta(str(path))


@st.cache_data
def _load_interactions(path):
    return load_biogrid_interactions(str(path))


fasta_data = _load_fasta(FASTA_PATH)
bg = _load_interactions(BIOGRID_PATH)

col1, col2 = st.columns(2)
col1.metric("Protein sequences (UniProt)", f"{len(fasta_data):,}")
col2.metric("Interaction records (BioGRID)", f"{len(bg):,}")

st.subheader("Sample protein sequences")
sample_rows = [
    {"Sequence ID": seq_id, "Length (aa)": len(seq)}
    for seq_id, seq in list(fasta_data.items())[:10]
]
st.dataframe(sample_rows)

st.subheader("Sample interaction records")
st.dataframe(bg.head(10))
