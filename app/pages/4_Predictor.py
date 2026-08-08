"""Predictor page: given two genes, predict interaction probability.

Expects `data/processed/model.pkl` (trained model) and a
`create_embedding(seq)` function in `src/ppi_utils.py` (added once the
notebook's Section 2 embedding pipeline exists). Once both are in place,
this page works without any further changes.
"""

import pickle
import sys
from pathlib import Path

import streamlit as st

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))

from src.ppi_utils import read_fasta, extract_gene_names  # noqa: E402

try:
    from src.ppi_utils import create_embedding
except ImportError:
    create_embedding = None

MODEL_PATH = REPO_ROOT / "data" / "processed" / "model.pkl"
FASTA_PATH = (
    REPO_ROOT / "data" / "raw" / "uniprotkb_proteome_UP000002311.fasta"
)

st.set_page_config(page_title="PPI Prediction | Predictor", page_icon="🧬")
st.title("Interaction Predictor")

if create_embedding is None or not MODEL_PATH.exists():
    missing = []
    if create_embedding is None:
        missing.append("`create_embedding()` in `src/ppi_utils.py`")
    if not MODEL_PATH.exists():
        missing.append(f"`{MODEL_PATH.relative_to(REPO_ROOT)}`")
    st.info(
        "Not available yet. This page needs the model pipeline from the "
        "notebook's Sections 2-3 to be finished first:\n\n"
        + "\n".join(f"- {m}" for m in missing)
    )
    st.stop()


@st.cache_data
def _load_fasta(path):
    return read_fasta(str(path))


@st.cache_resource
def _load_model(path):
    with open(path, "rb") as f:
        return pickle.load(f)


fasta_data = _load_fasta(FASTA_PATH)
seq_dic = extract_gene_names(fasta_data)
model = _load_model(MODEL_PATH)

gene_options = sorted(seq_dic.keys())
col1, col2 = st.columns(2)
gene_a = col1.selectbox("Protein A", gene_options)
gene_b = col2.selectbox("Protein B", gene_options)

if st.button("Predict"):
    embedding_a = create_embedding(seq_dic[gene_a])
    embedding_b = create_embedding(seq_dic[gene_b])
    probability = model.predict_proba([embedding_a + embedding_b])[0][1]
    st.metric("Predicted interaction probability", f"{probability:.1%}")
