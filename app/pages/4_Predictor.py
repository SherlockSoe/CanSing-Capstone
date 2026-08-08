"""Predictor page: given two genes, predict interaction probability.

Expects `data/processed/model.pkl`, written by the notebook's Section 5
(Evaluate the Model). Uses `read_fasta`/`extract_gene_names`/
`create_embedding` from `src/ppi_utils.py` — the same embedding function
used to build the model's training features (feature vector = the two
proteins' 1280-dim embeddings concatenated into 2560 dims).
"""

import pickle
import sys
from pathlib import Path

import numpy as np
import streamlit as st

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))

from src.ppi_utils import (  # noqa: E402
    read_fasta,
    extract_gene_names,
    create_embedding,
)

MODEL_PATH = REPO_ROOT / "data" / "processed" / "model.pkl"
FASTA_PATH = (
    REPO_ROOT / "data" / "raw" / "uniprotkb_proteome_UP000002311.fasta"
)

st.set_page_config(page_title="PPI Prediction | Predictor", page_icon="🧬")
st.title("Interaction Predictor")

if not MODEL_PATH.exists():
    st.info(
        "Not available yet. This page needs the notebook's Section 5 "
        "(Evaluate the Model) to run first, which writes "
        f"`{MODEL_PATH.relative_to(REPO_ROOT)}`."
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
    with st.spinner("Generating embeddings and predicting..."):
        embedding_a = create_embedding(seq_dic[gene_a])
        embedding_b = create_embedding(seq_dic[gene_b])
        features = np.concatenate([embedding_a, embedding_b])
        probability = model.predict_proba([features])[0][1]
    st.metric("Predicted interaction probability", f"{probability:.1%}")
