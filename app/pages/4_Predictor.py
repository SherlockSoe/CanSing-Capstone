"""Predictor page: given two genes, predict interaction probability.

Expects `data/processed/model.keras`, written by the notebook's Section 4
(Train the Model). Uses `read_fasta`/`extract_gene_names`/
`get_esm_embedding`/`get_pairwise_features` from `src/ppi_utils.py` — the
same embedding and feature-pairing functions used to build the model's
training features (feature vector = the two proteins' 1280-dim embeddings
paired via absolute-difference + element-wise-product into 2560 dims).
"""

import sys
from pathlib import Path

import streamlit as st
from tensorflow import keras

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))

from src.ppi_utils import (  # noqa: E402
    read_fasta,
    extract_gene_names,
    get_esm_embedding,
    get_pairwise_features,
)

MODEL_PATH = REPO_ROOT / "data" / "processed" / "model.keras"
FASTA_PATH = (
    REPO_ROOT / "data" / "raw" / "uniprotkb_proteome_UP000002311.fasta"
)

st.set_page_config(page_title="PPI Prediction | Predictor", page_icon="🧬")
st.title("Interaction Predictor")

if not MODEL_PATH.exists():
    st.info(
        "Not available yet. This page needs the notebook's Section 4 "
        "(Train the Model) to run first, which writes "
        f"`{MODEL_PATH.relative_to(REPO_ROOT)}`."
    )
    st.stop()


@st.cache_data
def _load_fasta(path):
    return read_fasta(str(path))


@st.cache_resource
def _load_model(path):
    return keras.models.load_model(path)


fasta_data = _load_fasta(FASTA_PATH)
seq_dic = extract_gene_names(fasta_data)
model = _load_model(MODEL_PATH)

gene_options = sorted(seq_dic.keys())
col1, col2 = st.columns(2)
gene_a = col1.selectbox("Protein A", gene_options)
gene_b = col2.selectbox("Protein B", gene_options)

if st.button("Predict"):
    with st.spinner("Generating embeddings and predicting..."):
        embedding_a = get_esm_embedding(seq_dic[gene_a])
        embedding_b = get_esm_embedding(seq_dic[gene_b])
        features = get_pairwise_features(embedding_a, embedding_b)
        probability = model.predict(features[None, :])[0][0]
    st.metric("Predicted interaction probability", f"{probability:.1%}")
