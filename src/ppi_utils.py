"""Shared data-loading helpers, used by both the notebook and the app."""

import re

import pandas as pd
import torch
from transformers import AutoTokenizer, AutoModel

_ESM_MODEL_NAME = "facebook/esm2_t33_650M_UR50D"
_esm_state = {}


def read_fasta(file_path):
    sequences = {}
    current_id = None
    current_seq = []

    with open(file_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            # Lines starting with '>' indicate a new sequence ID
            if line.startswith(">"):
                if current_id:
                    sequences[current_id] = "".join(current_seq)
                current_id = line[1:]  # Remove the '>' character
                current_seq = []
            else:
                current_seq.append(line)

        # Add the final sequence from the file
        if current_id:
            sequences[current_id] = "".join(current_seq)

    return sequences


def extract_locuslink(s):
    match = re.search(r"locuslink:([^|]+)", s)
    return match.group(1) if match else None


def extract_gene_names(fasta_data):
    """Map UniProt sequence IDs to gene symbols via the GN= field."""
    genes = {}
    for seq_id, sequence in fasta_data.items():
        match = re.search(r"GN=(\S+)", seq_id)
        if match:
            genes[match.group(1)] = sequence
    return genes


def load_biogrid_interactions(file_path):
    """Load a BioGRID mitab file into a cleaned interactions dataframe."""
    bg = pd.read_csv(file_path, sep="\t", header=None, dtype=str)
    bg.columns = bg.iloc[0]
    bg = bg[1:]
    bg.columns.name = None
    bg = bg.reset_index(drop=True)
    return bg


def get_interactors(bg, protein_id):
    # Protein appears in column A
    a_partners = bg.loc[bg["Gene_A"] == protein_id, "Gene_B"]

    # Protein appears in column B
    b_partners = bg.loc[bg["Gene_B"] == protein_id, "Gene_A"]

    # Combine, deduplicate, and return as a list
    return pd.concat([a_partners, b_partners]).drop_duplicates().tolist()


def _load_esm_model():
    """Lazily load the ESM2 model/tokenizer once, cached at module scope."""
    if not _esm_state:
        device = "mps" if torch.backends.mps.is_available() else "cpu"
        _esm_state["tokenizer"] = AutoTokenizer.from_pretrained(
            _ESM_MODEL_NAME
        )
        _esm_state["model"] = (
            AutoModel.from_pretrained(_ESM_MODEL_NAME).to(device).eval()
        )
        _esm_state["device"] = device
    return _esm_state["tokenizer"], _esm_state["model"], _esm_state["device"]


def _mean_pool_residues(last_hidden_state, attention_mask):
    """Mean-pool per-residue hidden states, excluding the CLS/EOS tokens."""
    mask = attention_mask.clone()
    seq_lens = attention_mask.sum(dim=1)
    mask[:, 0] = 0  # CLS token
    mask[torch.arange(mask.size(0)), seq_lens - 1] = 0  # EOS token
    mask = mask.unsqueeze(-1).float()
    return (last_hidden_state * mask).sum(1) / mask.sum(1)


def create_embedding(seq):
    """Generate the 1280-dim ESM2 embedding vector for one protein sequence."""
    tokenizer, model, device = _load_esm_model()
    inputs = tokenizer(
        seq, return_tensors="pt", truncation=True, max_length=1022
    ).to(device)
    with torch.no_grad():
        output = model(**inputs)
    pooled = _mean_pool_residues(
        output.last_hidden_state, inputs["attention_mask"]
    )
    return pooled[0].cpu().numpy()


def create_embeddings_batch(seqs, batch_size=8):
    """Generate ESM2 embeddings for many sequences, batched for throughput.

    Sequences are processed longest-first so same-batch padding is
    minimized (a handful of very long outlier sequences otherwise force
    heavy padding on every batch they land in). Returns embeddings in the
    same order as `seqs`.
    """
    tokenizer, model, device = _load_esm_model()
    order = sorted(range(len(seqs)), key=lambda i: len(seqs[i]), reverse=True)
    embeddings = [None] * len(seqs)

    for start in range(0, len(order), batch_size):
        batch_idx = order[start : start + batch_size]
        batch_seqs = [seqs[i] for i in batch_idx]
        inputs = tokenizer(
            batch_seqs,
            return_tensors="pt",
            padding=True,
            truncation=True,
            max_length=1022,
        ).to(device)
        with torch.no_grad():
            output = model(**inputs)
        pooled = _mean_pool_residues(
            output.last_hidden_state, inputs["attention_mask"]
        )
        pooled = pooled.cpu().numpy()
        for i, idx in enumerate(batch_idx):
            embeddings[idx] = pooled[i]

    return embeddings
