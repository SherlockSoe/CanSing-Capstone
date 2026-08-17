"""Shared data-loading helpers, used by both the notebook and the app."""

import re

import esm
import numpy as np
import pandas as pd
import torch

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
    """Lazily load the ESM2 model/alphabet once, cached at module scope.

    Model loading and embedding extraction (this function, plus
    `get_esm_embedding` / `get_esm_embeddings_batch` / `_mean_pool_residues`
    below) were adapted from the ESM2-Tutorial repository:
    https://github.com/ProteinVision/ESM2-Tutorial/blob/main/ESM2.ipynb
    """
    if not _esm_state:
        if torch.cuda.is_available():
            device = "cuda"
        elif torch.backends.mps.is_available():
            device = "mps"
        else:
            device = "cpu"

        model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
        model.eval()
        model.to(device)

        _esm_state["model"] = model
        _esm_state["alphabet"] = alphabet
        _esm_state["batch_converter"] = alphabet.get_batch_converter()
        _esm_state["last_layer"] = len(model.layers)
        _esm_state["device"] = device
    return (
        _esm_state["model"],
        _esm_state["alphabet"],
        _esm_state["batch_converter"],
        _esm_state["last_layer"],
        _esm_state["device"],
    )


def _mean_pool_residues(token_representations, batch_tokens, padding_idx):
    """Mean-pool per-residue hidden states, excluding the CLS/EOS/padding
    positions."""
    mask = batch_tokens != padding_idx
    seq_lens = mask.sum(dim=1)
    mask = mask.clone()
    mask[:, 0] = False  # CLS token
    mask[torch.arange(mask.size(0)), seq_lens - 1] = False  # EOS token
    mask = mask.unsqueeze(-1).float()
    return (token_representations * mask).sum(1) / mask.sum(1)


def get_esm_embedding(seq):
    """Generate the 1280-dim ESM2 embedding vector for one protein sequence."""
    model, alphabet, batch_converter, last_layer, device = _load_esm_model()
    _, _, batch_tokens = batch_converter([("P", seq)])
    batch_tokens = batch_tokens.to(device)
    with torch.no_grad():
        results = model(
            batch_tokens, repr_layers=[last_layer], return_contacts=False
        )
    token_representations = results["representations"][last_layer]
    pooled = _mean_pool_residues(
        token_representations, batch_tokens, alphabet.padding_idx
    )
    return pooled[0].cpu().numpy()


def get_esm_embeddings_batch(seqs, batch_size=8, max_tokens_per_batch=4096):
    """Generate ESM2 embeddings for many sequences, batched for throughput.

    Sequences are processed longest-first so same-batch padding is
    minimized (a handful of very long outlier sequences otherwise force
    heavy padding on every batch they land in). Batch size is additionally
    capped so that `batch_size * longest_seq_in_batch` stays under
    `max_tokens_per_batch` (mirrors fair-esm's own `extract.py` default) —
    attention memory grows with the *square* of sequence length, so a
    handful of long outliers (e.g. yeast's ~4900-residue MDN1) batched at
    a fixed size can blow past the MPS backend's ~4GB single-buffer limit.
    Returns embeddings in the same order as `seqs`.
    """
    model, alphabet, batch_converter, last_layer, device = _load_esm_model()
    order = sorted(range(len(seqs)), key=lambda i: len(seqs[i]), reverse=True)
    embeddings = [None] * len(seqs)

    start = 0
    while start < len(order):
        longest_in_batch = len(seqs[order[start]])
        max_count = max_tokens_per_batch // longest_in_batch
        count = max(1, min(batch_size, max_count))
        batch_idx = order[start : start + count]
        start += count
        batch_data = [("P", seqs[i]) for i in batch_idx]
        _, _, batch_tokens = batch_converter(batch_data)
        batch_tokens = batch_tokens.to(device)
        with torch.no_grad():
            results = model(
                batch_tokens, repr_layers=[last_layer], return_contacts=False
            )
        token_representations = results["representations"][last_layer]
        pooled = _mean_pool_residues(
            token_representations, batch_tokens, alphabet.padding_idx
        )
        pooled = pooled.cpu().numpy()
        for i, idx in enumerate(batch_idx):
            embeddings[idx] = pooled[i]

    return embeddings


def get_pairwise_features(emb1, emb2):
    """Symmetric pairwise feature vector for two embeddings: the absolute
    difference concatenated with the element-wise product."""
    return np.concatenate([np.abs(emb1 - emb2), emb1 * emb2])
