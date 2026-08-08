"""Shared data-loading helpers, used by both the notebook and the app."""

import re

import pandas as pd


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
