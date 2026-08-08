# ---
# jupyter:
#   jupytext:
#     formats: notebooks//ipynb,notebooks//py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.5
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Background
#
# The following may be helpful when reading or debugginfg this notebook:
#
#
# ### Standard Amino Acids
# FASTA represents each amino acid with a single character as follows:
#
# | 1-Letter Code | 3-Letter Code | Amino Acid Name |
# |---|---|---|
# | A | Ala | Alanine |
# | C | Cys | Cysteine |
# | D | Asp | Aspartic acid |
# | E | Glu | Glutamic acid |
# | F | Phe | Phenylalanine |
# | G | Gly | Glycine |
# | H | His | Histidine |
# | I | Ile | Isoleucine |
# | K | Lys | Lysine |
# | L | Leu | Leucine |
# | M | Met | Methionine |
# | N | Asn | Asparagine |
# | P | Pro | Proline |
# | Q | Gln | Glutamine |
# | R | Arg | Arginine |
# | S | Ser | Serine |
# | T | Thr | Threonine |
# | V | Val | Valine |
# | W | Trp | Tryptophan |
# | Y | Tyr | Tyrosine |
#
# ### Ambiguous & Special Characters
# FASTA also uses the following abbreviations for ambiguous amino acid
# identification and special characters:
#
# | 1-Letter Code | Description / Meaning |
# |---|---|
# | B | Aspartic acid (D) or Asparagine (N) |
# | J | Leucine (L) or Isoleucine (I) |
# | X | Unknown or any amino acid |
# | Z | Glutamic acid (E) or Glutamine (Q) |
# | * | Translation stop codon |
# | - | Gap of missing or unsequenced amino acid |

# %% [markdown]
# # Environment Setup
#
# In this section we're setting up the environment for the rest of the
# notebook.

# %%
# Core libraries
import pandas as pd
import numpy as np
import warnings
import sys

# File read / write
from pathlib import Path
import pickle

from Bio import SeqIO

# Environment Settings
pd.set_option("display.max_rows", None)
pd.set_option("display.max_columns", None)
np.random.seed(42)
warnings.filterwarnings("ignore")


# Filepaths for source data files
# Resolve the repo root by walking up from the working directory until a
# marker file is found, so this notebook runs the same way regardless of
# whether Jupyter/VSCode launches it from the repo root or from notebooks/
# (this differs by machine and editor between collaborators).
def find_repo_root(marker="requirements.txt"):
    for candidate in [Path.cwd().resolve(), *Path.cwd().resolve().parents]:
        if (candidate / marker).exists():
            return candidate
    raise FileNotFoundError(f"Could not locate repo root (missing {marker})")


REPO_ROOT = find_repo_root()
DATA_DIR = REPO_ROOT / "data" / "raw"
PROCESSED_DIR = REPO_ROOT / "data" / "processed"
PROCESSED_DIR.mkdir(parents=True, exist_ok=True)

# Make src/ importable so helper functions can be shared with the
# Streamlit app in app/ instead of being duplicated in both places.
sys.path.insert(0, str(REPO_ROOT))

path1 = DATA_DIR / "uniprotkb_proteome_UP000002311.fasta"
path2 = (
    DATA_DIR
    / "BIOGRID-ORGANISM-Saccharomyces_cerevisiae_S288c-5.0.259.mitab.txt"
)


# %%
# #%pip install biopython

# %% [markdown]
# # Helper Functions
#
# In this section we're defining some helper functions for use in the rest
# of the notebook.


# %%
# Shared with app/ (Streamlit GUI) via src/ppi_utils.py, so the data
# loading logic has a single source of truth instead of being duplicated.
from src.ppi_utils import (
    read_fasta,
    extract_locuslink,
    extract_gene_names,
    get_interactors,
    load_biogrid_interactions,
)


# %% [markdown]
# # Section 1 - Data Import
#
# In this section we are importing amino acid sequences for all 6,067 amino
# acids in S. cerevisiae.

# %%
# Import dataset
for record in SeqIO.parse(path1, "fasta"):
    print(f"ID: {record.id}")
    print(f"Sequence: {record.seq}")
    print(f"Length: {len(record.seq)}\n")

# %%
# Create the protein-sequence dictionary and a dictionary of protein
# residue length
fasta_data = read_fasta(path1)
seq_dic = extract_gene_names(fasta_data)
len_dic = {gene: len(seq) for gene, seq in seq_dic.items()}

# Create a list of all proteins
gene_list = list(seq_dic.keys())

# %% [markdown]
# Here we are importing the interaction data from BioGrid.

# %%
# Inspect the results
# seq_dic
# len_dic
gene_list
len(gene_list)

# %%
# Check if a particular gene is included in the Uniprot data
if "CDC73" in gene_list:
    print("Item found!")

# %%
# Now extract the Biogrid data and create the Biogrid dataframe
bg = load_biogrid_interactions(path2)

print(bg.shape)

# %%
bg.head()

# %%
# Check the Biogrid data at a particular location
print(bg.at[1, "Alt IDs Interactor A"])

# %%
# Extract the interactor A and B genes and create a new column for each
bg["Gene_A"] = bg["Alt IDs Interactor A"].apply(extract_locuslink)
bg["Gene_B"] = bg["Alt IDs Interactor B"].apply(extract_locuslink)

bg.head()

# %%
# Now we make the interaction dictionary
ia_dic = {}

for gene in gene_list:
    interactors = get_interactors(bg, gene)
    ia_dic[gene] = interactors

# %%
# Clean up the interaction dictionary

# Write the interaction dictionary to file
with open(PROCESSED_DIR / "interactions.pkl", "wb") as file:
    pickle.dump(ia_dic, file)

# %%
# Reading dictionary from the binary file
with open(PROCESSED_DIR / "interactions.pkl", "rb") as file:
    interaction_dictionary = pickle.load(file)

# %%
# Inspect the results
interaction_dictionary

# %% [markdown]
# # Section 2 - Create Embeddings
#
# In this section we are creating the embeddings for all 6,067 amino acids
# in S. cerevisiae.

# %%

# %% [markdown]
# # Section 3 - Train the Model
#
# In this section we are

# %%

# %% [markdown]
# # Section 4 - Evaluate the Model
#
# In this section we are

# %%

# %%
