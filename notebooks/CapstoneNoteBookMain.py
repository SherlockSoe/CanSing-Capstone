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
import json
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
# # Section 2 - Data Cleaning
#
# In this section we narrow the 6,067 candidate proteins down to a smaller,
# non-redundant set. We cluster near-identical sequences with CD-HIT and
# keep only one representative protein per cluster, which limits how much
# sequence overlap can leak across the train/test split in Section 4.

# %%
# Cluster protein sequences with CD-HIT at a 40% identity threshold (a
# standard low-redundancy threshold in PPI-prediction literature) and keep
# only each cluster's representative sequence. Requires the `cd-hit`
# binary (`brew install brewsci/bio/cd-hit`, see README).
import subprocess

CDHIT_INPUT = PROCESSED_DIR / "cdhit_input.fasta"
CDHIT_OUTPUT = PROCESSED_DIR / "cdhit_out"

with open(CDHIT_INPUT, "w") as f:
    for gene, seq in seq_dic.items():
        f.write(f">{gene}\n{seq}\n")

subprocess.run(
    [
        "cd-hit",
        "-i",
        str(CDHIT_INPUT),
        "-o",
        str(CDHIT_OUTPUT),
        "-c",
        "0.4",  # 40% sequence identity threshold
        "-n",
        "2",  # word length matching the 0.4-0.5 identity range
        "-M",
        "0",  # no memory limit
        "-d",
        "0",  # keep full sequence names in the .clstr file
    ],
    check=True,
)

# %%
# Parse the .clstr file: each cluster's representative is marked with "*".
cluster_of = {}
representative_of_cluster = {}
cluster_id = None

with open(f"{CDHIT_OUTPUT}.clstr") as f:
    for line in f:
        if line.startswith(">Cluster"):
            cluster_id = int(line.split()[-1])
        else:
            gene = line.split(">")[1].split("...")[0]
            cluster_of[gene] = cluster_id
            if line.strip().endswith("*"):
                representative_of_cluster[cluster_id] = gene

kept_genes = set(representative_of_cluster.values())

print(f"Proteins before cleaning: {len(seq_dic)}")
print(f"Clusters found: {len(representative_of_cluster)}")
print(f"Proteins kept (cluster representatives): {len(kept_genes)}")

# %%
# Write a Sankey-ready summary for the Streamlit app's Data Cleaning page:
# all proteins -> kept/dropped, then kept -> has/lacks any recorded BioGRID
# interaction (informational only; the physical-vs-genetic interaction
# filter used for modeling is applied later, in Section 4).
n_kept = len(kept_genes)
n_dropped = len(seq_dic) - n_kept

genes_with_any_interaction = set(bg["Gene_A"]) | set(bg["Gene_B"])
kept_with_interactions = len(kept_genes & genes_with_any_interaction)
kept_without_interactions = n_kept - kept_with_interactions

cleaning_summary = {
    "labels": [
        "All proteins",
        "Kept (non-redundant)",
        "Dropped (redundant)",
        "Has interaction data",
        "No interaction data",
    ],
    "source": [0, 0, 1, 1],
    "target": [1, 2, 3, 4],
    "value": [
        n_kept,
        n_dropped,
        kept_with_interactions,
        kept_without_interactions,
    ],
}

with open(PROCESSED_DIR / "cleaning_summary.json", "w") as f:
    json.dump(cleaning_summary, f)

# %%
# Keep only cluster-representative proteins and the interactions between
# them.
gene_list = [gene for gene in gene_list if gene in kept_genes]
seq_dic = {gene: seq for gene, seq in seq_dic.items() if gene in kept_genes}
len_dic = {gene: n for gene, n in len_dic.items() if gene in kept_genes}
bg = bg[bg["Gene_A"].isin(kept_genes) & bg["Gene_B"].isin(kept_genes)]
bg = bg.reset_index(drop=True)

print(f"Interaction rows after cleaning: {len(bg)}")

# %% [markdown]
# # Section 3 - Create Embeddings
#
# In this section we generate ESM-2 embeddings for each of the remaining
# (post-cleaning) proteins, using the `create_embeddings_batch` helper
# from `src/ppi_utils.py` (the same module the Streamlit predictor uses).

# %%
from src.ppi_utils import create_embeddings_batch

EMBEDDINGS_PATH = PROCESSED_DIR / "embeddings.pkl"

if EMBEDDINGS_PATH.exists():
    with open(EMBEDDINGS_PATH, "rb") as f:
        embedding_dic = pickle.load(f)
else:
    embedded_gene_order = list(seq_dic.keys())
    vectors = create_embeddings_batch(list(seq_dic.values()), batch_size=8)
    embedding_dic = dict(zip(embedded_gene_order, vectors))
    with open(EMBEDDINGS_PATH, "wb") as f:
        pickle.dump(embedding_dic, f)

print(f"Embedded {len(embedding_dic)} proteins")
print(f"Embedding dimension: {len(next(iter(embedding_dic.values())))}")

# %% [markdown]
# ## Validating the embeddings
#
# As a sanity check, we visualize the embedding space with PCA and confirm
# that known interacting protein pairs are, on average, more similar
# (higher cosine similarity) than random protein pairs -- a simple signal
# that the embeddings capture something biologically meaningful before we
# build a classifier on top of them.

# %%
from sklearn.decomposition import PCA
import plotly.express as px

genes_embedded = list(embedding_dic.keys())
embedding_matrix = np.array([embedding_dic[g] for g in genes_embedded])

coords = PCA(n_components=2, random_state=42).fit_transform(embedding_matrix)
pca_df = pd.DataFrame(
    {
        "PC1": coords[:, 0],
        "PC2": coords[:, 1],
        "length": [len_dic[g] for g in genes_embedded],
        "gene": genes_embedded,
    }
)
px.scatter(
    pca_df,
    x="PC1",
    y="PC2",
    color="length",
    hover_name="gene",
    title="ESM-2 embeddings (PCA), colored by sequence length",
)

# %%
# Compare cosine similarity for known-interacting pairs vs. random pairs.


def cosine_sim(a, b):
    return np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b))


interacting_pairs = list(
    bg[["Gene_A", "Gene_B"]]
    .drop_duplicates()
    .itertuples(index=False, name=None)
)
interacting_sims = [
    cosine_sim(embedding_dic[a], embedding_dic[b])
    for a, b in interacting_pairs
    if a in embedding_dic and b in embedding_dic
]

validation_rng = np.random.default_rng(42)
random_sims = []
while len(random_sims) < len(interacting_sims):
    a, b = validation_rng.choice(genes_embedded, size=2, replace=False)
    random_sims.append(cosine_sim(embedding_dic[a], embedding_dic[b]))

print(
    "Mean cosine similarity, interacting pairs: "
    f"{np.mean(interacting_sims):.4f}"
)
print(f"Mean cosine similarity, random pairs:      {np.mean(random_sims):.4f}")

# %% [markdown]
# # Section 4 - Train the Model
#
# In this section we build a labeled dataset of interacting (positive) and
# non-interacting (negative) protein pairs, split it into train/test sets
# without leaking sequence-similar proteins across the split, and train
# several candidate models for comparison.

# %% [markdown]
# ## Defining a positive interaction
#
# BioGRID's yeast data is dominated by large-scale *genetic* interaction
# screens (synthetic lethality, genetic suppression, dosage effects, etc.)
# rather than direct physical binding. Since this project defines PPI as a
# physical association (see Background), we restrict positive pairs to
# BioGRID's physical interaction types: physical association (MI:0915),
# direct interaction (MI:0407), and association (MI:0914) -- excluding
# colocalization and every genetic-interaction category.

# %%
PHYSICAL_INTERACTION_TYPES = ("MI:0915", "MI:0407", "MI:0914")
is_physical = bg["Interaction Types"].str.contains(
    "|".join(PHYSICAL_INTERACTION_TYPES)
)
bg_physical = bg[is_physical]

print(f"Interaction rows (all types): {len(bg)}")
print(f"Interaction rows (physical only): {len(bg_physical)}")

# %%
# Positive pairs: deduplicated physical interactions (unordered), excluding
# self-pairs and restricted to proteins we generated embeddings for.
positive_pairs = set()
for gene_a, gene_b in bg_physical[["Gene_A", "Gene_B"]].itertuples(
    index=False
):
    if gene_a == gene_b:
        continue
    if gene_a not in embedding_dic or gene_b not in embedding_dic:
        continue
    positive_pairs.add(tuple(sorted((gene_a, gene_b))))

positive_pairs = list(positive_pairs)
print(f"Positive pairs: {len(positive_pairs)}")

# %%
# Negative pairs: random protein pairs with no recorded BioGRID
# interaction (of any type), sampled 1:1 against the positive pairs. This
# is the standard approach for PPI datasets, since BioGRID only records
# observed interactions -- a randomly-sampled "negative" could in
# principle be a true, simply unrecorded interaction. We discuss this
# limitation in the written report.
all_known_pairs = {
    tuple(sorted((a, b)))
    for a, b in bg[["Gene_A", "Gene_B"]].itertuples(index=False)
    if a != b
}

sampling_rng = np.random.default_rng(42)
embedded_genes = np.array(list(embedding_dic.keys()))
negative_pairs = set()

while len(negative_pairs) < len(positive_pairs):
    gene_a, gene_b = sampling_rng.choice(embedded_genes, size=2, replace=False)
    pair = tuple(sorted((gene_a, gene_b)))
    if pair not in all_known_pairs and pair not in negative_pairs:
        negative_pairs.add(pair)

negative_pairs = list(negative_pairs)
print(f"Negative pairs: {len(negative_pairs)}")

# %%
# Split at the CD-HIT cluster level (Section 2), so no two proteins from
# the same cluster end up on both sides of the train/test split.
embedded_clusters = sorted({cluster_of[g] for g in embedded_genes})
split_rng = np.random.default_rng(42)
split_rng.shuffle(embedded_clusters)

split_point = int(len(embedded_clusters) * 0.8)
train_clusters = set(embedded_clusters[:split_point])
test_clusters = set(embedded_clusters[split_point:])
assert not (train_clusters & test_clusters)

gene_split = {
    gene: ("train" if cluster_of[gene] in train_clusters else "test")
    for gene in embedded_genes
}


def pair_split(gene_a, gene_b):
    side_a, side_b = gene_split[gene_a], gene_split[gene_b]
    return side_a if side_a == side_b else None


labeled_pairs = [(a, b, 1) for a, b in positive_pairs] + [
    (a, b, 0) for a, b in negative_pairs
]

X_train, y_train, X_test, y_test = [], [], [], []
for gene_a, gene_b, label in labeled_pairs:
    side = pair_split(gene_a, gene_b)
    if side is None:
        continue  # pair straddles the train/test boundary -- drop it
    features = np.concatenate([embedding_dic[gene_a], embedding_dic[gene_b]])
    if side == "train":
        X_train.append(features)
        y_train.append(label)
    else:
        X_test.append(features)
        y_test.append(label)

X_train, y_train = np.array(X_train), np.array(y_train)
X_test, y_test = np.array(X_test), np.array(y_test)

print(f"Train pairs: {len(X_train)}, positives: {y_train.sum()}")
print(f"Test pairs:  {len(X_test)}, positives: {y_test.sum()}")

# %%
# Train the candidate models: a linear baseline, two tree-ensemble
# methods, and the feed-forward neural network we ultimately select.
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.neural_network import MLPClassifier
from lightgbm import LGBMClassifier
from xgboost import XGBClassifier

candidate_models = {
    "Logistic Regression": LogisticRegression(max_iter=1000, random_state=42),
    "Random Forest": RandomForestClassifier(
        n_estimators=200, random_state=42, n_jobs=-1
    ),
    "LightGBM": LGBMClassifier(random_state=42),
    "XGBoost": XGBClassifier(
        random_state=42, eval_metric="logloss", n_jobs=-1
    ),
    "Feed-Forward Neural Network": MLPClassifier(
        hidden_layer_sizes=(512, 128),
        max_iter=200,
        early_stopping=True,
        random_state=42,
    ),
}

fitted_models = {}
for model_name, model in candidate_models.items():
    print(f"Training {model_name}...")
    model.fit(X_train, y_train)
    fitted_models[model_name] = model

# %% [markdown]
# # Section 5 - Evaluate the Model
#
# In this section we evaluate each candidate model on the held-out test
# set, compare their performance, and select the feed-forward neural
# network as the final model.

# %%
from sklearn.metrics import (
    accuracy_score,
    precision_score,
    recall_score,
    f1_score,
    roc_auc_score,
    roc_curve,
    precision_recall_curve,
    average_precision_score,
    confusion_matrix,
)

model_metrics = {}
for model_name, model in fitted_models.items():
    y_pred = model.predict(X_test)
    y_proba = model.predict_proba(X_test)[:, 1]

    fpr, tpr, _ = roc_curve(y_test, y_proba)
    pr_precision, pr_recall, _ = precision_recall_curve(y_test, y_proba)
    cm = confusion_matrix(y_test, y_pred)

    model_metrics[model_name] = {
        "accuracy": accuracy_score(y_test, y_pred),
        "precision": precision_score(y_test, y_pred),
        "recall": recall_score(y_test, y_pred),
        "f1": f1_score(y_test, y_pred),
        "auroc": roc_auc_score(y_test, y_proba),
        "fpr": fpr.tolist(),
        "tpr": tpr.tolist(),
        "auprc": average_precision_score(y_test, y_proba),
        "precision_curve": pr_precision.tolist(),
        "recall_curve": pr_recall.tolist(),
        "confusion_matrix": cm.tolist(),
    }

pd.DataFrame(model_metrics).T[
    ["accuracy", "precision", "recall", "f1", "auroc", "auprc"]
]

# %%
with open(PROCESSED_DIR / "model_metrics.json", "w") as f:
    json.dump({"models": model_metrics}, f)

# %%
# Save the selected model (the feed-forward neural network) for the
# Streamlit predictor page.
with open(PROCESSED_DIR / "model.pkl", "wb") as f:
    pickle.dump(fitted_models["Feed-Forward Neural Network"], f)
