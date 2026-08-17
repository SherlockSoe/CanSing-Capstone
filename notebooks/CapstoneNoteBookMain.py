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
#     display_name: Python 3 (ipykernel)
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
# The FASTA data format represents each amino acid with a single
# character as follows:
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
# The FASTA format also uses the following abbreviations for ambiguous
# amino acid identification and special characters:
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
import random
import sys

# File read / write
from pathlib import Path
import json
import pickle

# Machine Learning / Modeling
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
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

# Deep Learning
from tensorflow import keras
from tensorflow.keras import layers

# BioPython
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
# New installations required for this project

# #%pip install biopython
# pip install torch
# pip install fair-esm
# pip install tensorflow tensorflow-metal  # tensorflow-metal is macOS-only

# %% [markdown]
# # Helper Functions
#
# In this section we're defining some helper functions for use in the rest
# of the notebook.

# %%
# Shared with app/ (Streamlit GUI) via src/ppi_utils.py, so the data
# loading and embedding logic has a single source of truth instead of
# being duplicated.
from src.ppi_utils import (  # noqa: E402
    read_fasta,
    extract_locuslink,
    extract_gene_names,
    get_interactors,
    load_biogrid_interactions,
    get_esm_embedding,
    get_esm_embeddings_batch,
    get_pairwise_features,
)


# %%
# Function to generate a random amino acid sequence of fixed length, used
# below as an embedding sanity check (a random sequence should look
# "unrelated" to everything else).
def random_protein(length):
    return "".join(random.choices("ACDEFGHIKLMNPQRSTVWY", k=length))


# %% [markdown]
# # Section 1 - Data Import
#
# In this section we are importing amino acid sequences for all 6,067
# amino acids in S. cerevisiae, and interaction data from BioGrid.

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
sequence_dictionary = extract_gene_names(fasta_data)
length_dictionary = {
    gene: len(seq) for gene, seq in sequence_dictionary.items()
}

# Create a list of all proteins
gene_list = list(sequence_dictionary.keys())

# %%
# Inspect the results
len(gene_list)

# %%
# Check if a particular gene is included in the Uniprot data
if "CDC73" in gene_list:
    print("Item found!")

# %% [markdown]
# Here we are importing the interaction data from BioGrid.

# %%
# Now extract the Biogrid data and create the Biogrid dataframe
bg = load_biogrid_interactions(path2)

print(bg.shape)

# %%
bg.head()

# %%
# Check the Biogrid data at a particular dataframe location
print(bg.at[1, "Alt IDs Interactor A"])

# %%
# Extract the interactor A and B genes and create a new column for each
bg["Gene_A"] = bg["Alt IDs Interactor A"].apply(extract_locuslink)
bg["Gene_B"] = bg["Alt IDs Interactor B"].apply(extract_locuslink)

bg.head()

# %% [markdown]
# This cell takes quite a while to run locally; therefore writing the
# results to file for future use.

# %%
# Now we make the interaction dictionary and write to file for future use
ia_dic = {}

# Iterate through each gene and extract interactions
for gene in gene_list:
    interactors = get_interactors(bg, gene)
    ia_dic[gene] = interactors

# Write the interaction dictionary to file
with open(PROCESSED_DIR / "interactions.pkl", "wb") as file:
    pickle.dump(ia_dic, file)

# %%
# Read the interaction dictionary from the binary file
with open(PROCESSED_DIR / "interactions.pkl", "rb") as file:
    interaction_dictionary = pickle.load(file)

# %%
# Inspect the results
interaction_dictionary

# %% [markdown]
# # Section 2 - Data Cleaning & Pre-processing
#
# In this section we perform some data cleaning and pre-processing
# operations. First, we need to make sure we have sequence and interaction
# data for all proteins.

# %%
# Find entries in the interaction dictionary that are empty
empty_keys = [
    key for key, value in interaction_dictionary.items() if not value
]
print("Number of empty entries in the interaction dictionary:")
print(len(empty_keys))

# Remove them from the sequence and interaction dictionaries
for key in empty_keys:
    sequence_dictionary.pop(key, None)
    length_dictionary.pop(key, None)
    interaction_dictionary.pop(key, None)

print("")
print("Number of proteins for which we have sequence AND interaction data:")
print(len(sequence_dictionary))

new_gene_list = list(sequence_dictionary.keys())

# %%
# Write a summary for the Streamlit app's Data Cleaning page: all proteins
# -> has interaction data / no interaction data.
cleaning_summary = {
    "labels": [
        "All proteins",
        "Has interaction data",
        "No interaction data",
    ],
    "source": [0, 0],
    "target": [1, 2],
    "value": [len(new_gene_list), len(empty_keys)],
}

with open(PROCESSED_DIR / "cleaning_summary.json", "w") as f:
    json.dump(cleaning_summary, f)

# %%
# Next we need to remove any interactors for which we don't have sequence
# data
interaction_dictionary_filtered = {
    protein: [p for p in partners if p in new_gene_list]
    for protein, partners in interaction_dictionary.items()
}

# %%
interaction_dictionary_filtered

# %%
# Finally we create the overall interaction matrix using the interaction
# dictionary - a square matrix with 1's indicating interactions, 0's
# indicating no interaction. Note that this is a symmetric matrix.

# We will call this the 'Master Interaction Matrix' or MIM
MIM = pd.DataFrame(0, index=new_gene_list, columns=new_gene_list)

for protein, partners in interaction_dictionary_filtered.items():
    MIM.loc[protein, partners] = 1

# %%
# Perform some checks to make sure MIM is setup correctly, check sparsity
row = "CET1"
column = "STE2"

print(f"MIM entry at {row}, {column}:")
print(MIM.loc[row, column])
print("")

# See how sparse MIM is
count_ones = np.count_nonzero(MIM)
print("Number of non-zero entries in the MIM:")
print(count_ones)
print("")

print("Total number of entries in the MIM:")
num_entries = len(new_gene_list) ** 2
print(num_entries)
print("")

print("Sparsity of the MIM:")
print(count_ones / num_entries)

# %% [markdown]
# # Section 3 - Create Embeddings
#
# In this section we are creating the embeddings for all the proteins for
# which we have sequence AND interaction data. We will then perform some
# checks to make sure the embeddings are behaving as expected, by
# comparing the cosine similarity between various embeddings.

# %%
# Perform some basic checks of the embedding

# Define a sequence and confirm the embedding has the right dimensionality
seq1 = (
    "MSDAAPSLSNLFYDPTYNPGQSTINYTSIYGNGSTITFDELQGLVNSTVTQAIMFGVRCGAAALTLIVM"
    "WMTSRSRKTPIFIINQVSLFLIILHSALYFKYLLSNYSSVTYALTGFPQFISRGDVHVYGATNIIQVL"
    "LVASIETSLVFQIKVIFTGDNFKRIGLMLTSISFTLGIATVTMYFVSAVKGMIVTYNDVSATQDKYFN"
    "ASTILLASSINFMSFVLVVKLILAIRSRRFLGLKQFDSFHILLIMSCQSLLVPSIIFILAYSLKPNQG"
    "TDVLTTVATLLAVLSLPLSSMWATAANNASKTNTITSDFTTSTDRFYPGTLSSFQTDSINNDAKSSLR"
    "SRLYDLYPRRKETTSDKHSERTFVSETADDIEKNQFYQLPTPTSSKNTRIGPFADASYKEGEVEPVDM"
    "YTPDTAADEEARKFWTEDNNNL"
)
emb1 = get_esm_embedding(seq1)

# Confirm shape, content
print("Length of sequence embedding:")
print(len(emb1))
print("")
print("Values in embedding:")
print(emb1)

# %%
# Perform some more advanced checks

# Define two similar proteins - human alpha globin (HAG) and human beta
# globin (HBG)
HAG = (
    "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALT"
    "NAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTV"
    "LTSKYR"
)
HBG = (
    "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKV"
    "LGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVA"
    "GVANALAHKYH"
)

# Pick a random S. cerivisiae protein
CDC73 = sequence_dictionary["CDC73"]

# Create 2 random proteins of similar length
rand1 = random_protein(150)
rand2 = random_protein(150)

# Check the cosine similarity between the two
alpha = get_esm_embedding(HAG)
beta = get_esm_embedding(HBG)
delta = get_esm_embedding(CDC73)
gamma1 = get_esm_embedding(rand1)
gamma2 = get_esm_embedding(rand2)

# %%
# Calculate the cosine similarity between HAG and HBG
print(
    "Cosine similarity between human alpha globin (HAG) and human beta "
    "globin (HBG) embeddings:"
)
a = alpha
b = beta
similarity = np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b))
print(similarity)
print("")

# Calculate the cosine similarity between HAG and a protein sequence from
# S. cerivisiae
print("Cosine similarity between human alpha globin (HAG) and CDC73:")
a = alpha
b = delta
similarity = np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b))
print(similarity)
print("")

# Calculate the cosine similarity between HAG and a random protein
# sequence
print(
    "Cosine similarity between human alpha globin (HAG) and a random "
    "protein sequence:"
)
a = alpha
b = gamma1
similarity = np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b))
print(similarity)
print("")

# Calculate the cosine similarity between HAG and the other random
# protein sequence
print(
    "Cosine similarity between human alpha globin (HAG) and another "
    "random protein sequence:"
)
a = alpha
b = gamma2
similarity = np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b))
print(similarity)

# %% [markdown]
# This cell takes quite a while to run locally; therefore writing the
# results to file for future use.

# %%
# Create a dictionary of embeddings for all proteins, and write to file
EMBEDDINGS_PATH = PROCESSED_DIR / "embeddings.pkl"

if EMBEDDINGS_PATH.exists():
    with open(EMBEDDINGS_PATH, "rb") as file:
        embedding_dict = pickle.load(file)
else:
    embedded_gene_order = list(new_gene_list)
    vectors = get_esm_embeddings_batch(
        [sequence_dictionary[gene] for gene in embedded_gene_order],
        batch_size=8,
    )
    embedding_dict = dict(zip(embedded_gene_order, vectors))
    with open(EMBEDDINGS_PATH, "wb") as file:
        pickle.dump(embedding_dict, file)

# %%
embedding_dict

# %% [markdown]
# # Section 4 - Train the Model
#
# In this section we begin training our feed-forward neural network
# (a multi-layer perceptron).

# %%
# First we need to create a test/train split
train_proteins, test_proteins = train_test_split(
    new_gene_list,
    test_size=0.20,  # Train on 80% of the data, test on the other 20%
    random_state=42,
)

# Now form test and train interaction matrices
MIM_train = MIM.loc[train_proteins, train_proteins]
MIM_test = MIM.loc[test_proteins, test_proteins]


# %%
def create_balanced_dataset(
    interaction_matrix, embeddings, n=1000, random_state=42
):
    rng = np.random.default_rng(random_state)

    proteins = interaction_matrix.index.to_numpy()

    # Get upper-triangle indices, excluding diagonal
    i, j = np.triu_indices(len(proteins), k=1)

    # Get corresponding interaction values
    values = interaction_matrix.values[i, j]

    # Separate positive and negative pairs
    positive_idx = np.where(values == 1)[0]
    negative_idx = np.where(values == 0)[0]

    # Check that enough pairs exist
    if len(positive_idx) < n:
        raise ValueError(
            f"Only {len(positive_idx)} positive pairs available; "
            f"cannot sample {n}."
        )

    if len(negative_idx) < n:
        raise ValueError(
            f"Only {len(negative_idx)} negative pairs available; "
            f"cannot sample {n}."
        )

    # Randomly select n of each
    positive_sample = rng.choice(positive_idx, size=n, replace=False)
    negative_sample = rng.choice(negative_idx, size=n, replace=False)

    selected = np.concatenate([positive_sample, negative_sample])

    # Shuffle positive and negative examples together
    rng.shuffle(selected)

    X = []
    y = []
    pairs = []

    for idx in selected:
        protein_a = proteins[i[idx]]
        protein_b = proteins[j[idx]]

        emb_a = embeddings[protein_a]
        emb_b = embeddings[protein_b]

        # Symmetric pair representation
        pair_embedding = get_pairwise_features(emb_a, emb_b)

        X.append(pair_embedding)
        y.append(values[idx])
        pairs.append((protein_a, protein_b))

    return np.array(X), np.array(y), pairs


# %%
# Create test and train dataframes
X_train, y_train, train_pairs = create_balanced_dataset(
    MIM_train,
    embedding_dict,
    n=100000,  # Positive+negative interactions pulled from the train MIM
    random_state=42,
)

X_test, y_test, test_pairs = create_balanced_dataset(
    MIM_test,
    embedding_dict,
    n=10000,  # Positive+negative interactions pulled from the test MIM
    random_state=42,
)

# Scale the train and test sets
scaler = StandardScaler()

X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)

# %%
# Define a simple neural network to start with
model = keras.Sequential(
    [
        layers.Input(shape=(X_train.shape[1],)),
        layers.Dense(512, activation="relu"),
        layers.Dropout(0.3),
        layers.Dense(256, activation="relu"),
        layers.Dropout(0.3),
        layers.Dense(64, activation="relu"),
        layers.Dense(1, activation="sigmoid"),
    ]
)

model.compile(
    optimizer=keras.optimizers.Adam(learning_rate=0.001),
    loss="binary_crossentropy",
    metrics=[
        "accuracy",
        keras.metrics.AUC(name="auc"),
        keras.metrics.Precision(name="precision"),
        keras.metrics.Recall(name="recall"),
    ],
)

model.summary()

# %%
# Train the simple neural network
early_stopping = keras.callbacks.EarlyStopping(
    monitor="val_auc",
    mode="max",
    patience=10,
    restore_best_weights=True,
)

history = model.fit(
    X_train,
    y_train,
    validation_split=0.2,
    epochs=100,
    batch_size=32,
    callbacks=[early_stopping],
    verbose=1,
)

# %%
# Save it for future use
model.save(PROCESSED_DIR / "model.keras")

# %%
# Re-load the previously trained model (as opposed to re-training)
model = keras.models.load_model(PROCESSED_DIR / "model.keras")

# %% [markdown]
# # Section 5 - Evaluate the Model
#
# In this section we evaluate the trained model on the held-out test set.

# %%
# Evaluate the model
results = model.evaluate(X_test, y_test, verbose=0)

for name, value in zip(model.metrics_names, results):
    print(f"{name}: {value:.4f}")

results = model.evaluate(X_test, y_test, verbose=0, return_dict=True)

print(results)

# %%
y_prob = model.predict(X_test).ravel()
y_pred = (y_prob >= 0.5).astype(int)

print("ROC-AUC:", roc_auc_score(y_test, y_prob))
print("PR-AUC :", average_precision_score(y_test, y_prob))

# %%
# Write metrics to file for the Streamlit app's Model Results page
fpr, tpr, _ = roc_curve(y_test, y_prob)
pr_precision, pr_recall, _ = precision_recall_curve(y_test, y_prob)
cm = confusion_matrix(y_test, y_pred)

model_metrics = {
    "Feed-Forward Neural Network": {
        "accuracy": accuracy_score(y_test, y_pred),
        "precision": precision_score(y_test, y_pred),
        "recall": recall_score(y_test, y_pred),
        "f1": f1_score(y_test, y_pred),
        "auroc": roc_auc_score(y_test, y_prob),
        "fpr": fpr.tolist(),
        "tpr": tpr.tolist(),
        "auprc": average_precision_score(y_test, y_prob),
        "precision_curve": pr_precision.tolist(),
        "recall_curve": pr_recall.tolist(),
        "confusion_matrix": cm.tolist(),
    }
}

with open(PROCESSED_DIR / "model_metrics.json", "w") as f:
    json.dump({"models": model_metrics}, f)
