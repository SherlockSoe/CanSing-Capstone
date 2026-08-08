"""Landing page for the PPI prediction capstone GUI."""

import streamlit as st

st.set_page_config(page_title="PPI Prediction | Home", page_icon="🧬")

st.title("Predicting Protein-Protein Interactions in S. cerevisiae")

st.markdown(
    """
Use the pages in the sidebar to explore the dataset, the data-cleaning
step, model results, and an interactive predictor.

## What Are Protein-Protein Interactions (PPI)?

A protein-protein interaction is a specific, physical association between
two or more proteins that enables them to perform, regulate, or coordinate
biological functions within a cell. These interactions can be temporary,
such as those that occur during cell signaling or enzyme regulation, or
long-lasting, resulting in the formation of protein complexes such as those
involved in DNA replication, protein synthesis, or other metabolic
activity. PPI are essential for nearly every cellular process.

## Why Predict PPI?

Identifying and studying PPIs can help researchers understand cellular
pathways, identify disease mechanisms, and develop drugs that target or
disrupt specific PPIs. Historically, PPIs have been identified using a
combination of experimental approaches — no single method is sufficient,
since each has different strengths, limitations, and biases. Researchers
often validate an interaction using multiple complementary techniques,
including:

- Yeast two-hybrid (Y2H)
- Affinity purification-mass spectrometry (AP-MS)
- Cross-linking mass spectrometry (XL-MS)
- Fluorescence-based methods
- Structural methods, such as X-ray crystallography, cryogenic electron
  microscopy, and nuclear magnetic resonance spectroscopy

These techniques are expensive and time-consuming, and given the size and
complexity of the human proteome, identifying all PPI using these methods
is not practical. As a result, computational approaches for predicting PPI,
such as the one explored in this project, are becoming more widespread.

## Why *Saccharomyces cerevisiae*?

The organism chosen for this project is *S. cerevisiae*, commonly referred
to as baker's yeast. This single-celled organism is one of the most widely
used model organisms for studying PPI, as it is both:

- **Biologically relevant** — as a eukaryotic cell, many of the metabolic
  pathways in yeast are also found in more complex organisms, including
  humans.
- **Experimentally simple** — given how quickly yeast can reproduce at a
  relatively low cost, and their compatibility with high-throughput assays,
  a large amount of experimental data already exists.
"""
)

st.info(
    "This app tracks the notebook's progress: pages populate with real "
    "data and results as each notebook section is completed."
)
