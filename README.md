# CanSing Capstone

## Setup

1. Clone the repo and install dependencies:

   ```bash
   git clone <repo-url>
   cd CanSing-Capstone
   pip install -r requirements.txt
   ```

   `tensorflow` is pinned to `2.16.2` in `requirements.txt` — newer 2.20.x
   crashes on import in this environment (`mutex lock failed`), and
   `tensorflow-metal` (Apple GPU acceleration) only pairs reliably with
   this older release anyway. Don't bump it without checking both.

2. Download the raw data. These files are too large to commit to GitHub, so
   they're gitignored (`data/raw/`) and fetched on demand via the scripts in
   `scripts/`:

   ```bash
   python scripts/download_biogrid.py   # BioGRID S. cerevisiae interactions (~670 MB)
   python scripts/download_uniprot.py   # UniProt S. cerevisiae proteome FASTA (~4 MB)
   ```

   Both scripts save their output into `data/raw/` under fixed filenames, so
   the notebook finds them the same way on every machine.

3. Launch the notebook:

   ```bash
   jupyter notebook notebooks/CapstoneNoteBookMain.ipynb
   ```

   File paths (`path1`, `path2`) resolve automatically to `data/raw/` — no
   machine-specific editing needed, whether you launch Jupyter from the repo
   root or from `notebooks/`.

   **Expect Section 3 (Create Embeddings) to take a while the first time**:
   it downloads the ESM-2 650M-parameter model (~2.6 GB, one-time, cached
   by `fair-esm` under `~/.cache/torch/hub/checkpoints/` afterward) and
   then embeds every protein — roughly 30-60 minutes on a laptop GPU
   (Apple Silicon MPS or CUDA); much longer on CPU only. The resulting
   embeddings are cached to `data/processed/embeddings.pkl`, so re-running
   the notebook after that point is fast. Section 4 (Train the Model)
   trains the feed-forward network and typically takes a few more minutes.

## Long-running cells (`CELL_TAKE_TIME`)

Three cells in the notebook are marked with a `# CELL_TAKE_TIME` comment
because they're expensive to (re-)run — the interaction lookup takes a
few minutes, the embedding cells 30-60+ minutes each depending on
hardware (see step 3 above). Each writes its result to a pickle in
`data/processed/`, and the very next cell loads that pickle back with
`pickle.load` instead of recomputing:

| `CELL_TAKE_TIME` cell builds... | ...and writes | which is loaded by |
|---|---|---|
| the BioGRID gene-interaction dictionary (`get_interactors` over every gene in `gene_list`) | `data/processed/interactions.pkl` | the next cell |
| ESM2 embeddings for every protein using the full 650M-parameter model (`esm2_t33_650M_UR50D`) | `data/processed/embeddings.pkl` | the next cell |
| ESM2 embeddings for every protein using the smaller 8M-parameter model (`esm2_t6_8M_UR50D`), used for the model-size comparison in the results section | `data/processed/embeddings_small.pkl` | the next cell |

All three `.pkl` files are already committed to the repo, so on a fresh
clone none of these cells need to run — they're commented out by default
and the notebook just loads the existing files. Only uncomment and re-run
one if you've changed the upstream data (new BioGRID/UniProt download) or
the extraction logic (`get_interactors`/`get_esm_embedding` in
`src/ppi_utils.py`) and need to regenerate it — delete or rename the
corresponding `.pkl` first so a stale file can't get loaded by mistake.

## Running the GUI

The project includes a Streamlit app (`app/`) that presents the report's
narrative, dataset stats, data-cleaning Sankey diagram, model comparison
results, and an interactive predictor. Run it from the repo root:

```bash
streamlit run app/Home.py
```

- **Home** — project background.
- **Data Overview** — real numbers from the downloaded UniProt/BioGRID data
  (step 3 above); shows a setup reminder if you haven't run the download
  scripts yet.
- **Data Cleaning**, **Model Results**, **Predictor** — populate
  automatically once the notebook writes their backing artifacts
  (`data/processed/cleaning_summary.json` from Section 2,
  `data/processed/model.keras` from Section 4, and
  `data/processed/model_metrics.json` from Section 5); until then they
  show what's still needed. No app changes are required once those exist
  — see the docstring at the top of each page in `app/pages/` for the
  exact expected file/field names.

Shared data-loading and modeling logic (`read_fasta`, `extract_gene_names`,
`load_biogrid_interactions`, `extract_locuslink`, `get_interactors`,
`get_esm_embedding`, `get_esm_embeddings_batch`, `get_pairwise_features`)
lives in `src/ppi_utils.py` and is imported by both the notebook and the
app, so there's one implementation to keep correct — in particular, the
Predictor page builds its feature vector with the exact same
`get_esm_embedding`/`get_pairwise_features` functions and pairing scheme
(absolute difference + element-wise product) used to train the model in
the notebook.

## Collaborating on the notebook

- **Avoid parallel edits to the same cells** where possible — `.ipynb` merge
  conflicts are harder to resolve than plain-text ones, since the file
  format includes cell outputs and execution counts alongside the code.

## Code style

The project's Python code — `scripts/*.py`, `src/*.py`, `app/` — is
checked against PEP-8 with `flake8` (79-char line length, configured in
`.flake8`). Before submitting or opening a PR, verify it's clean:

```bash
flake8 scripts/ src/ app/
```

No output means no violations. If you need to fix something, `ruff` (already
in `requirements.txt`, configured in `ruff.toml`) auto-fixes most of it:

```bash
ruff format scripts/ src/ app/
ruff check --fix scripts/ src/ app/
```

## Data access and licensing

Raw data is not committed to this repo (see step 3 above) — it's fetched
on demand from the original sources:

- **UniProt** (`scripts/download_uniprot.py`) — S. cerevisiae S288c
  reference proteome (`UP000002311`), streamed from the UniProt REST API.
  Owned/maintained by the UniProt Consortium (EMBL-EBI, SIB, PIR) and
  distributed under the [Creative Commons Attribution 4.0 International
  (CC BY 4.0) License](https://www.uniprot.org/help/license) — reuse is
  permitted with attribution.
- **BioGRID** (`scripts/download_biogrid.py`) — S. cerevisiae S288c
  protein-protein interaction data, release `5.0.259`, from the BioGRID
  release archive. Copyright Mike Tyers / TyersLab.com, distributed under
  the [MIT License](https://wiki.thebiogrid.org/doku.php/terms_and_conditions)
  — the copyright and permission notice must be retained in any copy or
  substantial portion of the data.

Both sources' terms should be reviewed directly before any redistribution
beyond this project.

## Attributions

- ESM2 model loading and embedding extraction in `src/ppi_utils.py`
  (`_load_esm_model`, `get_esm_embedding`, `get_esm_embeddings_batch`,
  `_mean_pool_residues`) was adapted from the [ESM2-Tutorial
  repository](https://github.com/ProteinVision/ESM2-Tutorial/blob/main/ESM2.ipynb).

## Project structure

```
app/              Streamlit GUI (streamlit run app/Home.py)
data/raw/         Downloaded source data (gitignored, populate via scripts/)
data/processed/   Generated intermediate artifacts (interactions.pkl, embeddings*.pkl, etc.) — committed to the repo
notebooks/        Analysis and modeling notebook
scripts/          Data-download scripts
src/              Shared data-loading logic, imported by both the notebook and app/
```
