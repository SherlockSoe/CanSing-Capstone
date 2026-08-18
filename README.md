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

[`CapstoneNotebookTeamCanSing.pdf`](CapstoneNotebookTeamCanSing.pdf), at
the repo root, is a static export of the notebook with every cell already
run — it's the code that generates every result and figure in the report,
so you can read it straight through without installing anything or
waiting on the long-running cells below.

## Long-running cells (`CELL_TAKE_TIME`)

Five cells in the notebook are marked with a `# CELL_TAKE_TIME` comment
because they're expensive to (re-)run. Each writes its result to
`data/processed/`, and a later cell loads that file back instead of
recomputing it:

| `CELL_TAKE_TIME` cell builds... | ...and writes | which is loaded by |
|---|---|---|
| the BioGRID gene-interaction dictionary (`get_interactors` over every gene in `gene_list`) | `data/processed/interactions.pkl` | the next cell |
| ESM2 embeddings for every protein using the full 650M-parameter model (`esm2_t33_650M_UR50D`) | `data/processed/embeddings.pkl` | the next cell |
| ESM2 embeddings for every protein using the smaller 8M-parameter model (`esm2_t6_8M_UR50D`), used for the model-size comparison in the results section | `data/processed/embeddings_small.pkl` | the next cell |
| the feed-forward PPI model on the (scaled) full-size embeddings — `model.fit(...)`, 100 epochs | `data/processed/ppi_esm2_model.keras` | the model-reload cell (`model2 = keras.models.load_model(...)`) |
| the feed-forward PPI model on the small-embedding set — `model_small.fit(...)`, 100 epochs | `data/processed/ppi_esm2_model_small.keras` | the same reload cell (`model_small = keras.models.load_model(...)`) |

The first three (interaction lookup, then the two embedding cells) are
commented out by default — the code itself is inert until you uncomment
it. The two training cells work differently: `model.fit`/`model_small.fit`
run every time the cell executes (100 epochs each, several minutes), and
only the trailing comment tells you to comment them out yourself once
you've got a `.keras` file you're happy with — the notebook doesn't skip
them automatically the way it does for the pickle cells.

All five output files (`interactions.pkl`, `embeddings.pkl`,
`embeddings_small.pkl`, `ppi_esm2_model.keras`,
`ppi_esm2_model_small.keras` — plus `ppi_esm2_model_no_scaling.keras`,
the unscaled counterpart also loaded by the reload cell) are already
committed to `data/processed/`, so a fresh clone doesn't need to run any
of these cells. Only rerun one if you've changed the upstream data, the
extraction logic in `src/ppi_utils.py`, or the model architecture —
delete or rename the corresponding output file first so a stale one
can't get loaded by mistake, and for the training cells remember to
re-comment `model.fit`/`model_small.fit` afterward.

Shared data-loading and modeling logic (`read_fasta`, `extract_gene_names`,
`load_biogrid_interactions`, `extract_locuslink`, `get_interactors`,
`get_esm_embedding`, `get_esm_embeddings_batch`, `get_pairwise_features`)
lives in `src/ppi_utils.py` and is imported by the notebook.

## Collaborating on the notebook

- **Avoid parallel edits to the same cells** where possible — `.ipynb` merge
  conflicts are harder to resolve than plain-text ones, since the file
  format includes cell outputs and execution counts alongside the code.

## Code style

The project's Python code — `scripts/*.py`, `src/*.py` — is
checked against PEP-8 with `flake8` (79-char line length, configured in
`.flake8`). Before submitting or opening a PR, verify it's clean:

```bash
flake8 scripts/ src/
```

No output means no violations. If you need to fix something, `ruff` (already
in `requirements.txt`, configured in `ruff.toml`) auto-fixes most of it:

```bash
ruff format scripts/ src/
ruff check --fix scripts/ src/
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
CapstoneNotebookTeamCanSing.pdf   Static export of the notebook — code, results, and figures used in the report
data/raw/         Downloaded source data (gitignored, populate via scripts/)
data/processed/   Generated intermediate artifacts (interactions.pkl, embeddings*.pkl, etc.) — committed to the repo
notebooks/        Analysis and modeling notebook
scripts/          Data-download scripts
src/              Shared data-loading logic, imported by the notebook
```
