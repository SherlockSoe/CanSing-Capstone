# CanSing Capstone

## Setup

1. Clone the repo and install dependencies:

   ```bash
   git clone <repo-url>
   cd CanSing-Capstone
   pip install -r requirements.txt
   ```

   **macOS only:** `lightgbm` needs the OpenMP runtime, which isn't bundled
   with the pip package. If you hit `Library not loaded: @rpath/libomp.dylib`,
   run:

   ```bash
   brew install libomp
   ```

2. Install the `nbstripout` git filter (**one-time, per clone**). This strips
   notebook cell outputs/execution counts before they're committed, so two
   people editing the notebook don't generate noisy diffs or merge conflicts
   from stale output cells:

   ```bash
   nbstripout --install --attributes .gitattributes
   ```

3. Download the raw data. These files are too large to commit to GitHub, so
   they're gitignored (`data/raw/`) and fetched on demand via the scripts in
   `scripts/`:

   ```bash
   python scripts/download_biogrid.py   # BioGRID S. cerevisiae interactions (~670 MB)
   python scripts/download_uniprot.py   # UniProt S. cerevisiae proteome FASTA (~4 MB)
   ```

   Both scripts save their output into `data/raw/` under fixed filenames, so
   the notebook finds them the same way on every machine.

4. Launch the notebook:

   ```bash
   jupyter notebook notebooks/CapstoneNoteBookMain.ipynb
   ```

   File paths (`path1`, `path2`) resolve automatically to `data/raw/` — no
   machine-specific editing needed, whether you launch Jupyter from the repo
   root or from `notebooks/`.

## Running the GUI

The project includes a Streamlit app (`app/`) that presents the report's
narrative, dataset stats, and (once the model pipeline is built out) results
and an interactive predictor. Run it from the repo root:

```bash
streamlit run app/Home.py
```

- **Home** — project background.
- **Data Overview** — real numbers from the downloaded UniProt/BioGRID data
  (step 3 above); shows a setup reminder if you haven't run the download
  scripts yet.
- **Data Cleaning**, **Model Results**, **Predictor** — populate
  automatically once the notebook writes their backing artifacts
  (`data/processed/cleaning_summary.json`, `model_metrics.json`,
  `model.pkl`); until then they show what's still needed. No app changes
  are required when those become available — see the docstring at the top
  of each page in `app/pages/` for the exact expected file/field names.

Shared data-loading logic (`read_fasta`, `extract_gene_names`,
`load_biogrid_interactions`, `extract_locuslink`, `get_interactors`) lives
in `src/ppi_utils.py` and is imported by both the notebook and the app, so
there's one implementation to keep correct.

## Collaborating on the notebook

- **Outputs are stripped automatically.** `nbstripout` (step 2 above) removes
  cell outputs/execution counts at commit time via a git filter, so `git diff`
  and PR reviews only show real code changes, not re-run noise. Your working
  copy still shows outputs locally — only what gets committed is stripped.
- **Review changes as plain Python.** The notebook is paired with
  `notebooks/CapstoneNoteBookMain.py` via `jupytext` (in percent-cell
  format). Both files are kept in sync automatically whenever you save the
  notebook in Jupyter (or run `jupytext --sync notebooks/CapstoneNoteBookMain.ipynb`
  from the CLI). The `.py` file is what's actually easy to diff and merge on
  GitHub — prefer resolving merge conflicts there, then run `jupytext --sync`
  to bring the `.ipynb` back in line.
- **Avoid parallel edits to the same cells** where possible — `.ipynb` merge
  conflicts are still harder to resolve than plain-text ones even with the
  tooling above.

## Code style

The project's Python code — `scripts/*.py`, `src/*.py`, `app/`, and the
notebook's `jupytext` mirror, `notebooks/CapstoneNoteBookMain.py` — is
checked against PEP-8 with `flake8` (79-char line length, configured in
`.flake8`). Before submitting or opening a PR, verify it's clean:

```bash
flake8 scripts/ src/ app/ notebooks/CapstoneNoteBookMain.py
```

No output means no violations. If you need to fix something, `ruff` (already
in `requirements.txt`, configured in `ruff.toml`) auto-fixes most of it:

```bash
ruff format scripts/ src/ app/ notebooks/CapstoneNoteBookMain.py
ruff check --fix scripts/ src/ app/ notebooks/CapstoneNoteBookMain.py
```

Then run `jupytext --sync notebooks/CapstoneNoteBookMain.ipynb` to carry the
fixes from the `.py` mirror back into the notebook.

## Project structure

```
app/              Streamlit GUI (streamlit run app/Home.py)
data/raw/         Downloaded source data (gitignored, populate via scripts/)
data/processed/   Generated intermediate artifacts, e.g. interactions.pkl (gitignored)
notebooks/        Analysis and modeling notebook (paired .ipynb + .py via jupytext)
scripts/          Data-download scripts
src/              Shared data-loading logic, imported by both the notebook and app/
```
