# CanSing Capstone

## Setup

1. Clone the repo and install dependencies:

   ```bash
   git clone <repo-url>
   cd CanSing-Capstone
   pip install -r requirements.txt
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
   jupyter notebook "notebooks/Capstone Notebook - Main.ipynb"
   ```

   File paths (`path1`, `path2`) resolve automatically to `data/raw/` — no
   machine-specific editing needed, whether you launch Jupyter from the repo
   root or from `notebooks/`.

## Collaborating on the notebook

- **Outputs are stripped automatically.** `nbstripout` (step 2 above) removes
  cell outputs/execution counts at commit time via a git filter, so `git diff`
  and PR reviews only show real code changes, not re-run noise. Your working
  copy still shows outputs locally — only what gets committed is stripped.
- **Review changes as plain Python.** The notebook is paired with
  `notebooks/Capstone Notebook - Main.py` via `jupytext` (in percent-cell
  format). Both files are kept in sync automatically whenever you save the
  notebook in Jupyter (or run `jupytext --sync "notebooks/Capstone Notebook - Main.ipynb"`
  from the CLI). The `.py` file is what's actually easy to diff and merge on
  GitHub — prefer resolving merge conflicts there, then run `jupytext --sync`
  to bring the `.ipynb` back in line.
- **Avoid parallel edits to the same cells** where possible — `.ipynb` merge
  conflicts are still harder to resolve than plain-text ones even with the
  tooling above.

## Project structure

```
data/raw/         Downloaded source data (gitignored, populate via scripts/)
data/processed/   Generated intermediate artifacts, e.g. interactions.pkl (gitignored)
notebooks/        Analysis and modeling notebook (paired .ipynb + .py via jupytext)
scripts/          Data-download scripts
```
