# CanSing Capstone

## Setup

1. Clone the repo and install dependencies:

   ```bash
   git clone <repo-url>
   cd CanSing-Capstone
   pip install -r requirements.txt
   ```

2. Download the raw data. These files are too large to commit to GitHub, so
   they're gitignored (`data/raw/`) and fetched on demand via the scripts in
   `scripts/`:

   ```bash
   python scripts/download_biogrid.py   # BioGRID S. cerevisiae interactions (~670 MB)
   python scripts/download_uniprot.py   # UniProt S. cerevisiae proteome FASTA (~4 MB)
   ```

   Both scripts save their output into `data/raw/`.

3. Launch the notebook:

   ```bash
   jupyter notebook "notebooks/Capstone Notebook - Main.ipynb"
   ```

   Note: the notebook's `path1` / `path2` variables currently point at
   `Downloads/...`. Update them to point at the files under `data/raw/`
   before running.

## Project structure

```
data/raw/     Downloaded source data (gitignored, populate via scripts/)
notebooks/    Analysis and modeling notebook
scripts/      Data-download scripts
```
