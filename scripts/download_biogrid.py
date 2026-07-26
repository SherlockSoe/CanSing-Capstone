"""Download the BioGRID S. cerevisiae S288c interaction dataset (organism mitab archive).

Fetches BIOGRID-ORGANISM-5.0.259.mitab.zip from the BioGRID release archive,
extracts BIOGRID-ORGANISM-Saccharomyces_cerevisiae_S288c-5.0.259.mitab.txt,
and saves it to data/raw/.
"""

import io
import zipfile
from pathlib import Path

import requests

BIOGRID_VERSION = "5.0.259"
ZIP_URL = (
    f"https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/"
    f"BIOGRID-{BIOGRID_VERSION}/BIOGRID-ORGANISM-{BIOGRID_VERSION}.mitab.zip"
)
TARGET_MEMBER = f"BIOGRID-ORGANISM-Saccharomyces_cerevisiae_S288c-{BIOGRID_VERSION}.mitab.txt"

REPO_ROOT = Path(__file__).resolve().parent.parent
OUTPUT_DIR = REPO_ROOT / "data" / "raw"


def download_biogrid_mitab(output_dir: Path = OUTPUT_DIR) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / TARGET_MEMBER

    print(f"Downloading {ZIP_URL} ...")
    response = requests.get(ZIP_URL, timeout=60)
    response.raise_for_status()

    with zipfile.ZipFile(io.BytesIO(response.content)) as zf:
        if TARGET_MEMBER not in zf.namelist():
            raise FileNotFoundError(
                f"{TARGET_MEMBER} not found in archive. Available members: {zf.namelist()}"
            )
        with zf.open(TARGET_MEMBER) as src, open(output_path, "wb") as dst:
            dst.write(src.read())

    print(f"Saved {output_path}")
    return output_path


if __name__ == "__main__":
    download_biogrid_mitab()
