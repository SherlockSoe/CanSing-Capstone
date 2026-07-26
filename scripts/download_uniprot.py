"""Download the UniProt reference proteome FASTA for S. cerevisiae S288c (UP000002311).

Streams the proteome from the UniProt REST API and saves it to data/raw/.
"""

from pathlib import Path

import requests

PROTEOME_ID = "UP000002311"
STREAM_URL = (
    f"https://rest.uniprot.org/uniprotkb/stream"
    f"?query=proteome:{PROTEOME_ID}&format=fasta&compressed=false"
)

REPO_ROOT = Path(__file__).resolve().parent.parent
OUTPUT_DIR = REPO_ROOT / "data" / "raw"


def download_uniprot_fasta(output_dir: Path = OUTPUT_DIR) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / f"uniprotkb_proteome_{PROTEOME_ID}.fasta"

    print(f"Downloading {STREAM_URL} ...")
    response = requests.get(STREAM_URL, timeout=60)
    response.raise_for_status()

    output_path.write_bytes(response.content)

    print(f"Saved {output_path}")
    return output_path


if __name__ == "__main__":
    download_uniprot_fasta()
