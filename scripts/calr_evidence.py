"""Load the CALR public-evidence register into a deterministic TSV-like structure."""
from pathlib import Path
import csv

ROOT = Path(__file__).resolve().parents[1]
TABLE = ROOT / "results" / "tables" / "evidence_register.tsv"

def load_evidence():
    with TABLE.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))

if __name__ == "__main__":
    for row in load_evidence():
        print(f"{row['id']}: {row['reference_variant_or_topic']} — {row['claim']}")
