#!/usr/bin/env python3
"""Validate the structure of the public cohort registry."""
from pathlib import Path
import csv

REQUIRED = {"accession", "title", "platform", "design", "samples", "role"}

def validate(path: str) -> None:
    rows = list(csv.DictReader(Path(path).open(), delimiter="\t"))
    assert rows, "registry is empty"
    assert REQUIRED.issubset(rows[0]), "missing required columns"
    accessions = [r["accession"] for r in rows]
    assert len(accessions) == len(set(accessions)), "duplicate accession"
    for r in rows:
        assert r["accession"].startswith("GSE")
        assert int(r["samples"]) > 0
        assert r["role"]

if __name__ == "__main__":
    validate("data/public/validation/cohort_registry.tsv")
    print("cohort registry valid")
