#!/usr/bin/env python3
"""Retrieve and validate the CALR MANE Select reference sequence from NCBI.

The repository intentionally does not vendor the external RefSeq sequence.
This script records provenance and writes the downloaded FASTA locally.
"""
from __future__ import annotations
import argparse, hashlib, pathlib, urllib.request

ACCESSION = "NM_004343.4"
NCBI_EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

def fetch_reference(output: str) -> dict:
    params = f"?db=nuccore&id={ACCESSION}&rettype=fasta&retmode=text"
    url = NCBI_EFETCH + params
    data = urllib.request.urlopen(url, timeout=60).read()
    text = data.decode("utf-8")
    if not text.startswith(">") or ACCESSION not in text.splitlines()[0]:
        raise ValueError("NCBI response is not the expected CALR FASTA")
    path = pathlib.Path(output)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    sha256 = hashlib.sha256(data).hexdigest()
    sequence = "".join(x.strip() for x in text.splitlines()[1:])
    return {"accession": ACCESSION, "length": len(sequence), "sha256": sha256, "url": url, "output": str(path)}

def main():
    p = argparse.ArgumentParser()
    p.add_argument("-o", "--output", default="data/raw/NM_004343.4.fasta")
    args = p.parse_args()
    print(fetch_reference(args.output))

if __name__ == "__main__":
    main()
