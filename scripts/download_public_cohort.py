#!/usr/bin/env python3
"""Download GSE156336 processed/raw counts with provenance and integrity checks."""
from __future__ import annotations
import argparse, hashlib, json, pathlib, urllib.request

URLS = {
    "GSE156336_CPM": "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE156336&file=GSE156336_Feature_Counts_RNA_seq_CPM.txt.gz&format=file",
    "GSE156336_RAW": "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE156336&file=GSE156336_Raw_Feature_Counts_RNA_seq.txt.gz&format=file",
}

def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", choices=sorted(URLS), default="GSE156336_CPM")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()
    out = pathlib.Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    url = URLS[args.dataset]
    try:
        urllib.request.urlretrieve(url, out)
    except Exception as exc:
        raise SystemExit(f"DOWNLOAD_FAILED: {exc}\nURL={url}")
    if out.stat().st_size == 0:
        raise SystemExit("DOWNLOAD_FAILED: empty file")
    meta = {"dataset": args.dataset, "accession": "GSE156336", "source": "NCBI GEO", "url": url,
            "filename": out.name, "bytes": out.stat().st_size, "sha256": sha256(out)}
    meta_path = pathlib.Path(str(out) + ".provenance.json")
    meta_path.write_text(json.dumps(meta, indent=2) + "\n")
    print(json.dumps(meta, indent=2))

if __name__ == "__main__":
    main()
