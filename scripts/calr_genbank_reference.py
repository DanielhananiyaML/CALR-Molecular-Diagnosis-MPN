#!/usr/bin/env python3
"""Fetch and parse the CALR MANE Select GenBank record.

The script deliberately retrieves the reference at runtime. It extracts the
CDS and exon feature coordinates from the GenBank flat file and records a
SHA-256 checksum for provenance. It does not perform clinical interpretation.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import pathlib
import re
import urllib.request
from typing import Iterable

ACCESSION = "NM_004343.4"
URL = (
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    "?db=nuccore&id=NM_004343.4&rettype=gb&retmode=text"
)


def fetch_genbank(output: str) -> tuple[str, bytes]:
    request = urllib.request.Request(
        URL,
        headers={"User-Agent": "CALR-Molecular-Diagnosis-MPN/0.3 (+https://github.com/DanielhananiyaML)"},
    )
    data = urllib.request.urlopen(request, timeout=60).read()
    text = data.decode("utf-8")
    if "LOCUS" not in text or ACCESSION not in text:
        raise ValueError("NCBI response is not the expected CALR GenBank record")
    path = pathlib.Path(output)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(data)
    return text, data


def _feature_blocks(text: str) -> Iterable[tuple[str, str, list[str]]]:
    in_features = False
    current_key = None
    current_loc = None
    qualifiers: list[str] = []
    for line in text.splitlines():
        if line.startswith("FEATURES"):
            in_features = True
            continue
        if in_features and line.startswith("ORIGIN"):
            if current_key:
                yield current_key, current_loc or "", qualifiers
            return
        if not in_features:
            continue
        m = re.match(r"^\s{5}(\S+)\s+(.+)$", line)
        if m:
            if current_key:
                yield current_key, current_loc or "", qualifiers
            current_key, current_loc, qualifiers = m.group(1), m.group(2).strip(), []
            continue
        q = re.match(r'^\s{21}(/.+)$', line)
        if q and current_key:
            qualifiers.append(q.group(1))


def parse_location(location: str) -> dict:
    clean = location.replace(" ", "")
    # This parser records simple contiguous locations; compound/join locations
    # are retained verbatim and flagged rather than silently mis-parsed.
    m = re.fullmatch(r"<?(\d+)\.\.>?(\d+)", clean)
    if not m:
        return {"raw": location, "start": None, "end": None, "compound": True}
    return {"raw": location, "start": int(m.group(1)), "end": int(m.group(2)), "compound": False}


def _qualifier_value(qualifiers: list[str], name: str) -> str | None:
    prefix = f'/{name}='
    for q in qualifiers:
        if q.startswith(prefix):
            return q[len(prefix):].strip('"')
    return None


def parse_record(text: str) -> dict:
    features = []
    for key, loc, qualifiers in _feature_blocks(text):
        parsed = parse_location(loc)
        features.append({
            "type": key,
            "location": parsed,
            "gene": _qualifier_value(qualifiers, "gene"),
            "product": _qualifier_value(qualifiers, "product"),
            "number": _qualifier_value(qualifiers, "number"),
            "note": _qualifier_value(qualifiers, "note"),
        })

    cds = next((f for f in features if f["type"] == "CDS" and f["location"]["start"] is not None), None)
    exons = [f for f in features if f["type"] == "exon"]
    return {"accession": ACCESSION, "cds": cds, "exons": exons, "features": features}


def extract_origin(text: str) -> str:
    if "ORIGIN" not in text:
        raise ValueError("GenBank ORIGIN section not found")
    body = text.split("ORIGIN", 1)[1].split("//", 1)[0]
    seq = re.sub(r"[^A-Za-z]", "", body).upper()
    if not seq or not re.fullmatch(r"[ACGTN]+", seq):
        raise ValueError("Unexpected bases in GenBank ORIGIN")
    return seq


def cds_sequence(text: str, cds: dict) -> str:
    loc = cds["location"]
    if loc["compound"] or loc["start"] is None:
        raise ValueError("CDS location is compound/unsupported; manual review required")
    seq = extract_origin(text)
    return seq[loc["start"] - 1 : loc["end"]]


def build_report(text: str, raw: bytes) -> dict:
    parsed = parse_record(text)
    cds_seq = cds_sequence(text, parsed["cds"])
    return {
        "accession": ACCESSION,
        "sha256": hashlib.sha256(raw).hexdigest(),
        "source_url": URL,
        "cds": parsed["cds"],
        "exons": parsed["exons"],
        "cds_length": len(cds_seq),
        "cds_sha256": hashlib.sha256(cds_seq.encode()).hexdigest(),
        "canonical_events": {
            "type_1": "NM_004343.4:c.1099_1150del",
            "type_2": "NM_004343.4:c.1154_1155insTTGTC",
        },
    }


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--genbank", default="data/raw/NM_004343.4.gb")
    p.add_argument("--report", default="results/tables/calr_genbank_reference.json")
    args = p.parse_args()
    text, raw = fetch_genbank(args.genbank)
    report = build_report(text, raw)
    out = pathlib.Path(args.report)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
