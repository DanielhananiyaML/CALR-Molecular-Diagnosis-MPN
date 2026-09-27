#!/usr/bin/env python3
"""Reference-sequence source adapters for CALR.

NCBI RefSeq is the primary source. Ensembl MANE Select is an explicit fallback
for environments where NCBI E-utilities are unreachable. The fallback is only
accepted when the returned transcript/CDS is validated for the expected CALR
MANE transcript and CDS length. No sequence is hard-coded here.
"""
from __future__ import annotations

import hashlib
import urllib.request

NCBI_GB = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NM_004343.4&rettype=gb&retmode=text"
ENSEMBL_CDS = "https://rest.ensembl.org/sequence/id/ENST00000316448.10?type=cds"
EXPECTED_TRANSCRIPT = "ENST00000316448.10"
EXPECTED_CDS_LENGTH = 1254


def fetch(url: str, headers: dict[str, str] | None = None, timeout: int = 60) -> bytes:
    req = urllib.request.Request(url, headers=headers or {})
    return urllib.request.urlopen(req, timeout=timeout).read()


def fetch_ncbi() -> tuple[bytes, dict]:
    data = fetch(NCBI_GB, {"User-Agent": "CALR-Molecular-Diagnosis-MPN/0.3"})
    text = data.decode("utf-8", errors="strict")
    if "LOCUS" not in text or "NM_004343.4" not in text:
        raise ValueError("NCBI response is not the expected CALR RefSeq record")
    return data, {"provider": "NCBI RefSeq", "accession": "NM_004343.4", "url": NCBI_GB,
                  "sha256": hashlib.sha256(data).hexdigest()}


def fetch_ensembl_cds() -> tuple[str, dict]:
    data = fetch(ENSEMBL_CDS, {"Content-Type": "text/x-fasta"})
    text = data.decode("utf-8", errors="strict")
    lines = [x.strip() for x in text.splitlines() if x.strip()]
    if not lines or not lines[0].startswith(">" + EXPECTED_TRANSCRIPT):
        raise ValueError("Ensembl response is not the expected CALR MANE transcript")
    seq = "".join(lines[1:]).upper()
    if len(seq) != EXPECTED_CDS_LENGTH or any(b not in "ACGTN" for b in seq):
        raise ValueError(f"Unexpected Ensembl CDS: length={len(seq)}")
    return seq, {"provider": "Ensembl", "transcript": EXPECTED_TRANSCRIPT, "url": ENSEMBL_CDS,
                 "sha256": hashlib.sha256(seq.encode()).hexdigest(), "length": len(seq)}
