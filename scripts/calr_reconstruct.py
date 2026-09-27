#!/usr/bin/env python3
"""Reconstruct canonical CALR Type 1/Type 2 coding consequences.

Primary source: NCBI RefSeq NM_004343.4 GenBank.
Fallback: Ensembl MANE Select ENST00000316448.10 CDS when NCBI is unreachable.
The source and SHA-256 checksum are recorded in the output table.
"""
from __future__ import annotations
import argparse
import csv
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS))
from calr_consequence import consequence
from calr_genbank_reference import build_report, fetch_genbank, cds_sequence
from calr_sequence_sources import fetch_ensembl_cds

VARIANTS = [
    ("Type 1", "NM_004343.4:c.1099_1150del", "p.Leu367fs / commonly p.L367fs*46"),
    ("Type 2", "NM_004343.4:c.1154_1155insTTGTC", "p.Lys385fs / commonly p.K385fs*47"),
]


def get_reference(source: str):
    if source in ("auto", "ncbi"):
        try:
            gb_path = ROOT / "data/raw/NM_004343.4.gb"
            text, raw = fetch_genbank(str(gb_path))
            report = build_report(text, raw)
            return cds_sequence(text, report["cds"]), report["provider"] if "provider" in report else "NCBI RefSeq", report["sha256"]
        except Exception:
            if source == "ncbi":
                raise
    seq, meta = fetch_ensembl_cds()
    return seq, meta["provider"], meta["sha256"]


def run(output: str = "results/tables/calr_reconstruction.tsv", source: str = "auto") -> pathlib.Path:
    cds, provider, ref_sha = get_reference(source)
    out = ROOT / output
    out.parent.mkdir(parents=True, exist_ok=True)
    fields = ["class", "hgvs", "expected_protein_label", "reference_cds_length", "altered_cds_length",
              "net_change_bp", "frameshift", "reference_protein_length", "altered_protein_length",
              "reference_provider", "reference_sha256"]
    with out.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for label, hgvs, expected in VARIANTS:
            result = consequence(cds, hgvs)
            writer.writerow({
                "class": label, "hgvs": hgvs, "expected_protein_label": expected,
                "reference_cds_length": result["reference_cds_length"],
                "altered_cds_length": result["altered_cds_length"],
                "net_change_bp": result["net_change_bp"], "frameshift": result["frameshift"],
                "reference_protein_length": len(result["reference_protein"]),
                "altered_protein_length": len(result["altered_protein"]),
                "reference_provider": provider, "reference_sha256": ref_sha,
            })
    return out


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--source", choices=["auto", "ncbi", "ensembl"], default="auto")
    p.add_argument("--output", default="results/tables/calr_reconstruction.tsv")
    args = p.parse_args()
    print(run(args.output, args.source))
