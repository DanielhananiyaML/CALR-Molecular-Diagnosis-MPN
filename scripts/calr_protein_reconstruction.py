#!/usr/bin/env python3
"""Protein-level reconstruction of canonical CALR Type 1 and Type 2 consequences.

This layer is deliberately distinct from nucleotide-level HGVS reconstruction.
It starts from the authoritative NP_004334.1 reference protein and applies the
published novel C-terminal peptide consequences for the two canonical CALR exon 9
frameshifts. It is a validation bridge, not a substitute for nucleotide-level
reconstruction from NM_004343.4.
"""
from __future__ import annotations
import argparse, csv, hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REFERENCE = ROOT / "data/raw/NP_004334.1.fasta"
VARIANTS = {
    "Type 1": {
        "hgvs": "NM_004343.4:c.1099_1150del",
        "protein_change": "p.Leu367ThrfsTer46 (p.L367fs*46)",
        "anchor": 366,
        "tail": "TRRMMRTKMRMRRMRRTRRKMRRKMSPARPRTSCREACLQGWTEA",
        "expected_length": 411,
    },
    "Type 2": {
        "hgvs": "NM_004343.4:c.1154_1155insTTGTC",
        "protein_change": "p.Lys385AsnfsTer47 (p.K385fs*47)",
        "anchor": 384,
        "tail": "NCRRMMRTKMRMRRMRRTRRKMRRKMSPARPRTSCREACLQGWTEA",
        "expected_length": 430,
    },
}

def read_fasta(path: Path) -> str:
    return ''.join(line.strip() for line in path.read_text().splitlines() if line.strip() and not line.startswith('>')).upper()

def reconstruct(reference: str, spec: dict) -> str:
    return reference[:spec['anchor']] + spec['tail']

def run(output='results/tables/calr_protein_reconstruction.tsv'):
    reference = read_fasta(REFERENCE)
    ref_sha = hashlib.sha256(REFERENCE.read_bytes()).hexdigest()
    rows=[]
    for label,spec in VARIANTS.items():
        mutant=reconstruct(reference,spec)
        rows.append({
            'class': label,
            'hgvs': spec['hgvs'],
            'protein_change': spec['protein_change'],
            'reference_accession': 'NP_004334.1',
            'reference_length': len(reference),
            'mutant_length': len(mutant),
            'expected_length': spec['expected_length'],
            'length_match': len(mutant)==spec['expected_length'],
            'novel_c_terminal_length': len(spec['tail']),
            'c_terminal_sequence': spec['tail'],
            'reference_fasta_sha256': ref_sha,
        })
    out=ROOT/output
    out.parent.mkdir(parents=True,exist_ok=True)
    with out.open('w',newline='') as fh:
        w=csv.DictWriter(fh,fieldnames=rows[0].keys(),delimiter='\t'); w.writeheader(); w.writerows(rows)
    return out

if __name__=='__main__':
    ap=argparse.ArgumentParser(); ap.add_argument('--output',default='results/tables/calr_protein_reconstruction.tsv'); args=ap.parse_args(); print(run(args.output))
