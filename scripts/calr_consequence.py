#!/usr/bin/env python3
"""Sequence-aware CALR exon 9 indel consequence engine.

This engine operates only when a trusted reference sequence is supplied.
It validates the requested nucleotide event, calculates frame status, and
translates the altered coding sequence. It does not perform clinical variant
classification.
"""
from __future__ import annotations
import re
from dataclasses import dataclass

CODON_TABLE = {
'ATA':'I','ATC':'I','ATT':'I','ATG':'M','ACA':'T','ACC':'T','ACG':'T','ACT':'T','AAC':'N','AAT':'N','AAA':'K','AAG':'K',
'AGC':'S','AGT':'S','AGA':'R','AGG':'R','CTA':'L','CTC':'L','CTG':'L','CTT':'L','CCA':'P','CCC':'P','CCG':'P','CCT':'P',
'CAC':'H','CAT':'H','CAA':'Q','CAG':'Q','CGA':'R','CGC':'R','CGG':'R','CGT':'R','GTA':'V','GTC':'V','GTG':'V','GTT':'V',
'GCA':'A','GCC':'A','GCG':'A','GCT':'A','GAC':'D','GAT':'D','GAA':'E','GAG':'E','GGA':'G','GGC':'G','GGG':'G','GGT':'G',
'TCA':'S','TCC':'S','TCG':'S','TCT':'S','TTC':'F','TTT':'F','TTA':'L','TTG':'L','TAC':'Y','TAT':'Y','TAA':'*','TAG':'*',
'TGC':'C','TGT':'C','TGA':'*','TGG':'W'}

@dataclass(frozen=True)
class Event:
    kind: str
    start: int
    end: int
    sequence: str = ""

def parse(hgvs: str) -> Event:
    m = re.fullmatch(r"NM_004343\.4:c\.(\d+)(?:_(\d+))?del(\d+)?", hgvs)
    if m:
        s,e,n = int(m.group(1)), int(m.group(2) or m.group(1)), m.group(3)
        length = e-s+1
        if n and int(n) != length: raise ValueError('Deletion length mismatch')
        return Event('del',s,e,'')
    m = re.fullmatch(r"NM_004343\.4:c\.(\d+)(?:_(\d+))?ins([ACGT]+)", hgvs)
    if m:
        return Event('ins',int(m.group(1)),int(m.group(2) or m.group(1)),m.group(3))
    raise ValueError(f'Unsupported HGVS: {hgvs}')

def apply_coding_event(cds: str, event: Event) -> str:
    if not re.fullmatch('[ACGTacgt]+', cds): raise ValueError('CDS must contain DNA bases only')
    s = cds.upper()
    if event.kind == 'del':
        if not 1 <= event.start <= event.end <= len(s): raise IndexError('Deletion outside CDS')
        return s[:event.start-1] + s[event.end:]
    # HGVS c. insertion is between c.start and c.end; for c.1154_1155ins this is after 1154.
    if not 1 <= event.start <= event.end <= len(s): raise IndexError('Insertion boundary outside CDS')
    if event.end != event.start + 1:
        raise ValueError('This implementation expects adjacent insertion boundaries')
    return s[:event.start] + event.sequence.upper() + s[event.start:]

def translate(cds: str) -> str:
    usable = len(cds) - (len(cds) % 3)
    return ''.join(CODON_TABLE.get(cds[i:i+3], 'X') for i in range(0, usable, 3))

def consequence(cds: str, hgvs: str) -> dict:
    event = parse(hgvs)
    altered = apply_coding_event(cds, event)
    return {
        'hgvs': hgvs,
        'reference_cds_length': len(cds),
        'altered_cds_length': len(altered),
        'net_change_bp': len(altered)-len(cds),
        'frameshift': (len(altered)-len(cds)) % 3 != 0,
        'reference_protein': translate(cds),
        'altered_protein': translate(altered),
    }
