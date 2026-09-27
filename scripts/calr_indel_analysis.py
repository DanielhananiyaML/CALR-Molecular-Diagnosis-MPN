#!/usr/bin/env python3
"""Transparent CALR exon 9 indel parser/classifier.

This module intentionally performs sequence-independent HGVS event analysis.
It does not claim clinical pathogenicity and does not replace sequence-level
variant validation.
"""
from __future__ import annotations
import re
from dataclasses import dataclass, asdict

TYPE1_LEGACY = "NM_004343.3:c.1092_1143del52"
TYPE1_MANE = "NM_004343.4:c.1099_1150del"
TYPE2 = "NM_004343.4:c.1154_1155insTTGTC"

@dataclass(frozen=True)
class Indel:
    hgvs: str
    kind: str
    start: int
    end: int
    sequence: str

    @property
    def length(self) -> int:
        return len(self.sequence) if self.kind == "ins" else self.end - self.start + 1

    @property
    def net_change(self) -> int:
        return self.length if self.kind == "ins" else -self.length

    @property
    def frameshift(self) -> bool:
        return self.net_change % 3 != 0


def parse_exon9_hgvs(hgvs: str) -> Indel:
    s = hgvs.strip()
    m = re.fullmatch(r"NM_004343\.(?:3|4):c\.(\d+)(?:_(\d+))?del(\d+)?", s)
    if m:
        start = int(m.group(1)); end = int(m.group(2) or start)
        declared = int(m.group(3) or (end - start + 1))
        if declared != end - start + 1:
            raise ValueError(f"Deletion length mismatch in {s}")
        return Indel(s, "del", start, end, "")
    m = re.fullmatch(r"NM_004343\.(?:3|4):c\.(\d+)(?:_(\d+))?ins([ACGT]+)", s)
    if m:
        start = int(m.group(1)); end = int(m.group(2) or start)
        return Indel(s, "ins", start, end, m.group(3))
    raise ValueError(f"Unsupported CALR HGVS indel: {s}")


def classify(indel: Indel) -> str:
    if indel.hgvs == TYPE1_LEGACY or indel.hgvs == TYPE1_MANE:
        return "Type 1"
    if indel.hgvs == TYPE2:
        return "Type 2"
    if not indel.frameshift:
        return "Non-classical in-frame indel"
    return "Non-classical frameshift — review"


def analyse(hgvs: str) -> dict:
    indel = parse_exon9_hgvs(hgvs)
    cls = classify(indel)
    return {
        "hgvs": indel.hgvs,
        "kind": indel.kind,
        "start": indel.start,
        "end": indel.end,
        "length": indel.length,
        "net_change": indel.net_change,
        "frameshift": indel.frameshift,
        "classification": cls,
        "clinical_interpretation": "Research classification only; sequence and clinical evidence required."
    }

if __name__ == "__main__":
    examples = [TYPE1_LEGACY, TYPE1_MANE, TYPE2, "NM_004343.4:c.1120_1132del"]
    for x in examples:
        print(asdict(parse_exon9_hgvs(x)))
        print(analyse(x))
