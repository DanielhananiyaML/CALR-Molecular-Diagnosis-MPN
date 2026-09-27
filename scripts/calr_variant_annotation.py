#!/usr/bin/env python3
"""Create a structured CALR canonical-variant evidence table."""
import csv
from pathlib import Path

ROWS = [
    {
        "class": "Type 1",
        "legacy_hgvs": "NM_004343.3:c.1092_1143del52",
        "mane_hgvs": "NM_004343.4:c.1099_1150del",
        "protein": "p.Leu367fs (commonly p.L367fs*46)",
        "net_change_bp": -52,
        "consequence": "frameshift",
        "context": "CALR exon 9; canonical MPN-associated variant",
    },
    {
        "class": "Type 2",
        "legacy_hgvs": "",
        "mane_hgvs": "NM_004343.4:c.1154_1155insTTGTC",
        "protein": "p.Lys385fs (commonly p.K385fs*47)",
        "net_change_bp": 5,
        "consequence": "frameshift",
        "context": "CALR exon 9; canonical MPN-associated variant",
    },
]

out = Path(__file__).resolve().parents[1] / "results" / "tables" / "calr_canonical_variants.tsv"
out.parent.mkdir(parents=True, exist_ok=True)
with out.open("w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=ROWS[0].keys(), delimiter="\t")
    w.writeheader(); w.writerows(ROWS)
print(out)
