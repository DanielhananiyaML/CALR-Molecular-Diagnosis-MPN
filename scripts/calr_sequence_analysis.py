#!/usr/bin/env python3
"""Reference metadata and molecular interpretation report for CALR."""
from pathlib import Path
import json

REPORT = {
    "gene": "CALR",
    "transcript": "NM_004343.4",
    "ensembl": "ENST00000316448.10",
    "protein": "NP_004334.1 / UniProt P27797",
    "chromosome": "19",
    "assembly": "GRCh38",
    "focus": "exon 9 indels",
    "canonical_type1_current_mane": "NM_004343.4:c.1099_1150del",
    "canonical_type1_legacy": "NM_004343.3:c.1092_1143del52",
    "canonical_type2": "NM_004343.4:c.1154_1155insTTGTC",
    "note": "Historical and current transcript versions must be retained together for literature/database reconciliation."
}
out = Path(__file__).resolve().parents[1] / "results" / "tables" / "calr_reference_metadata.json"
out.parent.mkdir(parents=True, exist_ok=True)
out.write_text(json.dumps(REPORT, indent=2) + "\n")
print(out)
