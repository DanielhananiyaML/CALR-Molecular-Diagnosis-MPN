import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts.calr_evidence import load_evidence

def test_evidence_register_has_core_references():
    rows = load_evidence()
    ids = {r["id"] for r in rows}
    assert {"NCBI-CALR-REFSEQ", "CALR-TYPE1", "CALR-TYPE2", "CALR-MPN"}.issubset(ids)

def test_evidence_register_is_complete():
    rows = load_evidence()
    assert all(r["claim"].strip() for r in rows)
    assert all(r["interpretation_guardrail"].strip() for r in rows)
