import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from calr_consequence import parse, apply_coding_event

def test_parse_type2():
    e=parse('NM_004343.4:c.1154_1155insTTGTC')
    assert e.kind=='ins' and e.sequence=='TTGTC'

def test_deletion_changes_length_by_52():
    e=parse('NM_004343.4:c.1099_1150del')
    seq='A'*1200
    assert len(apply_coding_event(seq,e))==1148

def test_insertion_changes_length_by_5():
    e=parse('NM_004343.4:c.1154_1155insTTGTC')
    seq='A'*1200
    assert len(apply_coding_event(seq,e))==1205
