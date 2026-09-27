import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from calr_sequence_sources import EXPECTED_CDS_LENGTH

def test_expected_mane_cds_length():
    assert EXPECTED_CDS_LENGTH == 1254
