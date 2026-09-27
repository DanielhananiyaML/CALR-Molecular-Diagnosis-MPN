import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from calr_genbank_reference import extract_origin, parse_record, cds_sequence


FIXTURE = '''LOCUS       NM_TEST                 30 bp    mRNA    linear   PRI 01-JAN-2026
FEATURES             Location/Qualifiers
     source          1..30
                     /organism="Homo sapiens"
     exon            1..10
                     /number="1"
     exon            11..30
                     /number="2"
     CDS             4..24
                     /gene="CALR"
                     /product="calreticulin"
ORIGIN
        1 aaaATGAAACCCGGGTTTAAACCCggg
//
'''


def test_origin_extraction():
    assert extract_origin(FIXTURE) == 'AAAATGAAACCCGGGTTTAAACCCGGG'


def test_feature_and_cds_parsing():
    parsed = parse_record(FIXTURE)
    assert parsed['cds']['location']['start'] == 4
    assert parsed['cds']['location']['end'] == 24
    assert len(parsed['exons']) == 2
    assert parsed['exons'][1]['number'] == '2'


def test_cds_sequence_extraction():
    parsed = parse_record(FIXTURE)
    assert cds_sequence(FIXTURE, parsed['cds']) == 'ATGAAACCCGGGTTTAAACCC'
