from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from calr_protein_reconstruction import read_fasta, reconstruct, VARIANTS, REFERENCE

def test_reference_length():
    assert len(read_fasta(REFERENCE)) == 417

def test_type1_reconstruction_length():
    seq=read_fasta(REFERENCE); assert len(reconstruct(seq, VARIANTS['Type 1'])) == 411

def test_type2_reconstruction_length():
    seq=read_fasta(REFERENCE); assert len(reconstruct(seq, VARIANTS['Type 2'])) == 430

def test_type1_tail():
    seq=read_fasta(REFERENCE); assert reconstruct(seq, VARIANTS['Type 1']).endswith(VARIANTS['Type 1']['tail'])

def test_type2_tail():
    seq=read_fasta(REFERENCE); assert reconstruct(seq, VARIANTS['Type 2']).endswith(VARIANTS['Type 2']['tail'])
