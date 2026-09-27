import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.calr_protein_benchmark import BENCHMARKS, validate_benchmark


def test_type1_benchmark():
    b = BENCHMARKS["Type_1"]
    seq = "X" * (b["expected_length"] - len(b["c_terminal"])) + b["c_terminal"]
    result = validate_benchmark("Type_1", seq)
    assert result["length_match"]
    assert result["c_terminal_match"]


def test_type2_benchmark():
    b = BENCHMARKS["Type_2"]
    seq = "X" * (b["expected_length"] - len(b["c_terminal"])) + b["c_terminal"]
    result = validate_benchmark("Type_2", seq)
    assert result["length_match"]
    assert result["c_terminal_match"]
