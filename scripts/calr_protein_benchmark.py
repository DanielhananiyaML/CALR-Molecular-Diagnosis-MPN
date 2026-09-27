"""Benchmark CALR reconstructed protein consequences against published sequence features.

This module deliberately separates published benchmark evidence from runtime RefSeq
reconstruction. It does not claim clinical interpretation or substitute for sequence
retrieval from NM_004343.4.
"""

BENCHMARKS = {
    "Type_1": {
        "protein_change": "p.L367fs*46",
        "expected_length": 411,
        "c_terminal": "RMRRMRRTRRKMRRKMSPARPRTSCREACLQGWTEA",
    },
    "Type_2": {
        "protein_change": "p.K385fs*47",
        "expected_length": 430,
        "c_terminal": "RRMMRTKMRMRRMRRTRRKMRRKMSPARPRTSCREACLQGWTEA",
    },
}


def validate_benchmark(variant, protein_sequence):
    """Return simple reproducibility checks against published mutant benchmarks."""
    if variant not in BENCHMARKS:
        raise ValueError(f"Unknown benchmark variant: {variant}")
    expected = BENCHMARKS[variant]
    seq = protein_sequence.strip().upper()
    return {
        "variant": variant,
        "protein_length": len(seq),
        "expected_length": expected["expected_length"],
        "length_match": len(seq) == expected["expected_length"],
        "c_terminal_match": seq.endswith(expected["c_terminal"]),
        "protein_change": expected["protein_change"],
    }
