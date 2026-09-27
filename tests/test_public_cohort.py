import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("analyze_public_cohort", ROOT / "scripts" / "analyze_public_cohort.py")
mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)

def test_calr_identifier_resolution():
    assert mod.resolve_calr(["ENSG00000179218", "ENSG00000123456"]) == ["ENSG00000179218"]

def test_manifest_has_expected_gse156336_groups():
    import pandas as pd
    m = pd.read_csv(ROOT / "data/public/GSE156336/manifest.tsv", sep="\t")
    assert len(m) == 22
    assert m.group.value_counts().to_dict() == {"Triple_negative": 12, "Healthy": 4, "CALR": 3, "JAK2": 3}


def test_mean_difference_ci_direction_and_order():
    difference, low, high = mod.mean_difference_ci(
        [100, 110, 120],
        [200, 210, 220],
    )

    assert difference == -100.0
    assert low < difference < high


def test_cohens_d_direction():
    d = mod.cohens_d(
        [100, 110, 120],
        [200, 210, 220],
    )

    assert d < 0


def test_exploratory_comparisons_groups():
    import pandas as pd

    result = pd.DataFrame(
        {
            "group": [
                "CALR", "CALR", "CALR",
                "Healthy", "Healthy", "Healthy",
                "JAK2", "JAK2", "JAK2",
            ],
            "CALR_expression": [
                100, 110, 120,
                200, 210, 220,
                150, 160, 170,
            ],
        }
    )

    comparisons = mod.exploratory_comparisons(result)

    assert set(comparisons["comparison_group"]) == {
        "Healthy",
        "JAK2",
    }

    assert all(
        comparisons["reference_group"] == "CALR"
    )
