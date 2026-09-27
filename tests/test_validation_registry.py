from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from validate_cohort_registry import validate

def test_registry_valid():
    validate("data/public/validation/cohort_registry.tsv")
