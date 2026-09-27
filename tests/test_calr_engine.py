import sys
from pathlib import Path
import unittest
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from calr_indel_analysis import analyse

class TestCALREngine(unittest.TestCase):
    def test_type1_legacy(self):
        r = analyse("NM_004343.3:c.1092_1143del52")
        self.assertEqual(r["length"], 52)
        self.assertEqual(r["net_change"], -52)
        self.assertTrue(r["frameshift"])
        self.assertEqual(r["classification"], "Type 1")

    def test_type1_mane_alias(self):
        r = analyse("NM_004343.4:c.1099_1150del")
        self.assertEqual(r["classification"], "Type 1")

    def test_type2(self):
        r = analyse("NM_004343.4:c.1154_1155insTTGTC")
        self.assertEqual(r["length"], 5)
        self.assertTrue(r["frameshift"])
        self.assertEqual(r["classification"], "Type 2")

    def test_in_frame_is_not_type1_or_type2(self):
        r = analyse("NM_004343.4:c.1132_1134del")
        self.assertEqual(r["length"], 3)
        self.assertFalse(r["frameshift"])
        self.assertEqual(r["classification"], "Non-classical in-frame indel")

if __name__ == "__main__":
    unittest.main()
