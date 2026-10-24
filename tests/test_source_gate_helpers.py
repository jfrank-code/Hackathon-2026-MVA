import importlib.util
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / 'src' / 'pivot_mva_source_gate.py'
spec = importlib.util.spec_from_file_location('pivot_mva_source_gate', SCRIPT)
mod = importlib.util.module_from_spec(spec)
assert spec and spec.loader
spec.loader.exec_module(mod)

class SourceGateHelpersTest(unittest.TestCase):
    def test_contig_aliases(self):
        self.assertEqual(mod.contig_aliases('chr15', {'15'}), '15')
        self.assertEqual(mod.contig_aliases('15', {'chr15'}), 'chr15')
        self.assertEqual(mod.contig_aliases('chr1', {'chr1'}), 'chr1')
        self.assertIsNone(mod.contig_aliases('chr20', {'15', 'chr1'}))

    def test_genotype_alt_state(self):
        self.assertEqual(mod.genotype_alt_state('0/1', 1), 'ALT_CALLED')
        self.assertEqual(mod.genotype_alt_state('0|2', 2), 'ALT_CALLED')
        self.assertEqual(mod.genotype_alt_state('0/0', 1), 'REF_ONLY')
        self.assertEqual(mod.genotype_alt_state('./.', 1), 'MISSING_GT')
        self.assertEqual(mod.genotype_alt_state('0/1', 2), 'OTHER_ALT_CALLED')

    def test_ad_for_alt(self):
        self.assertEqual(mod.ad_for_alt('10,8,3', 2), (10, 3))
        self.assertEqual(mod.ad_for_alt('15,12', 1), (15, 12))
        self.assertEqual(mod.ad_for_alt('.', 1), (None, None))

if __name__ == '__main__':
    unittest.main()
